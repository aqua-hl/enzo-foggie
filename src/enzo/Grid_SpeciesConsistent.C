/***********************************************************************
/
/  GRID CLASS (METAL_DENSITY IS THE SUM OF THE SIX GAS-PHASE METAL RESERVOIRS)
/
/  written by: Harrison Lo
/  date:       October 2026
/  modified1:  November 2026. Metal_Density is derived from the six
/              reservoirs; they are no longer rescaled to it.
/
/  PURPOSE: With UseDustSpeciesTrack the gas-phase metals are carried
/    twice: as Metal_Density and as the five element fields plus
/    Metal_Density_Other (N, Ne, S, ...). Grackle's dust update infers the
/    untracked metals as max(M - E, 0), so a cell whose elements E exceed
/    M = Metal_Density would create metal. The six reservoirs are the
/    record: before each Grackle call M is set to their sum, so the
/    elements can never exceed it and no reservoir is rescaled. M is then
/    capped at 90% of the gas density, as GrackleWrapper caps every metal
/    field, and where that cap bites the six are scaled down with M.
/
/    The PPM sweeps keep M equal to the sum themselves
/    (ConsistentMetalColourFluxes.C), so this call covers what acts
/    between the last sweep and the chemistry and should change little.
/    A large change means some other operation moves Metal_Density
/    without the six, and that operation is what needs fixing.
/
/    A negative reservoir is zeroed and counted; a cell with a non-finite
/    metal field is left alone and counted.
/
/    The changes are accumulated per process (density x cell volume in
/    code units, over every grid including cells covered by subgrids) and
/    written to SpeciesConsistent.out by SpeciesConsistentReport at each
/    data dump, with the counters of the sweeps and of the flux
/    correction (grid::SettleGrainDeficit).
/
************************************************************************/

#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int FindField(int field, int farray[], int numfields);
int CommunicationSumValues(double *Values, int Number);
void ConsistentMetalCounters(double values[6]);
void GrainDeficitCounters(double values[5]);

/* Since the last report: 0 cells seen; 1 change of Metal_Density by the
   sum; 2 negative inputs zeroed; 3 non-finite cells skipped; 4 cells
   capped at 90% of the gas density; 5 change of Metal_Density by the cap. */
static double SCAccum[6] = {0, 0, 0, 0, 0, 0};
static int SCWarned = FALSE;

int grid::SpeciesConsistent(float *MetalPointer)
{
  if (ProcessorNumber != MyProcessorNumber || !UseDustSpeciesTrack ||
      MetalPointer == NULL)
    return SUCCESS;

  int types[6] = {MetalDensityCarbon, MetalDensityOxygen,
                  MetalDensityMagnesium, MetalDensitySilicon,
                  MetalDensityIron, MetalDensityOther};
  float *f[6];
  int q, num;
  for (q = 0; q < 6; q++) {
    num = FindField(types[q], FieldType, NumberOfBaryonFields);
    if (num < 0)
      ENZO_FAIL("SpeciesConsistent: a tracked metal field is missing.\n");
    f[q] = BaryonField[num];
  }
  int DensNum = FindField(Density, FieldType, NumberOfBaryonFields);
  if (DensNum < 0)
    ENZO_FAIL("SpeciesConsistent: the density field is missing.\n");

  double vol = 1.0;
  for (int dim = 0; dim < GridRank; dim++)
    vol *= CellWidth[dim][0];

  int i, j, k, n;
  double M, S, cap, scale, v[6];
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
        n = (k*GridDimension[1] + j)*GridDimension[0] + i;
        M = MetalPointer[n];
        int bad = !isfinite(M);
        for (q = 0; q < 6; q++) {
          v[q] = f[q][n];
          if (!isfinite(v[q])) bad = TRUE;
        }
        if (bad) {
          SCAccum[3] += 1;
          if (!SCWarned) {
            fprintf(stderr, "SpeciesConsistent: non-finite metal field in "
                    "grid %"ISYM" cell (%"ISYM",%"ISYM",%"ISYM"), left alone\n",
                    ID, i, j, k);
            SCWarned = TRUE;
          }
          continue;
        }
        S = 0;
        for (q = 0; q < 6; q++) {
          if (v[q] < 0) {
            v[q] = 0;
            SCAccum[2] += 1;
          }
          S += v[q];
        }

        MetalPointer[n] = (float) S;
        SCAccum[1] += (MetalPointer[n] - M) * vol;
        for (q = 0; q < 6; q++)
          f[q][n] = (float) v[q];

        cap = 0.9 * BaryonField[DensNum][n];
        if (S > cap) {
          scale = cap / S;
          for (q = 0; q < 6; q++)
            f[q][n] *= scale;
          MetalPointer[n] = (float) cap;
          SCAccum[4] += 1;
          SCAccum[5] += (cap - S) * vol;
        }
        SCAccum[0] += 1;
      }

  return SUCCESS;
}

/* Restart from a dump written before Metal_Density_Other existed: add the
   field on every processor (AddFields allocates only on the owner, and the
   boundary gets it in Group_ReadAllData) and seed it on the owner with
   max(metal - elements, 0), so the six gas-metal reservoirs partition
   Metal_Density from the restart on. Any mismatch already in the dump is
   frozen into it here. Both Group_ReadGrid readers call this. */

int grid::AddOtherMetalField()
{
  if (!UseDustSpeciesTrack || NumberOfBaryonFields <= 0 ||
      FindField(MetalDensityOther, FieldType, NumberOfBaryonFields) >= 0)
    return SUCCESS;

  int ElType[5] = {MetalDensityCarbon, MetalDensityOxygen,
                   MetalDensityMagnesium, MetalDensitySilicon,
                   MetalDensityIron};
  int ElNum[5], q;
  int MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
  if (MetalNum < 0)
    return SUCCESS;
  for (q = 0; q < 5; q++) {
    ElNum[q] = FindField(ElType[q], FieldType, NumberOfBaryonFields);
    if (ElNum[q] < 0)
      return SUCCESS;
  }

  char *OtherName = "Metal_Density_Other";
  int OtherNum = NumberOfBaryonFields;
  int OtherToAdd = MetalDensityOther;
  this->AddFields(&OtherToAdd, 1);
  if (DataLabel[OtherNum] == NULL) {
    if (debug)
      fprintf(stdout, "Adding Metal_Density_Other field.\n");
    DataLabel[OtherNum] = OtherName;
  }

  if (MyProcessorNumber == ProcessorNumber &&
      BaryonField[OtherNum] != NULL && BaryonField[MetalNum] != NULL) {
    int size = 1, n;
    for (int dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
    for (n = 0; n < size; n++) {
      float e = 0;
      for (q = 0; q < 5; q++)
        e += BaryonField[ElNum[q]][n];
      BaryonField[OtherNum][n] = max(BaryonField[MetalNum][n] - e,
                                     (float) 0.0);
    }
  }
  return SUCCESS;
}

/* Sum the counters over all processors, append them to
   SpeciesConsistent.out on the root processor and reset them. Collective:
   every processor must call it (done at each data dump). */

int SpeciesConsistentReport(FLOAT Time, int CycleNumber)
{
  if (!UseDustSpeciesTrack)
    return SUCCESS;

  double total[17];
  int q;
  for (q = 0; q < 6; q++)
    total[q] = SCAccum[q];
  ConsistentMetalCounters(total + 6);
  GrainDeficitCounters(total + 12);
  CommunicationSumValues(total, 17);

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *fp = fopen("SpeciesConsistent.out", "a");
    if (fp != NULL) {
      if (ftell(fp) == 0)
        fprintf(fp, "# Metal_Density = sum of its six components, since the "
                "previous line. Masses are code units (density x cell "
                "volume) over all grids, including covered cells.\n"
                "# Before Grackle (grid::SpeciesConsistent): n_cells, "
                "dMetal by the sum, negative components zeroed, non-finite "
                "cells skipped, cells capped at 0.9 of the gas, dMetal by "
                "the cap.\n"
                "# In the sweeps (ConsistentMetalColourFluxes.C): net and "
                "absolute dMetal by the sum, fraction profiles built, "
                "profiles steepened, negative face fractions, stencil "
                "cells without metal.\n"
                "# In the flux correction (grid::SettleGrainDeficit): grain "
                "fields found at or below zero, their deficit, mass the "
                "gas phase paid, fields left at the floor, mass created "
                "in metals plus dust.\n"
                "# cycle time n_cells dMetal n_negative n_nonfinite "
                "n_capped dMetal_cap sweep_dMetal sweep_abs_dMetal "
                "n_profiles n_steepened n_negative_faces n_zero_metal "
                "reflux_n_deficit reflux_deficit reflux_paid "
                "reflux_n_floored reflux_created\n");
      fprintf(fp, "%"ISYM" %.10"GSYM" %.0f %.6e %.0f %.0f %.0f %.6e "
              "%.6e %.6e %.0f %.0f %.0f %.0f %.0f %.6e %.6e %.0f %.6e\n",
              CycleNumber, Time,
              total[0], total[1], total[2], total[3], total[4], total[5],
              total[6], total[7], total[8], total[9], total[10], total[11],
              total[12], total[13], total[14], total[15], total[16]);
      fclose(fp);
    }
  }

  for (q = 0; q < 6; q++)
    SCAccum[q] = 0;
  return SUCCESS;
}
