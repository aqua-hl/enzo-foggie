/***********************************************************************
/
/  GRID CLASS (A GRAIN DEFICIT LEFT BY THE FLUX CORRECTION IS PAID BY THE
/              GAS-PHASE METALS THE GRAINS ARE MADE OF)
/
/  written by: Harrison Lo
/  date:       November 2026
/
/  PURPOSE: With FluxCorrection = 2 each colour is corrected with its own
/    flux difference, and the correction can take more of a field out of a
/    coarse cell than the cell holds: the coarse step put mass into the
/    cell that the fine steps say never crossed the face, and by the time
/    of the correction part of it has moved on to other cells or been
/    changed by the chemistry. Setting such a field to the colour floor
/    adds the missing mass.
/
/    For the grains that mass need not be created. Each element is
/    conserved as gas plus grains, and Grackle moves mass between a grain
/    field and the gas-phase fields of its elements (and with them
/    Metal_Density and the gas density, which holds the gas-phase metals
/    but not the dust). So a cell that owes more grains than it holds pays
/    the rest from the gas phase:
/
/      1. A grain field left at or below zero is set to the colour floor.
/         The mass that takes is its deficit.
/      2. Each gas-phase element gives its share of the deficit, by the
/         composition of the grain: carbonaceous from C; Mg-silicate from
/         Mg, Si and O; Fe-silicate from Fe, Si and O. An element that
/         cannot give its share is left at the floor, and only that
/         remainder is created.
/      3. Metal_Density is set to the sum of its six components, and the
/         gas density loses the mass the elements gave.
/
/    This is exact where the grains are missing because the chemistry
/    returned them to the gas after the hydro step. Elsewhere it moves
/    the missing mass from the gas phase to the grains instead of
/    creating it.
/
/    A component that its own correction left at or below zero is floored
/    in the same pass. With the bulk dust model (one Dust_Density field,
/    no species) the deficit is paid by Metal_Density.
/
/    Called by grid::CorrectForRefinedFluxes on the cells it corrected,
/    after the faces are back in primitive form, so the velocity and the
/    specific energy of the cell are kept, as in Grackle's dust update.
/
/    Counters are kept per process (density x cell volume in code units)
/    and reported in SpeciesConsistent.out by SpeciesConsistentReport.
/    They include cells covered by subgrids, which the projection then
/    overwrites, so they are upper limits for the leaf cells.
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
void DustSilicateComposition(float MgSilicate[4], float FeSilicate[4]);

/* Since the last report: 0 grain fields found at or below zero; 1 their
   deficit; 2 mass paid by the gas phase; 3 fields left at the floor
   because the gas phase could not pay (or their own correction emptied
   them); 4 mass created in metals plus dust. */
static double SGAccum[5] = {0, 0, 0, 0, 0};

void GrainDeficitCounters(double values[5])
{
  for (int q = 0; q < 5; q++) {
    values[q] = SGAccum[q];
    SGAccum[q] = 0;
  }
}

int grid::SettleGrainDeficit(int Start[], int End[], int Offset,
                             int CorrectLeft, int CorrectRight, float Floor)
{
  if (ProcessorNumber != MyProcessorNumber || NumberOfBaryonFields <= 0 ||
      (!CorrectLeft && !CorrectRight))
    return SUCCESS;

  /* The grain fields, the gas-phase fields that pay for them, and the
     mass fraction of each grain that each payer gives. */

  int GrainNum[3], PayerNum[6], NumberOfGrains, NumberOfPayers, g, p;
  float Share[3][6];
  for (g = 0; g < 3; g++)
    for (p = 0; p < 6; p++)
      Share[g][p] = 0;

  int DensNum = FindField(Density, FieldType, NumberOfBaryonFields);
  int MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);

  if (UseDustSpeciesTrack) {
    int GrainType[3] = {DustDensityMgSilicate, DustDensityFeSilicate,
                        DustDensityCarbonaceous};
    int PayerType[6] = {MetalDensityCarbon, MetalDensityOxygen,
                        MetalDensityMagnesium, MetalDensitySilicon,
                        MetalDensityIron, MetalDensityOther};
    NumberOfGrains = 3;
    NumberOfPayers = 6;
    for (g = 0; g < 3; g++)
      GrainNum[g] = FindField(GrainType[g], FieldType, NumberOfBaryonFields);
    for (p = 0; p < 6; p++)
      PayerNum[p] = FindField(PayerType[p], FieldType, NumberOfBaryonFields);
    for (g = 0; g < 3; g++)
      if (GrainNum[g] < 0)
        ENZO_FAIL("SettleGrainDeficit: a dust species field is missing.\n");
    for (p = 0; p < 6; p++)
      if (PayerNum[p] < 0)
        ENZO_FAIL("SettleGrainDeficit: a tracked metal field is missing.\n");
    if (MetalNum < 0 || DensNum < 0)
      ENZO_FAIL("SettleGrainDeficit: the metal or density field is missing.\n");

    /* Payers are C, O, Mg, Si, Fe, other; the silicate fractions come as
       Mg, Fe, Si, O. The other metals are in no grain. */
    float MgSilicate[4], FeSilicate[4];
    DustSilicateComposition(MgSilicate, FeSilicate);
    Share[0][2] = MgSilicate[0];  Share[1][2] = FeSilicate[0];
    Share[0][4] = MgSilicate[1];  Share[1][4] = FeSilicate[1];
    Share[0][3] = MgSilicate[2];  Share[1][3] = FeSilicate[2];
    Share[0][1] = MgSilicate[3];  Share[1][1] = FeSilicate[3];
    Share[2][0] = 1;
  } else {
    NumberOfGrains = 1;
    NumberOfPayers = 1;
    GrainNum[0] = FindField(DustDensity, FieldType, NumberOfBaryonFields);
    PayerNum[0] = MetalNum;
    if (GrainNum[0] < 0 || MetalNum < 0 || DensNum < 0)
      return SUCCESS;
    Share[0][0] = 1;
  }

  double vol = 1.0;
  for (int dim = 0; dim < GridRank; dim++)
    vol *= CellWidth[dim][0];

  int i, j, k, n, side, index, settle;
  double deficit, need[6], value, target, paid, created;

  for (k = Start[2]; k <= End[2]; k++)
    for (j = Start[1]; j <= End[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0] + Start[0];
      for (i = Start[0]; i <= End[0]; i++, index++)
        for (side = 0; side < 2; side++) {
          if ((side == 0 && !CorrectLeft) || (side == 1 && !CorrectRight))
            continue;
          n = (side == 0) ? index : index + Offset;

          settle = (UseDustSpeciesTrack && BaryonField[MetalNum][n] <= 0);
          for (g = 0; g < NumberOfGrains; g++)
            if (BaryonField[GrainNum[g]][n] <= 0)
              settle = TRUE;
          for (p = 0; p < NumberOfPayers; p++)
            if (BaryonField[PayerNum[p]][n] <= 0)
              settle = TRUE;
          if (!settle)
            continue;

          /* 1. Grains to the floor; what that takes is owed by the gas. */

          created = 0;
          for (p = 0; p < NumberOfPayers; p++)
            need[p] = 0;
          for (g = 0; g < NumberOfGrains; g++)
            if (BaryonField[GrainNum[g]][n] <= 0) {
              deficit = Floor - BaryonField[GrainNum[g]][n];
              BaryonField[GrainNum[g]][n] = Floor;
              for (p = 0; p < NumberOfPayers; p++)
                need[p] += Share[g][p] * deficit;
              created += deficit;
              SGAccum[0] += 1;
              SGAccum[1] += deficit * vol;
            }

          /* 2. Each payer gives its share, down to the floor. */

          paid = 0;
          for (p = 0; p < NumberOfPayers; p++) {
            value = BaryonField[PayerNum[p]][n];
            target = value - need[p];
            if (target <= 0) {
              target = Floor;
              SGAccum[3] += 1;
            }
            if (target < value)
              paid += value - target;
            created += target - value;
            BaryonField[PayerNum[p]][n] = (float) target;
          }

          /* 3. The metal density follows its components, and the gas
                density loses what the gas phase paid. */

          if (UseDustSpeciesTrack) {
            value = 0;
            for (p = 0; p < NumberOfPayers; p++)
              value += BaryonField[PayerNum[p]][n];
            BaryonField[MetalNum][n] = (float) value;
          }
          if (paid > 0 && paid < BaryonField[DensNum][n])
            BaryonField[DensNum][n] -= (float) paid;

          SGAccum[2] += paid * vol;
          SGAccum[4] += created * vol;
        }
    }

  return SUCCESS;
}
