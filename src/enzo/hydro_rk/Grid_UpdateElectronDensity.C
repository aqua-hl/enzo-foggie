/***********************************************************************
/
/  GRID CLASS (CALCULATE ELECTRON DENSITY FROM SPECIES)
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"

int FindField(int field, int farray[], int numfields);

int grid::UpdateElectronDensity(void)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  if (MultiSpecies == 0)
    return SUCCESS;

  int i, n, dim, size, nfield, n0;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
			      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

  for (dim = 0, size = 1; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  for (i = 0; i < size; i++)
    BaryonField[DeNum][i] = BaryonField[HIINum][i] + 0.25*BaryonField[HeIINum][i] +
      0.5*BaryonField[HeIIINum][i]; 

  if (MultiSpecies > 1)
    for (i = 0; i < size; i++)
      BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] - BaryonField[HMNum][i];

  if (MultiSpecies > 2)
    for (i = 0; i < size; i++)
      BaryonField[DeNum][i] += 0.5*BaryonField[DIINum][i];

  /* Grackle metal chemistry: singly ionized metal species, weighted by
     1/mass number as in Grackle's make_consistent. */
  if (UseMetalChemistry) {
    int mc[NUM_METAL_CHEMISTRY_SPECIES];
    if (IdentifyMetalChemistrySpeciesFields(mc) == FAIL)
      ENZO_FAIL("Error in grid->IdentifyMetalChemistrySpeciesFields.\n");
    for (i = 0; i < size; i++)
      BaryonField[DeNum][i] +=
        BaryonField[mc[MC_CII]][i]   / 12.0 +
        BaryonField[mc[MC_COII]][i]  / 28.0 +
        BaryonField[mc[MC_OII]][i]   / 16.0 +
        BaryonField[mc[MC_OHII]][i]  / 17.0 +
        BaryonField[mc[MC_H2OII]][i] / 18.0 +
        BaryonField[mc[MC_H3OII]][i] / 19.0 +
        BaryonField[mc[MC_O2II]][i]  / 32.0;
  }

  return SUCCESS;

}
