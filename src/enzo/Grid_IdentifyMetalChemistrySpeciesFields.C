/***********************************************************************
/
/  GRID CLASS (IDENTIFY THE GRACKLE METAL CHEMISTRY SPECIES FIELDS)
/
/  PURPOSE: Locate the 19 metal species fields evolved by Grackle's
/           non-equilibrium metal network (metal_chemistry = 1). The
/           species order is fixed by MetalChemistrySpecies.h. Fails if
/           UseMetalChemistry is set and any species field is missing;
/           otherwise missing fields are reported as -1.
/
/  RETURNS:
/    SUCCESS or FAIL
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

/* function prototypes */

int FindField(int f, int farray[], int n);

int grid::IdentifyMetalChemistrySpeciesFields(int SpeciesNum[NUM_METAL_CHEMISTRY_SPECIES])
{

  int isp, missing = FALSE;

  for (isp = 0; isp < NUM_METAL_CHEMISTRY_SPECIES; isp++) {
    SpeciesNum[isp] = FindField(MetalChemistrySpeciesTypes[isp], FieldType,
                                NumberOfBaryonFields);
    if (SpeciesNum[isp] < 0)
      missing = TRUE;
  }

  if (UseMetalChemistry && missing) {
    fprintf(stderr, "UseMetalChemistry = 1 but metal chemistry species "
            "field(s) missing:");
    for (isp = 0; isp < NUM_METAL_CHEMISTRY_SPECIES; isp++)
      if (SpeciesNum[isp] < 0)
        fprintf(stderr, " %s", MetalChemistrySpeciesLabels[isp]);
    fprintf(stderr, "\n");
    ENZO_FAIL("Missing metal chemistry species field.\n");
  }

  return SUCCESS;
}
