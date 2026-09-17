/***********************************************************************
/
/  GRID CLASS (CHECK THAT MORE BARYON FIELDS CAN BE ADDED)
/
/  PURPOSE: Fail with a clear message before FieldType[] / BaryonField[]
/           would be written past MAX_NUMBER_OF_BARYON_FIELDS. One slot
/           past NumberOfBaryonFields must stay free because it is used
/           as scratch space (e.g. the refinement flag field in
/           StarParticleHandler); this is the rule InitializeNew enforces
/           after the fact, applied here before the writes happen.
/
/  RETURNS:
/    SUCCESS (or throws via ENZO_VFAIL)
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

int grid::CheckBaryonFieldCapacity(int NumberToAdd, const char *Description)
{

  if (NumberOfBaryonFields + NumberToAdd >= MAX_NUMBER_OF_BARYON_FIELDS) {
    ENZO_VFAIL("Cannot add %"ISYM" baryon field(s) (%s): %"ISYM" already "
               "present plus one scratch slot exceeds "
               "MAX_NUMBER_OF_BARYON_FIELDS = %"ISYM". Rebuild with "
               "make max-baryons-###.\n",
               NumberToAdd, Description, NumberOfBaryonFields,
               (int) MAX_NUMBER_OF_BARYON_FIELDS)
  }

  return SUCCESS;
}
