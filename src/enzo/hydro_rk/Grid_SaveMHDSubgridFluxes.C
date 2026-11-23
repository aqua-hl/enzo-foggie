/***********************************************************************
/
/  GRID CLASS (SAVE SUBGRID BOUNDARY FLUXES)
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::SaveMHDSubgridFluxes(fluxes *SubgridFluxes[], int NumberOfSubgrids,
			       float *Flux3D[], int flux, float fluxcoef, float dt)
{

  return this->SaveSubgridFluxes(SubgridFluxes, NumberOfSubgrids, Flux3D,
                                flux, fluxcoef, dt);
}
