/***********************************************************************
/
/  ELEMENT MASS FRACTIONS OF THE SILICATE GRAINS
/
/  written by: Harrison Lo
/  date:       November 2026
/
/  PURPOSE: Returns the mass fractions of Mg, Fe, Si and O, in that order,
/    in Mg-silicate (Mg2SiO4) and Fe-silicate (Fe2SiO4) grains. With the
/    species-resolved dust model they are the ones Grackle's grain growth
/    and destruction use, so every Enzo routine that moves mass between
/    the grain fields and the gas-phase element fields agrees with Grackle
/    element by element. Carbonaceous grains are pure carbon.
/
************************************************************************/

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

void DustSilicateComposition(float MgSilicate[4], float FeSilicate[4])
{
  MgSilicate[0] = 0.345521;  MgSilicate[1] = 0.0;
  MgSilicate[2] = 0.199630;  MgSilicate[3] = 0.454849;
  FeSilicate[0] = 0.0;       FeSilicate[1] = 0.548111;
  FeSilicate[2] = 0.137825;  FeSilicate[3] = 0.314064;
#if defined(USE_GRACKLE) && defined(GRACKLE_NEW_DUST_MODEL)
  if (grackle_data != NULL) {
    MgSilicate[0] = (float) grackle_data->dust_mg_silicate_f_Mg;
    MgSilicate[1] = (float) grackle_data->dust_mg_silicate_f_Fe;
    MgSilicate[2] = (float) grackle_data->dust_mg_silicate_f_Si;
    MgSilicate[3] = (float) grackle_data->dust_mg_silicate_f_O;
    FeSilicate[0] = (float) grackle_data->dust_fe_silicate_f_Mg;
    FeSilicate[1] = (float) grackle_data->dust_fe_silicate_f_Fe;
    FeSilicate[2] = (float) grackle_data->dust_fe_silicate_f_Si;
    FeSilicate[3] = (float) grackle_data->dust_fe_silicate_f_O;
  }
#endif
}
