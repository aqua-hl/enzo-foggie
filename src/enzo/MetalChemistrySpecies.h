/***********************************************************************
/
/  GRACKLE METAL CHEMISTRY SPECIES (metal_chemistry = 1)
/
/  PURPOSE: Ordered list of the 19 metal species evolved by Grackle's
/           non-equilibrium C/O/Si network, in Grackle's own order
/           (src/clib/field_data_evolved_species.def). Eleven species
/           reuse the field types originally defined for Simon Glover's
/           network in typedefs.h; the remaining eight were added for
/           Grackle. Enabled in Enzo with UseMetalChemistry = 1.
/
/           The output labels follow the existing Glover convention for
/           the neutral species (e.g. COI_Density) so that the names
/           already known to yt are kept.
/
************************************************************************/

#ifndef METAL_CHEMISTRY_SPECIES_DEFINED__
#define METAL_CHEMISTRY_SPECIES_DEFINED__

#define NUM_METAL_CHEMISTRY_SPECIES 19

/* Index of each species in the arrays below and in the field-number
   arrays filled by grid::IdentifyMetalChemistrySpeciesFields. */
enum metal_chemistry_species {
  MC_CI = 0, MC_CII, MC_CO, MC_CO2, MC_OI, MC_OH, MC_H2O, MC_O2,
  MC_SiI, MC_SiOI, MC_SiO2I, MC_CH, MC_CH2, MC_COII, MC_OII, MC_OHII,
  MC_H2OII, MC_H3OII, MC_O2II
};

static const int MetalChemistrySpeciesTypes[NUM_METAL_CHEMISTRY_SPECIES] = {
  CIDensity,    CIIDensity,   COIDensity,   CO2IDensity,  OIDensity,
  OHIDensity,   H2OIDensity,  O2IDensity,   SiIDensity,   SiOIDensity,
  SiO2IDensity, CHIDensity,   CH2IDensity,  COIIDensity,  OIIDensity,
  OHIIDensity,  H2OIIDensity, H3OIIDensity, O2IIDensity
};

static const char *const MetalChemistrySpeciesLabels[NUM_METAL_CHEMISTRY_SPECIES] = {
  "CI_Density",    "CII_Density",   "COI_Density",   "CO2I_Density",  "OI_Density",
  "OHI_Density",   "H2OI_Density",  "O2I_Density",   "SiI_Density",   "SiOI_Density",
  "SiO2I_Density", "CHI_Density",   "CH2I_Density",  "COII_Density",  "OII_Density",
  "OHII_Density",  "H2OII_Density", "H3OII_Density", "O2II_Density"
};

#endif /* METAL_CHEMISTRY_SPECIES_DEFINED__ */
