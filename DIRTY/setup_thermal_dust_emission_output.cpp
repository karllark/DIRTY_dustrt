// ======================================================================
//   Procedure for the setup of the output of the thermal dust emission
// when we are keeping track of the different grain/emission types.
//
// 2008 Mar/KDG - written
// ======================================================================

#include "setup_thermal_dust_emission_output.h"

void setup_thermal_dust_emission_output(runinfo_struct &runinfo,
                                        output_struct &de_output,
                                        output_struct &output,
                                        photon_data &photon)

{
  int i = 0;
  for (i = 0; i < 2; i++)
    de_output.image_size[i] = output.image_size[i];
  de_output.file_base = output.file_base;
  de_output.emission_type =
      "_de"; // setup the emission_type string for dust emission
  de_output.type = output.type;
  de_output.num_outputs = runinfo.n_emission_grain_types;
  de_output.arrays_allocated = 0;

  // now allocate the space needed in the photon structure to save the
  // probabilities that the emitted photon is of a particular grain/emission
  // type.
  for (i = 0; i < de_output.num_outputs; i++)
    photon.birth_photon_type_prob.push_back(0.0);
}
