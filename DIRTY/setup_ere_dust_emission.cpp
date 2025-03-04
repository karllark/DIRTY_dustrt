// ======================================================================
//   Procedure for the setup of the output of the ere dust emission.
//
// 2008 Sep/KDG - written
// ======================================================================

#include "setup_ere_dust_emission_output.h"

void setup_ere_dust_emission_output(output_struct &ere_output,
                                    output_struct &output)

{
  int i = 0;
  for (i = 0; i < 2; i++)
    ere_output.image_size[i] = output.image_size[i];
  ere_output.file_base = output.file_base;
  ere_output.emission_type = "_ere"; // setup the emission_type string for ere
  ere_output.type = output.type;
  ere_output.num_outputs = 1;
  ere_output.arrays_allocated = 0;
}
