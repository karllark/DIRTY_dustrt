// ======================================================================
//   Procedure to allocate and initialize the output stucture
//
// 2005 May/KDG - written
// 2008 Jun/KDG - moved computation of the transformation(rotation) matrix
//                to a subroutine (need to use it elsewhere)
// ======================================================================
#include "initialize_output.h"

void initialize_output (output_struct& output,
			geometry_struct& geometry)

{
  int i;

  // allocate outputs if not already allocated
  if (!output.arrays_allocated) {
    one_output single_output;
    single_output.num_photons_xy.MSize(output.image_size[0],output.image_size[1]);
    single_output.stellar_weight_xy.MSize(output.image_size[0],output.image_size[1]);
    single_output.stellar_weight_xy_x2.MSize(output.image_size[0],output.image_size[1]);
    single_output.scattered_weight_xy.MSize(output.image_size[0],output.image_size[1]);
    single_output.scattered_weight_xy_x2.MSize(output.image_size[0],output.image_size[1]);
    for (i = 0; i < output.num_outputs; i++) {
      output.outputs.push_back(single_output);
    }
    output.arrays_allocated = 1;
  }

  for (i = 0; i < output.num_outputs; i++) {
    // compute the transformation matrix for rotating positions
    // so the observer is along the z-axis
    compute_observer_trans_matrix(output, geometry, i);

    // zero output quantities
    output.outputs[i].total_num_photons = 0.0;
    output.outputs[i].total_num_scattered_photons = 0.0;
    output.outputs[i].total_stellar_weight = 0.0;
    output.outputs[i].total_stellar_weight_x2 = 0.0;
    output.outputs[i].total_scattered_weight = 0.0;
    output.outputs[i].total_scattered_weight_x2 = 0.0;
    output.outputs[i].ave_first_tau = 0.0;
    output.outputs[i].ave_first_tau_x2 = 0.0;

    int j,k;
    for (j = 0; j < output.image_size[0]; j++)
      for (k = 0; k < output.image_size[1]; k++) {
	output.outputs[i].num_photons_xy(j,k) = 0.0;
	output.outputs[i].stellar_weight_xy(j,k) = 0.0;
	output.outputs[i].stellar_weight_xy_x2(j,k) = 0.0;
	output.outputs[i].scattered_weight_xy(j,k) = 0.0;
	output.outputs[i].scattered_weight_xy_x2(j,k) = 0.0;
      }
    
  }

}
