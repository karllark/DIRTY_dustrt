// ======================================================================
//   Header file with output definitions.
//
// 2005 May/KDG - written
// ======================================================================

#ifndef _DIRTY_OUTPUT_DEF_
#define _DIRTY_OUTPUT_DEF_

#include "NumUtils.h"
//#include "vect_utils.h"

// structure with the output for one line-of-sight, albedo, etc.
struct one_output {
  // number of photons
  float total_num_photons; // total
  NumUtils::Matrix<float> num_photons_xy; // photons per image pixel

  // number of scattered photons
  float total_num_scattered_photons; // total

  // stellar photons
  double total_stellar_weight; // total
  double total_stellar_weight_x2; // total squared sum (for unc calc)
  NumUtils::Matrix<float> stellar_weight_xy; // weight per image pixel
  NumUtils::Matrix<float> stellar_weight_xy_x2; // squared sum weight per image pixel

  // scattered photons
  double total_scattered_weight; // total
  double total_scattered_weight_x2; // total squared sum (for unc calc)
  NumUtils::Matrix<float> scattered_weight_xy; // weight per image pixel
  NumUtils::Matrix<float> scattered_weight_xy_x2; // squared sum weight per image pixel

  // transfrom used to put the observer on the z-axis
  // used by rotate_zaxis_to_observer routine
  float rotate_transform[3][3];

  // position of observer in 3D
  float observer_position[3];
};

// structure with all the outputs for one run
struct output_struct {
  int num_outputs;  // number of outputs
  int arrays_allocated; // 1=all arrays have been allocated
  long image_size[2]; // size of output images in pixels
  string file_base; // base name of output files
  string type;  // type of output (ratio/flux current options)

  vector<one_output> outputs; // vector of outputs
};

#endif
