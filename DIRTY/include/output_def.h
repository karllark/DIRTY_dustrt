// ======================================================================
//   Header file with output definitions.
//
// 2005 May/KDG - written
// ======================================================================

#ifndef _DIRTY_OUTPUT_DEF_
#define _DIRTY_OUTPUT_DEF_

#include "NumUtils.h"
// #include "vect_utils.h"

// structure with the output for one line-of-sight, albedo, etc.
struct one_output
{
    // number of photons
    double total_num_photons;                        // total
    NumUtils::Matrix<double> num_stellar_photons_xy; // photons per image pixel
    NumUtils::Matrix<double> num_photons_xy;         // scattered photons per image pixel

    // number of scattered photons
    double total_num_scattered_photons; // total

    // stellar photons
    double total_stellar_weight;                   // total
    double total_stellar_weight_x2;                // total squared sum (for unc calc)
    NumUtils::Matrix<double> stellar_weight_xy;    // weight per image pixel
    NumUtils::Matrix<double> stellar_weight_xy_x2; // squared sum weight per image pixel

    // scattered photons
    double total_scattered_weight;                   // total
    double total_scattered_weight_x2;                // total squared sum (for unc calc)
    NumUtils::Matrix<double> scattered_weight_xy;    // weight per image pixel
    NumUtils::Matrix<double> scattered_weight_xy_x2; // squared sum weight per image pixel

    // other
    double ave_first_tau;    // average of tau (starts as a sum)
    double ave_first_tau_x2; // total squared sum (for unc calc)

    // transfrom used to put the observer on the z-axis
    // used by rotate_zaxis_to_observer routine
    float rotate_transform[3][3];

    // position of observer in 3D
    float observer_position[3];
};

// structure with all the outputs for one run
struct output_struct
{
    int num_outputs;          // number of outputs
    int arrays_allocated;     // 1=all arrays have been allocated
    long image_size[2];       // size of output images in pixels
    string file_base;         // base name of output files
    string emission_type;     // string to identify the output type ("" if stellar)
    string type;              // type of output (ratio/flux current options)
    int do_output_model_grid; // output the model grid?

    vector<one_output> outputs; // vector of outputs
};

#endif
