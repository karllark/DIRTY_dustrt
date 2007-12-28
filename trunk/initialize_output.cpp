// ======================================================================
//   Procedure to allocate and initialize the output stucture
//
// 2005 May/KDG - written
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
    // **note added 3 Jul 2006**
    //   this transformation is clearly not the standard rotation matrix in 3D
    //   while I can't remember fully, I believe that I derived this transformation
    //   matrix specifically to rotate so the observer position is along the z-axis
    //   if the observer is put at the standard theta,phi location
    //   (need to check this)
    float theta = geometry.observer_angles[0][i];
    float phi = geometry.observer_angles[1][i];

#ifdef DEBUG_INITOUT
    cout << "observer (theta,phi) = (" << theta*180./M_PI << "," << phi*180./M_PI << ")" << endl;
#endif

    output.outputs[i].rotate_transform[0][0] = cos(theta)*cos(phi);
    output.outputs[i].rotate_transform[0][1] = cos(theta)*sin(phi);
    output.outputs[i].rotate_transform[0][2] = -sin(theta);
    
    output.outputs[i].rotate_transform[1][0] = -sin(phi);
    output.outputs[i].rotate_transform[1][1] = cos(phi);
    output.outputs[i].rotate_transform[1][2] = 0.0;
    
    output.outputs[i].rotate_transform[2][0] = sin(theta)*cos(phi);
    output.outputs[i].rotate_transform[2][1] = sin(theta)*sin(phi);
    output.outputs[i].rotate_transform[2][2] = cos(theta);
    

    // attempting to have everything referenced to a source on the x-axis
    // all angles negative as this transform is setup to rotate the source
    // onto the z-axis so that the yz postions can be used to find the location
    // of each photon in the output image
    // ***doesn't work*** probably need to do the full spherical transformation
    //   which is what I remember doing previously and is above (which does work!)
    
//     output.outputs[i].rotate_transform[0][0] = cos(-phi)*cos(-theta);
//     output.outputs[i].rotate_transform[0][1] = sin(-phi);
//     output.outputs[i].rotate_transform[0][2] = cos(-phi)*sin(-theta);
    
//     output.outputs[i].rotate_transform[1][0] = -sin(-phi)*cos(-theta);
//     output.outputs[i].rotate_transform[1][1] = cos(-phi);
//     output.outputs[i].rotate_transform[1][2] = -sin(-phi)*sin(-theta);
    
//     output.outputs[i].rotate_transform[2][0] = -sin(-theta);
//     output.outputs[i].rotate_transform[2][1] = 0.0;
//     output.outputs[i].rotate_transform[2][2] = cos(-theta);

    int j,k;
#ifdef DEBUG_INITOUT
    cout << "transformation matrix for (theta,phi) = "; 
    cout << "(" << theta << "," << phi << ")" << endl;
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++)
	cout << output.outputs[i].rotate_transform[j][k] << " ";
      cout << endl;
    }
#endif

    output.outputs[i].observer_position[0] = geometry.distance*sin(theta)*cos(phi);
    output.outputs[i].observer_position[1] = geometry.distance*sin(theta)*sin(phi);
    output.outputs[i].observer_position[2] = geometry.distance*cos(theta);

//     cout << theta << " " << phi << endl;
//     cout << output.outputs[i].observer_position[0] << " ";
//     cout << output.outputs[i].observer_position[1] << " ";
//     cout << output.outputs[i].observer_position[2] << " ";
//     cout << endl;

    // zero output quantities
    output.outputs[i].total_num_photons = 0.0;
    output.outputs[i].total_num_scattered_photons = 0.0;
    output.outputs[i].total_stellar_weight = 0.0;
    output.outputs[i].total_stellar_weight_x2 = 0.0;
    output.outputs[i].total_scattered_weight = 0.0;
    output.outputs[i].total_scattered_weight_x2 = 0.0;

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
