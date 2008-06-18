// ======================================================================
//   Procedure to compute the transformation/rotation matrix for the 
// observer.
//
// 2008 Jun/KDG - written (taken from initialize_output.cpp)
// ======================================================================
#include "compute_observer_trans_matrix.h"

void compute_observer_trans_matrix (output_struct& output,
				    geometry_struct& geometry,
				    int i)

{
  // **note added 3 Jul 2006**
  //   this transformation is clearly not the standard rotation matrix in 3D
  //   while I can't remember fully, I believe that I derived this transformation
  //   matrix specifically to rotate so the observer position is along the z-axis
  //   if the observer is put at the standard theta,phi location
  //   (need to check this)
  float theta; 
  float phi;
  
  if (geometry.num_observers > 1) {
    theta = geometry.observer_angles[0][i];
    phi = geometry.observer_angles[1][i];
  } else { // handle the case when the number of outputs is set by something other than multiple lines-of-sight
    theta = geometry.observer_angles[0][0];
    phi = geometry.observer_angles[1][0];
  }
  
#ifdef DEBUG_COTM
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
  
#ifdef DEBUG_COTM
  int j,k;
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

}
