// ======================================================================
//   Procedure to generate a photon and set up the stats on it.          
// This procedure is used for a single star located somewhere in the 
// dust distribution
//
// 2006 Nov/KDG - written
// 8 Mar 2007/KDG - fixed setup of initial position (rotation matrix)
//                  correct decission between isotropic/file radiation field
// ======================================================================
#include "new_photon_diffuse_source.h"
//#define DEBUG_NPDS

void new_photon_diffuse_source (photon_data& photon,
				geometry_struct& geometry,
				random_dirty& random_obj)


{
  // setup the weights
  photon.stellar_weight = 1.0;
  photon.scat_weight = 0.0;
  
  // initialize statistics variables
  photon.num_scat = 0;
  
  // save the photon info for later
  //   photon_data save_photon = photon;
  
  double phi = 0.0;
  int i = 0;
  if (geometry.new_photon_source_type == NEW_PHOTON_DIFFUSE_ISOTROPIC) {
    // direction of photon; assuming an isotropic source
    // in direction cosines...
    phi = M_PI*(2.0*random_obj.random_num() - 1.0);
    photon.dir_cosines[2] = 2.0*random_obj.random_num() - 1.0;
    // testing of single illumination directions (should probably make this an option)
//     phi = 0.*(M_PI/180.);
//     photon.dir_cosines[2] = 0.0;
    double temp = sqrt(1.0 - pow(photon.dir_cosines[2],2));
    photon.dir_cosines[0] = cos(phi)*temp;
    photon.dir_cosines[1] = sin(phi)*temp;
  } else {
    // determine which bin the photon emerges from
    double ran_num = random_obj.random_num();
    while ((i < int(geometry.diffuse_source_sum_intensity.size())) && (ran_num > geometry.diffuse_source_sum_intensity[i])) i++; 
    int pos_index = i;
//     cout << geometry.diffuse_source_theta[pos_index] << " ";
//     cout << geometry.diffuse_source_phi[pos_index] << endl;
#ifdef DEBUG_NPDS
    cout << "theta = " << geometry.diffuse_source_theta[pos_index] << endl;
    cout << "phi = " << geometry.diffuse_source_phi[pos_index] << endl;
#endif
    
    // direction of photon 
    // from diffuse source location
    phi = geometry.diffuse_source_phi[pos_index];
    photon.dir_cosines[2] = cos(geometry.diffuse_source_theta[pos_index]+(M_PI/2.));
    double temp = sqrt(1.0 - pow(photon.dir_cosines[2],2));
    photon.dir_cosines[0] = cos(phi)*temp;
    photon.dir_cosines[1] = sin(phi)*temp; 
  }
  
#ifdef DEBUG_NPDS
  cout << "photon.dir_cosines[0..2] = ";
  cout << photon.dir_cosines[0] << " ";
  cout << photon.dir_cosines[1] << " ";
  cout << photon.dir_cosines[2] << endl;
#endif
  // start the photon in a plane which is perpendicular to the direction
  // first setup the photon in the xy plane
  double pos_phi = M_PI*(2.0*random_obj.random_num() - 1.0);
  // ensure a uniform distribution in a circle
  //    r_1 = sqrt(random*r_max^2)
  double pos_radius = sqrt(pow(0.95*geometry.radius,2.)*random_obj.random_num());
  photon.position[0] = 0.0;
  photon.position[1] = pos_radius*cos(pos_phi);
  photon.position[2] = pos_radius*sin(pos_phi);
  // now rotate the photon position so that it is perpendicular to the direction
  //  use theta & phi determined above
  float rotate_transform[3][3];
  // adjust for theta_dir_cosine = 90 - theta
  double sin_theta = photon.dir_cosines[2];
  double cos_theta = sqrt(1.0 - pow(sin_theta,2));
  // derived as the matrix multiplication between rotation in the xy(phi) and xz(theta) planes
  rotate_transform[0][0] = cos_theta*cos(phi);
  rotate_transform[0][1] = -sin(phi);
  rotate_transform[0][2] = -sin_theta*cos(phi);

  rotate_transform[1][0] = cos_theta*sin(phi);
  rotate_transform[1][1] = cos(phi);
  rotate_transform[1][2] = -sin_theta*sin(phi);

  rotate_transform[2][0] = sin_theta;
  rotate_transform[2][1] = 0.0;
  rotate_transform[2][2] = cos_theta;

#ifdef DEBUG_NPDS
  cout << "np; start rotate" << " ";
  cout.flush();
#endif
  // apply the rotation (use the existing routine)
  rotate_zaxis_for_observer(rotate_transform, photon);
#ifdef DEBUG_NPDS
  cout << "np; done rotate" << endl;
  cout << "photon # = " << photon.number << endl;
  for (i = 0; i < 3; i++)
    cout << i << " " << photon.position[i] << " " << photon.dir_cosines[i] << endl;
  cout.flush();
#endif

  // now determine the position indexes of the photon
  determine_photon_position_index_initial(geometry, photon);
#ifdef DEBUG_NPDS
  cout << "np; done index pos" << " ";
  cout.flush();
#endif

  // now move the photon to the edge
  double target_tau = 1e20;
  int escape = 0;
  double distance_traveled = 0.0;
  double tau_to_surface = 0.0;
  photon.path_cur_cells = -1;  // set to -1 *not* to save cells tranversed
#ifdef DEBUG_NPDS
  cout << "np; start distance traveled" << " ";
  cout.flush();
#endif
  distance_traveled = calc_photon_trajectory(photon, geometry, target_tau, escape, tau_to_surface);
#ifdef DEBUG_NPDS
  cout << "np; done distance traveled" << " ";
  cout.flush();
#endif

  // now flip dir_cosines so the photon is heading into the nebula, not out
  for (i = 0; i < 3; i++) {
//     if (photon.number == 453) {
//       cout << i << " ";
//       cout << photon.dir_cosines[i] << " ";
//       cout << photon.position[i] << endl;
//     }
    photon.dir_cosines[i] *= -1.;
    photon.birth_position[i] = photon.position[i];
  }
  //if (photon.number == 453) exit(8);

#ifdef DEBUG_NPDS
  cout << "1: photon position and dir_cosine at edge" << endl;
  for (i = 0; i < 3; i++)
    cout << i << " " << photon.position[i] << " " << photon.dir_cosines[i] << endl;
#endif

  // now determine the position indexes of the photon (again to handle the edges)
  determine_photon_position_index_initial(geometry, photon);

#ifdef DEBUG_NPDS
  cout << "2: photon position and dir_cosine at edge" << endl;
  cout << "photon # = " << photon.number << endl;
  for (i = 0; i < 3; i++)
    cout << i << " " << photon.position[i] << " " << photon.dir_cosines[i] << endl;
#endif

  //if (photon.number == 453) exit(8);

}
