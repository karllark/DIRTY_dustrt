// ======================================================================
//   Procedure to generate a photon and set up the stats on it.          
// This procedure is used for stars distributed in a double 
// exponential disk.
//
// 2008 Aug/KDG - written
// ======================================================================
#include "new_photon_dexp_disk.h"

void new_photon_dexp_disk (photon_data& photon,
			   geometry_struct& geometry,
			   random_dirty random_obj)


{
  // setup the weights
  photon.stellar_weight = 1.0;
  photon.scat_weight = 0.0;

  // initialize statistics variables
  photon.num_scat = 0;

  // direction of photon; assuming an isotropic source
  // in direction cosines...
  double phi = M_PI*(2.0*random_obj.random_num() - 1.0);
  photon.dir_cosines[2] = 2.0*random_obj.random_num() - 1.0;
  double temp = sqrt(1.0 - pow(photon.dir_cosines[2],2));
  photon.dir_cosines[0] = cos(phi)*temp;
  photon.dir_cosines[1] = sin(phi)*temp;

  int i = 0;
#ifdef DEBUG_NPDD
  for (i = 0; i < 3; i++)
    cout << photon.dir_cosines[i] << " ";
  cout << "starting dir cosines" << endl;
#endif

  // determine which star the photon will be emitted from
  // assuming a double exponential disk
  double ran_num = random_obj.random_num();
  // z first as it is analytic
  photon.position[2] = -1.0*geometry.stellar_scaleheight*log(1.0 - ran_num/geometry.stellar_emit_constant_z);
  // now determine if this is above or below the plane (using a random number
  if (random_obj.random_num() > 0.5) photon.position[2] *= -1.0;

  // xy is analytic, but not easily in a simple form
  ran_num = random_obj.random_num();
  double min_xy = 0.0;
  double max_xy = geometry.radius;
  double test_xy = geometry.radius/2.0;
  float right_side = log(1.0 - ran_num/geometry.stellar_emit_constant_xy);
  float left_side = -1.0*test_xy/geometry.stellar_scalelength + log(test_xy/geometry.stellar_scalelength + 1.0);
//   cout << -1.0*min_xy/geometry.stellar_scalelength + log(min_xy/geometry.stellar_scalelength + 1.0) << " ";
//   cout << -1.0*max_xy/geometry.stellar_scalelength + log(max_xy/geometry.stellar_scalelength + 1.0) << endl;
//   cout << left_side << " " << right_side << " " << test_xy << endl;
//   float tans = 0.0;
  while (fabs((left_side - right_side)/left_side) > 0.01) {
//     cout << left_side << " " << right_side << " " << test_xy << endl;
    if (left_side > right_side) min_xy = test_xy; else max_xy = test_xy;
//     cout << min_xy << " " << test_xy << " " << max_xy << endl;
    test_xy = (max_xy + min_xy)/2.0;
    left_side = -1.0*test_xy/geometry.stellar_scalelength + log(test_xy/geometry.stellar_scalelength + 1.0);
//     cin >> tans;
  }
  // now that we have rho, get xy (ran number for phi)
  phi = M_PI*(2.0*random_obj.random_num() - 1.0);
  photon.position[0] = cos(phi)*test_xy;
  photon.position[1] = sin(phi)*test_xy;

  // save the birth position
  for (i = 0; i < 3; i++) {
//     cout << photon.position[i] << " ";
    photon.birth_position[i] = photon.position[i];
  }
//   cout << endl;

  // now determine the position indexes of the photon
  determine_photon_position_index_initial(geometry, photon);

}
