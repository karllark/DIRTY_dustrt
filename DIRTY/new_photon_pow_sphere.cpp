// ======================================================================
//   Procedure to generate a photon and set up the stats on it.
// This procedure is used for a continuous spherical stellar distribution
// with a power-law density
//
// 2008 Dec/KHL - created from "new_photon_discrete_stars.cpp"
// ======================================================================
#include "new_photon_pow_sphere.h"

void
new_photon_pow_sphere (photon_data &photon, geometry_struct &geometry, random_dirty &random_obj)

{
  // setup the weights
  photon.stellar_weight = 1.0;
  photon.scat_weight = 0.0;

  // initialize statistics variables
  photon.num_scat = 0;

  // direction of photon; assuming an isotropic source
  // in direction cosines...
  double phi = M_PI * (2.0 * random_obj.random_num () - 1.0);
  photon.dir_cosines[2] = 2.0 * random_obj.random_num () - 1.0;
  double temp = sqrt (1.0 - pow (photon.dir_cosines[2], 2));
  photon.dir_cosines[0] = cos (phi) * temp;
  photon.dir_cosines[1] = sin (phi) * temp;

#ifdef DEBUG_NPDS
  for (int i = 0; i < 3; i++)
    cout << photon.dir_cosines[i] << " ";
  cout << "starting dir cosines" << endl;
#endif

  // determine where the photon will be emitted from
  double r = pow (geometry.pow_sphere_constant1 * random_obj.random_num () + geometry.pow_sphere_constant2,
                  geometry.pow_sphere_constant3);
  double cos_theta = 2.0 * random_obj.random_num () - 1.0; // theta = 0 on x-y plane, not the polar angle
  double r_sin_theta = r * sqrt (1.0 - cos_theta * cos_theta);
  phi = M_PI * (2.0 * random_obj.random_num () - 1.0);
  photon.position[0] = r_sin_theta * cos (phi);
  photon.position[1] = r_sin_theta * sin (phi);
  photon.position[2] = r * cos_theta;

  // save the birth photon position
  for (int i = 0; i < 3; i++)
    photon.birth_position[i] = photon.position[i];

  // now determine the position indexes of the photon
  determine_photon_position_index_initial (geometry, photon);
}
