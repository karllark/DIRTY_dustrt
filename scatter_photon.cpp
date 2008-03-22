// ======================================================================
//   Procedure to scatter a photon into a new direction.  
//
// 2004 Dec/KDG - written
// ======================================================================
#include "scatter_photon.h"

void scatter_photon (geometry_struct& geometry,
		     photon_data& photon,
		     random_dirty& random_obj)

{

  double sqr_g = 0.0;
  double cos_alpha = 0.0;

  // determine the cosine of the scattering angle
  //   if g is to near zero (isotropic), just randomly determine angle
  //   equations tested to work as long as g > 1e-15  (KDG 28 Dec 2004)
  if (geometry.g > ROUNDOFF_ERR_TRIG) {
    sqr_g = pow(geometry.g,2);
    cos_alpha = (1.0 + sqr_g) - 
      pow((1.0 - sqr_g)/(1.0 - geometry.g + 2.0*geometry.g*random_obj.random_num()),2);
    cos_alpha /= (2.0*geometry.g);
  } else
    cos_alpha = 2.0*random_obj.random_num() - 1.0;

  // calculate sine of scattering angle
  double sin_alpha = 1.0 - pow(cos_alpha,2);
  if (sin_alpha != 0.0) sin_alpha = sqrt(sin_alpha);

  // determine the angle perpendicular to the photon direction
  //   assuming unaligned grains -> random
  double phi = M_PI*(2.0*random_obj.random_num() - 1.0);
  // get the sine and cosine of this angle
  double cos_phi = cos(phi);
  double sin_phi = sin(phi);

  // calculate the new direction cosines
  //   equations 22 & 23 of Witt (1977, ApJS, 35, 1)
  //   roundoff error good to approx. 1e-16 (checked KDG 28 Dec 2004)
  double new_dir_cosines[3];
  if ((1.0 - fabs(photon.dir_cosines[2])) > ROUNDOFF_ERR_TRIG) {
    double tmp = sqrt(1.0 - pow(photon.dir_cosines[2],2));
    new_dir_cosines[0] = sin_alpha*(cos_phi*photon.dir_cosines[2]*photon.dir_cosines[0] -
			    sin_phi*photon.dir_cosines[1]);
    new_dir_cosines[0] = (new_dir_cosines[0]/tmp) + cos_alpha*photon.dir_cosines[0];
    new_dir_cosines[1] = sin_alpha*(cos_phi*photon.dir_cosines[2]*photon.dir_cosines[1] +
			    sin_phi*photon.dir_cosines[0]);
    new_dir_cosines[1] = (new_dir_cosines[1]/tmp) + cos_alpha*photon.dir_cosines[1];
    new_dir_cosines[2] = -sin_alpha*cos_phi*tmp + cos_alpha*photon.dir_cosines[2];
  } else {
    new_dir_cosines[0] = sin_alpha*cos_phi;
    new_dir_cosines[1] = sin_alpha*sin_phi;
    new_dir_cosines[2] = cos_alpha*photon.dir_cosines[2];
  }

  // update the direction cosines of the photon
  int i;
  for (i = 0; i < 3; i++) 
    photon.dir_cosines[i] = new_dir_cosines[i];

  // update the absorbed energy in the grid
  int k = photon.current_grid_num;
  long grid_val = photon.grid_number[k];
  geometry.grids[grid_val].grid(photon.position_index[k][0],photon.position_index[k][1],photon.position_index[k][2]).absorbed_energy[geometry.abs_energy_wave_index] +=
    (1. - geometry.albedo)*photon.scat_weight;
  
  // update the scattered weight
  photon.scat_weight *= geometry.albedo;

  // update the number of scatterings
  photon.num_scat++;
}