// ======================================================================
//   Procedure to scatter a photon into a new direction.  
//
// 2004 Dec/KDG - written
// 2008 Aug/KDG - added continous absorption
// 2009 Mar/KDG - added limited polychromaticism
// ======================================================================
#include "scatter_photon.h"
//#define DEBUG_SP

void scatter_photon (geometry_struct& geometry,
		     photon_data& photon,
		     random_dirty& random_obj)

{

  double sqr_g = 0.0;
  double cos_alpha = 0.0;

  // determine the cosine of the scattering angle
  //   if g is to near zero (isotropic), just randomly determine angle
  //   equations tested to work as long as g > 1e-15  (KDG 28 Dec 2004)
  if (geometry.g == -2) { // model phase function
    double rannum = random_obj.random_num();
    uint i1,i2,i3;
    i1 = 0;
    i3 = geometry.phi.size();
    i2 = (i1 + i3)/2;
    while ((i3 - i1) > 1) {
      if (rannum > geometry.phi_sum[i2]) i1 = i2; else i3 = i2;
//       cout << rannum << " ";
//       cout << geometry.phi[i1] << " " << geometry.phi[i3] << " ";
//       cout << i1 << " " << i3 << endl;
      i2 = (i1 + i3)/2;
    }
    cos_alpha = geometry.phi_angle[i1] - 
      ((rannum - geometry.phi_sum[i3])/(geometry.phi_sum[i1] - geometry.phi_sum[i3]))*
      (geometry.phi_angle[i1] - geometry.phi_angle[i3]);
//     cout << float(i1)/float(geometry.phi.size()) << " ";
//     cout << (rannum - geometry.phi_sum[i3])/(geometry.phi_sum[i1] - geometry.phi_sum[i3]) << endl;
//     cout << "cos_angle = " << cos_alpha << endl;
//     exit(8);
//     cout << geometry.phi_angle[i1] << " ";
//     cout << cos_alpha << " ";
//     cout << geometry.phi_angle[i3] << endl;
//     cout.flush();
  } else if (geometry.g > ROUNDOFF_ERR_TRIG) {
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
//   int k = photon.current_grid_num;
//   long grid_val = photon.grid_number[k];
//   geometry.grids[grid_val].grid(photon.position_index[k][0],photon.position_index[k][1],photon.position_index[k][2]).absorbed_energy[geometry.abs_energy_wave_index] +=
//     (1. - geometry.albedo)*photon.scat_weight;

  // continuous absorption (instead of point absorption - should help convergence)
#ifdef DEGUG_SP
  cout << photon.num_scat << endl;
  float tmp_sum = 0.0;
#endif
  float abs_weight = (1. - geometry.albedo)*photon.scat_weight;

//   if (runinfo.limited_polychromatic) {
//     // now setup the weights for limited polychromaticism
//     int k = 0;
//     for (i = 0; i < photon.n_polyc; i++) {
//       k = i+photon.polyc_min_index;
//       photon.polyc_weights[i] = exp((1.0 - runinfo.tau_to_tau_ref[k])*photon.target_tau)*
// 	runinfo.tau_to_tau_ref[k]*
// 	(((1.0 - pow(runinfo.g[k],2))/(4.0*M_PI*pow(1.0 + pow(runinfo.g[k],2) - 2.0*runinfo.g[k]*cos_alpha,1.5)))/
// 	 ((1.0 - pow(geometry.g,2))/(4.0*M_PI*pow(1.0 + pow(geometry.g,2) - 2.0*geometry.g*cos_alpha,1.5))))*
// 	(1. - runinfo.albedo[k)*photon.scat_weight;
//     }
//   }

  for (i = 0; i < photon.path_cur_cells; i++) {
    geometry.grids[photon.path_pos_index[0][i]].grid(photon.path_pos_index[1][i],photon.path_pos_index[2][i],photon.path_pos_index[3][i]).absorbed_energy[geometry.abs_energy_wave_index] += abs_weight*(photon.path_tau[i]/photon.target_tau);
    geometry.grids[photon.path_pos_index[0][i]].grid(photon.path_pos_index[1][i],photon.path_pos_index[2][i],photon.path_pos_index[3][i]).absorbed_energy_x2[geometry.abs_energy_wave_index] += pow(abs_weight*(photon.path_tau[i]/photon.target_tau),2.0);
    geometry.grids[photon.path_pos_index[0][i]].grid(photon.path_pos_index[1][i],photon.path_pos_index[2][i],photon.path_pos_index[3][i]).absorbed_energy_num_photons[geometry.abs_energy_wave_index]++;
#ifdef DEBUG_SP
    cout << photon.number << " ";
    cout << "path: ";
    cout << i << " ";
    cout << photon.path_pos_index[0][i] << " ";
    cout << photon.path_pos_index[1][i] << " ";
    cout << photon.path_pos_index[2][i] << " ";
    cout << photon.path_pos_index[3][i] << " ";
    cout << photon.path_tau[i] << " ";
    cout << geometry.grids[photon.path_pos_index[0][i]].grid(photon.path_pos_index[1][i],photon.path_pos_index[2][i],photon.path_pos_index[3][i]).dust_tau_per_pc*geometry.tau_to_tau_ref*geometry.grids[photon.path_pos_index[0][i]].phys_cube_size[0] << " ";
    cout << abs_weight*(photon.path_tau[i]/photon.target_tau) << " ";
    cout << abs_weight << " ";
    cout << endl;
    //    tmp_sum += photon.path_tau[i];
#endif
  }
#ifdef DEGUG_SP
  cout << "total tau = " << photon.target_tau << " " << tmp_sum << endl;
#endif  
  // update the scattered weight
  photon.scat_weight *= geometry.albedo;

  // update the number of scatterings
  photon.num_scat++;
}
