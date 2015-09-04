// ======================================================================
//   Procedure to scatter a photon into a new direction.  
//
// 2004 Dec/KDG - written
// 2008 Aug/KDG - added continous absorption
// 2009 Mar/KDG - added limited polychromaticism
// 2013 Oct/KDG - updated to work with subgrids
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

#ifdef DEBUG_SP
  cout << "# scat = " << photon.num_scat << endl;
  cout << "total tau = " << photon.target_tau << endl;
#endif  

  // update the scattered weight
  photon.scat_weight *= geometry.albedo;

  // updated to include the weight reduction due to continuous absorption - KDG 17 jun 15
  // does not seem to work - KDG 17 jun 15
  //photon.scat_weight *= geometry.albedo*(1.0 - exp(-1.0*photon.prev_tau_surface));

  // update the number of scatterings
  photon.num_scat++;

  // new continuous absorption (written by Ka-Hei Law)
  // appears after adjusting the photon.scat_weight as this absorption
  // is for the *next* scattering/interaction
  {
    // find path_tau[]
    photon_data dummy_photon = photon;
    dummy_photon.current_grid_num = 0;  // set to the base grid to start tarjectory correctly
    dummy_photon.path_cur_cells = 0; // set to 0 to save cells traversed

    double target_tau = 1e20;
    int escape = 0;
    double tau_to_surface = 0.0;
    calc_photon_trajectory(dummy_photon, geometry, target_tau, escape, tau_to_surface);

    photon.prev_tau_surface = tau_to_surface;

    /*
     * Continuous absorption
     * tau_entering/tau_leaving:
     *   total tau traveled when the photon packet enters/leaves the grid cell
     * prob_entering/prob_leaving:
     *   probability that the photon packet enters/leaves the grid cell
     */
    const double abs_weight_init = (1. - geometry.albedo)*dummy_photon.scat_weight;

    double tau_entering = 0.;
    double prob_entering = 1.;
    
#ifdef DEBUG_SP
    double tot_abs_temp = 0.;
#endif

    double tau_leaving = 0.0;
    double prob_leaving = 0.0;
    double abs_weight = 0.0;

    for (int i = 0; i < dummy_photon.path_cur_cells; i++) {
      // find the absorbed weight
      tau_leaving = tau_entering + dummy_photon.path_tau[i];
      prob_leaving = exp(-tau_leaving);
      abs_weight = abs_weight_init*(prob_entering - prob_leaving);

#ifdef DEBUG_SP
      cout << i << " ";
      cout << tau_entering << " ";
      cout << tau_leaving << " ";
      cout << prob_entering << " ";
      cout << prob_leaving << " ";
      cout << abs_weight << endl;
      
      tot_abs_temp += abs_weight;
#endif  

      // deposit the energy
      grid_cell& this_cell = geometry.grids[dummy_photon.path_pos_index[0][i]].grid(dummy_photon.path_pos_index[1][i],dummy_photon.path_pos_index[2][i],dummy_photon.path_pos_index[3][i]);
      this_cell.absorbed_energy[geometry.abs_energy_wave_index] += abs_weight;
      this_cell.absorbed_energy_x2[geometry.abs_energy_wave_index] += abs_weight*abs_weight;
      this_cell.absorbed_energy_num_photons[geometry.abs_energy_wave_index]++;

      // move to the next grid cell
      tau_entering = tau_leaving;
      prob_entering = prob_leaving;
    }

#ifdef DEBUG_SP
  cout << abs_weight_init << " ";
  cout << tot_abs_temp << endl;
  if (photon.target_tau > 0.5) exit(8);
#endif

  }
}
