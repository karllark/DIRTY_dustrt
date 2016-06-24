// ======================================================================
//   Procedure to compute the first scattering of a photon, using the 
// technique of forcing the scattering to happen inside the dust.
//
// 2004 Dec/KDG - written
// 28 Feb 2007/KDG - updated to move calculation of stellar weight to 
//                   classify stellar photon
// 22 Oct 2013/KDG - added option to deal with tau_to_surface = 0
// 03 Sep 2015/KDG - added sampling from uniform distribution in tau
//                   1/2 of the time to better sample high optical depth scattering
// ======================================================================
#include "forced_first_scatter.h"
#include <unistd.h>
//#define DEBUG_FFS
//#define OUTNUM 0

int forced_first_scatter (geometry_struct& geometry,
			  photon_data& photon,
			  random_dirty& random_obj)

{
#ifdef DEBUG_FFS
  if (photon.number == OUTNUM) {
    cout << "ffs started; "; cout.flush();
    cout << "ffs1: photon.path_cur_cells = " << photon.path_cur_cells << endl;
  }

#endif

  // variables
  double target_tau = 1e20;
  int ffs_escape = 0;
  double distance_traveled = 0.0;
  double tau_traveled = 0.0;
  double tau_to_surface = 0.0;
  double ran_num = 0.0;
  
  photon_data dummy_photon; // copy of photon for surface tau calculation

  // get the optical depth to the edge of the dust in the current direction
  dummy_photon = photon;
  //dummy_photon.path_cur_cells = -1;  // set to -1 *not* to save cells tranversed
  dummy_photon.path_cur_cells = 0;  // set to 0 to save cells transversed

  distance_traveled = calc_photon_trajectory(dummy_photon, geometry, target_tau, ffs_escape, tau_to_surface);
  photon.first_tau = tau_to_surface;
  photon.prev_tau_surface = tau_to_surface;

  // cout << photon.number << " ";
  // cout << tau_to_surface << " ";
  // cout << ffs_escape << " ";
  // cout << distance_traveled << endl;

#ifdef DEBUG_FFS
  if (photon.number == OUTNUM) cout << "ffs2: photon.path_cur_cells = " << photon.path_cur_cells << endl;
  if (photon.number == OUTNUM) cout << "surface tau done; "; cout.flush();
#endif

  //  cout << photon.number << " " << tau_to_surface << endl;

  ffs_escape = 0;
  if (tau_to_surface == 0.0) {
    ffs_escape = 1;
#ifdef DEBUG_FFS
    if (photon.number == OUTNUM) cout << "photon escaping w/o interacting" << endl;
#endif
  } else {

    // do the weights
    // calculate the stellar weight (part of photon which escapes)
    //    calculation good to very high taus (700+), zero afterwards (KDG 28 Dec 2004)
    double new_stellar_weight = exp(-photon.first_tau);

    // determine optical depth to first scattering
    // unlike 2nd, 3rd, etc. scatterings, this scattering is forced to be 
    // in the dust distribution and is why the weights are not unity

    photon.scat_weight = photon.stellar_weight*(1. - new_stellar_weight);

    // calculate the weight which is scattered
    if (random_obj.random_num() >= geometry.scat_bias_fraction) { // classical forced scattering

      target_tau = -log(1.0 - random_obj.random_num()*(1.0 - new_stellar_weight));

    } else {// uniformaly sampled in optical depth

      ran_num = random_obj.random_num();
      target_tau = ran_num*photon.first_tau;
      
    }

    // calculate the biased weight factor
    double biased_weight_factor = 0.0;
    biased_weight_factor = (1.0 - geometry.scat_bias_fraction) +
      geometry.scat_bias_fraction*(1.0 - exp(-photon.first_tau))*exp(target_tau)/photon.first_tau;

    photon.scat_weight /= biased_weight_factor;

    photon.target_tau = target_tau;

    // update the stellar weight
    // *not done* this is done correctly in the classify_stellar_photon routine
    // KDG - 28 Feb 2007
    // photon.stellar_weight *= new_stellar_weight;


#ifdef DEBUG_FFS
    if (photon.number == OUTNUM) cout << "target_tau = " << target_tau << endl;
    if (photon.number == OUTNUM) cout << "tau_to_surface = " << tau_to_surface << endl;
#endif

    // check if the optical depth to the forced scattering is effectively
    // the same as the optical depth to the surface - adjust as necessary
    if ((1.0 - (target_tau/tau_to_surface)) <= ROUNDOFF_ERR_TRIG) {
#ifdef DEBUG_FFS
      cout << "tau_to_surface and target_tau the same ";
      cout << tau_to_surface << " " << target_tau << " ";
      cout << (1.0 - (target_tau/tau_to_surface)) << endl;
#endif
      target_tau = (1.0 - ROUNDOFF_ERR_TRIG)*tau_to_surface;
    }

  // move photon to location determined by target_tau and get the distance traveled
    int escape = 0;  // reset escape
    //photon.path_cur_cells = 0;  // set to 0 to save cells tranversed
#ifdef DEBUG_FFS
    if (photon.number == OUTNUM) cout << "tau_traveled in = " << tau_traveled << endl;
#endif
    // old code (23 Jun 2016)
    //distance_traveled = calc_photon_trajectory(photon, geometry, target_tau, escape, tau_traveled);

    // use the already computed photon track in dummy_photon instead of recalculating
    distance_traveled = calc_photon_trajectory_from_track(photon, dummy_photon, target_tau,
							  escape, tau_traveled);

#ifdef DEBUG_FFS

    // test if the two routines gave the same answer
    int i;
    cout << "ori: ";
    for (i = 0; i < 3; i++)
      cout << photon.position[i] << " ";
    cout << endl;
    cout << "new: ";
    for (i = 0; i < 3; i++)
      cout << save_photon.position[i] << " ";
    cout << endl;

    cout << "ori: ";
    for (i = 0; i < 3; i++)
      cout << photon.position_index[0][i] << " ";
    cout << endl;
    cout << "new: ";
    for (i = 0; i < 3; i++)
      cout << save_photon.position_index[0][i] << " ";
    cout << endl;
    
    if (photon.number == OUTNUM) {
      cout << "target_tau = " << target_tau << endl;
      cout << "tau_traveled = " << tau_traveled << endl;
      cout << "diff = " << target_tau - tau_traveled << endl;
      cout << "distance_traveled = " << distance_traveled << endl;
      cout << "escape = " << escape << endl;
    }
#endif
//     if (fabs((target_tau - tau_traveled)/target_tau) > ROUNDOFF_ERR_TRIG) {
//       cout << "*****error*****forced_first_scatter****" << endl;
//       cout << "target_tau = " << target_tau << endl;
//       cout << "tau_traveled = " << tau_traveled << endl;
//       cout << "diff = " << target_tau - tau_traveled << endl;
//       cout << "diff/target_tau = " << (target_tau - tau_traveled)/target_tau << endl;
//       cout << ROUNDOFF_ERR_TRIG << endl;
//       cout << "photon # = " << photon.number << endl;
//       cout << "tau_to_surface = " << photon.first_tau << endl;
//       exit(8);
//   } else {
//     // This is a special case where the photon has traveled the correct distance
//     // but has also just exited a subgrid.  This means that the photon index
//     // for the subgrid is below or above the actual subgrid size for one dimension.
//     // Easy to fix here by just redetermining the photon indexes [and does not 
//     // happen often at all (very rare)].
//     int m = 0;
//     int k = photon.current_grid_num;
//     int redetermine
//     for (m = 0; m < 3; m++) {
//       if (photon.position_index[k][m] < 0) {
// 	cout << "m = " << m << " ";
// 	cout << "tmp_photon.position_index = " << photon.position_index[k][m] << endl;
//       }
//     }

//     determine_photon_position_index_initial(geometry, photon);
//     }
    
    {
      /*
       * Continuous absorption
       * tau_entering/tau_leaving:
       *   total tau traveled when the photon packet enters/leaves the grid cell
       * prob_entering/prob_leaving:
       *   probability that the photon packet enters/leaves the grid cell
       */
      const double abs_weight_init = (1. - geometry.albedo)*dummy_photon.stellar_weight; // not the reduced scat_weight
      double tau_entering = 0.;
      double prob_entering = 1.;
      
      for (int i = 0; i < dummy_photon.path_cur_cells; i++) {
	// find the absorbed weight
	double tau_leaving = tau_entering + dummy_photon.path_tau[i];
	double prob_leaving = exp(-tau_leaving);
	double abs_weight = abs_weight_init*(prob_entering - prob_leaving);
	
	// deposit the energy
	grid_cell& this_cell = geometry.grids[dummy_photon.path_pos_index[0][i]].grid(dummy_photon.path_pos_index[1][i],dummy_photon.path_pos_index[2][i],dummy_photon.path_pos_index[3][i]);
	this_cell.absorbed_energy[geometry.abs_energy_wave_index] += abs_weight;
	this_cell.absorbed_energy_x2[geometry.abs_energy_wave_index] += abs_weight*abs_weight;
	this_cell.absorbed_energy_num_photons[geometry.abs_energy_wave_index]++;
	
	// move to the next grid cell
	tau_entering = tau_leaving;
	prob_entering = prob_leaving;
      }
      
    }
    
  }
  
  return(ffs_escape);
}
