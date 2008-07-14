// ======================================================================
//   Procedure to compute the first scattering of a photon, using the 
// technique of forcing the scattering to happen inside the dust.
//
// 2004 Dec/KDG - written
// 28 Feb 2007/KDG - updated to move calculation of stellar weight to 
//                   classify stellar photon
// ======================================================================
#include "forced_first_scatter.h"
#include <unistd.h>
//#define DEBUG_FFS
//#define OUTNUM 222074

void forced_first_scatter (geometry_struct& geometry,
			   photon_data& photon,
			   random_dirty& random_obj)

{
#ifdef DEBUG_FFS
  if (photon.number == OUTNUM) cout << "ffs started; "; cout.flush();
#endif

  // variables
  double target_tau = 1e20;
  int escape = 0;
  double distance_traveled = 0.0;
  double tau_traveled = 0.0;
  double tau_to_surface = 0.0;

  photon_data dummy_photon; // copy of photon for surface tau calculation

  // get the optical depth to the edge of the dust in the current direction
  dummy_photon = photon;
  distance_traveled = calc_photon_trajectory(dummy_photon, geometry, target_tau, escape, tau_to_surface);
  photon.first_tau = tau_to_surface;
#ifdef DEBUG_FFS
  if (photon.number == OUTNUM) cout << "surface tau done; "; cout.flush();
#endif

  // do the weights
  // calculate the stellar weight (part of photon which escapes)
  //    calculation good to very high taus (700+), zero afterwards (KDG 28 Dec 2004)
  double new_stellar_weight = exp(-photon.first_tau);
  // calculate the weight which is scattered
  photon.scat_weight = photon.stellar_weight - new_stellar_weight;
  // update the stellar weight
  // *not done* this is done correctly in the classify_stellar_photon routine
  // KDG - 28 Feb 2007
  // photon.stellar_weight *= new_stellar_weight;

  // determine optical depth to first scattering
  // unlike 2nd, 3rd, etc. scatterings, this scattering is forced to be 
  // in the dust distribution and is why the weights are done as above
  target_tau = -log(1.0 - random_obj.random_num()*(1.0 - new_stellar_weight));

//   cout << "tau_to_surf = " << tau_to_surface << endl;
#ifdef DEBUG_FFS
  if (photon.number == OUTNUM) cout << "target_tau = " << target_tau << endl;
  if (photon.number == OUTNUM) cout << "tau_to_surface = " << tau_to_surface << endl;
#endif
		   
  // move photon to location determined by target_tau and get the distance traveled
  escape = 0;  // reset escape
#ifdef DEBUG_FFS
  if (photon.number == OUTNUM) cout << "tau_traveled in = " << tau_traveled << endl;
#endif
  distance_traveled = calc_photon_trajectory(photon, geometry, target_tau, escape, tau_traveled);

#ifdef DEBUG_FFS
  if (photon.number == OUTNUM) {
    cout << "target_tau = " << target_tau << endl;
    cout << "tau_traveled = " << tau_traveled << endl;
    cout << "diff = " << target_tau - tau_traveled << endl;
    cout << "distance_traveled = " << distance_traveled << endl;
    cout << "escape = " << escape << endl;
  }
#endif
  if (fabs(target_tau - tau_traveled) > ROUNDOFF_ERR_TRIG) {
    cout << "*****error*****forced_first_scatter****" << endl;
    cout << "target_tau = " << target_tau << endl;
    cout << "tau_traveled = " << tau_traveled << endl;
    cout << "diff = " << target_tau - tau_traveled << endl;
    exit(8);
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
  }

}
