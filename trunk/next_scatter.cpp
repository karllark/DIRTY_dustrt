// ======================================================================
//   Procedure to scatter a photon into a new direction.  This procedure
// returns a 0 if the photon scatters inside the dust and a 1 if it 
// escapes.
//
// 2004 Dec/KDG - written
// ======================================================================
#include "next_scatter.h"
//#define DEBUG_NS

int next_scatter (geometry_struct& geometry,
		  photon_data& photon,
		  random_dirty& random_obj)

{
  // determine the optical depth to the next scattering
  double target_tau = 0.0;
  target_tau = -log(random_obj.random_num());
  
  // check to see if we will start in a subgrid
  if (photon.current_grid_num > 0) {
#ifdef DEBUG_NS
    cout << "starting in a subgrid" << endl;
#endif
    photon.current_grid_num = 0;
  }

  // determine the site of the next scattering
  int escape = 0;
  double distance_traveled = 0.0;
  double tau_traveled = 0.0;
  distance_traveled = calc_photon_trajectory(photon, geometry, target_tau, escape, tau_traveled);
#ifdef DEBUG_NS
  cout << "ns cpt done; ";
  cout << "distance_traveled = " << distance_traveled << endl;
  cout << "target_tau = " << target_tau << endl;
  cout << "photon.scat_weight = " << photon.scat_weight << endl;
#endif

//   int j = 0;
//   for (j = 0; j < 3; j++) 
//     cout << photon.position[j] << " ";
//   cout << endl;

  // check if the photon has left the dust
  if ((target_tau - tau_traveled) > ROUNDOFF_ERR_TRIG)
    escape = 1;
  else
    escape = 0;

#ifdef DEBUG_NS
  cout << "ns escape = " << escape << endl;
#endif
#ifdef DEBUG_NS
  if ((target_tau - tau_traveled) < -ROUNDOFF_ERR_TRIG) {
    cout << "*****error*****next_scatter*****" << endl;
    cout << "target_tau = " << target_tau << endl;
    cout << "tau_traveled = " << tau_traveled << endl;
    cout << "diff = " << target_tau - tau_traveled << endl;
  }
#endif
  
  // return escape (1 = yes, 0 = no)
  return(escape);
}
