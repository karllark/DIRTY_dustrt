// ======================================================================
//   Function to calculate the trajectory of a photon.  The stopping
// point is either reaching the target optical depth (target_tau) or 
// reaching the edge of the distribution (marked by negative optical
// depth between (0 and, but not including -1: integer negative 
// optical depths denote subgrid indexes.
//
// This function returns the distance traveled.
//
// 2003 Jun/KDG - written
// ======================================================================
#include "calc_photon_trajectory.h"
// #define PHOTON_POS
//#define DEBUG_CPT

double calc_photon_trajectory (photon_data& photon,
			       geometry_struct& geometry,
			       double target_tau,
			       int& escape,
			       double& tau_traveled)

{
  double distance_traveled = 0.0;
  double tau_left = target_tau;  // reduce till zero = done
  double delta_tau = 0.0;
  double delta_dist = 0.0;

  // move through the grid until escaping or reaching the target tau
  //   need to also handle the case where the photon starts several subgrids down
  //   multiple calls to calc_photon_trajectory needed
  while ((tau_left > 0.0) && (!escape)) {

#ifdef DEBUG_CPT
    if (photon.number > OUTNUM) {
      cout << "before delta distance & tau_left & delta_tau & escape = ";
      cout << delta_dist << " ";
      cout << tau_left << " ";
      cout << delta_tau << " ";
      cout << escape << endl;
      cout << "current_grid_num, num_current grids = ";
      cout << photon.current_grid_num << " ";
      cout << photon.num_current_grids << endl;
    }
#endif
    delta_dist = calc_delta_dist(photon, geometry, tau_left, escape, delta_tau);
#ifdef PHOTON_POS
    int i = 0;
    for (i = 0; i < 3; i++)
      cout << photon.position[i] << " ";
    cout << endl;
#endif    
    tau_traveled += delta_tau;
    tau_left -= delta_tau;
    distance_traveled += delta_dist;
#ifdef DEBUG_CPT
    if (photon.number > OUTNUM) {
      cout << "delta distance & tau_left & delta_tau & escape = ";
      cout << delta_dist << " ";
      cout << tau_left << " ";
      cout << delta_tau << " ";
      cout << escape << endl;
      cout << "current_grid_num, num_current grids = ";
      cout << photon.current_grid_num << " ";
      cout << photon.num_current_grids << endl;
    }
#endif
  }

  return(distance_traveled);
}
