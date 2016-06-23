// ======================================================================
//   Function to calculate the trajectory of a photon.  The stopping
// point is either reaching the target optical depth (target_tau) or 
// reaching the edge of the distribution (marked by negative optical
// depth between (0 and, but not including -1: integer negative 
// optical depths denote subgrid indexes.
//
// This function returns the distance traveled.
//
// This function uses a previously calculated photon track to save
// computations.  Usually, such a track has been calculated to the
// surface for the forced scattering or continuous absorption.
//
// 2016 Jun/KDG - written
// ======================================================================
#include "calc_photon_trajectory_from_track.h"

double calc_photon_trajectory_from_track (photon_data& photon,
					  photon_data& photon_track,
					  double target_tau)

{
  double tau_left = target_tau;  // reduce till zero = done
  double distance_traveled = 0.0;

  // move through the grid until escaping or reaching the target tau
  // use the photon record with the previously calculated track
  int i = 0;
  int fin_i = -1;
  while (tau_left > ROUNDOFF_ERR_TRIG) {

    if (photon_track.path_tau[i] > tau_left) { // photon interacts inside the cell
      distance_traveled += (tau_left/photon_track.path_tau[i])*photon_track.path_distance[i];
      tau_left = 0.0;
      fin_i = i;
    } else {  // transverse entire cell
      tau_left -= photon_track.path_tau[i];
      distance_traveled += photon_track.path_distance[i];
    }
    
    if (i >= photon_track.path_cur_cells) {
      cout << "calc_photon_trajectory_from_track: photon traveled to " << endl;
      cout << "the end of track before reaching the desired tau" << endl;
      exit(8);
    } else
      i++;
  }

  if (fin_i >= 0) {
    // update the photon position for the distance traveled
    photon.current_grid_num = photon_track.path_pos_index[0][fin_i];
    photon.num_current_grids = photon_track.path_num_current_grids[fin_i];
    int k = photon_track.path_current_grid_num[fin_i];
    for (i = 0; i < 3; i++) {
      photon.position[i] += distance_traveled*photon.dir_cosines[i];
      photon.position_index[k][i] = photon_track.path_pos_index[1+i][fin_i];
    }
  } else {
    cout << "calc_photon_trajectory_from_track: fin_i not set" << endl;
    exit(8);
  }

  return(distance_traveled);
}
