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
					  double target_tau,
					  int& escape,
					  double& tau_traveled)

{
  double tau_left = target_tau;  // reduce till zero = done
  double distance_traveled = 0.0;
  escape = 0;
  
  // move through the grid until escaping or reaching the target tau
  // use the photon record with the previously calculated track
  int i = 0;
  int fin_i = -1;
  while ((tau_left > ROUNDOFF_ERR_TRIG) && (!escape)) {

    if (photon_track.path_tau[i] > tau_left) { // photon interacts inside the cell
      distance_traveled += (tau_left/photon_track.path_tau[i])*photon_track.path_distance[i];
      tau_left = 0.0;
      fin_i = i;
    } else {  // transverse entire cell
      tau_left -= photon_track.path_tau[i];
      distance_traveled += photon_track.path_distance[i];
    }

    cout << photon_track.path_cur_cells << endl;
    cout << i << " " << tau_left << " " << escape << endl;
    
    if (i++ >= photon_track.path_cur_cells) {
      // reached the end of the track - escaped from model
      escape = 1;
      fin_i = i;
    } else
      i++;
  }

  tau_traveled = target_tau - tau_left;

  if (fin_i >= 0) {
    // update the photon position for the distance traveled
    photon.current_grid_num = photon_track.path_pos_index[0][fin_i];
    photon.num_current_grids = photon_track.path_num_current_grids[fin_i];
    int k = photon_track.path_current_grid_num[fin_i];
    for (i = 0; i < 3; i++) {
      photon.position[i] += distance_traveled*photon.dir_cosines[i];
      photon.position_index[k][i] = photon_track.path_pos_index[1+i][fin_i];
    }

    // length of new track
    photon.path_cur_cells = fin_i + 1;

    // transfer the photon portion of the photon track information
    if (photon.path_cur_cells >= photon.path_max_cells)
      for (i = photon.path_max_cells - 1; i < photon.path_cur_cells; i++) {
	// lengthen the photon.path variables
	photon.path_tau.push_back(0.0);
	photon.path_distance.push_back(0.0);
	photon.path_num_current_grids.push_back(0);
	photon.path_current_grid_num.push_back(0);
	photon.path_pos_index[0].push_back(0);
	photon.path_pos_index[1].push_back(0);
	photon.path_pos_index[2].push_back(0);
	photon.path_pos_index[3].push_back(0);
	photon.path_max_cells++;
      } 

    // now fill the photon track with the cached version up to the point of interaction
    for (i = 0; i <= fin_i; i++) {
      for (int j = 0; j < 4; j++)
	photon.path_pos_index[j][i] = photon_track.path_pos_index[j][i];
      photon.path_num_current_grids[i] = photon_track.path_num_current_grids[i];
      photon.path_current_grid_num[i] = photon_track.path_current_grid_num[i];
      photon.path_tau[i] = photon_track.path_tau[i];
      photon.path_distance[i] = photon_track.path_distance[i];
    }

    
  } else {
    cout << "calc_photon_trajectory_from_track: fin_i not set" << endl;
    exit(8);
  }

  return(distance_traveled);
}
