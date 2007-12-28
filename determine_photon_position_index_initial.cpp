// ======================================================================
//   Procedure to determine the position indexes of a photon given
// it's physical position and the active grid.  This involves determining
// the grid and subgrid indexes.
//
// 2003 Jun/KDG - written
// ======================================================================
#include "determine_photon_position_index_initial.h"

void determine_photon_position_index_initial (geometry_struct& geometry,
					      photon_data& photon)

{
  // find the position in the grid (assuming grid is made of equal cubes)

  int done = 0;
  int k = 0;
  int cur_grid_num = 0;
  while (not done) {
    int i = 0;
    for (i = 0; i < 3; i++) {
      photon.position_index[k][i] = int((photon.position[i] - geometry.grids[cur_grid_num].positions[i][0])/
					geometry.grids[cur_grid_num].phys_cube_size[i]);
      // make sure that photon is headed into the cell indexed (edges need special treatment)
      if ((photon.position[i] == geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]]) &&
	  (photon.dir_cosines[i] < 0.0))
	photon.position_index[k][i]--;
      else if ((photon.position[i] == geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]+1]) &&
	       (photon.dir_cosines[i] > 0.0))
	photon.position_index[k][i]++;

//       cout << photon.position_index[k][i] << " ";
//       cout << geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]] << " ";
//       cout << photon.position[i] << " ";
//       cout << geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]+1] << endl;
    }

    if (geometry.grids[cur_grid_num].grid(photon.position_index[k][0],photon.position_index[k][1],photon.position_index[k][2]).dust_tau_per_pc <= -1.0) {
      cur_grid_num = -int(geometry.grids[cur_grid_num].grid(photon.position_index[k][0],photon.position_index[k][1],photon.position_index[k][2]).dust_tau_per_pc);
      photon.grid_number[k+1] = cur_grid_num;
//       cout << "subgrid! =" << " ";
//       cout << cur_grid_num << endl;
    } else {
      done = 1;
    }

    // increment counter
    k = k + 1;

  }

  // record the number of current grids the photon is nested in
  photon.num_current_grids = k;
  photon.current_grid_num = 0;

}
