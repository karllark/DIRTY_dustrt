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
  int tmp_save_pindex = -1;
  while (not done) {
    int i = 0;
    for (i = 0; i < 3; i++) {
      photon.position_index[k][i] = int((photon.position[i] - geometry.grids[cur_grid_num].positions[i][0])/
					geometry.grids[cur_grid_num].phys_cube_size[i]);
      tmp_save_pindex = photon.position_index[k][i];
      // make sure that photon is headed into the cell indexed (edges need special treatment)
      if ((photon.position[i] == geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]]) &&
	  (photon.dir_cosines[i] < 0.0))
	photon.position_index[k][i]--;
      else if ((photon.position[i] == geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]+1]) &&
	       (photon.dir_cosines[i] > 0.0))
	photon.position_index[k][i]++;

      // check that the photon is in the grid cell just determined
      if ((photon.position[i] < geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]]) ||
	  (photon.position[i] > geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]+1])) {
	// now make changes if this is a roundoff error issue (very near the grid min/max)
	if ((photon.position[i] < geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]]) &&
	    ((geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]] - photon.position[i]) < ROUNDOFF_ERR_INDEX)) {
	  photon.position_index[k][i]--;
	} else if ((photon.position[i] > geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]+1]) &&
		   ((photon.position[i] - geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]+1]) < ROUNDOFF_ERR_INDEX)) {
	  photon.position_index[k][i]++;
	} else {
	  cout << "photon not in designated grid cell (determine_photon.position_index_initial.cpp)." << endl;
	  cout << "out of bounds of correction." << endl;
	  cout << "photon # = " << photon.number << endl;
	  cout << "i = " << i << endl;
	  cout << "photon.position[i] = " << photon.position[i] << endl;
	  cout << "photon.birth_position[i] = " << photon.birth_position[i] << endl;
	  cout << "photon.dir_consines[i] = " << photon.dir_cosines[i] << endl;
	  cout << "photon.num_scat = " << photon.num_scat << endl;
	  cout << "photon.position_index[k][i] = " << photon.position_index[k][i] << endl;
	  cout << "tmp_save_pindex = " << tmp_save_pindex << endl;
	  cout << "real value of pindex = " << (photon.position[i] - geometry.grids[cur_grid_num].positions[i][0])/
	    geometry.grids[cur_grid_num].phys_cube_size[i] << endl;
	  cout << "cell min position = " << geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]] << endl;
	  cout << "cell max position = " << geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]+1] << endl;
	  cout << "geometry.grids[cur_grid_num].positions[i][0] = " << geometry.grids[cur_grid_num].positions[i][0] << endl;
	  cout << "geometry.grids[cur_grid_num].phys_cube_size[i] = " << geometry.grids[cur_grid_num].phys_cube_size[i] << endl;
	  exit(8);
	}
      }
      
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
