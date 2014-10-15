// ======================================================================
//   Procedure to determine the position indexes of a photon given
// it's physical position and the active grid.  This involves determining
// the grid and subgrid indexes.
//
// 2003 Jun/KDG - written
// ======================================================================
#include "determine_photon_position_index_initial.h"
//#define DEBUG_DPPII
//#define OUTNUM 0

void determine_photon_position_index_initial (geometry_struct& geometry,
					      photon_data& photon)

{
  // find the position in the grid (assuming grid is made of equal cubes)

  int done = 0;
  int k = 0;
  photon.grid_number[0] = 0;
  int cur_grid_num = 0;
  int tmp_save_pindex = -1;
  while (not done) {
    int i = 0;
    for (i = 0; i < 3; i++) {
#ifdef DEBUG_DPPII
      if (photon.number == OUTNUM) {
	cout << "i = " << i << endl;
	cout << "photon.position[i] = " << photon.position[i] << endl;
	cout << "grid min = " << geometry.grids[cur_grid_num].positions[i][0] << endl;
	cout << "grid size = " << geometry.grids[cur_grid_num].phys_cube_size[i] << endl;
      }
#endif      

      photon.position_index[k][i] = int((photon.position[i] - geometry.grids[cur_grid_num].positions[i][0])/
					geometry.grids[cur_grid_num].phys_cube_size[i]);
      tmp_save_pindex = photon.position_index[k][i];
      // make sure that photon is headed into the cell indexed (edges need special treatment)
      if ((photon.position[i] == geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]]) &&
	  (photon.dir_cosines[i] <= 0.0))
	photon.position_index[k][i]--;
      else if ((photon.position[i] == geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]+1]) &&
	       (photon.dir_cosines[i] > 0.0))
	photon.position_index[k][i]++;

//       // debug stuff
//       if (photon.position_index[k][i] >= geometry.grids[cur_grid_num].index_dim[i]) {
// 	cout << tmp_save_pindex << endl;
// 	cout << photon.position[i] << endl;
// 	cout << "photon.position[i] = " << photon.position[i] << endl;
// 	cout << "photon.birth_position[i] = " << photon.birth_position[i] << endl;
// 	cout << "photon.dir_consines[i] = " << photon.dir_cosines[i] << endl;
// 	cout << "photon.num_scat = " << photon.num_scat << endl;
// 	cout << "photon.position_index[k][i] = " << photon.position_index[k][i] << endl;
// 	cout << "real value of pindex = " << (photon.position[i] - geometry.grids[cur_grid_num].positions[i][0])/
// 	  geometry.grids[cur_grid_num].phys_cube_size[i] << endl;
// 	cout << "geometry.grids[cur_grid_num].positions[i][0] = " << geometry.grids[cur_grid_num].positions[i][0] << endl;
// 	cout << "geometry.grids[cur_grid_num].phys_cube_size[i] = " << geometry.grids[cur_grid_num].phys_cube_size[i] << endl;
// 	cout << "found the bastard ! " << endl;
// 	exit(8);
//       }
//       // back to regular code

#ifdef DEBUG_DPPII
      if (photon.number == OUTNUM) {
	cout << "pre check (debug)" << endl;
// 	cout << "photon not in designated grid cell (determine_photon.position_index_initial.cpp)." << endl;
// 	cout << "out of bounds of correction." << endl;
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
	cout << "min edge diff = " << geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]] - photon.position[i] << endl;
	cout << "max edge diff = " << photon.position[i] - geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]+1] << endl;
	cout << "ROUNDOFF_ERR_INDEX = " << ROUNDOFF_ERR_INDEX << endl;
      }
#endif

      // check that the photon is in the grid cell just determined
      if ((photon.position[i] < geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]]) ||
	  (photon.position[i] > geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]+1])) {
	// now make changes if this is a roundoff error issue (very near the grid min/max)
	int not_ok = 0;
	if (photon.position_index[k][i] == geometry.grids[cur_grid_num].index_dim[i]) {
#ifdef DEBUG_DPPII
	  cout << "photon # = " << photon.number << " ";
	  cout << "frac_miss--(exact max boundary)" << endl;
#endif
	  photon.position_index[k][i]--;
	  //	  exit(8);
	} else if (photon.position[i] < geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]]) {
	  float frac_miss = (geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]] - photon.position[i])/
	    geometry.grids[cur_grid_num].phys_cube_size[k];
	  cout << "frac_miss--(initial) = " << frac_miss << endl;
	  if (frac_miss < ROUNDOFF_ERR_INDEX)
	    photon.position_index[k][i]--;
	  else
	    not_ok = 1;
	} else if (photon.position[i] > geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]+1]) {
	  float frac_miss = (photon.position[i] - geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]+1])/
	    geometry.grids[cur_grid_num].phys_cube_size[k];
	  cout << "frac_miss++(initial) = " << frac_miss << endl;
          cout << "ROUNDOFF_ERR_INDEX = " << ROUNDOFF_ERR_INDEX << endl;
	  if (frac_miss < ROUNDOFF_ERR_INDEX)
	    photon.position_index[k][i]++;
	  else
	    not_ok = 1;
	}

	if (not_ok) {
	  cout << "photon not in designated grid cell (determine_photon.position_index_initial.cpp)." << endl;
	  cout << "out of bounds of correction." << endl;
	  cout << "photon # = " << photon.number << endl;
	  cout << "i = " << i << endl;
	  cout << "cur_grid_num = " << cur_grid_num << endl;
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
	  cout << "min edge diff/cube size = " << (geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]] - photon.position[i])/
	       geometry.grids[cur_grid_num].phys_cube_size[k] << endl;
	  cout << "max edge diff/cub size = " << (photon.position[i] - geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]+1])/
	       geometry.grids[cur_grid_num].phys_cube_size[k] << endl;
          cout << "ROUNDOFF_ERR_INDEX = " << ROUNDOFF_ERR_INDEX << endl;
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

#ifdef DEBUG_DPPII
  if (photon.number == OUTNUM) {
    cout << "made it to the end.." << endl;
  }
#endif

  // record the number of current grids the photon is nested in
  photon.num_current_grids = k;
  photon.current_grid_num = 0;

}
