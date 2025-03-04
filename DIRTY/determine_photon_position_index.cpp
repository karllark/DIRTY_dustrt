// ======================================================================
//   Procedure to determine the position indexes of a photon given
// it's physical position and the active grid.  This involves determining
// the grid and subgrid indexes.
//
// This procedure assumes that the position indexes in the current grid
// are know and only the subgrid indexes need determining.
//
// See determine_photon_position_index_initial if the indexes in the
// current grid are not known.
//
// 2003 Jun/KDG - written
// ======================================================================
#include "determine_photon_position_index.h"
// #define OUTNUM 2567
// #define DEBUG_DPPI

void
determine_photon_position_index (geometry_struct &geometry, photon_data &photon)

{
  // find the position in the grid (assuming grid is made of equal cubes)

  int done = 0;
  int k = photon.current_grid_num;
  int cur_grid_num = photon.grid_number[k];

  int k1, k2, mid_k;
  while (not done)
    {
      int i = 0;
#ifdef DEBUG_DPPI
      if (photon.number == OUTNUM)
        {
          cout << endl;
          cout << "in dppi" << endl;
        }
#endif
      for (i = 0; i < 3; i++)
        {
          // now find the index
          // binary search for index of position in grid
          k1 = 0;
          k2 = geometry.grids[cur_grid_num].index_dim[i];
          // cout << k1 << " " << k2 << endl;
          while ((k2 - k1) > 1)
            {
              mid_k = int (0.5 * (k1 + k2));
              if (photon.position[i] < geometry.grids[cur_grid_num].positions[i][mid_k])
                k2 = mid_k;
              else
                k1 = mid_k;
              // cout << k1 << " " << k2 << endl;
            }
          photon.position_index[k][i] = k1;
          // photon.position_index[k][i] = int((photon.position[i] -
          // geometry.grids[cur_grid_num].positions[i][0])/
          //                                  geometry.grids[cur_grid_num].phys_cube_size[i]);

          // take care of case when we are at the edge of the grid in the
          // maximum sense
          if (photon.position_index[k][i] == geometry.grids[cur_grid_num].index_dim[i])
            photon.position_index[k][i]--;
#ifdef DEBUG_DPPI
          if (photon.number == OUTNUM)
            {
              cout << "before adjust" << endl;
              cout << "i = " << i << endl;
              cout << photon.number << endl;
              cout << photon.position_index[k][i] << " ";
              cout << geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]] << " ";
              cout << photon.position[i] << " ";
              cout << geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i] + 1] << " ";
              cout << photon.dir_cosines[i] << endl;
            }
#endif
          // make sure that photon is headed into the cell indexed (edges need
          // special treatment)
          if ((photon.position[i] == geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]])
              && (photon.dir_cosines[i] < 0.0) && (photon.position_index[k][i] > 0))
            photon.position_index[k][i]--;
          else if ((photon.position[i] == geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i] + 1])
                   && (photon.dir_cosines[i] > 0.0))
            photon.position_index[k][i]++;
          // take care of case when we are at the edge of the grid in the
          // maximum sense
          if (photon.position_index[k][i] == geometry.grids[cur_grid_num].index_dim[i])
            photon.position_index[k][i]--;
#ifdef DEBUG_DPPI
          if (photon.number == OUTNUM)
            {
              cout << "after adjust" << endl;
              cout << "i = " << i << endl;
              cout << photon.position_index[k][i] << " ";
              cout << geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i]] << " ";
              cout << photon.position[i] << " ";
              cout << geometry.grids[cur_grid_num].positions[i][photon.position_index[k][i] + 1] << " ";
              cout << photon.dir_cosines[i] << endl;
            }
#endif
          // make sure the photon is in this grid
          if ((photon.position[i] < geometry.grids[cur_grid_num].positions[i][0])
              || (photon.position[i]
                  > geometry.grids[cur_grid_num].positions[i][geometry.grids[cur_grid_num].index_dim[i]]))
            {
              // now make changes if this is a roundoff error issue (very near
              // the grid min/max)
              int not_ok = 0;
              if (photon.position[i] < geometry.grids[cur_grid_num].positions[i][0])
                {
                  float frac_miss
                      = (geometry.grids[cur_grid_num].positions[i][0] - photon.position[i])
                        / (geometry.grids[cur_grid_num].positions[i][1] - geometry.grids[cur_grid_num].positions[i][0]);
                  cout << "frac_miss-- = " << frac_miss << endl;
                  if (frac_miss < ROUNDOFF_ERR_INDEX)
                    photon.position[i] = geometry.grids[cur_grid_num].positions[i][0];
                  // 	    photon.position_index[k][i]--;
                  else
                    not_ok = 1;
                }
              else if (photon.position[i]
                       > geometry.grids[cur_grid_num].positions[i][geometry.grids[cur_grid_num].index_dim[i]])
                {
                  float frac_miss
                      = (photon.position[i]
                         - geometry.grids[cur_grid_num].positions[i][geometry.grids[cur_grid_num].index_dim[i]])
                        / (geometry.grids[cur_grid_num].positions[i][geometry.grids[cur_grid_num].index_dim[i]]
                           - geometry.grids[cur_grid_num].positions[i][geometry.grids[cur_grid_num].index_dim[i] - 1]);
                  cout << "frac_miss++ = " << frac_miss << endl;
                  if (frac_miss < ROUNDOFF_ERR_INDEX)
                    photon.position[i]
                        = geometry.grids[cur_grid_num].positions[i][geometry.grids[cur_grid_num].index_dim[i]];
                  // 	    photon.position_index[k][i]++;
                  else
                    not_ok = 1;
                }

              if (not_ok)
                {
                  cout << "outside of bounds of current grid (in "
                          "determine_photon_position_index.cpp)"
                       << endl;
                  cout << "This should not happen." << endl;
                  cout << "photon # = " << photon.number << endl;
                  cout << "i = " << i << endl;
                  cout << "photon.position[i] = " << photon.position[i] << endl;
                  cout << "photon.birth_position[i] = " << photon.birth_position[i] << endl;
                  cout << "photon.dir_consines[i] = " << photon.dir_cosines[i] << endl;
                  cout << "photon.num_scat = " << photon.num_scat << endl;
                  cout << "cur_grid_num = " << cur_grid_num << endl;
                  cout << "min/max of grid = " << geometry.grids[cur_grid_num].positions[i][0] << " ";
                  cout << geometry.grids[cur_grid_num].positions[i][geometry.grids[cur_grid_num].index_dim[i]] << endl;
                  //	  cout << "size of cube [pc] = " <<
                  // geometry.grids[cur_grid_num].phys_cube_size[i] << endl;
                  // cout << "proposed index = " << int((photon.position[i] -
                  // geometry.grids[cur_grid_num].positions[i][0])/
                  // 				     geometry.grids[cur_grid_num].phys_cube_size[i])
                  // << endl; cout << "real value of pindex = " <<
                  // (photon.position[i] -
                  // geometry.grids[cur_grid_num].positions[i][0])/
                  //   geometry.grids[cur_grid_num].phys_cube_size[i] << endl;
                  cout << "out of a possible " << geometry.grids[cur_grid_num].index_dim[i] << endl;
                  cout << "out of a possible2 " << geometry.grids[cur_grid_num].positions[i].size () << endl;
                  int parent_grid_num = geometry.grids[cur_grid_num].parent_grid_num;
                  cout << "parent grid cell min = "
                       << geometry.grids[parent_grid_num].positions[i][photon.position_index[k - 1][i]] << endl;
                  cout << "parent grid cell max = "
                       << geometry.grids[parent_grid_num].positions[i][photon.position_index[k - 1][i] + 1] << endl;
                  exit (8);
                }
            }
        }

      if (geometry.grids[cur_grid_num]
              .grid (photon.position_index[k][0], photon.position_index[k][1], photon.position_index[k][2])
              .dust_tau_per_pc
          <= -1.0)
        {
          cur_grid_num
              = -int (geometry.grids[cur_grid_num]
                          .grid (photon.position_index[k][0], photon.position_index[k][1], photon.position_index[k][2])
                          .dust_tau_per_pc);
          photon.grid_number[k + 1] = cur_grid_num;
          //       cout << "subgrid! =" << " ";
          //       cout << cur_grid_num << endl;
        }
      else
        {
          done = 1;
        }

      // increment counter
      k = k + 1;
    }

  // record the number of current grids the photon is nested in
  photon.num_current_grids = k;

#ifdef DEBUG_DPPI
  if (photon.number == OUTNUM)
    {
      cout << "end of dppi" << endl;
    }
#endif
}
