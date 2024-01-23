// ======================================================================
// Check a grid to make sure that all subdivided cells have an assocaited
// subgrid.  Use recursion to search subgrids of subgrids, etc.
//
// KDG Oct 2009 - written
// ======================================================================
#include "setup_dust_grid_check_grid.h"
// #define DEBUG_SDGCG

void setup_dust_grid_check_grid(geometry_struct& geometry, int cur_grid,
                                int par_grid, vector<int>& par_idim)

{
  int i, j, k = 0;
  float dust_tau_per_pc = 0.0;
  bool errors_found = false;
  float min_tau = 0.0;
  float max_tau = 0.0;

  //   cout << "checking grid # = " << cur_grid << " ";
  //   cout << "out of " << geometry.grids.size() << " grids" << endl;

  // check the parent grid is correctly indicated
  if (geometry.grids[cur_grid].parent_grid_num != par_grid) {
    cout << "parent grid number set to "
         << geometry.grids[cur_grid].parent_grid_num;
    cout << " but should be " << par_grid << endl;
    cout << "for current grid index = " << cur_grid << endl;
    errors_found = true;
  }

  // check if the min/max positions of the subgrid match the positions of the
  // parent grid cell
  if (par_grid >= 0) {
    for (i = 0; i < 3; i++) {
      // check the min values
      if (geometry.grids[cur_grid].positions[i][0] !=
          geometry.grids[par_grid].positions[i][par_idim[i]]) {
        cout << "min subgrid i = " << i
            << " position does not match parent cell position" << endl;
        cout << "subgrid min = " << geometry.grids[cur_grid].positions[i][0];
        cout << "; parent cell min = "
            << geometry.grids[par_grid].positions[i][par_idim[i]] << endl;
        cout << "subgrid number = " << cur_grid
            << "; parent number = " << par_grid << endl;
        errors_found = true;
      }
      // check the max values
      if (geometry.grids[cur_grid].positions[i][geometry.grids[cur_grid].index_dim[i]] !=
          geometry.grids[par_grid].positions[i][par_idim[i]+1]) {
        cout << "max subgrid i = " << i
            << " position does not match parent cell position" << endl;
        cout << "subgrid max = " << geometry.grids[cur_grid].positions[i][geometry.grids[cur_grid].index_dim[i]];
        cout << "; parent cell max = "
            << geometry.grids[par_grid].positions[i][par_idim[i]+1] << endl;
        cout << "subgrid number = " << cur_grid
            << "; parent number = " << par_grid << endl;
        errors_found = true;
      }
    }
  }

  for (k = 0; k < geometry.grids[cur_grid].index_dim[2]; k++)
    for (j = 0; j < geometry.grids[cur_grid].index_dim[1]; j++) {
      for (i = 0; i < geometry.grids[cur_grid].index_dim[0]; i++) {
        dust_tau_per_pc =
            geometry.grids[cur_grid].grid(i, j, k).dust_tau_per_pc;
        if (!isfinite(dust_tau_per_pc)) {
          cout << "Non-finite value of " << dust_tau_per_pc << " detected in grid # = " << cur_grid;
          cout << " at cell = (" << i << "," << j << "," << k << ")" << endl;
          errors_found = true;
        } else if (dust_tau_per_pc < -0.5) {
          if (abs(dust_tau_per_pc) <= geometry.grids.size()) {
            par_idim[0] = i;
            par_idim[1] = j;
            par_idim[2] = k;
            setup_dust_grid_check_grid(geometry, abs(dust_tau_per_pc), cur_grid,
                                       par_idim);
          } else {
            cout << "subgrid # = " << abs(dust_tau_per_pc) << " does not exist."
                 << endl;
            cout << "in grid # = " << cur_grid << endl;
            errors_found = true;
          }
        } else {
          if (dust_tau_per_pc < min_tau)
            min_tau = dust_tau_per_pc;
          if (dust_tau_per_pc > max_tau)
            max_tau = dust_tau_per_pc;
        }
        // 	cout << endl;
      }
    }

  // check that the grid is not all filled with -0.5 values
  // -0.5 and subgrids allowed
    if ((max_tau == -0.5) && (min_tau == -0.5)) {
      cout << "grid # = " << cur_grid << " is filled with only -0.5 values, not allowed" << endl;
      errors_found = true;
    }

  // stop if any errors found
  if (errors_found == true)
    exit(8);
}
