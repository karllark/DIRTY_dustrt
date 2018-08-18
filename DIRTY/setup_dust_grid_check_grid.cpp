// ======================================================================
// Check a grid to make sure that all subdivided cells have an assocaited
// subgrid.  Use recursion to search subgrids of subgrids, etc.
//
// KDG Oct 2009 - written
// ======================================================================
#include "setup_dust_grid_check_grid.h"
//#define DEBUG_SDGCG

void setup_dust_grid_check_grid (geometry_struct& geometry,
				 int cur_grid)

{
  int i,j,k = 0;
  float dust_tau_per_pc = 0.0;

//   cout << "checking grid # = " << cur_grid << " ";
//   cout << "out of " << geometry.grids.size() << " grids" << endl;

  for (k = 0; k < geometry.grids[cur_grid].index_dim[2]; k++)
    for (j = 0; j < geometry.grids[cur_grid].index_dim[1]; j++) {
      for (i = 0; i < geometry.grids[cur_grid].index_dim[0]; i++) {
	dust_tau_per_pc = geometry.grids[cur_grid].grid(i,j,k).dust_tau_per_pc;
// 	cout << dust_tau_per_pc << " ";
	if (dust_tau_per_pc < -0.5) {
	  if (abs(dust_tau_per_pc) <= geometry.grids.size())
	    setup_dust_grid_check_grid(geometry, abs(dust_tau_per_pc));
	  else {
	    cout << "subgrid # = " << abs(dust_tau_per_pc) << " does not exist." << endl;
	    cout << "in grid # = " << cur_grid << endl;
	    exit(8);
	  }
	}
// 	cout << endl;
      }
	
      }

}
