// ======================================================================
//   Procedure to verify the integrity of the dust grid.  Used as a debugging
// tool for now.  Will use to verify user input grids later.
//
// 2008 Jun/KDG - written
// ======================================================================
#include "verify_dust_grid.h"
//#define DEBUG_VDG

void verify_dust_grid (geometry_struct& geometry)

{
  // loop through the dust density grid and make sure the subgrids are all 
  // fully contained in a single cell of the parent gird
  int i,j,k,m = 0;
  int subgrid_index = 0;
  if (geometry.grids.size() > 1) {
    for (m = 0; m < int(geometry.grids.size()); m++) {
      
      // loop of the cells in this grid
      for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
	for (j = 0; j < geometry.grids[m].index_dim[1]; j++)
	  for (i = 0; i < geometry.grids[m].index_dim[0]; i++) {
	    if (geometry.grids[m].grid(i,j,k).dust_tau_per_pc <= -1.0) {
	      subgrid_index = int(-1.*geometry.grids[m].grid(i,j,k).dust_tau_per_pc);
	      if ((geometry.grids[m].positions[0][i] != geometry.grids[subgrid_index].positions[0][0]) ||
		  (geometry.grids[m].positions[0][i+1] != geometry.grids[subgrid_index].positions[0][geometry.grids[subgrid_index].index_dim[0]])) {
		cout << "bad subgrid beginning x postion (does not match parent grid cell min)" << endl;
		cout << "dust_tau_per_pc = " << geometry.grids[m].grid(i,j,k).dust_tau_per_pc << endl;
		cout << "n_grids = " << geometry.grids.size() << endl;
		cout << "parent cell (m,i,j,k) = " << m << " " << i << " " << j << " " << k << endl;
		cout << "parent min/max = " << geometry.grids[m].positions[0][i] << " " << geometry.grids[m].positions[0][i+1] << endl;
		cout << "subgrid min/max = " << geometry.grids[subgrid_index].positions[0][0] << " ";
		cout << geometry.grids[subgrid_index].positions[0][geometry.grids[subgrid_index].index_dim[0]] << endl;
		exit(8);
	      }
	    }
	  }
    }
  }

}
