#include "check_absorbed_energy_grid.h"
//#define DEBUG_SAEG

void check_absorbed_energy_grid (geometry_struct& geometry,
				 runinfo_struct& runinfo)

{
  // loop through the dust density grid and convert to radiation field density
  int i,j,k,m = 0;
  uint x;
  for (x = 0; x < runinfo.wavelength.size(); x++) {
    // loop over all the defined grids
    double total_abs_energy = 0.0;
    total_abs_energy = 0.0;

    for (m = 0; m < int(geometry.grids.size()); m++) {
      
      // loop of the cells in this grid
      for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
	for (j = 0; j < geometry.grids[m].index_dim[1]; j++)
	  for (i = 0; i < geometry.grids[m].index_dim[0]; i++) {
	    total_abs_energy += geometry.grids[m].grid(i,j,k).absorbed_energy[x]*
	      Constant::FPI*runinfo.ave_C_abs[x]*geometry.grids[m].grid(i,j,k).num_H;
	    if (x == 1) cout << geometry.grids[m].grid(i,j,k).num_H << " ";
	  }
    }
    cout << "caeg; x = " << x << " ";
    cout << runinfo.wavelength[x] << " ";
    cout << total_abs_energy << " ";
    cout << runinfo.sed_lum[x] << " ";
    cout << geometry.abs_energy_wave_index << " ";
    cout << endl;
  }

}
