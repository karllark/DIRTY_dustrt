// ======================================================================
//   Procedure to setup the absorbed energy grid for the current wavelength.
// Should allow either memory (speed) or disk (space) to be used.
//
// 2007 Jun/KDG - written
// ======================================================================
#include "setup_absorbed_energy_grid.h"

void setup_absorbed_energy_grid (geometry_struct& geometry,
				 int wave_index,
				 int doing_dust_emission)

{
  if (geometry.abs_energy_storage_type == 0)
    geometry.abs_energy_wave_index = wave_index;
  else
    geometry.abs_energy_wave_index = 0;

  int i,j,k,m = 0;
  // loop over all the defined grids
  for (m = 0; m < int(geometry.grids.size()); m++) {
    // loop of the cells in this grid
    for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
      for (j = 0; j < geometry.grids[m].index_dim[1]; j++)
	for (i = 0; i < geometry.grids[m].index_dim[0]; i++) {
	  if (doing_dust_emission) {
	    // save the existing value of the radiation field (needed for dust emission)
	    geometry.grids[m].grid(i,j,k).save_radiation_field_density =
	      geometry.grids[m].grid(i,j,k).absorbed_energy[geometry.abs_energy_wave_index];
	    // zero out the current absorbed energy value
	    geometry.grids[m].grid(i,j,k).absorbed_energy[geometry.abs_energy_wave_index] = 0.0;
	  } else {
	    // if this is the first time or if memory storage requested
	    // then push zero into the absorbed_energy variable in each grid cell
	    if ((!geometry.abs_energy_grid_initialized) || (geometry.abs_energy_storage_type == 0))
	      geometry.grids[m].grid(i,j,k).absorbed_energy.push_back(0.0);
	    else // otherwise set the first element in the absorbed_energy vector to zero
	      geometry.grids[m].grid(i,j,k).absorbed_energy[0] = 0.0;
	  }
	  geometry.grids[m].grid(i,j,k).save_radiation_field_density = 0.0;
	}
  }
  // set to know that the absorbed energy grid has been initialized the first time
  // each cell has an absorbed_energy vector of 1 element
  if (!geometry.abs_energy_grid_initialized) geometry.abs_energy_grid_initialized = 1;
}
