// ======================================================================
//   Procedure to get the dust thermal emission spectra at each point in the
// model.  Bascially, loops over the grid(s) and calls Karl M.'s routine
// to get the dust emission from each grid cell.
//
// 2007 Oct/KDG - written
// ======================================================================
#include "get_dust_thermal_emission.h"

void get_dust_thermal_emission (geometry_struct& geometry,
				runinfo_struct& runinfo, 
				GrainModel& CurGrainModel)

{
  int i,j,k,m,z = 0;
  uint x = 0;

  // get the number of emission components
  int n_emit_components = 1;
  if (runinfo.do_emission_grain) n_emit_components += 2*CurGrainModel.getNComp();
  cout << "n_emit_components = " << n_emit_components << endl;

  // initialize the total emitted sum (by component and wavelength)
  if (!geometry.emitted_energy_grid_initialized) {
    runinfo.emitted_lum.resize(n_emit_components);
    for (z = 0; z < n_emit_components; z++)
      runinfo.emitted_lum[z].resize(runinfo.wavelength.size(),0.0);
  } else {
    for (z = 0; z < n_emit_components; z++)
      for (x = 0; x < runinfo.wavelength.size(); x++)
	runinfo.emitted_lum[z][x] = 0.0;
  }
	    
  // loop through the dust density grid and convert to radiation field density
  // loop over all the defined grids
  for (m = 0; m < int(geometry.grids.size()); m++) {
    // loop of the cells in this grid
    for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
      for (j = 0; j < geometry.grids[m].index_dim[1]; j++)
	for (i = 0; i < geometry.grids[m].index_dim[0]; i++) {

	  // setup emitted energy array
	  if (!geometry.emitted_energy_grid_initialized) {
	    geometry.grids[m].grid(i,j,k).emitted_energy.resize(n_emit_components);
	    for (z = 0; z < n_emit_components; z++)
	      geometry.grids[m].grid(i,j,k).emitted_energy[z].resize(runinfo.wavelength.size(),0.0);
	  } else
	    for (z = 0; z < n_emit_components; z++)
	      for (x = 0; x < runinfo.wavelength.size(); x++)
		geometry.grids[m].grid(i,j,k).emitted_energy[z][x] = 0.0;

	  // determine if there is any energy absorbed needing to be remitted
	  double tot_abs_energy = 0.0;
	  for (x = 0; x < runinfo.wavelength.size(); x++)
	    tot_abs_energy += geometry.grids[m].grid(i,j,k).absorbed_energy[x];
		
	  if (tot_abs_energy > 0.) {
	  // get the dust emission spectrum given the input wavlength vector and radiation field vector
// 	  dust_thermal_emission(CurGrainModel, runinfo.wavelength,
// 				geometry.grids[m].grid(i,j,k).absorbed_energy,
// 				runinfo.do_emission_grain,
// 				geometry.grids[m].grid(i,j,k).emitted_energy);

	  // add to the emitted energy sums
	  // **probably need to think about units**
	    for (z = 0; z < n_emit_components; z++)
	      for (x = 0; x < runinfo.wavelength.size(); x++) {
		// temp stuff
		geometry.grids[m].grid(i,j,k).emitted_energy[z][x] = 1.0;
		// perm stuff
		runinfo.emitted_lum[z][x] += geometry.grids[m].grid(i,j,k).emitted_energy[z][x];
	      }
	  }
	}
    
  }
  geometry.emitted_energy_grid_initialized = 1;
}
