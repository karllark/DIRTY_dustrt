// ======================================================================
//   Procedure to get the dust thermal emission spectra at each point in the
// model.  Bascially, loops over the grid(s) and calls Karl M.'s routine
// to get the dust emission from each grid cell.
//
// 2007 Oct/KDG - written
// 2008 Mar/KDG - updated calculation to do the units of the emitted energy correctly
// ======================================================================
#include "get_dust_thermal_emission.h"

void get_dust_thermal_emission (geometry_struct& geometry,
				runinfo_struct& runinfo, 
				GrainModel& CurGrainModel)

{
  int i,j,k,m,z = 0;
  uint x = 0;

  bool DoStochastic=false; 

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
	    // emitted energy returned is in units of ergs s^-1 HI atom^-1

	    ComputeDustEmission(geometry.grids[m].grid(i,j,k).absorbed_energy,
				CurGrainModel, 
				geometry.grids[m].grid(i,j,k).emitted_energy,
				DoStochastic); 

// 	    dust_thermal_emission(CurGrainModel, runinfo.wavelength,
// 				  geometry.grids[m].grid(i,j,k).absorbed_energy,
// 				  runinfo.do_emission_grain,
// 				  geometry.grids[m].grid(i,j,k).emitted_energy);

	    // add to the emitted energy sums
	    for (z = 0; z < n_emit_components; z++)
	      for (x = 0; x < runinfo.wavelength.size(); x++) {
		// temp stuff
		geometry.grids[m].grid(i,j,k).emitted_energy[z][x] = 1.0;

		// perm stuff
		// need to multiply the emitted energy passed back by dust_thermal_emission
		// by the number of HI atoms in the cell to get the total emitted energy
		geometry.grids[m].grid(i,j,k).emitted_energy[z][x] *= geometry.grids[m].grid(i,j,k).num_H;
		// add up the emitted energy to the total emitted (per wavelength & component)
		runinfo.emitted_lum[z][x] += geometry.grids[m].grid(i,j,k).emitted_energy[z][x];
	      }
	  }
	}
    
  }
  geometry.emitted_energy_grid_initialized = 1;
}