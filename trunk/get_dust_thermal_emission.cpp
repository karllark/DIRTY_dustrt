// ======================================================================
//   Procedure to get the dust thermal emission spectra at each point in the
// model.  Bascially, loops over the grid(s) and calls Karl M.'s routine
// to get the dust emission from each grid cell.
//
// 2007 Oct/KDG - written
// 2008 Mar/KDG - updated calculation to do the units of the emitted energy correctly
// 2008 Sept 2 - added explicit cast to NumUtils::integrate
// ======================================================================
#include "get_dust_thermal_emission.h"
//#define DEBUG_GDTE

void get_dust_thermal_emission (geometry_struct& geometry,
				runinfo_struct& runinfo, 
				GrainModel& CurGrainModel)

{
  int i,j,k,m,z = 0;
  uint x = 0;

  bool DoStochastic=false; 

  double global_total_emitted = 0.;
  double global_total_absorbed = 0.;

  // get the number of emission components
  int n_emit_components = 1;
  if (runinfo.do_emission_grain) n_emit_components += 2*CurGrainModel.getNComp();
  if (runinfo.verbose >= 2) 
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

  // calculate the total absorbed energy from the previously saved total abosrbed energy by wavelength
  vector<double> tmp_wave;
  tmp_wave.resize(runinfo.wavelength.size());
  vector<double> tmp_abs_energy;
  tmp_abs_energy.resize(runinfo.wavelength.size());
  for (x = 0; x < runinfo.wavelength.size(); x++) {
    tmp_wave[x] = double(runinfo.wavelength[x]);
    tmp_abs_energy[x] = runinfo.absorbed_energy[x]*4.*Constant::PI*runinfo.ave_C_abs[x];
  }
  global_total_absorbed = NumUtils::integrate<double>(tmp_wave,tmp_abs_energy);

  double min_enough_energy = 0.1*runinfo.energy_conserve_target*global_total_absorbed/geometry.num_cells;
  long num_cells_enough = 0;
  long num_cells_not_enough = 0;
  long num_cells_zero = 0;
  long num_cells_too_few_waves = 0;

#ifdef DEBUG_GDTE 
  cout << "min enough energy = " << min_enough_energy << endl;
#endif

  // define here to save doing it in the loop
  double tot_abs_energy = 0.0;
  double max_abs_energy = 0.0;
  int tot_nonzero = 0;

  // loop over all the defined grids
  for (m = 0; m < int(geometry.grids.size()); m++) {
    // loop of the cells in this grid
    for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
      for (j = 0; j < geometry.grids[m].index_dim[1]; j++) {

	if (runinfo.verbose >= 1)
	  cout << "working on dust emission grid (m,k,j) = " << m << " " << k << " " << j << endl;

	for (i = 0; i < geometry.grids[m].index_dim[0]; i++) {

	  // setup emitted energy array
	  if (!geometry.emitted_energy_grid_initialized) {
	    geometry.grids[m].grid(i,j,k).emitted_energy.resize(n_emit_components);
	    for (z = 0; z < n_emit_components; z++)
	      geometry.grids[m].grid(i,j,k).emitted_energy[z].resize(runinfo.wavelength.size(),0.0);
	  } else
	    // will need code that does not set the nonequilibrium thermal emission to zero
	    // when we are doing the next iteration for the case where we donot want to 
	    // recompute the nonequilibrium case
	    // also make sure to set the total emission to the sum of the nonequilibrium emission
	    for (z = 0; z < n_emit_components; z++)
	      for (x = 0; x < runinfo.wavelength.size(); x++)
		geometry.grids[m].grid(i,j,k).emitted_energy[z][x] = 0.0;

	  // determine if there is any energy absorbed needing to be remitted
	  tot_abs_energy = 0.0;
	  max_abs_energy = 0.0;
	  tot_nonzero = 0;

	  max_abs_energy = 0.0;
	  for (x = 0; x < runinfo.wavelength.size(); x++) {
// 	    tmp_wave[x] = double(runinfo.wavelength[x]);
	    tmp_abs_energy[x] = geometry.grids[m].grid(i,j,k).absorbed_energy[x]*4.*Constant::PI*runinfo.ave_C_abs[x];
	    if (tmp_abs_energy[x] > max_abs_energy) max_abs_energy = tmp_abs_energy[x];
	    if (geometry.grids[m].grid(i,j,k).absorbed_energy[x] > 0.) tot_nonzero++;
	  }
 	  tot_abs_energy = NumUtils::integrate<double>(tmp_wave,tmp_abs_energy)*geometry.grids[m].grid(i,j,k).num_H;
		
#ifdef DEBUG_GDTE
	  cout << "total nonzero = " << tot_nonzero << endl;
	  cout << "total absorbed energy = " << tot_abs_energy << endl;
	  cout << "max absorbed energy = " << max_abs_energy << endl;
#endif

	  if ((tot_abs_energy > 0.0) && (tot_nonzero < int(0.5*runinfo.wavelength.size()))) {

	    num_cells_too_few_waves++;

	  } else if (tot_abs_energy >= min_enough_energy) {

	    num_cells_enough++;
	    
// 	    cout << "total min = " << tot_abs_energy << " " << min_enough_energy << endl;
// 	  if ((tot_abs_energy > 0.) && (tot_nonzero > int(0.75*runinfo.wavelength.size())) && (tot_abs_energy > 1e-35)) {
#ifdef DEBUG_GDTE
	    // output the J
	    for (x = 0; x < runinfo.wavelength.size(); x++) {
	      cout << geometry.grids[m].grid(i,j,k).absorbed_energy[x] << " ";
	    }
	    cout << endl;
#endif

	    // get the dust emission spectrum given the input wavlength vector and radiation field vector
	    // emitted energy returned is in units of ergs s^-1 HI atom^-1
	    ComputeDustEmission(geometry.grids[m].grid(i,j,k).absorbed_energy,
				CurGrainModel, 
				geometry.grids[m].grid(i,j,k).emitted_energy,
				DoStochastic); 

	    double total_emit_energy = 0.0;
	    total_emit_energy = NumUtils::integrate<double>(tmp_wave,geometry.grids[m].grid(i,j,k).emitted_energy[0])*geometry.grids[m].grid(i,j,k).num_H;
	    
	    global_total_emitted += total_emit_energy;

	    if (fabs(1.0 - (total_emit_energy/tot_abs_energy)) > 1e-3) {
	      cout << "energy conservations worse than 1e-3" << endl;
	      cout << "total_abs/H atom = " << tot_abs_energy << endl;
	      cout << "total_emit_energy/H atom = " << total_emit_energy << endl;
	      cout << "ratio emit/abs = " << total_emit_energy/tot_abs_energy << endl;
// 	      exit(8);
	    }

	    // add to the emitted energy sums
	    for (z = 0; z < n_emit_components; z++)
	      for (x = 0; x < runinfo.wavelength.size(); x++) {
		// temp stuff
// 		geometry.grids[m].grid(i,j,k).emitted_energy[z][x] = 1.0;

		// perm stuff
		// need to multiply the emitted energy passed back by dust_thermal_emission
		// by the number of HI atoms in the cell to get the total emitted energy
		if (!finite(geometry.grids[m].grid(i,j,k).emitted_energy[z][x])) {
		  cout << z << " " << x << endl;
		  cout << geometry.grids[m].grid(i,j,k).emitted_energy[z][x] << endl;
		  cout << " emitted energy is not finite (before num_H calc) " << endl;
		  exit(8);
		}
		geometry.grids[m].grid(i,j,k).emitted_energy[z][x] *= geometry.grids[m].grid(i,j,k).num_H;
		// add up the emitted energy to the total emitted (per wavelength & component)
		if (!finite(geometry.grids[m].grid(i,j,k).emitted_energy[z][x])) {
		  cout << z << " " << x << endl;
		  cout << geometry.grids[m].grid(i,j,k).emitted_energy[z][x] << endl;
		  cout << " emitted energy is not finite (after num_H calc) " << endl;
		  exit(8);
		}
		runinfo.emitted_lum[z][x] += geometry.grids[m].grid(i,j,k).emitted_energy[z][x];
	      }
	  } else {

	    num_cells_not_enough++;
	    if (tot_abs_energy == 0) num_cells_zero++;

// 	    cout << "not enough energy (min) = " << tot_abs_energy << " (" << min_enough_energy << ")" << endl;
	  }
	}
      }
    
  }

  runinfo.total_absorbed_energy = global_total_absorbed;

  runinfo.num_cells_enough = num_cells_enough;
  runinfo.num_cells_not_enough = num_cells_not_enough;
  runinfo.num_cells_zero = num_cells_zero;
  runinfo.num_cells_too_few_waves = num_cells_too_few_waves;
  if (fabs(1.0 - (global_total_emitted/global_total_absorbed)) > runinfo.energy_conserve_target) {
    cout << "energy conservation will not be meet..." << endl;
    cout << "need to change the min_enough_energy target in the code" << endl;
    cout << "or" << endl;
    cout << "(more likely) need to add more photons to get better radiation fields." << endl;

    cout << "global energy conservations check" << endl;
    cout << "total_abs = " << global_total_absorbed << endl;
    cout << "total_emit_energy = " << global_total_emitted << endl;
    cout << "ratio emit/abs = " << global_total_emitted/global_total_absorbed << endl;
    
    cout << "num cells = " << geometry.num_cells << endl;
    cout << "num cells w/ enough energy = " << num_cells_enough << endl;
    cout << "num cells w/ not enough energy (includes zero) = " << num_cells_not_enough << endl;
    cout << "num cells w/ zero energy = " << num_cells_zero << endl;
    cout << "num cells w/ too few wavelengths = " << num_cells_too_few_waves << endl;

  } else if (runinfo.verbose >= 1) {
    cout << "global energy conservations check" << endl;
    cout << "total_abs = " << global_total_absorbed << endl;
    cout << "total_emit_energy = " << global_total_emitted << endl;
    cout << "ratio emit/abs = " << global_total_emitted/global_total_absorbed << endl;
    
    cout << "num cells w/ enough energy = " << num_cells_enough << endl;
    cout << "num cells w/ not enough energy (includes zero) = " << num_cells_not_enough << endl;
    cout << "num cells w/ zero energy = " << num_cells_zero << endl;
    cout << "num cells w/ too few wavelengths = " << num_cells_too_few_waves << endl;
  }

  geometry.emitted_energy_grid_initialized = 1;
}
