// ======================================================================
//   Procedure to get the dust thermal emission spectra at each point in the
// model.  Bascially, loops over the grid(s) and calls Karl M.'s routine
// to get the dust emission from each grid cell.
//
// 2007 Oct/KDG - written
// 2008 Mar/KDG - updated calculation to do the units of the emitted energy correctly
// 2008 Sept 2 - added explicit cast to NumUtils::integrate
// 2009 Jun 3-4 - Modified handling of ComputeDustEmission returns to be failure 
//                mode aware. Code will populate Failure instance with information 
//                if dust heating code fails for any bin.
//                If more failure modes (not related to dust emission) need to be added,
//                it should be straitforward...
// ======================================================================
#include "get_dust_thermal_emission.h"
//#define DEBUG_GDTE

void get_dust_thermal_emission (geometry_struct& geometry,
				runinfo_struct& runinfo, 
				GrainModel& CurGrainModel, 
				DirtyFailure * Failure)

{
  float _FailureSz; 
  int _FailureComp;

  int i,j,k,m,z = 0;
  uint x = 0;

  int status; 

  bool DoStochastic = false;
  if (runinfo.do_stochastic_dust_emission) DoStochastic = true;
  // For now, turn off stochastic heating when using effective grain heating.  
  if (runinfo.effective_grain_heating) DoStochastic = false; 
  // define the cutoffs for passing to the dust emission code
//   float cutoff_frac = 0.05;
//   float min_energy_frac = 0.001;
  int good_enough_photons = 10;
  float cutoff_frac = 0.5;
  float min_energy_frac = 0.1;
  if (DoStochastic) {
    cutoff_frac = 0.5;
    min_energy_frac = 0.1;
  }

  double global_total_emitted = 0.;
  double global_total_absorbed = 0.;
  double total_emit_energy = 0.0;
  // get the number of emission components
  int n_emit_components = runinfo.n_emission_grain_types;
  //if (runinfo.do_emission_grain) n_emit_components += 2*CurGrainModel.getNComp();
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

  cout << "global_total_absorbed = " << global_total_absorbed << endl;
  if (global_total_absorbed == 0.0) exit(8);

  double min_enough_energy = min_energy_frac*runinfo.energy_conserve_target*global_total_absorbed/geometry.num_cells;
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

  vector<double> tmp_wave_for_interpol;
  vector<double> tmp_j_for_interpol;
  vector<double> tmp_new_j_from_interpol;

  int cur_plane_good = 0;

  // Required energy conservation per bin - input parameter? 
  float econs_tolerance = 5.0e-3; 
  //float econs_tolerance_throwaway = 0.5; 
  float econs_tolerance_throwaway = 10.; 

  // loop over all the defined grids
  for (m = 0; m < int(geometry.grids.size()); m++) {

    // loop of the cells in this grid
    for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
      for (j = 0; j < geometry.grids[m].index_dim[1]; j++) {

	if (runinfo.verbose >= 1) {
	  cout << "working on dust emission grid (m,k,j) = " << m << " " << k << " " << j;
	  cout.flush();
	}

	cur_plane_good = 0;
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
	    if (geometry.grids[m].grid(i,j,k).absorbed_energy_num_photons[x] >= good_enough_photons) {
	      tot_nonzero++;
	      tmp_abs_energy[x] = geometry.grids[m].grid(i,j,k).absorbed_energy[x]*4.*Constant::PI*runinfo.ave_C_abs[x];
	      if (tmp_abs_energy[x] > max_abs_energy) max_abs_energy = tmp_abs_energy[x];
	    } else tmp_abs_energy[x] = 0.0;
	  }
 	  tot_abs_energy = NumUtils::integrate<double>(tmp_wave,tmp_abs_energy)*geometry.grids[m].grid(i,j,k).num_H;
		
#ifdef DEBUG_GDTE
	  cout << "total nonzero = " << tot_nonzero << endl;
	  cout << "total absorbed energy = " << tot_abs_energy << endl;
	  cout << "max absorbed energy = " << max_abs_energy << endl;
#endif
	  if ((tot_abs_energy > 0.0) && (tot_nonzero < int(cutoff_frac*runinfo.wavelength.size()))) {

	    num_cells_too_few_waves++;

	  } else if (tot_abs_energy >= min_enough_energy) {

	    num_cells_enough++;
	    
	    // interoplate the radiative field to fill all the wavelength points
	    // if it has any nonzero points
	    if ((tot_nonzero != int(runinfo.wavelength.size())) && (DoStochastic)) {
	      // 	    if (tot_nonzero != int(runinfo.wavelength.size())) {
	      for (x = 0; x < runinfo.wavelength.size(); x++) {
#ifdef DEBUG_GDTE
 		cout << geometry.grids[m].grid(i,j,k).absorbed_energy[x] << " ";
#endif
// 		if (geometry.grids[m].grid(i,j,k).absorbed_energy[x] > 0) {
		if (geometry.grids[m].grid(i,j,k).absorbed_energy_num_photons[x] >= good_enough_photons) {
		  tmp_wave_for_interpol.push_back(double(runinfo.wavelength[x]));
		  tmp_j_for_interpol.push_back(double(geometry.grids[m].grid(i,j,k).absorbed_energy[x]));
		}
	      }
#ifdef DEBUG_GDTE
	      cout << endl;
#endif
	      
#ifdef DEBUG_GDTE
	      cout << "entering j interpol..." << endl;
#endif
	      tmp_new_j_from_interpol = NumUtils::interpol(tmp_j_for_interpol,tmp_wave_for_interpol,tmp_wave);
#ifdef DEBUG_GDTE
	      cout << "...leaving j interpol" << endl;
#endif
	      tmp_wave_for_interpol.resize(0);
	      tmp_j_for_interpol.resize(0);

	      for (x = 0; x < runinfo.wavelength.size(); x++) {
//  		cout << tmp_wave[x] << " ";
//  		cout << runinfo.wavelength[x] << " ";
// 		cout << tmp_new_j_from_interpol[x] << " ";
// 		cout << geometry.grids[m].grid(i,j,k).absorbed_energy[x] << " ";
// 		cout << endl;
		geometry.grids[m].grid(i,j,k).absorbed_energy[x] = tmp_new_j_from_interpol[x];
	      }

	      tmp_new_j_from_interpol.resize(0);

// 	      exit(1);
	    }

// 	    cout << "total min = " << tot_abs_energy << " " << min_enough_energy << endl;
// 	  if ((tot_abs_energy > 0.) && (tot_nonzero > int(0.75*runinfo.wavelength.size())) && (tot_abs_energy > 1e-35)) {
#ifdef DEBUG_GDTE
	    //output the J
	    float rad_unc = 0.0;
	    for (x = 0; x < runinfo.wavelength.size(); x++) {
	      cout << geometry.grids[m].grid(i,j,k).absorbed_energy_num_photons[x] << " ";
	      cout << geometry.grids[m].grid(i,j,k).absorbed_energy[x] << " ";
	      cout << geometry.grids[m].grid(i,j,k).absorbed_energy_x2[x] << " ";
	      rad_unc = (geometry.grids[m].grid(i,j,k).absorbed_energy_x2[x]/geometry.grids[m].grid(i,j,k).absorbed_energy_num_photons[x]) -
		pow(double(geometry.grids[m].grid(i,j,k).absorbed_energy[x]/geometry.grids[m].grid(i,j,k).absorbed_energy_num_photons[x]),double(2.0));
	      cout << rad_unc << " ";
	      if (rad_unc > 0.0)
		rad_unc = sqrt(rad_unc);///geometry.grids[m].grid(i,j,k).absorbed_energy_num_photons[x];
	      else
		rad_unc = 0.0;
	      cout << rad_unc << " ";
// 	      rad_unc *= geometry.grids[m].grid(i,j,k).absorbed_energy[x];
// 	      cout << rad_unc << " ";
	      // the S/N reduces to the below equation (modulo a factor of (N-1)/N)
	      cout << geometry.grids[m].grid(i,j,k).absorbed_energy[x]/rad_unc << " ";
	      cout << endl;
	    }
	    cout << endl;
#endif
	    // get the dust emission spectrum given the input wavlength vector and radiation field vector
	    // emitted energy returned is in units of ergs s^-1 HI atom^-1
#ifdef DEBUG_GDTE
	    cout << "entering ComputeDustEmission..." << endl;
#endif
	    cur_plane_good++;
	    status = ComputeDustEmission(geometry.grids[m].grid(i,j,k).absorbed_energy,
					 CurGrainModel, 
					 geometry.grids[m].grid(i,j,k).emitted_energy,
					 DoStochastic,runinfo.effective_grain_heating,
					 _FailureSz,_FailureComp); 
#ifdef DEBUG_GDTE
	    cout << "...leaving ComputeDustEmission" << endl;
#endif
	    if (status != Flags::FSUCCESS) { //  Returned a failure from ComputeDustEmission
	      Failure->AddFailure(status); 
	      Failure->AddCellBook(m,i,j,k); 
	      Failure->AddGrainInfo(CurGrainModel.getModelName(),_FailureSz,_FailureComp); 
	      Failure->AddEnergyInfo(tot_abs_energy,geometry.grids[m].grid(i,j,k).absorbed_energy); 
	    }  else {
	      total_emit_energy = 0.0;
	      total_emit_energy = NumUtils::integrate<double>(tmp_wave,geometry.grids[m].grid(i,j,k).emitted_energy[0])*geometry.grids[m].grid(i,j,k).num_H;

	      float energy_diff = fabs(1.0 - (total_emit_energy/tot_abs_energy));
	      if (energy_diff > econs_tolerance) {
		// Add a failure status for this?  - maybe
		// right now - testing by not using this cell
		cout << "energy conservations worse than " << econs_tolerance << endl;
		cout << "total_abs/H atom = " << tot_abs_energy << endl;
		cout << "total_emit_energy/H atom = " << total_emit_energy << endl;
		cout << "ratio emit/abs = " << total_emit_energy/tot_abs_energy << " ";
		cout << total_emit_energy << endl;
		// 	      exit(8)
		// zeroing out the emitted energy - test
		if (energy_diff > econs_tolerance_throwaway) {
		  cout << "not using the energy emitted in this cell" << endl;
		  for (z = 0; z < n_emit_components; z++)
		    for (x = 0; x < runinfo.wavelength.size(); x++) {
		      geometry.grids[m].grid(i,j,k).emitted_energy[z][x] = 0.0;
		    }
		} else {
		  global_total_emitted += total_emit_energy;
		}
              } else {
                global_total_emitted += total_emit_energy;
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
		    // Add failure status for this? 
		    cout << endl;
		    cout << i << " ";
		    cout << j << " ";
		    cout << k << endl;
		    uint xx = 0;
		    for (xx = 0; xx < runinfo.wavelength.size(); xx++) {
		      cout << xx << " " << geometry.grids[m].grid(i,j,k).absorbed_energy[xx] << " ";
		      cout << geometry.grids[m].grid(i,j,k).absorbed_energy_num_photons[xx] << " ";
		      cout << geometry.grids[m].grid(i,j,k).save_radiation_field_density[xx] << " ";
		      cout << geometry.grids[m].grid(i,j,k).emitted_energy[0][xx] << " ";
		      cout << geometry.grids[m].grid(i,j,k).emitted_energy[1][xx] << " ";
		      cout << geometry.grids[m].grid(i,j,k).emitted_energy[2][xx] << " ";
		      cout << geometry.grids[m].grid(i,j,k).emitted_energy[3][xx] << " ";
		      cout << endl;
		    }
		    cout << geometry.grids[m].grid(i,j,k).emitted_energy.size() << endl;
		    cout << geometry.grids[m].grid(i,j,k).emitted_energy[z].size() << endl;
		    cout << z << " " << x << endl;
		    cout << geometry.grids[m].grid(i,j,k).emitted_energy[z][x] << endl;
		    cout << " emitted energy is not finite (before num_H calc) " << endl;
		    exit(8);
		  }
		  geometry.grids[m].grid(i,j,k).emitted_energy[z][x] *= geometry.grids[m].grid(i,j,k).num_H;
		  // add up the emitted energy to the total emitted (per wavelength & component)
		  if (!finite(geometry.grids[m].grid(i,j,k).emitted_energy[z][x])) {
		    // Add failure status for this? 
		    cout << z << " " << x << endl;
		    cout << geometry.grids[m].grid(i,j,k).emitted_energy[z][x] << endl;
		    cout << " emitted energy is not finite (after num_H calc) " << endl;
		    exit(8);
		  }
		  runinfo.emitted_lum[z][x] += geometry.grids[m].grid(i,j,k).emitted_energy[z][x];
		}
	    }
	  } else {
	    
	    num_cells_not_enough++;
	    if (tot_abs_energy == 0) num_cells_zero++;

// 	    cout << "not enough energy (min) = " << tot_abs_energy << " (" << min_enough_energy << ")" << endl;
	  }
	}
	cout << "; " << cur_plane_good << endl;
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
