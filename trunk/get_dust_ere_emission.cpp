// ======================================================================
//   Procedure to get the dust ERE emission spectra at each point in the
// model.  Bascially, loops over the grid(s) and uses a simple model of
// the ERE carrier to compute the ERE emission.
//
// 2008 Sep/KDG - created from get_dust_thermal_emission
// 2009 Jan/KDG - finished (fixed the ERE scaling problem - do in absorbed
//                          absorbed energy, not J units)
// ======================================================================
#include "get_dust_ere_emission.h"
//#define DEBUG_GDEE

void get_dust_ere_emission (geometry_struct& geometry,
			    runinfo_struct& runinfo)

{
  int i,j,k,m,z = 0;
  uint x = 0;

  int n_emit_components = 1;

  double global_total_emitted = 0.;
  double global_total_absorbed = 0.;
  double total_emit_energy = 0.0;

#ifdef DEBUG_GDEE
  cout << "start ere initialize" << endl;
#endif

  // initialize the total emitted sum (by component and wavelength)
  if (!runinfo.emitted_ere_energy_grid_initialized) {
    runinfo.emitted_ere_lum.resize(n_emit_components);
    for (z = 0; z < n_emit_components; z++)
      runinfo.emitted_ere_lum[z].resize(runinfo.wavelength.size(),0.0);
  } else {
    for (z = 0; z < n_emit_components; z++)
      for (x = 0; x < runinfo.wavelength.size(); x++)
	runinfo.emitted_ere_lum[z][x] = 0.0;
  }


#ifdef DEBUG_GDEE
  cout << "end ere initialize" << endl;
#endif

  long num_cells_enough = 0;
  long num_cells_not_enough = 0;
  long num_cells_zero = 0;
  long num_cells_too_few_waves = 0;

  // define here to save doing it in the loop
  double tot_abs_energy = 0.0;
  double max_abs_energy = 0.0;
  int tot_nonzero = 0;

  vector<double> tmp_wave_for_interpol;
  vector<double> tmp_j_for_interpol;
  vector<double> tmp_new_j_from_interpol;

  double total_ere_photons = 0.0;
  double total_ere_energy = 0.0;
  double ere_sigma_ratio = 0.0;

  // calculate the total absorbed energy from the previously saved total abosrbed energy by wavelength
  vector<double> tmp_wave;
  tmp_wave.resize(runinfo.wavelength.size());
  vector<double> tmp_abs_energy;
  tmp_abs_energy.resize(runinfo.wavelength.size());

  // loop over all the defined grids
  for (m = 0; m < int(geometry.grids.size()); m++) {

    // loop of the cells in this grid
    for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
      for (j = 0; j < geometry.grids[m].index_dim[1]; j++) {

	if (runinfo.verbose >= 1)
	  cout << "working on dust ere emission (m,k,j) = " << m << " " << k << " " << j << endl;

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
	
	  if (geometry.grids[m].grid(i,j,k).dust_tau_per_pc > 0.0) {
	    // determine if there is any energy absorbed needing to be remitted
	    tot_abs_energy = 0.0;
	    max_abs_energy = 0.0;
	    tot_nonzero = 0;
	    
	    max_abs_energy = 0.0;
	    for (x = 0; x < runinfo.wavelength.size(); x++) {
	      tmp_wave[x] = double(runinfo.wavelength[x]);
	      tmp_abs_energy[x] = geometry.grids[m].grid(i,j,k).absorbed_energy[x]*4.*Constant::PI*runinfo.ave_C_abs[x];
	      if (tmp_abs_energy[x] > max_abs_energy) max_abs_energy = tmp_abs_energy[x];
	      if (geometry.grids[m].grid(i,j,k).absorbed_energy[x] > 0.) tot_nonzero++;
	    }
	    if (tot_nonzero >= int(0.5*runinfo.wavelength.size()))
	      tot_abs_energy = NumUtils::integrate<double>(tmp_wave,tmp_abs_energy)*geometry.grids[m].grid(i,j,k).num_H;
	    else
	      tot_abs_energy = 0.0;
	    
#ifdef DEBUG_GDEE
	    cout << "total nonzero = " << tot_nonzero << endl;
	    cout << "total absorbed energy = " << tot_abs_energy << endl;
	    cout << "max absorbed energy = " << max_abs_energy << endl;
	    cout << "num_H = " << geometry.grids[m].grid(i,j,k).num_H << endl;
#endif
	    
	    if ((tot_abs_energy > 0.0) && (tot_nonzero < int(0.5*runinfo.wavelength.size()))) {
	      
	      num_cells_too_few_waves++;
	      
	    } else if (tot_abs_energy == 0.0) {
	      
	      num_cells_zero++;
	      
	    } else {
	      
	      num_cells_enough++;
	      
	      // interoplate the radiative field to fill all the wavelength points
	      // if it has any nonzero points
	      if (tot_nonzero != int(runinfo.wavelength.size())) {
		for (x = 0; x < runinfo.wavelength.size(); x++) {
		  if (geometry.grids[m].grid(i,j,k).absorbed_energy[x] > 0) {
		    tmp_wave_for_interpol.push_back(runinfo.wavelength[x]);
		    tmp_j_for_interpol.push_back(geometry.grids[m].grid(i,j,k).absorbed_energy[x]);
		  }
		}
		
		tmp_new_j_from_interpol = NumUtils::interpol(tmp_j_for_interpol,tmp_wave_for_interpol,tmp_wave);
		tmp_wave_for_interpol.resize(0);
		tmp_j_for_interpol.resize(0);
		
		for (x = 0; x < runinfo.wavelength.size(); x++) {
// 		  cout << tmp_wave[x] << " ";
// 		  cout << runinfo.wavelength[x] << " ";
// 		  cout << tmp_new_j_from_interpol[x] << " ";
// 		  cout << geometry.grids[m].grid(i,j,k).absorbed_energy[x] << " ";
// 		  cout << endl;
		  geometry.grids[m].grid(i,j,k).absorbed_energy[x] = tmp_new_j_from_interpol[x];
		}
		
		tmp_new_j_from_interpol.resize(0);
		
		// 	      exit(1);
	      }
		
	      // 	    cout << "total min = " << tot_abs_energy << " " << min_enough_energy << endl;
	      // 	  if ((tot_abs_energy > 0.) && (tot_nonzero > int(0.75*runinfo.wavelength.size())) && (tot_abs_energy > 1e-35)) {
#ifdef DEBUG_GDEE
	      // output the J
	      for (x = 0; x < runinfo.wavelength.size(); x++) {
		cout << geometry.grids[m].grid(i,j,k).absorbed_energy[x] << " ";
	      }
	      cout << endl;
#endif
	      
	      // get the ERE emission for this cell
	      // simple model currently
	      //   photon conversion efficiency, max wavelength for excitation
	      //   central wavelength, FWHM
	      
	      // first integrate below the max wavelength in photon units
	      //   resuse the j vectors for the integrations
	      x = 0;
	      while ((runinfo.wavelength[x] <= runinfo.ere_excitation_wavelength) && 
		     (x < runinfo.wavelength.size())) {
		tmp_wave_for_interpol.push_back(runinfo.wavelength[x]);
		tmp_j_for_interpol.push_back(geometry.grids[m].grid(i,j,k).absorbed_energy[x]*
					     Constant::IPLANCKLIGHT*runinfo.wavelength[x]);
		// remove the energy just absorbed from J (it is now in ERE) (only ere peak photon energy)
 		geometry.grids[m].grid(i,j,k).absorbed_energy[x] -= tmp_j_for_interpol[x]*runinfo.ere_efficiency*
 		  Constant::PLANCKLIGHT/runinfo.ere_peak_wavelength;
		// convert from rad field density to the absorbed energy density
		tmp_j_for_interpol[x] *= 4.*Constant::PI*runinfo.ave_C_abs[x];
		x++;
	      }
	      total_ere_photons = NumUtils::integrate<double>(tmp_wave_for_interpol,tmp_j_for_interpol);
	      
	      // convert to the number of photons to emit in ERE
	      total_ere_photons *= runinfo.ere_efficiency;
	      
#ifdef DEBUG_GDEE
	      cout << "total_ere_photons/H atom = " << total_ere_photons << endl;
#endif
	      
	      // convert to ERE energy
	      total_ere_photons *= Constant::PLANCKLIGHT/runinfo.ere_peak_wavelength;
	      
#ifdef DEBUG_GDEE	
	      cout << "total_ere_energy/H atom = " << total_ere_photons << endl;
#endif

	      // convert to peak ERE energy
	      total_ere_photons /= (sqrt(2.0)*(runinfo.ere_fwhm/2.35)*pow(M_PI,0.5));

	      // now create the Gaussian ERE emission
	      //   need to check the conversion between fwhm and sigma - on a plane and no internet)
 	      // adjust to account for non-zero width of emission
	      //   this is a factor that works for a specific set of values
	      //   need to work out the correct correction factor (takes into account width)
	      //   KDG - 9 Jan 2009
	      //   this is the 1.42 factor
	      for (x = 0; x < runinfo.wavelength.size(); x++) {
		ere_sigma_ratio = fabs(runinfo.wavelength[x] - runinfo.ere_peak_wavelength)/(runinfo.ere_fwhm/2.35);
		// 	      cout << ere_sigma_ratio << " ";
		if (ere_sigma_ratio < 4.0) {
		  geometry.grids[m].grid(i,j,k).emitted_energy[0][x] = (total_ere_photons/1.42)*exp(-(0.5*ere_sigma_ratio*ere_sigma_ratio));
		  //  		cout << runinfo.wavelength[x] << " " << geometry.grids[m].grid(i,j,k).emitted_energy[0][x];
		} else
		  geometry.grids[m].grid(i,j,k).emitted_energy[0][x] = 0.0;
		// 	      cout << endl;
	      }
	      
	      tmp_wave_for_interpol.resize(0);
	      tmp_j_for_interpol.resize(0);
	      
	      total_emit_energy = 0.0;
	      total_emit_energy = NumUtils::integrate<double>(tmp_wave,geometry.grids[m].grid(i,j,k).emitted_energy[0])*geometry.grids[m].grid(i,j,k).num_H;
	      
	      global_total_emitted += total_emit_energy;
	      
	      total_ere_energy = (total_ere_photons*(runinfo.ere_fwhm/2.35)*pow(M_PI,0.5)*geometry.grids[m].grid(i,j,k).num_H);
	      if (fabs(1.0 - (total_emit_energy/total_ere_energy) > 1e-1)) {
		cout << "energy conservations worse than 1e-1" << endl;
		cout << "total_abs/H atom = " << tot_abs_energy << endl;
		cout << "total_emit_energy/H atom = " << total_emit_energy << endl;
		cout << "ratio emit/abs = " << total_emit_energy/tot_abs_energy << endl;
		cout << "ratio emit/ere = " << total_emit_energy/(total_ere_photons*(runinfo.ere_fwhm/2.35)*pow(M_PI,0.5)*geometry.grids[m].grid(i,j,k).num_H) << endl;
		cout << "ere/H atom = " << (total_ere_photons*(runinfo.ere_fwhm/2.35)*pow(M_PI,0.5))*geometry.grids[m].grid(i,j,k).num_H << endl;
		// 	      exit(8);
	      }
	      
	      // add to the emitted energy sums
	      for (z = 0; z < n_emit_components; z++)
		for (x = 0; x < runinfo.wavelength.size(); x++) {
		  // need to multiply the emitted energy passed back by dust_thermal_emission
		  // by the number of HI atoms in the cell to get the total emitted energy
		  if (!finite(geometry.grids[m].grid(i,j,k).emitted_energy[z][x])) {
		    cout << z << " " << x << endl;
		    cout << geometry.grids[m].grid(i,j,k).emitted_energy[z][x] << endl;
		    cout << " emitted ere energy is not finite (before num_H calc) " << endl;
		    exit(8);
		  }
		  geometry.grids[m].grid(i,j,k).emitted_energy[z][x] *= geometry.grids[m].grid(i,j,k).num_H;
		  // add up the emitted energy to the total emitted (per wavelength & component)
		  if (!finite(geometry.grids[m].grid(i,j,k).emitted_energy[z][x])) {
		    cout << z << " " << x << endl;
		    cout << geometry.grids[m].grid(i,j,k).emitted_energy[z][x] << endl;
		    cout << " emitted ere energy is not finite (after num_H calc) " << endl;
		    exit(8);
		  }
		  runinfo.emitted_ere_lum[z][x] += geometry.grids[m].grid(i,j,k).emitted_energy[z][x];
		}
	    }
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
    cout << "energy conservation will not be meet (ERE)..." << endl;
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
    cout << "global energy conservations check (ERE)" << endl;
    cout << "total_abs = " << global_total_absorbed << endl;
    cout << "total_emit_energy = " << global_total_emitted << endl;
    cout << "ratio emit/abs = " << global_total_emitted/global_total_absorbed << endl;
    
    cout << "num cells w/ enough energy = " << num_cells_enough << endl;
    cout << "num cells w/ not enough energy (includes zero) = " << num_cells_not_enough << endl;
    cout << "num cells w/ zero energy = " << num_cells_zero << endl;
    cout << "num cells w/ too few wavelengths = " << num_cells_too_few_waves << endl;
  }

  runinfo.emitted_ere_energy_grid_initialized = 0;
}
