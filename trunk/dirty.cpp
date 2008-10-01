// ======================================================================
//   Program to compute the radiative transfer and dust emission
// for any type of astrophysical object.
// This code is extensivly based on the original version of DIRTY
// written by Karl Gordon and Karl Misselt.
// This new version was written to add an active/adaptive mesh to 
// the code.  In addition, the code will be much cleaner and better
// documented.  It should be possible to release this version to
// others.
//
// 2003 Apr/KDG - started (active grid setup)
// 2003 Jun/KDG - added photon tracking through active grid
// 2004 Dec/KDG - added photon scattering
// 2006 Apr/KDG - added output results
// 2006 Jun/KDG - added parameter file
// 2006 Nov/KDG - added multiple wavelengths (empirical dust)
// 2007 Apr/KDG - added dustgrains (from KM)
// 2007 Dec/KDG - added arbitary grid emission (for dust emission)
// 2008 Jan/KDG - added global output
// 2008 Mar/KDG - fixed global output to use a FITS ASCII table
// 2008 Aug/KDG - added continous absorption
// 2008 Sep/KDG - added ERE
// ======================================================================
#include "dirty.h"
//#define DEBUG_DIRTY

int main(int argc, char* argv[]) 

{
  // parse the command line
  string param_filename;
  if (argc == 2) {
    param_filename = argv[1];
    // check that the file exists
    ifstream test_file(param_filename.c_str());
    if (test_file.fail()) {
      cout << "Parameter file (" << param_filename << ") does not exist." << endl;
      exit(8);
    }
    test_file.close();
  } else {
    cout << "Usage: dirty parameter_file" << endl;
    exit(8);
  }

  // read parameter file
  ConfigFile param_data(param_filename);

  geometry_struct geometry;  // structure with geometry info (dust grid, etc.)
  output_struct output;  // stucture with the output info (images, etc.)
  photon_data photon;   // structure with the photon info (position, direction, weight, etc.)
  runinfo_struct runinfo; // structure with information about the run

  random_dirty random_obj;  // object for random number generator
  random_obj.random_num(long(-987654321));  // initialize

  GrainModel CurGrainModel;  // object for grain model

  // setup the dust grid with dust density (tau/pc), positions, etc.
  setup_dust_grid(param_data, geometry, photon, random_obj);
#ifdef DEBUG_DIRTY
  cout << "sdg done; ";
  cout.flush();
#endif

  // get the run parameters (basic info)
  get_run_parameters(param_data, output, geometry, runinfo);
#ifdef DEBUG_DIRTY
  cout << "grp done; ";
  cout.flush();
#endif
  
  // output model_grid info
  if (output.do_output_model_grid)
    output_model_grid(geometry, output);

  // get the dust grain parameters
  get_dust_parameters(param_data, CurGrainModel, runinfo);
#ifdef DEBUG_DIRTY
  cout << "gdp done; ";
  cout.flush();
#endif
  
  // read SED parameters
  get_sed_parameters(param_data, runinfo, CurGrainModel);
#ifdef DEBUG_DIRTY
  cout << "gsp done; ";
  cout.flush();
#endif

  // loop over all wavelengths
  int i;
  for (i = 0; i < runinfo.n_waves; i++) {
    if (runinfo.verbose >= 1) {
      cout << "working on wavelength [micron] = " << runinfo.wavelength[i]*Constant::CM_UM << endl;
      cout << "tau = " << runinfo.tau_to_tau_ref[i]*geometry.tau << " ";
      cout << "albedo = " << runinfo.albedo[i] << " ";
      cout << "g = " << runinfo.g[i] << " ";
      cout << endl;
    }

    // setup/(re)initialize absorbed energy grid
    setup_absorbed_energy_grid(geometry, runinfo, i, 0);

//     // check the absorbed energy grid (temp needed as energy not conserved)
//     // KDG - 23 Mar 2008
//     check_absorbed_energy_grid(geometry, runinfo);

    // setup dust grains for this wavelength
    get_dust_scat_parameters(i, runinfo, geometry);

    // do RT part
#ifdef DEBUG_DIRTY
    cout << "rt: rt begin; ";
    cout.flush();
#endif
    radiative_transfer(geometry, runinfo, output, photon, random_obj);
#ifdef DEBUG_DIRTY
    cout << "rt: rt end; ";
    cout.flush();
#endif

    // output RT results
    output_results(output, geometry, runinfo, i);

//     // check the absorbed energy grid (temp needed as energy not conserved)
//     // KDG - 23 Mar 2008
//     cout << "**test2**" << endl;
//     check_absorbed_energy_grid(geometry, runinfo);

    // store the result (either in memory or on disk)
    // remember to zero out the absorbed energy grid
    if ((runinfo.do_dust_emission) || (runinfo.do_ere_emission)) {
#ifdef DEBUG_DIRTY
      cout << endl;
      cout << "wave [cm] = " << geometry.wavelength << endl;
      cout << "te: saeg start; ";
      cout.flush();
#endif
      store_absorbed_energy_grid(geometry, runinfo, output, i, 0);
#ifdef DEBUG_DIRTY
      cout << "te: saeg end; ";
      cout.flush();
#endif
    }

//     // check the absorbed energy grid (temp needed as energy not conserved)
//     // KDG - 23 Mar 2008
// 	check_absorbed_energy_grid(geometry, runinfo);
//     if (i == 1) exit(8);
//  	exit(8);
  }

  // check the absorbed energy grid (temp needed as energy not conserved)
  // KDG - 23 Mar 2008
//   check_absorbed_energy_grid(geometry, runinfo);

  // do gas emission (iteration needed as the gas provides 
  //    an additional opacity source)

#ifdef DEBUG_DIRTY
  cout << "RT done; ";
  cout.flush();
#endif

  // set this here to allow the ERE part to do the initialization
  geometry.emitted_energy_grid_initialized = 0;

  // do ERE part (no iteration needed)
  if (runinfo.do_ere_emission) {
    runinfo.emitted_ere_energy_grid_initialized = 0;
    // setup a new output stucture to handle the ere
    output_struct ere_output;  // stucture with the output info (images, etc.)
#ifdef DEBUG_DIRTY
  cout << "setup ere output ";
  cout.flush();
#endif
    setup_ere_dust_emission_output(ere_output, output);
    runinfo.out_sed_lum_offset += 2;  // increment to save ERE direct/scattered luminosity

    // get ERE emission
#ifdef DEBUG_DIRTY
  cout << "get ere ";
  cout.flush();
#endif
    get_dust_ere_emission(geometry, runinfo);

    // get ready for RT
#ifdef DEBUG_DIRTY
  cout << "setup grid for MC ";
  cout.flush();
#endif
    setup_emitted_grid_for_montecarlo(geometry, runinfo, CurGrainModel);
      
    // do ERE RT: loop over all wavelengths
    for (i = 0; i < runinfo.n_waves; i++) {

      if (runinfo.verbose >= 1)
	cout << "working on ere wavelength [micron] = " << runinfo.wavelength[i]*Constant::CM_UM << endl;

#ifdef DEBUG_DIRTY
      cout << endl;
      cout << "i = " << i << endl;
      cout << "emitted_lum = " << runinfo.emitted_ere_lum[0][i] << endl;
#endif

      if (runinfo.emitted_ere_lum[0][i] > 1e-3*runinfo.out_sed_lum[0][i]*runinfo.sed_lum[i]) {

#ifdef DEBUG_DIRTY
	cout << "doing ERE RT" << endl;

	cout << "ee: saeg start; ";
	cout.flush();
#endif
	// setup/(re)initialize absorbed energy grid
	setup_absorbed_energy_grid(geometry, runinfo, i, 1);
#ifdef DEBUG_DIRTY
	cout << "ee: saeg done; ";
	cout.flush();
#endif
      
#ifdef DEBUG_DIRTY
	cout << "ee: gdsp start; ";
	cout.flush();
#endif
	// setup dust grains for this wavelength
	get_dust_scat_parameters(i, runinfo, geometry);
#ifdef DEBUG_DIRTY
	cout << "ee: gdsp end; ";
	cout.flush();
#endif
	
#ifdef DEBUG_DIRTY
	cout << "ee: rt start; ";
	cout.flush();
#endif
	// do RT part
	radiative_transfer(geometry, runinfo, ere_output, photon, random_obj);
#ifdef DEBUG_DIRTY
	cout << "ee: rt end; ";
	cout.flush();
#endif
	
	// output RT results
	output_results(ere_output, geometry, runinfo, i);
	
	// store the result (either in memory or on disk)
	// remember to zero out the absorbed energy grid
	store_absorbed_energy_grid(geometry, runinfo, ere_output, i, 1);
      }
    }
  }
//   emit_ere_photons();

  // start RT+DE iteration (only if DE flag set)
  //    do DE part
  //    do DE/RT part
  // iterate until converged for max iterations reached

  int iter_num = 1;
  int iter_max = 10;
  int iter_done = 0;

#ifdef DEBUG_DIRTY
  cout << "RT done2; ";
  cout.flush();
#endif

  // setup a new output stucture to handle the different grain types
  output_struct de_output;  // stucture with the output info (images, etc.)

#ifdef DEBUG_DIRTY
  cout << "RT done3; ";
  cout.flush();
#endif

  if (!runinfo.do_dust_emission) {
    iter_done = 1;
  } else {
    if (runinfo.do_emission_grain && (geometry.num_observers > 1)) {
      cout << "not possible to do dust emission output of emission by grain type and" << endl;
      cout << "multiple-lines-of-sight without new code." << endl;
      cout << "Set do_emission_type=0 to get multiple-lines-of-sight and total IR emission." << endl;
      exit(8);
    } else if (runinfo.do_emission_grain) {
#ifdef DEBUG_DIRTY
      cout << "te: stdeo start";
      cout.flush();
#endif
      setup_thermal_dust_emission_output(runinfo, de_output, output, photon);
#ifdef DEBUG_DIRTY
      cout << "te: stdeo done";
      cout.flush();
#endif
    }

    // increment to save dust emission direct/scattered luminosity
    runinfo.out_sed_lum_offset += 2;
    runinfo.dust_thermal_emission = 1;
  }

  while (!iter_done) {

#ifdef DEBUG_DIRTY
    cout << "te: gdt start; ";
    cout.flush();
#endif
    // get the dust emission at each point in the model
    get_dust_thermal_emission(geometry, runinfo, CurGrainModel);
#ifdef DEBUG_DIRTY
    cout << "te: gdt done; ";
    cout.flush();
#endif

    // setup the grid emission
    setup_emitted_grid_for_montecarlo(geometry, runinfo, CurGrainModel);
      
    // loop over all wavelengths
    for (i = 0; i < runinfo.n_waves; i++) {

      if (runinfo.verbose >= 1)
	cout << "working on wavelength [micron] = " << runinfo.wavelength[i]*Constant::CM_UM << endl;

#ifdef DEBUG_DIRTY
      cout << endl;
      cout << "i = " << i << endl;
      cout << "emitted_lum = " << runinfo.emitted_lum[0][i] << endl;
#endif

      if (runinfo.emitted_lum[0][i] > 1e-3*runinfo.out_sed_lum[0][i]*runinfo.sed_lum[i]) {

#ifdef DEBUG_DIRTY
	cout << "doing RT" << endl;

	cout << "te: saeg start; ";
	cout.flush();
#endif
	// setup/(re)initialize absorbed energy grid
	setup_absorbed_energy_grid(geometry, runinfo, i, iter_num);
#ifdef DEBUG_DIRTY
	cout << "te: saeg done; ";
	cout.flush();
#endif
      
#ifdef DEBUG_DIRTY
	cout << "te: gdsp start; ";
	cout.flush();
#endif
	// setup dust grains for this wavelength
	get_dust_scat_parameters(i, runinfo, geometry);
#ifdef DEBUG_DIRTY
	cout << "te: gdsp end; ";
	cout.flush();
#endif
	
#ifdef DEBUG_DIRTY
	cout << "te: rt start; ";
	cout.flush();
#endif
	// do RT part
	radiative_transfer(geometry, runinfo, de_output, photon, random_obj);
#ifdef DEBUG_DIRTY
	cout << "te: rt end; ";
	cout.flush();
#endif
	
	// output RT results
	output_results(de_output, geometry, runinfo, i);
	
	// store the result (either in memory or on disk)
	// remember to zero out the absorbed energy grid
	store_absorbed_energy_grid(geometry, runinfo, de_output, i, iter_num);
      }
    }

    // determine if energy is conserved well enough or if another
    // iteration is needed due to dust self-absorption
    check_de_energy_conservation(runinfo, iter_done);

    // limit the max iterations
    if (iter_num >= iter_max) iter_done = 1; else iter_num++;
  }

  // output global, multiwavelength luminosities
  if (runinfo.do_global_output) { 
    if (runinfo.do_dust_emission)
      output_global_results(runinfo, de_output, geometry);
    else
      output_global_results(runinfo, output, geometry);
  }
  return 0;

}

