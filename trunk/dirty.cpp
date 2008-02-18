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
// ======================================================================
#include "dirty.h"

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
  
  // get the dust grain parameters
  get_dust_parameters(param_data, CurGrainModel, runinfo);
  
  // read SED parameters
  get_sed_parameters(param_data, runinfo);

  // loop over all wavelengths
  int i;
  for (i = 0; i < runinfo.n_waves; i++) {
    // setup/(re)initialize absorbed energy grid
    setup_absorbed_energy_grid(geometry, i, 0);

    // setup dust grains for this wavelength
    get_dust_scat_parameters(i, runinfo, geometry);

    // do RT part
    radiative_transfer(geometry, runinfo, output, photon, random_obj);

    // output RT results
    output_results(output, geometry, runinfo, i);

    // store the result (either in memory or on disk)
    // remember to zero out the absorbed energy grid
    store_absorbed_energy_grid(geometry, runinfo, output, i, 0);
  }

  // do gas emission (iteration needed as the gas provides 
  //    an additional opacity source)

  // set this here to allow the ERE part to do the initialization
  geometry.emitted_energy_grid_initialized = 0;

  // do ERE part (no iteration needed)
  //runinfo.out_sed_lum_offset += 2;  // increment to save ERE direct/scattered luminosity
  //output.emission_type = "_ere"; // setup the emission_type string for dust emission
  //emit_ere_photons();

  // start RT+DE iteration (only if DE flag set)
  //    do DE part
  //    do DE/RT part
  // iterate until converged for max iterations reached

  int iter_num = 1;
  int iter_max = 1;
  int iter_done = 0;
  if (!runinfo.do_dust_emission) {
    iter_done = 1;
  } else {
    runinfo.out_sed_lum_offset += 2;  // increment to save dust emission direct/scattered luminosity
    output.emission_type = "_de"; // setup the emission_type string for dust emission
  }

  while (!iter_done) {
    // get the dust emission at each point in the model
    get_dust_thermal_emission(geometry, runinfo, CurGrainModel);
    
    // setup the grid emission
    setup_emitted_grid_for_montecarlo(geometry, runinfo, CurGrainModel);
      
    // loop over all wavelengths
    for (i = 0; i < runinfo.n_waves; i++) {

      // setup/(re)initialize absorbed energy grid
      setup_absorbed_energy_grid(geometry, i, 1);
      
      // setup dust grains for this wavelength
      get_dust_scat_parameters(i, runinfo, geometry);

      // do RT part
      radiative_transfer(geometry, runinfo, output, photon, random_obj);
      
      // output RT results
      output_results(output, geometry, runinfo, i);
      
      // store the result (either in memory or on disk)
      // remember to zero out the absorbed energy grid
      store_absorbed_energy_grid(geometry, runinfo, output, i, 1);
    }

    iter_num++;
    // limit the max iterations
    if (iter_num >= iter_max) iter_done = 1;
  }

  // output final RT+DE images/totals [TBD]
  // output global, multiwavelength luminosities
  if (runinfo.do_global_output)
    output_global_results(runinfo);

  return 0;

}

