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
// 2009 Jun/KAM - create DirtyFailure instance on each iteration and 
//                output iteration specific file containing dust emission 
//                failure info if instructed by the parameter file. 
// ======================================================================
#include "dirty.h"
// #define DEBUG_DIRTY

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

  // Are we outputing a failure log - do it here so I don't have to mess with KDGs structures. 
  // Output failure will only be set if we find a 'yes' entry in the param file. 
  bool OutputFailure=param_data.BValue("Run","Output Failure Log"); 
 
  geometry_struct geometry;  // structure with geometry info (dust grid, etc.)
  output_struct output;  // stucture with the output info (images, etc.)
  photon_data photon;   // structure with the photon info (position, direction, weight, etc.)
  runinfo_struct runinfo; // structure with information about the run

  runinfo.param_filename = param_filename; // save the filename

  random_dirty random_obj(long(-987654321));  // object for random number generator

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
  
  // get the dust grain parameters
  get_dust_parameters(param_data, CurGrainModel, geometry, runinfo);

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

  // do the radiative transfer over all the wavelengths
  radiative_transfer_many_waves(geometry, runinfo, output, photon, random_obj, REG_RT, 0);

#ifdef DEBUG_DIRTY
  cout << "RT done; ";
  cout.flush();
#endif

  // do gas emission (iteration needed as the gas provides 
  //    an additional opacity source)

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
      
    // do the ERE radiative transfer over all the wavelengths
    radiative_transfer_many_waves(geometry, runinfo, ere_output, photon, random_obj, ERE_RT, 1);
  }

  // output model_grid info
  if (output.do_output_model_grid)
    output_model_grid(geometry, output, runinfo);

  // start RT+DE iteration (only if DE flag set)
  int iter_num = 1;
//   int iter_max = 5;
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

      // set this here to allow the DE part to do the initialization (or reinitialize)
      geometry.emitted_energy_grid_initialized = 0;

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

  // hold output failure file name. 
  string fFailureFilename; 

  while (!iter_done) {

    // Create a DirtyFailure Object
    // Will remain empty if we don't do an output.  Do all the allocations regardless 
    // if whether OutputFailure is true or not to avoid compile  warnings. 
    fFailureFilename=output.file_base+"_iter"+StringManip::vtos(iter_num)+"_failure.log"; 
    DirtyFailure * Failure = new DirtyFailure(fFailureFilename,runinfo.n_waves);

#ifdef DEBUG_DIRTY
    cout << "te: gdt start; ";
    cout.flush();
#endif
    // get the dust emission at each point in the model
    get_dust_thermal_emission(geometry, runinfo, CurGrainModel, Failure);
#ifdef DEBUG_DIRTY
    cout << "te: gdt done; ";
    cout.flush();
#endif

    // setup the grid emission
    setup_emitted_grid_for_montecarlo(geometry, runinfo, CurGrainModel);
      
    // do the DE radiative transfer over all the wavelengths
    radiative_transfer_many_waves(geometry, runinfo, de_output, photon, random_obj, DE_RT, iter_num);

    // determine if energy is conserved well enough or if another
    // iteration is needed due to dust self-absorption
    check_de_energy_conservation(runinfo, iter_done);

    // limit the max iterations
    if (iter_num >= runinfo.iter_max) iter_done = 1; else iter_num++;
    
    // If output puting failure log, then do it. 
    if (OutputFailure) Failure->WriteFailureLog();
    // Destroy the Failure Object. 
    delete Failure; 
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

