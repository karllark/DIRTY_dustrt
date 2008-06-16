// ======================================================================
//   Procedure to get the basic run parameters from the ConfigFile.  If
// not present, print error and exit.
//
// 2006 Jun/KDG - written
// ======================================================================
#include "get_run_parameters.h"

void get_run_parameters (ConfigFile& param_data,
			 output_struct& output,
			 geometry_struct& geometry,
			 runinfo_struct& runinfo)
{
  geometry.n_photons = long(param_data.DValue("Run","num_photons"));
  check_input_param("num_photons",geometry.n_photons,1,1000000000);

  // get type of energy grid storing (memory/disk) desired
  geometry.abs_energy_storage_type = param_data.IValue("Run","abs_energy_storage");
  check_input_param("abs_energy_storage_type",geometry.abs_energy_storage_type,0,1);

  // indicate that the absorbed energy grid has not been initialized
  geometry.abs_energy_grid_initialized = 0;

  geometry.n_photons = long(param_data.DValue("Run","num_photons"));
  check_input_param("num_photons",geometry.n_photons,1,1000000000);

  int image_size = param_data.IValue("Run","output_image_size");
  check_input_param("output_image_size",image_size,1,50000);
  output.image_size[0] = image_size;
  output.image_size[1] = image_size;

  output.file_base = param_data.SValue("Run","output_filebase");
  check_input_param("output_filebase",output.file_base,"dirty_test");
  
  // setup the emission_type string (added to file_base later) to be "" for stellar
  output.emission_type = "";

  output.type = param_data.SValue("Run","output_type");
  check_input_param("output_filebase",output.type,"ratio");

  runinfo.verbose = param_data.IValue("Run","verbose");
  if (runinfo.verbose == -99) runinfo.verbose = 0;  // set to no output if not initially set
  check_input_param("verbose",runinfo.verbose,0,2);

  runinfo.do_dust_emission = param_data.IValue("Run","do_dust_emission");
  if (runinfo.do_dust_emission == -99) runinfo.do_dust_emission = 0;  // set to no if not initially set
  check_input_param("do_dust_emission",runinfo.do_dust_emission,0,1);

  runinfo.do_emission_grain = param_data.IValue("Run","do_emission_grain");
  if (runinfo.do_emission_grain == -99) runinfo.do_emission_grain = 0;  // set to no if not initially set
  check_input_param("do_emission_grain",runinfo.do_emission_grain,0,1);

  if (runinfo.do_dust_emission) {
    // now see what energy conservation is required
    runinfo.energy_conserve_target = param_data.IValue("Run","energy_conserve_target");
    check_input_param("energy_conserve_target",runinfo.energy_conserve_target,0.,1.);
  }

  // output info

  runinfo.do_image_output = param_data.IValue("Run","do_image_output");
  if (runinfo.do_image_output == -99) runinfo.do_image_output = 0;  // set to no if not initially set
  check_input_param("do_image_output",runinfo.do_image_output,0,1);

  runinfo.do_global_output = param_data.IValue("Run","do_global_output");
  if (runinfo.do_global_output == -99) runinfo.do_global_output = 0;  // set to no if not initially set
  check_input_param("do_global_output",runinfo.do_global_output,0,1);

//   runinfo.global_output_fits = param_data.IValue("Run","global_output_fits");
//   if (runinfo.global_output_fits == -99) runinfo.global_output_fits = 0;  // set to no if not initially set
//   check_input_param("global_output_fits",runinfo.global_output_fits,0,1);

  // check that at least one type of output is picked
  if (!runinfo.do_image_output && !runinfo.do_global_output) {
    cout << "Neither images nor global,multiwavelength output picked = no output!" << endl;
    cout << "One or the other is required." << endl;
    cout << "Add 'do_image_output=1' or 'do_global_output=1' to the 'Run' section of the parameter file." << endl;
    exit(8);
  }

  output.arrays_allocated = 0;
  output.num_outputs = geometry.num_observers;
  runinfo.dust_thermal_emission = 1;

}
