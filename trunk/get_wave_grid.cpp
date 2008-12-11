// ======================================================================
//   Procedure to define the wavelength grid given inputs in the configfile.
//
// 2007 Apr/KDG - written
// ======================================================================
#include "get_wave_grid.h"
//#define DEBUG_GWG

void get_wave_grid (ConfigFile& param_data,
		    runinfo_struct& runinfo)
{
  // get info on grid type
  string grid_type = param_data.SValue("Dust Grains","wave_type");
  check_input_param("wave_type",grid_type,"res");
  if (grid_type == "res") {
    float wave_min = param_data.FValue("Dust Grains","wave_min");
    check_input_param("dust_model: wave_min",wave_min,0.05,1e4);
    float wave_max = param_data.FValue("Dust Grains","wave_max");
    check_input_param("dust_model: wave_min",wave_min,0.05,1e4);
    float wave_res = param_data.FValue("Dust Grains","wave_resolution");
    check_input_param("dust_model: wave_resolution",wave_res,0.05,1e4);
 
    // determine the number of wavelength points needed
    runinfo.n_waves = int(log10(wave_max/wave_min)/
			  log10((1.0 + 2.0*wave_res)/(2.0*wave_res - 1.0)) + 1.0);
#ifdef DEBUG_GWG    
    cout << "n_waves = " << runinfo.n_waves << endl;
#endif

    // construct the wavelength grid
    float log_wave_min = log10(wave_min);
    float delta_log_wave = (log10(wave_max) - log_wave_min)/(runinfo.n_waves - 1);
    float cur_wave = 0.0;
    int i = 0;
    for (i = 0; i < runinfo.n_waves; i++) {
      cur_wave = log_wave_min + delta_log_wave*float(i);
      cur_wave = pow(double(10.0),double(cur_wave));
      // convert wavelengths from micron to cm
      cur_wave *= Constant::UM_CM;
      runinfo.wavelength.push_back(cur_wave);
#ifdef DEBUG_GWG    
      cout << "wave (" << i+1 << ") = " << cur_wave << endl;
#endif
    }
  } else if (grid_type == "file") {
    // get the filename
    string wave_filename = param_data.SValue("Dust Grains","wave_file");
    // check that the file exists
    ifstream wave_file(wave_filename.c_str());
    if (wave_file.fail()) {
      cout << "Empirical wave file (" << wave_filename << ") does not exist." << endl;
      exit(8);
    }
    wave_file.close();

    // get the wavelength grid
    vector<double> wavelength;
    DataFile(wave_filename, wavelength);

    // go through the values and make sure they are within bounds
    int i;
    double check_wave = 0.005;
    runinfo.n_waves = wavelength.size();
    for (i = 0; i < runinfo.n_waves; i++) {
      check_input_param("wavelength from file",wavelength[i],check_wave,1e4);
      check_wave = wavelength[i];
      runinfo.wavelength.push_back(wavelength[i]*Constant::UM_CM);
    }
    
  } else {
    cout << "wavelength grid wave_type = " << grid_type << " not known" << endl;
    cout << "allowed types = (res)" << endl;
    exit(8);
  }
}
