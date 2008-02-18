// ======================================================================
//   Procedure to get the sed paramters from the Config file and
// read in the necessary secondary files.
//
// 2007 Sep/KDG - written
// ======================================================================
#include "get_sed_parameters.h"

void get_sed_parameters (ConfigFile& param_data,
			 runinfo_struct& runinfo)
{

  using namespace NumUtils;

  runinfo.sed_type = param_data.SValue("SED","type");
  if (runinfo.sed_type == "NULL") {
    // no SED, just fill sed_lum with ones
        
    unsigned int i = 0;
    for (i = 0; i < runinfo.wavelength.size(); i++) 
      runinfo.sed_lum.push_back(1.);

  } else if (runinfo.sed_type == "ssp_file") {

    // get the filename
    string sed_filename = param_data.SValue("SED","sed_file");
    // check that the file exists
    ifstream sed_file(sed_filename.c_str());
    if (sed_file.fail()) {
      cout << "SED file (" << sed_filename << ") does not exist." << endl;
      exit(8);
    }
    sed_file.close();
    
    // get the sed parameters
    //   assuming units are in ergs/s/Hz
    vector<double> wavelength, luminosity, stellar, nebular, total;
    DataFile(sed_filename, wavelength, luminosity, stellar, nebular, total);

    // get the sfr(constant) or mass(burst)
    float sfr_or_mass = param_data.FValue("SED","sfr_or_mass");
    check_input_param("SED sfr or mass",sfr_or_mass,0.,1e20);

    // convert input SED from log(a) to a and scale by sfr/mass and change from Hz^-1 to cm^-1
    // convert the input wavelength to microns
    //   also check that the input SED values are within allowed ranges
    unsigned int i = 0;
    for (i = 0; i < wavelength.size(); i++) {
      // wavelengths input as Angstroms
      check_input_param("SED wavelength",wavelength[i],0,1e7);
      // luminosity input as log(a)
      check_input_param("SED wavelength",luminosity[i],0,40);
	
      wavelength[i] *= Constant::ANG_CM;
      luminosity[i] = sfr_or_mass*pow(10.,luminosity[i])*(Constant::LIGHT)*pow(wavelength[i],2.0);
    }


    // get the runinfo.wavelength in double form for interpol function
    vector<double> run_wave;
    for (i = 0; i < runinfo.wavelength.size(); i++) 
      run_wave.push_back(double(runinfo.wavelength[i]));

    // interpolate SED onto wavelength grid
    runinfo.sed_lum = interpol(luminosity, wavelength, run_wave);

  } else {
    cout << "Type of SED [" << runinfo.sed_type << "] not recognized." << endl;
    exit(8);
  }

  if (runinfo.do_global_output) {
    runinfo.out_sed_lum_offset = 0;  // doing stellar direct/scattered storage

    // allocate the vector to store the transmitted stellar/dust luminosity
    int n_out_sed = 2;  // maybe allow this to be 4 if we are doing ERE emission
    runinfo.out_sed_lum.resize(n_out_sed);
    runinfo.out_sed_lum_unc.resize(n_out_sed);
    int k = 0;
    for (k = 0; k < n_out_sed; k++) {
      runinfo.out_sed_lum[k].resize(runinfo.wavelength.size(),0.0);
      runinfo.out_sed_lum_unc[k].resize(runinfo.wavelength.size(),0.0);
    }
  }

}
