// ======================================================================
//   Procedure to get the sed paramters from the Config file and
// read in the necessary secondary files.
//
// 2007 Sep/KDG - written
// 2008 Mar/KDG - udpated to handle global sed output setup
// 2008 Nov/KDG - updated to bin then interpolate to conserve all the input energy
// ======================================================================
#include "get_sed_parameters.h"
//#define DEBUG_GSP

void get_sed_parameters (ConfigFile& param_data,
			 runinfo_struct& runinfo,
			 GrainModel& CurGrainModel)

{
  using namespace NumUtils;

  runinfo.sed_type = param_data.SValue("SED","type");
  if (runinfo.sed_type == "NULL") {
    // no SED, just fill sed_lum with ones
        
    unsigned int i = 0;
    for (i = 0; i < runinfo.wavelength.size(); i++) 
      runinfo.sed_lum.push_back(1.);

  } else {
    // get the filename
    string sed_filename = param_data.SValue("SED","sed_file");
    // check that the file exists
    ifstream sed_file(sed_filename.c_str());
    if (sed_file.fail()) {
      cout << "SED file (" << sed_filename << ") does not exist." << endl;
      exit(8);
    }
    sed_file.close();
    
    unsigned int i = 0;
    vector<double> wavelength, luminosity;

    if (runinfo.sed_type == "ssp_file") {

      // get the sed parameters
      //   assuming units are in ergs/s/Hz
      vector<double> stellar, nebular, total;
      DataFile(sed_filename, wavelength, luminosity, stellar, nebular, total);

      // get the sfr(constant) or mass(burst)
      float sfr_or_mass = param_data.FValue("SED","sfr_or_mass");
      check_input_param("SED sfr or mass",sfr_or_mass,0.,1e20);

      // convert input SED from log(a) to a and scale by sfr/mass and change from Hz^-1 to cm^-1
      // convert the input wavelength to microns
      //   also check that the input SED values are within allowed ranges
      for (i = 0; i < wavelength.size(); i++) {
	// wavelengths input as Angstroms
	check_input_param("SED wavelength",wavelength[i],0,1e7);
	// luminosity input as log(a)
	check_input_param("SED luminosity",luminosity[i],0,40);
	
	wavelength[i] *= Constant::ANG_CM;
	luminosity[i] = sfr_or_mass*pow(10.,luminosity[i])*(Constant::LIGHT)/pow(wavelength[i],2.0);
      }

    } else if (runinfo.sed_type == "bb_file") {
      // get the sed parameters
      //   assuming units are in ergs/s/Hz
      //   and wavelengths are in microns
      DataFile(sed_filename, wavelength, luminosity);

      // change from Hz^-1 to cm^-1
      //   also check that the input SED values are within allowed ranges
      for (i = 0; i < wavelength.size(); i++) {
	// wavelengths input as microns
	check_input_param("SED wavelength",wavelength[i],0,1e9);
	// luminosity input
	check_input_param("SED wavelength",luminosity[i],0,1e40);
	
	wavelength[i] *= Constant::UM_CM;
	luminosity[i] *= (Constant::LIGHT)/pow(wavelength[i],2.0);

#ifdef DEBUG_GSP
	cout << i << " ";
	cout << wavelength[i] << " ";
	cout << luminosity[i] << endl;
#endif
      }
    } else {
      cout << "Type of SED [" << runinfo.sed_type << "] not recognized." << endl;
      exit(8);
    }

    // bin the input SED onto the model wavelength grid
    // to ensure all the SED energy is used in the model
    vector<float>::iterator closeIter;
    vector<double> temp_sed;
    vector<int> temp_sed_npts;
    temp_sed.resize(runinfo.wavelength.size(),0.0);
    temp_sed_npts.resize(runinfo.wavelength.size(),0);
    for (i = 0; i < luminosity.size(); i++) {
      if ((wavelength[i] >= runinfo.wavelength[0]) && (wavelength[i] <= runinfo.wavelength[runinfo.wavelength.size()-1])) {
	closeIter = find_if(runinfo.wavelength.begin(),runinfo.wavelength.end(),bind2nd(greater_equal<float>(),static_cast<float>(wavelength[i])));
	uint closeIndex = distance(runinfo.wavelength.begin(),closeIter);
#ifdef DEBUG_GSP
	cout << i << " ";
	cout << wavelength[i] << " ";
	cout << runinfo.wavelength[closeIndex] << " "; 
	cout << (wavelength[i] - runinfo.wavelength[closeIndex-1]) << " ";
	cout << (runinfo.wavelength[closeIndex] - wavelength[i]) << " ";
	cout << endl;
#endif
	if (closeIndex > 0) {
	  if ((runinfo.wavelength[closeIndex] - wavelength[i]) > (wavelength[i] - runinfo.wavelength[closeIndex-1])) closeIndex--;
	}
#ifdef DEBUG_GSP
	cout << i << " ";
	cout << wavelength[i] << " ";
	cout << runinfo.wavelength[closeIndex] << " "; 
// 	cout << (wavelength[i] - runinfo.wavelength[closeIndex-1]) << " ";
// 	cout << (runinfo.wavelength[closeIndex] - wavelength[i]) << " ";
	cout << endl;
#endif
	temp_sed[closeIndex] += luminosity[i];
	temp_sed_npts[closeIndex]++;
      }
    }

    // setup for the interpolation to get ride of zeros
    // get the runinfo.wavelength in double form for interpol function
    // put the luminosity in log space to make the interpolation better for rapidly changing seds
    vector<double> run_wave;
    vector<double> temp_sed_nozero;
    vector<double> temp_sed_nozero_wave;
    for (i = 0; i < runinfo.wavelength.size(); i++) {
      run_wave.push_back(log10l(double(runinfo.wavelength[i])));
//      run_wave.push_back(double(runinfo.wavelength[i]));
      if (temp_sed_npts[i] > 0) {
	temp_sed_nozero.push_back(log10l(double(temp_sed[i]/temp_sed_npts[i])));
//	temp_sed_nozero.push_back(temp_sed[i]/temp_sed_npts[i]);
	temp_sed_nozero_wave.push_back(log10l(runinfo.wavelength[i]));
// 	temp_sed_nozero_wave.push_back(runinfo.wavelength[i]);
#ifdef DEBUG_GSP
	cout << i << " ";
	cout << log10l(double(temp_sed[i]/temp_sed_npts[i])) << " ";
	cout << runinfo.wavelength[i] << " ";
	cout << temp_sed_npts[i] << " ";
	cout << endl;
#endif
      }
    }

    // extrapolate to the longest wavelength
    if (runinfo.wavelength.back() > temp_sed_nozero_wave.back()) {
      vector<double>::iterator y_iter = temp_sed_nozero.end() - 1;
      vector<double>::iterator x_iter = temp_sed_nozero_wave.end() - 1;
      //cout << *y_iter << "\t" << *(y_iter - 1) << endl;
      //cout << *x_iter << "\t" << *(x_iter - 1) << endl;
      //double x = runinfo.wavelength.back();
      //double y = (*y_iter - *(y_iter - 1))/(log10l(*x_iter) - log10l(*(x_iter - 1)))*(log10l(x) - log10l(*x_iter)) + *y_iter;
      double x = log10l(runinfo.wavelength.back());
      double y = (*y_iter - *(y_iter - 1))/(*x_iter - *(x_iter - 1))*(x - *x_iter) + *y_iter;
      //cout << x << "\t" << y << endl;
      
      temp_sed_nozero.push_back(y);
      temp_sed_nozero_wave.push_back(x);
    }

    // interpolate SED onto wavelength grid
    // numbers are power law coefficents for extrapolation on the left & right sides
    runinfo.sed_lum = interpol(temp_sed_nozero, temp_sed_nozero_wave, run_wave,0,0);
    //    runinfo.sed_lum = interpol(luminosity, wavelength, run_wave,2,-2);
    
    // now un-log10 the luminosity 
    for (i = 0; i < runinfo.wavelength.size(); i++) {
//      cout << i << " ";
//      cout << runinfo.wavelength[i] << " ";
//      cout << runinfo.sed_lum[i] << " ";
      runinfo.sed_lum[i] = pow(10.0,runinfo.sed_lum[i]);
//      cout << runinfo.sed_lum[i] << endl;
    }

    //    exit(8);
  }    

  if (runinfo.do_global_output) {
    runinfo.out_sed_lum_offset = 0;  // doing stellar direct/scattered storage

    // allocate the vector to store the transmitted stellar/dust luminosity
    int n_out_sed = 2;
    if (runinfo.do_ere_emission)
      n_out_sed += 1;
    if (runinfo.do_dust_emission) 
      n_out_sed += 2*CurGrainModel.getNComp();
    n_out_sed *= 2;  // multiply by 2 to allow for direct and scattered
    
    runinfo.out_sed_lum.resize(n_out_sed);
    runinfo.out_sed_lum_unc.resize(n_out_sed);
    int k = 0;
    for (k = 0; k < n_out_sed; k++) {
      runinfo.out_sed_lum[k].resize(runinfo.wavelength.size(),0.0);
      runinfo.out_sed_lum_unc[k].resize(runinfo.wavelength.size(),0.0);
    }
  }

}
