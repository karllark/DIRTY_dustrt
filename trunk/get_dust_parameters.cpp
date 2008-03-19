// ======================================================================
//   Procedure to get the dust parameters from the ConfigFile.  If
// not present, print error and exit.
//
// 2006 Nov/KDG - written
// 2007 Apr/KDG - added dustgrains (from KM)
// ======================================================================
#include "get_dust_parameters.h"
//#define DEBUG_GDP

void get_dust_parameters (ConfigFile& param_data,
			  GrainModel& CurGrainModel,
			  runinfo_struct& runinfo)
{
  // TBD: interpolate the dust scattering parameters to the SED grid
  runinfo.empir_dust = 0;
  runinfo.model_dust = 0;
  
  string dust_type = param_data.SValue("Dust Grains","type");
  if (strcmp(dust_type.c_str(),"single_wavelength") == 0) {
    runinfo.empir_dust = 1;
    runinfo.n_waves = 1;

    float albedo = param_data.FValue("Dust Grains","albedo");
    check_input_param("albedo",albedo,0.,1.);
    runinfo.albedo.push_back(albedo);

    float g = param_data.FValue("Dust Grains","g");
    check_input_param("g",g,-1.,1.);
    runinfo.g.push_back(g);

    float wavelength = param_data.FValue("Dust Grains","wavelength");
    check_input_param("wavelength",wavelength,0.001,1e5);
    runinfo.wavelength.push_back(wavelength);

    // set the tau_to_tau_ref value to 1.
    runinfo.tau_to_tau_ref.push_back(1.);

    // set the grain info we won't have or really use
    runinfo.tau_to_h.push_back(1.0);
    runinfo.ave_C_abs.push_back(1.);

  } else if (strcmp(dust_type.c_str(),"multi_wavelength") == 0) {
    runinfo.empir_dust = 1;

    // get the filename
    string empir_dust_filename = param_data.SValue("Dust Grains","file");
    // check that the file exists
    ifstream empir_dust_file(empir_dust_filename.c_str());
    if (empir_dust_file.fail()) {
      cout << "Empirical dust file (" << empir_dust_filename << ") does not exist." << endl;
      exit(8);
    }
    empir_dust_file.close();

    // get the dust scattering parameters
    vector<double> wavelength, tau, albedo, g;
    DataFile(empir_dust_filename, wavelength, tau, albedo, g);

    // go through the values and make sure they are within bounds
    int i;
    runinfo.n_waves = wavelength.size();
    for (i = 0; i < runinfo.n_waves; i++) {
      // wavelengths assumed to be in microns (10 A to 10 cm)
      check_input_param("empirical dust wavelength",wavelength[i],0.001,1e5);
      check_input_param("empirical dust tau",tau[i],0.0,1e6);
      check_input_param("empirical dust albedo",albedo[i],0.0,1.0);
      check_input_param("empirical dust g",g[i],-1.0,1.0);

      runinfo.wavelength.push_back(wavelength[i]);
      runinfo.tau_to_tau_ref.push_back(tau[i]);
      runinfo.albedo.push_back(albedo[i]);
      runinfo.g.push_back(g[i]);

      // set the grain info we won't have or really use
      runinfo.tau_to_h.push_back(1.0);
      runinfo.ave_C_abs.push_back(1.);
    }

  } else if (strcmp(dust_type.c_str(),"dust_model") == 0) {
    runinfo.model_dust = 1;
    get_wave_grid(param_data,runinfo);
    // get dust grain info
    CurGrainModel.MakeGrainModel(param_data,runinfo.wavelength);
    // now get the tau, albedo, and g values
    runinfo.albedo = CurGrainModel.getAlbedo();
    runinfo.g = CurGrainModel.getphFuncEff();
    runinfo.tau_to_h = CurGrainModel.getTau(); // getTau returns tau/H I atom
    runinfo.tau_to_tau_ref = runinfo.tau_to_h;
    runinfo.ave_C_abs = CurGrainModel.getCAbsEffNorm();  // getCAbsEffNorm returns cm^2/H I atom
    float norm_tau = CurGrainModel.getTau(0.55*(Constant::UM_CM));
    int i;
#ifdef DEBUG_GDP
      cout << "wavelength" << " ";
      cout << "tau/tau_norm" << " ";
      cout << "albedo" << " ";
      cout << "g" << endl;
#endif
    for (i = 0; i < runinfo.n_waves; i++) {
      runinfo.tau_to_tau_ref[i] /= norm_tau;
#ifdef DEBUG_GDP
      cout << runinfo.wavelength[i]*Constant::CM_UM << " ";
      cout << runinfo.tau_to_tau_ref[i] << " ";
      cout << runinfo.albedo[i] << " ";
      cout << runinfo.g[i] << endl;
#endif
    }
#ifdef DEBUG_GDP
    cout << "norm tau = " << norm_tau << endl;
#endif

  } else {
    cout << "Type of dust grains [" << dust_type << "] not recognized." << endl;
    exit(8);
  }

}
