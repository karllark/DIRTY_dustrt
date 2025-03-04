// ======================================================================
//   Procedure to get the dust parameters from the ConfigFile.  If
// not present, print error and exit.
//
// 2006 Nov/KDG - written
// 2007 Apr/KDG - added dustgrains (from KM)
// 2009 May/KDG - added model scattering phase function
// ======================================================================
#include "get_dust_parameters.h"
// #define DEBUG_GDP

void
get_dust_parameters (ConfigFile &param_data, GrainModel &CurGrainModel, geometry_struct &geometry,
                     runinfo_struct &runinfo)

{
  // TBD: interpolate the dust scattering parameters to the SED grid
  runinfo.empir_dust = 0;
  runinfo.model_dust = 0;

  string dust_type = param_data.SValue ("Dust Grains", "type");
  if (strcmp (dust_type.c_str (), "single_wavelength") == 0)
    {
      runinfo.empir_dust = 1;
      runinfo.n_waves = 1;

      float albedo = param_data.FValue ("Dust Grains", "albedo");
      check_input_param ("albedo", albedo, 0., 1.);
      runinfo.albedo.push_back (albedo);

      float g = param_data.FValue ("Dust Grains", "g");
      check_input_param ("g", g, -1., 1.);
      runinfo.g.push_back (g);

      float wavelength = param_data.FValue ("Dust Grains", "wavelength");
      check_input_param ("wavelength", wavelength, 0.001, 1e5);
      runinfo.wavelength.push_back (wavelength);

      // set the tau_to_tau_ref value to 1.
      runinfo.tau_to_tau_ref.push_back (1.);

      // set the grain info to have have the correction scaling
      // here we have effectively set Cext = 1.0,
      //    hence Cabs = (1-albedo)Cext = (1-albedo)
      runinfo.tau_to_h.push_back (1.0);
      runinfo.ave_C_abs.push_back (1. - albedo);

      // check that if a single wavlength run is picked, then do_global_output
      // is not
      if (runinfo.do_global_output)
        {
          cout << "Cannot do_global_output when only for a single wavelength "
                  "run"
               << endl;
          exit (8);
        }

      // check that if a single wavlength run is picked, then do_dust_emission
      // is not
      if (runinfo.do_dust_emission)
        {
          cout << "Cannot do_dust_emission when only for a single wavelength "
                  "run"
               << endl;
          exit (8);
        }
    }
  else if (strcmp (dust_type.c_str (), "single_wavelength_modelg") == 0)
    {
      runinfo.empir_dust = 1;
      runinfo.n_waves = 1;

      float albedo = param_data.FValue ("Dust Grains", "albedo");
      check_input_param ("albedo", albedo, 0., 1.);
      runinfo.albedo.push_back (albedo);

      // get the filename with the scattering phase function
      string phi_filename = param_data.SValue ("Dust Grains", "phi_file");
      // check that the file exists
      ifstream phi_file (phi_filename.c_str ());
      if (phi_file.fail ())
        {
          cout << "scattering phase function phi file (" << phi_filename << ") does not exist." << endl;
          exit (8);
        }
      phi_file.close ();

      // get the dust scattering parameters
      vector<double> angle, phi;
      DataFile (phi_filename, angle, phi);

      // go through the values and make sure they are within bounds
      uint i;

      double cur_weight = 0.;
      double weight = 0.;

      check_input_param ("model scattering phase function (phi)", phi[0], 0., 1e10);
      angle[0] = cos (angle[0] * Constant::PI / 180.0);
      geometry.phi.push_back (phi[0]);
      geometry.phi_sum.push_back (0.0);
      geometry.phi_angle.push_back (angle[0]);
      for (i = 1; i < angle.size (); i++)
        {
          // phi is between 0 and 1
          check_input_param ("model scattering phase function (phi)", phi[i], 0., 1e10);
          angle[i] = cos (angle[i] * Constant::PI / 180.0);

          cur_weight = 0.5 * (phi[i] + phi[i - 1]) * (angle[i - 1] - angle[i]);
          weight += cur_weight;
          geometry.phi_angle.push_back (angle[i]);
          geometry.phi.push_back (phi[i]);
          geometry.phi_sum.push_back (weight);
        }
      cout << "integral = " << geometry.phi_sum[angle.size () - 1] << endl;

      // now normalized so that it runs from 0 to 1
      double save_sum = geometry.phi_sum[angle.size () - 1];
      for (i = 1; i < angle.size (); i++)
        {
          geometry.phi[i] /= save_sum * (2. * Constant::PI);
          geometry.phi_sum[i] /= save_sum;
        }
      runinfo.g.push_back (-2); // set to -2 to trigger use of model phase function

      float wavelength = param_data.FValue ("Dust Grains", "wavelength");
      check_input_param ("wavelength", wavelength, 0.001, 1e5);
      runinfo.wavelength.push_back (wavelength);

      // set the tau_to_tau_ref value to 1.
      runinfo.tau_to_tau_ref.push_back (1.);

      // set the grain info we won't have or really use
      runinfo.tau_to_h.push_back (1.0);
      runinfo.ave_C_abs.push_back (1.);
    }
  else if (strcmp (dust_type.c_str (), "multi_wavelength") == 0)
    {
      runinfo.empir_dust = 1;

      // get the filename
      string empir_dust_filename = param_data.SValue ("Dust Grains", "file");
      // check that the file exists
      ifstream empir_dust_file (empir_dust_filename.c_str ());
      if (empir_dust_file.fail ())
        {
          cout << "Empirical dust file (" << empir_dust_filename << ") does not exist." << endl;
          exit (8);
        }
      empir_dust_file.close ();

      // get the dust scattering parameters
      vector<double> wavelength, tau, albedo, g;
      DataFile (empir_dust_filename, wavelength, tau, albedo, g);

      // go through the values and make sure they are within bounds
      int i;
      runinfo.n_waves = wavelength.size ();
      for (i = 0; i < runinfo.n_waves; i++)
        {
          // wavelengths assumed to be in microns (10 A to 10 cm)
          check_input_param ("empirical dust wavelength", wavelength[i], 0.001, 1e5);
          check_input_param ("empirical dust tau", tau[i], 0.0, 1e6);
          check_input_param ("empirical dust albedo", albedo[i], 0.0, 1.0);
          check_input_param ("empirical dust g", g[i], -1.0, 1.0);

          runinfo.wavelength.push_back (wavelength[i]);
          runinfo.tau_to_tau_ref.push_back (tau[i]);
          runinfo.albedo.push_back (albedo[i]);
          runinfo.g.push_back (g[i]);

          // set the grain info we won't have or really use
          runinfo.tau_to_h.push_back (1.0);
          runinfo.ave_C_abs.push_back (1.);
        }
    }
  else if (strcmp (dust_type.c_str (), "dust_model") == 0)
    {
      runinfo.model_dust = 1;
      get_wave_grid (param_data, runinfo);
      // save the dust grain filename (full path)
      string grainpath = param_data.SValue ("Model Book Keeping", "Path to Dust Properties");
      string grainsubdir = param_data.SValue ("Model Book Keeping", "Model SubDir");
      runinfo.dust_grain_filename = grainpath + grainsubdir + param_data.SValue ("Model Book Keeping", "Model Name");
      // get dust grain info
      CurGrainModel.MakeGrainModel (param_data, runinfo.wavelength);
      // now get the tau, albedo, and g values
      runinfo.effective_grain_heating = param_data.BValue ("Model Book Keeping", "Effective Grain for Heating");
#ifdef DEBUG_GDP
      cout << "Effective heating is " << runinfo.effective_grain_heating << endl;
#endif
      runinfo.albedo = CurGrainModel.getAlbedo ();
      runinfo.g = CurGrainModel.getphFuncEff ();
      runinfo.tau_to_h = CurGrainModel.getTau (); // getTau returns tau/H I atom
      runinfo.tau_to_tau_ref = runinfo.tau_to_h;
      runinfo.ave_C_abs = CurGrainModel.getCAbsEffNorm (); // getCAbsEffNorm returns cm^2/H I atom
      if (runinfo.effective_grain_heating)
        runinfo.n_emission_grain_types = 3;
      else
        runinfo.n_emission_grain_types = 1 + 2 * CurGrainModel.getNComp ();

#ifdef DEBUG_GDP
      cout << "getting tau_wave next..." << endl;
#endif

      runinfo.norm_tau_wave = param_data.FValue ("Geometry", "tau_wave");
      if (isnan (runinfo.norm_tau_wave))
        runinfo.norm_tau_wave = 0.55;

#ifdef DEBUG_GDP
      cout << "tau_wave = " << runinfo.norm_tau_wave << endl;
#endif

      check_input_param ("wavelength to normalize tau", runinfo.norm_tau_wave * (Constant::UM_CM),
                         *min_element (runinfo.wavelength.begin (), runinfo.wavelength.end ()),
                         *max_element (runinfo.wavelength.begin (), runinfo.wavelength.end ()));

#ifdef DEBUG_GDP
      cout << "just checked in the input parameter" << endl;
#endif

      float norm_tau = CurGrainModel.getTau (runinfo.norm_tau_wave * (Constant::UM_CM));
      //     cout << CurGrainModel.getTau(0.1*(Constant::UM_CM))/norm_tau <<
      //     endl; exit(8); cout << runinfo.norm_tau_wave << " ";
      // cout << norm_tau << endl;

      // norm_tau = 4.307e-22;  // according to Karl M. via phone 23 Jun 2015

      // norm_tau = 4.31968e-22;  // according to Simone

      int i;
#ifdef DEBUG_GDP
      cout << "wavelength" << " ";
      cout << "tau/tau_norm" << " ";
      cout << "albedo" << " ";
      cout << "g" << endl;
#endif
      for (i = 0; i < runinfo.n_waves; i++)
        {
          //   cout << runinfo.wavelength[i]*Constant::CM_UM << " ";
          //   cout << runinfo.tau_to_tau_ref[i] << " ";
          runinfo.tau_to_tau_ref[i] /= norm_tau;
          //   cout << runinfo.tau_to_tau_ref[i] << " ";
          //   cout << runinfo.albedo[i] << " ";
          //   cout << runinfo.g[i] << endl;
#ifdef DEBUG_GDP
          cout << runinfo.wavelength[i] * Constant::CM_UM << " ";
          cout << runinfo.tau_to_tau_ref[i] << " ";
          cout << runinfo.albedo[i] << " ";
          cout << runinfo.g[i] << endl;
#endif
        }
#ifdef DEBUG_GDP
      cout << "norm tau = " << norm_tau << endl;
#endif
    }
  else
    {
      cout << "Type of dust grains [" << dust_type << "] not recognized." << endl;
      exit (8);
    }
}
