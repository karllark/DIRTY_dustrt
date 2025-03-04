// ======================================================================
//   Procedure to get the sed paramters from the Config file and
// read in the necessary secondary files.
//
// 2007 Sep/KDG - written
// 2008 Mar/KDG - udpated to handle global sed output setup
// 2008 Nov/KDG - updated to bin then interpolate to conserve all the input
// energy
// ======================================================================
#include "get_sed_parameters.h"
// #define DEBUG_GSP

void
get_sed_parameters (ConfigFile &param_data, runinfo_struct &runinfo, GrainModel &CurGrainModel)

{
#ifdef DEBUG_GSP
  cout << "entering get_sed_parameters.." << endl;
#endif

  using namespace NumUtils;

  runinfo.sed_type = param_data.SValue ("SED", "type");
  if (runinfo.sed_type == "NULL")
    {
      // no SED, just fill sed_lum with ones

      unsigned int i = 0;
      for (i = 0; i < runinfo.wavelength.size (); i++)
        // reasonable for a star, needed to get radiation field correct for
        // single wavelength case
        runinfo.sed_lum.push_back (1.e45);
    }
  else
    {
#ifdef DEBUG_GSP
      cout << "getting the sed_filename and reading..." << endl;
#endif
      // get the filename
      string sed_filename = param_data.SValue ("SED", "sed_file");
      // check that the file exists
      ifstream sed_file (sed_filename.c_str ());
      if (sed_file.fail ())
        {
          cout << "SED file (" << sed_filename << ") does not exist." << endl;
          exit (8);
        }
      sed_file.close ();

      unsigned int i = 0;
      vector<double> wavelength, luminosity;

      if (runinfo.sed_type == "ssp_file")
        {
          // get the sed parameters
          //   assuming units are in ergs/s/Hz
          vector<double> stellar, nebular, total;
          DataFile (sed_filename, wavelength, luminosity, stellar, nebular, total);

          // get the sfr(constant) or mass(burst)
          float sfr_or_mass = param_data.FValue ("SED", "sfr_or_mass");
          check_input_param ("SED sfr or mass", sfr_or_mass, 0., 1e20);

          // convert input SED from log(a) to a and scale by sfr/mass and
          // change from Hz^-1 to cm^-1 convert the input wavelength to microns
          //   also check that the input SED values are within allowed ranges
          for (i = 0; i < wavelength.size (); i++)
            {
              // wavelengths input as Angstroms
              check_input_param ("SED wavelength", wavelength[i], 0, 1e7);
              // luminosity input as log(a)
              // check_input_param("SED luminosity",luminosity[i],0,40);
              // to use total luminosity
              check_input_param ("SED luminosity", total[i], -15, 40);
              luminosity[i] = total[i];

              wavelength[i] *= Constant::ANG_CM;
              luminosity[i] = sfr_or_mass * pow (10., luminosity[i]) * (Constant::LIGHT)
                              / pow (wavelength[i], 2.0);
            }

#ifdef DEBUG_GSP
          cout << "done reading ssp_file." << endl;
#endif
        }
      else if (runinfo.sed_type == "bb_file")
        {
          // get the sed parameters
          //   assuming units are in ergs/s/Hz
          //   and wavelengths are in microns
          DataFile (sed_filename, wavelength, luminosity);

          // change from Hz^-1 to cm^-1
          //   also check that the input SED values are within allowed ranges
          for (i = 0; i < wavelength.size (); i++)
            {
              // wavelengths input as microns
              check_input_param ("SED wavelength", wavelength[i], 0, 1e9);
              // luminosity input
              check_input_param ("SED wavelength", luminosity[i], 0, 1e40);

              wavelength[i] *= Constant::UM_CM;
              luminosity[i] *= (Constant::LIGHT) / pow (wavelength[i], 2.0);

#ifdef DEBUG_GSP
              cout << i << " ";
              cout << wavelength[i] << " ";
              cout << luminosity[i] << endl;
#endif
            }
        }
      else
        {
          cout << "Type of SED [" << runinfo.sed_type << "] not recognized." << endl;
          exit (8);
        }

      // determine if the mapping the input SED to the requested wavelength
      // grid will be done by binning or interpolation
      int mapping_bin;
      mapping_bin = param_data.IValue ("SED", "sed_bin");
      if (mapping_bin == -99)
        mapping_bin = 0;
      check_input_param ("SED mapping_bin", mapping_bin, 0, 1);

      // needed for both mapping methods
      vector<double> run_wave;
      for (i = 0; i < runinfo.wavelength.size (); i++)
        {
          run_wave.push_back (log10l (double (runinfo.wavelength[i])));
        }

      // now bin or interpolate
      if (mapping_bin == 1)
        {
          cout << "Binning in put SED" << endl;
#ifdef DEBUG_GSP
          cout << "starting SED binning..." << endl;
#endif
          // bin the input SED onto the model wavelength grid
          // to ensure all the SED energy is used in the model
          vector<float>::iterator closeIter;
          vector<double> temp_sed;
          vector<int> temp_sed_npts;
          temp_sed.resize (runinfo.wavelength.size (), 0.0);
          temp_sed_npts.resize (runinfo.wavelength.size (), 0);
#ifdef DEBUG_GSP
          cout << "starting SED binning loop..." << endl;
#endif
          for (i = 0; i < luminosity.size (); i++)
            {
#ifdef DEBUG_GSP
              cout << "i = " << i;
              cout << " " << runinfo.wavelength.size () << endl;
              cout << wavelength[i] << " ";
              cout << runinfo.wavelength[0] << " ";
              cout << runinfo.wavelength[runinfo.wavelength.size () - 1] << endl;
#endif
              if ((wavelength[i] >= runinfo.wavelength[0])
                  && (wavelength[i] <= runinfo.wavelength[runinfo.wavelength.size () - 1]))
                {
                  closeIter = find_if (
                      runinfo.wavelength.begin (), runinfo.wavelength.end (),
                      bind2nd (greater_equal<float> (), static_cast<float> (wavelength[i])));
#ifdef DEBUG_GSP
                  cout << "closeIter done" << endl;
#endif
                  uint closeIndex = distance (runinfo.wavelength.begin (), closeIter);
#ifdef DEBUG_GSP
                  cout << i << " ";
                  cout << wavelength[i] << " ";
                  cout << runinfo.wavelength[closeIndex] << " ";
                  if (closeIndex > 0)
                    {
                      cout << (wavelength[i] - runinfo.wavelength[closeIndex - 1]) << " ";
                    }
                  else
                    {
                      cout << wavelength[i] << " ";
                    }
                  cout << (runinfo.wavelength[closeIndex] - wavelength[i]) << " ";
                  cout << endl;
#endif
                  if (closeIndex > 0)
                    {
                      if ((runinfo.wavelength[closeIndex] - wavelength[i])
                          > (wavelength[i] - runinfo.wavelength[closeIndex - 1]))
                        closeIndex--;
                    }
#ifdef DEBUG_GSP
                  cout << i << " ";
                  cout << wavelength[i] << " ";
                  cout << runinfo.wavelength[closeIndex] << " ";
                  // 	cout << (wavelength[i] -
                  // runinfo.wavelength[closeIndex-1]) << "
                  // "; 	cout << (runinfo.wavelength[closeIndex] -
                  // wavelength[i])
                  // << " ";
                  cout << endl;
#endif
                  temp_sed[closeIndex] += luminosity[i];
                  temp_sed_npts[closeIndex]++;
                }
            }

          // setup for the interpolation to get ride of zeros
          // get the runinfo.wavelength in double form for interpol function
          // put the luminosity in log space to make the interpolation better
          // for rapidly changing seds
          vector<double> temp_sed_nozero;
          vector<double> temp_sed_nozero_wave;
          for (i = 0; i < runinfo.wavelength.size (); i++)
            {
              if (temp_sed_npts[i] > 0)
                {
                  temp_sed_nozero.push_back (log10l (double (temp_sed[i] / temp_sed_npts[i])));
                  temp_sed_nozero_wave.push_back (log10l (runinfo.wavelength[i]));
#ifdef DEBUG_GSP
                  cout << i << " ";
                  cout << log10l (double (temp_sed[i] / temp_sed_npts[i])) << " ";
                  cout << runinfo.wavelength[i] << " ";
                  cout << temp_sed_npts[i] << " ";
                  cout << endl;
#endif
                }
            }

          // extrapolate to the longest wavelength
          if (runinfo.wavelength.back () > temp_sed_nozero_wave.back ())
            {
              vector<double>::iterator y_iter = temp_sed_nozero.end () - 1;
              vector<double>::iterator x_iter = temp_sed_nozero_wave.end () - 1;
              // cout << *y_iter << "\t" << *(y_iter - 1) << endl;
              // cout << *x_iter << "\t" << *(x_iter - 1) << endl;
              // double x = runinfo.wavelength.back();
              // double y = (*y_iter - *(y_iter - 1))/(log10l(*x_iter) -
              // log10l(*(x_iter - 1)))*(log10l(x) - log10l(*x_iter)) +
              // *y_iter;
              double x = log10l (runinfo.wavelength.back ());
              double y
                  = (*y_iter - *(y_iter - 1)) / (*x_iter - *(x_iter - 1)) * (x - *x_iter) + *y_iter;
              // cout << x << "\t" << y << endl;

              temp_sed_nozero.push_back (y);
              temp_sed_nozero_wave.push_back (x);
            }

          // interpolate SED onto wavelength grid
          // numbers are power law coefficents for extrapolation on the left &
          // right sides
          runinfo.sed_lum = interpol (temp_sed_nozero, temp_sed_nozero_wave, run_wave, 0, 0);
          //    runinfo.sed_lum = interpol(luminosity, wavelength,
          //    run_wave,2,-2);
        }
      else
        {
          // code to interpolate instead of rebin

          vector<double> iwave_log;
          vector<double> ised_log;
          ised_log.resize (wavelength.size (), 0.0);
          iwave_log.resize (wavelength.size (), 0);
          for (i = 0; i < luminosity.size (); i++)
            {
              iwave_log[i] = log10l (wavelength[i]);
              ised_log[i] = log10l (luminosity[i]);
            }

          runinfo.sed_lum = interpol (ised_log, iwave_log, run_wave, 0, 0);
        }

      // now un-log10 the luminosity
      for (i = 0; i < runinfo.wavelength.size (); i++)
        {
          runinfo.sed_lum[i] = pow (10.0, runinfo.sed_lum[i]);
        }
    }

  if (runinfo.do_global_output)
    {
      runinfo.out_sed_lum_offset = 0; // doing stellar direct/scattered storage

      // allocate the vector to store the transmitted stellar/dust luminosity
      int n_out_sed = 1 + runinfo.n_emission_grain_types;
      //     int n_out_sed = 2;
      //     if (runinfo.do_dust_emission)
      //       n_out_sed += 2*CurGrainModel.getNComp();
      if (runinfo.do_ere_emission)
        n_out_sed += 1;
      n_out_sed *= 2; // multiply by 2 to allow for direct and scattered

      runinfo.out_sed_lum.resize (n_out_sed);
      runinfo.out_sed_lum_unc.resize (n_out_sed);
      int k = 0;
      for (k = 0; k < n_out_sed; k++)
        {
          runinfo.out_sed_lum[k].resize (runinfo.wavelength.size (), 0.0);
          runinfo.out_sed_lum_unc[k].resize (runinfo.wavelength.size (), 0.0);
        }
    }
}
