// ======================================================================
//   Procedure to output the results of dirty.
//
// 2006 Apr/KDG - written
// 2007 Mar/KDG - added comments about how to compute radiation field
//                in ergs s^-1 cm^-3
// ======================================================================
#include "output_results.h"
//#define DEBUG_OUTR

// helper routine
void output_2d_info(fitsfile *out_ptr,
		    char extname[100],
		    char itype[100])

{

  int status = 0;

  fits_write_key(out_ptr, TSTRING, "EXTNAME", extname, "Name of Extension", &status);
  fits_write_key(out_ptr, TSTRING, "ITYPE", itype, "Image Type", &status);
  check_fits_io(status,"fits_write_comment : output results");

}

// main routine
//
void output_results (output_struct& output,
		     geometry_struct& geometry,
		     runinfo_struct& runinfo,
		     int index)

{
  int i;
  
#ifdef DEBUG_OUTR
  cout << "# outputs = " << output.num_outputs << endl;
#endif

  for (i = 0; i < output.num_outputs; i++) {

    if (runinfo.verbose >= 2) {
      cout << "Output for observer = " << (i+1) << endl;
      cout << "# photons = " << output.outputs[i].total_num_photons << endl;
    }

    // compute total stellar weight uncertainty
    double total_stellar_weight_err = 0.0;  
    // first compute the uncertainty on the average weight of a single stellar photon
    if (output.outputs[i].total_num_photons > 0)
      total_stellar_weight_err = 
	output.outputs[i].total_stellar_weight_x2/output.outputs[i].total_num_photons -
	pow(output.outputs[i].total_stellar_weight/output.outputs[i].total_num_photons,2);
    if (total_stellar_weight_err > 0.0)
      total_stellar_weight_err = sqrt(total_stellar_weight_err/output.outputs[i].total_num_photons);
    else
      total_stellar_weight_err = 0.0;
    // this converts to a fractional uncertainty which is the same as the
    // fractional uncertainty on the total stellar weight
    total_stellar_weight_err /= output.outputs[i].total_stellar_weight/output.outputs[i].total_num_photons;

    if (runinfo.verbose >= 2) {
      cout << "Stellar weight = " << output.outputs[i].total_stellar_weight << endl;
    }

    double stellar_sl = output.outputs[i].total_stellar_weight/output.outputs[i].total_num_photons;
    if (runinfo.verbose >= 2) {
      cout << "Stellar S/L = " << stellar_sl;
      cout << " +/- " << stellar_sl*total_stellar_weight_err;
      cout << " (" << 100.0*total_stellar_weight_err << "%)" << endl;
    }

    // determine conversion from lum to flux 4*pi*d^2 in cm^2
    double lum_to_flux = 0.0;
    lum_to_flux = 1.0/(4.*M_PI*pow(geometry.distance*(Constant::PC_CM),2.0));
//     cout << "lum_to_flux = " << lum_to_flux << endl;

    // determine stellar flux
    double stellar_flux = 0.0;
    stellar_flux = stellar_sl*geometry.total_source_luminosity*lum_to_flux;
    if (runinfo.verbose >= 2) {
      cout << "Stellar flux [ergs cm^-2 s^-1 A^-1] = " << stellar_flux << endl;

      cout << "# scattered photons = " << output.outputs[i].total_num_scattered_photons << endl;
    }
    // save the result
    if (runinfo.do_global_output) {
//       cout << output.num_outputs << " ";
//       cout << runinfo.out_sed_lum.size() << " ";
//       cout << runinfo.out_sed_lum_offset+(2*i) << endl;
      runinfo.out_sed_lum[0+runinfo.out_sed_lum_offset+(2*i)][geometry.wave_index] = stellar_sl;
      runinfo.out_sed_lum_unc[0+runinfo.out_sed_lum_offset+(2*i)][geometry.wave_index] = stellar_sl*total_stellar_weight_err;
    }

    // compute total scattered weight uncertainty
    double total_scattered_weight_err = 0.0;  
    // first compute the uncertainty on the average weight of a single scattered photon
    if (output.outputs[i].total_num_photons > 0)
      total_scattered_weight_err = 
	output.outputs[i].total_scattered_weight_x2/output.outputs[i].total_num_scattered_photons -
	pow(output.outputs[i].total_scattered_weight/output.outputs[i].total_num_scattered_photons,2);
    if (total_scattered_weight_err > 0.0)
      total_scattered_weight_err = sqrt(total_scattered_weight_err/output.outputs[i].total_num_scattered_photons);
    else
      total_scattered_weight_err = 0.0;
    // this converts to a fractional uncertainty which is the same as the
    // fractional uncertainty on the total scattered weight
    total_scattered_weight_err /= output.outputs[i].total_scattered_weight/output.outputs[i].total_num_scattered_photons;

    double scattered_sl = output.outputs[i].total_scattered_weight/output.outputs[i].total_num_photons;
    if (runinfo.verbose >= 2) {
      cout << "Scattered weight = " << output.outputs[i].total_scattered_weight << endl;
      cout << "Scattered S/L = " << scattered_sl;
      cout << " +/- " << scattered_sl*total_scattered_weight_err;
      cout << " (" << 100.0*total_scattered_weight_err << "%)" << endl;
    }

    // save the result
    if (runinfo.do_global_output) {
//       cout << runinfo.out_sed_lum.size() << " ";
//       cout << 1+runinfo.out_sed_lum_offset+(2*i) << endl;
      runinfo.out_sed_lum[1+runinfo.out_sed_lum_offset+(2*i)][geometry.wave_index] = scattered_sl;
      runinfo.out_sed_lum_unc[1+runinfo.out_sed_lum_offset+(2*i)][geometry.wave_index] =  scattered_sl*total_scattered_weight_err;
//       cout << "done." << endl;
    }

//     cout << "weights = ";
//     cout << stellar_sl << " ";
//     cout << scattered_sl << " ";
//     cout << 1.0 - (stellar_sl + scattered_sl) << endl;

    // determine scattered flux
    double scattered_flux = 0.0;
    scattered_flux = scattered_sl*geometry.total_source_luminosity*lum_to_flux;
    if (runinfo.verbose >= 2) {
      cout << "Scattered flux [ergs cm^-2 s^-1 A^-1] = " << scattered_flux << endl;
    }

    if (runinfo.do_image_output) {

      NumUtils::Matrix<float> total_weight_xy;
      NumUtils::Matrix<float> total_weight_xy_x2;
      total_weight_xy.MSize(output.image_size[0],output.image_size[1]);
      total_weight_xy_x2.MSize(output.image_size[0],output.image_size[1]);

      // divide the weight images by the input luminosity
      // and calculate the uncertainty image
      int j,k;
      for (j = 0; j < output.image_size[0]; j++)
	for (k = 0; k < output.image_size[1]; k++) {
	  // calculate the stellar uncertainty image (reuse x2 location for uncertainty)
	  // first compute the uncertainty on the average weight of a single stellar photon
	  output.outputs[i].stellar_weight_xy_x2(j,k) = 
	    output.outputs[i].stellar_weight_xy_x2(j,k)/output.outputs[i].num_stellar_photons_xy(j,k) -
	    pow(output.outputs[i].stellar_weight_xy(j,k)/output.outputs[i].num_stellar_photons_xy(j,k),2);
	  if (output.outputs[i].stellar_weight_xy_x2(j,k) > 0.0)
	    output.outputs[i].stellar_weight_xy_x2(j,k) =
	      sqrt(output.outputs[i].stellar_weight_xy_x2(j,k)/output.outputs[i].num_stellar_photons_xy(j,k));
	  else
	    output.outputs[i].stellar_weight_xy_x2(j,k) = 0.0;

	  // this converts to a fractional uncertainty which is the same as the
	  // fractional uncertainty on the total scattered weight
	  output.outputs[i].stellar_weight_xy_x2(j,k) /=
	    output.outputs[i].stellar_weight_xy(j,k)/output.outputs[i].num_stellar_photons_xy(j,k);
	  
	  // calculate the scattered uncertainty image (reuse x2 location for uncertainty)
	  // first compute the uncertainty on the average weight of a single scattered photon
	  output.outputs[i].scattered_weight_xy_x2(j,k) = 
	    output.outputs[i].scattered_weight_xy_x2(j,k)/output.outputs[i].num_photons_xy(j,k) -
	    pow(output.outputs[i].scattered_weight_xy(j,k)/output.outputs[i].num_photons_xy(j,k),2);
	  if (output.outputs[i].scattered_weight_xy_x2(j,k) > 0.0)
	    output.outputs[i].scattered_weight_xy_x2(j,k) =
	      sqrt(output.outputs[i].scattered_weight_xy_x2(j,k)/output.outputs[i].num_photons_xy(j,k));
	  else
	    output.outputs[i].scattered_weight_xy_x2(j,k) = 0.0;
	  
	  // this converts to a fractional uncertainty which is the same as the
	  // fractional uncertainty on the total scattered weight
	  output.outputs[i].scattered_weight_xy_x2(j,k) /=
	    output.outputs[i].scattered_weight_xy(j,k)/output.outputs[i].num_photons_xy(j,k);

	  // divide by luminosity
	  output.outputs[i].scattered_weight_xy(j,k) /= output.outputs[i].total_num_photons;
	  output.outputs[i].stellar_weight_xy(j,k) /= output.outputs[i].total_num_photons;
	  
	  // divide by area of each pixel in sr to get a surface brightness
	  // needs to be done!
	  
	  // now convert the fractional uncertainty to a absolute uncertainty
	  output.outputs[i].stellar_weight_xy_x2(j,k) *= output.outputs[i].stellar_weight_xy(j,k);
	  output.outputs[i].scattered_weight_xy_x2(j,k) *= output.outputs[i].scattered_weight_xy(j,k);

	  // add the two together to get get the "observed" image
	  total_weight_xy(j,k) = output.outputs[i].stellar_weight_xy(j,k) + output.outputs[i].scattered_weight_xy(j,k);
	  total_weight_xy_x2(j,k) = sqrt(output.outputs[i].stellar_weight_xy_x2(j,k)*output.outputs[i].stellar_weight_xy_x2(j,k) +
					 output.outputs[i].scattered_weight_xy_x2(j,k)*output.outputs[i].scattered_weight_xy_x2(j,k));
	}
   

// #ifdef DEBUG_OUTR
//       int m,n;
//       cout << "Image of scattered weight" << endl;
//       for (m = 0; m < output.image_size[0]; m++) {
// 	for (n = 0; n < output.image_size[1]; n++)
// 	  cout << setw(8) << output.outputs[i].scattered_weight_xy(m,n) << " ";
// 	cout << endl;
//       }
      
//       cout << "Image of stellar weight" << endl;
//       for (m = 0; m < output.image_size[0]; m++) {
// 	for (n = 0; n < output.image_size[1]; n++)
// 	  cout << setw(6) << output.outputs[i].stellar_weight_xy(m,n) << " ";
// 	cout << endl;
//       }
      
//       cout << "Image of number of scattered photons" << endl;
//       for (m = 0; m < output.image_size[0]; m++) {
// 	for (n = 0; n < output.image_size[1]; n++)
// 	  cout << setw(6) << output.outputs[i].num_photons_xy(m,n) << " ";
// 	cout << endl;
//       }
// #endif

      // filename of the current output file
      string filename = "!" + output.file_base;
      // add the extra string (allows for ere/dust emission images)
      filename += output.emission_type;
      if (geometry.num_observers > 1) {
	// convert current integer line-of-sight index to a string
	// should be a more elegant way to do this!
	stringstream ss;
	string los_index;
	ss << (i+1);
	ss >> los_index;
	filename += "_los" + los_index;
      }
      if (runinfo.dust_thermal_emission && (runinfo.n_emission_grain_types > 1)) {
	// convert current integer grain_type index to a string
	// should be a more elegant way to do this!
	stringstream ss;
	string grain_index;
	ss << (i+1);
	ss >> grain_index;
	filename += "_ge" + grain_index;
      }
      if (runinfo.n_waves > 1) {
	// convert current integer wavelength index to a string
	// should be a more elegant way to do this!
	stringstream ss;
	string wave_index;
	ss << (index+1);
	ss >> wave_index;
	filename += "_w" + wave_index;

	// convert the wavelength in microns to a string
	filename += "_" + StringManip::vtos(runinfo.wavelength[index]*Constant::CM_UM) + "um";
      }
      filename += ".fits";
      
#ifdef DEBUG_OUTR
      cout << "filename for output = " << filename << endl;
#endif
      
      // create a FITS file with extensions to fill with the output of the model
      fitsfile *out_ptr;   // pointer to the output fits file
      int status = 0;
      fits_create_file(&out_ptr,filename.c_str(), &status);
      check_fits_io(status, "fits_create_file : output_results");
      
      // the primary header is a blank image (very small)
      fits_create_img(out_ptr, 8, 0, 0, &status);
      check_fits_io(status,"fits_create_img : output_results");
      
      // populate the primary header with the details of the run
      // TBD
      
      // final stuff for primary header
      fits_write_comment(out_ptr, "**---------------------------------**",&status);
      fits_write_comment(out_ptr, "Output of the DIRTY model",&status);
      fits_write_comment(out_ptr, "Karl D. Gordon & Karl A. Misselt", &status);
      fits_write_comment(out_ptr, "version v2.0prealpha (Jun 2007)", &status);
      fits_write_comment(out_ptr, "**---------------------------------**",&status);
      check_fits_io(status,"fits_write_comment : output results");
      
      fits_write_key(out_ptr, TFLOAT, "RADIUS", &geometry.radius, "model radius [pc]", &status);
      fits_write_key(out_ptr, TFLOAT, "DIST", &geometry.distance, "distance to model [pc]", &status);
      
      float tdusta = float(geometry.albedo);
      fits_write_key(out_ptr, TFLOAT, "DUST_A", &tdusta, "dust grain albedo", &status);
      float tdustg = float(geometry.g);
      fits_write_key(out_ptr, TFLOAT, "DUST_G", &tdustg, "dust grain g", &status);
      float cur_rad_tau = geometry.tau*geometry.tau_to_tau_ref;
      fits_write_key(out_ptr, TFLOAT, "DUST_TAU", &cur_rad_tau, "radial dust optical depth at current wavelength", &status); 
      fits_write_key(out_ptr, TFLOAT, "DUST_FF", &geometry.filling_factor, "dust filling factor", &status);
      fits_write_key(out_ptr, TFLOAT, "DUST_DR", &geometry.density_ratio, "dust density ratio (k2/k1)", &status);
      fits_write_key(out_ptr, TFLOAT, "DUST_K1", &geometry.clump_densities[0], "k1 optical depth per pc", &status);
      fits_write_key(out_ptr, TFLOAT, "DUST_K2", &geometry.clump_densities[1], "k2 optical depth per pc", &status);
      
      check_fits_io(status,"fits_write_key : output results (geometry details)");
      
      fits_write_key(out_ptr, TINT, "N_OBS", &geometry.num_observers, "number of observer line-of-sights", &status);
      int cur_obs_num = i + 1;
      fits_write_key(out_ptr, TINT, "CUR_OBS", &cur_obs_num, "current observer line-of-sight", &status);
      float theta_deg = geometry.observer_angles[0][i]*180.0/M_PI;
      fits_write_key(out_ptr, TFLOAT, "LOSTHETA", &theta_deg, "theta of line-of-sight [degrees]", &status);
      float phi_deg = geometry.observer_angles[1][i]*180.0/M_PI;
      fits_write_key(out_ptr, TFLOAT, "LOSPHI", &phi_deg, "phi of line-of-sight [degrees]", &status);
      
      check_fits_io(status,"fits_write_key : output results (run details 0.9)");

      if (!finite(stellar_sl)) stellar_sl = 0.0;
      fits_write_key(out_ptr, TDOUBLE, "STEL_SL", &stellar_sl, "stellar/emitted luminosity", &status);
      double stel_sle = stellar_sl*total_stellar_weight_err;
      if (!finite(stel_sle)) stel_sle = 0.0;
      fits_write_key(out_ptr, TDOUBLE, "STEL_SLE", &stel_sle, "STEL_SL uncertainty", &status);
      
      if (!finite(scattered_sl)) scattered_sl = 0.0;
      fits_write_key(out_ptr, TDOUBLE, "SCAT_SL", &scattered_sl, "scattered/emitted luminosity", &status);
      check_fits_io(status,"fits_write_key : output results (run details: 0.95)");
      double scat_sle = scattered_sl*total_scattered_weight_err;
      if (!finite(scat_sle)) scat_sle = 0.0;
      fits_write_key(out_ptr, TDOUBLE, "SCAT_SLE", &scat_sle, "SCAT_SL uncertainty", &status);
      
      check_fits_io(status,"fits_write_key : output results (run details 1.0)");
      
      // create and output the total intensity image
      fits_create_img(out_ptr, -32, 2, output.image_size, &status);
      check_fits_io(status,"fits_create_image : output_results (total intensity/luminosity)");
      output_2d_info(out_ptr, "TOT_SBoI", "Total I/L");
      
      fits_write_img(out_ptr, TFLOAT, 1, output.image_size[0]*output.image_size[1], 
		     &total_weight_xy[0], &status);
      
      // create and output the total intensity uncertainty image
      fits_create_img(out_ptr, -32, 2, output.image_size, &status);
      check_fits_io(status,"fits_create_image : output_results (total intensity unc/luminosity)");
      output_2d_info(out_ptr, "TOT_SBoI_unc", "Total unc I/L");
      
      fits_write_img(out_ptr, TFLOAT, 1, output.image_size[0]*output.image_size[1], 
		     &total_weight_xy_x2[0], &status);
      
      // create and output the scattered intensity image
      fits_create_img(out_ptr, -32, 2, output.image_size, &status);
      check_fits_io(status,"fits_create_image : output_results (scattered intensity/luminosity)");
      output_2d_info(out_ptr, "SCAT_SBoI", "Scattered I/L");
      
      fits_write_img(out_ptr, TFLOAT, 1, output.image_size[0]*output.image_size[1], 
		     &output.outputs[i].scattered_weight_xy[0], &status);
      
      // create and output the scattered intensity uncertainty image
      fits_create_img(out_ptr, -32, 2, output.image_size, &status);
      check_fits_io(status,"fits_create_image : output_results (scattered intensity unc/luminosity)");
      output_2d_info(out_ptr, "SCAT_SBoI_unc", "Scattered unc I/L");
      
      fits_write_img(out_ptr, TFLOAT, 1, output.image_size[0]*output.image_size[1], 
		     &output.outputs[i].scattered_weight_xy_x2[0], &status);
      
      // create and output the stellar intensity image
      fits_create_img(out_ptr, -32, 2, output.image_size, &status);
      check_fits_io(status,"fits_create_image : output_results (stellar intensity/luminosity)");
      output_2d_info(out_ptr, "STEL_SBoI", "Stellar I/L");
      
      fits_write_img(out_ptr, TFLOAT, 1, output.image_size[0]*output.image_size[1], 
		     &output.outputs[i].stellar_weight_xy[0], &status);
      
      // create and output the scattered intensity uncertainty image
      fits_create_img(out_ptr, -32, 2, output.image_size, &status);
      check_fits_io(status,"fits_create_image : output_results (stellar intensity unc/luminosity)");
      output_2d_info(out_ptr, "STEL_SBoI_unc", "Stellar unc I/L");
      
      fits_write_img(out_ptr, TFLOAT, 1, output.image_size[0]*output.image_size[1], 
		     &output.outputs[i].stellar_weight_xy_x2[0], &status);
      
      // create and output the number of stellar photons image
      fits_create_img(out_ptr, -32, 2, output.image_size, &status);
      check_fits_io(status,"fits_create_image : output_results (number of stellar photons)");
      output_2d_info(out_ptr, "N_STEL_PHOT", "# stellar photons");

      fits_write_img(out_ptr, TFLOAT, 1, output.image_size[0]*output.image_size[1], 
		     &output.outputs[i].num_stellar_photons_xy[0], &status);

      // create and output the number of scattered photons image
      fits_create_img(out_ptr, -32, 2, output.image_size, &status);
      check_fits_io(status,"fits_create_image : output_results (number of scattered photons)");
      output_2d_info(out_ptr, "N_SCAT_PHOT", "# scattered photons");

      fits_write_img(out_ptr, TFLOAT, 1, output.image_size[0]*output.image_size[1], 
		     &output.outputs[i].num_photons_xy[0], &status);

      // close FITS File
      fits_close_file(out_ptr, &status);
      check_fits_io(status,"fits_close_file : output results");
    }
  }
}
