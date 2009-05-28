// ======================================================================
//   Procedure to output the global results of dirty.
//
// 2008 Jan/KDG - written
// 2008 Aug/KDG - output more info into the header
// 2008 sept 2 - added explicit casts in NumUtils::integrate KAM
// ======================================================================
#include "output_global_results.h"
//#define DEBUG_OGR

void output_global_results (runinfo_struct& runinfo,
			    output_struct& output,
			    geometry_struct& geometry)

{
  int i,j = 0;

  // setup the columns of the output table
  int tfields = 9; // min number of fields
  if (runinfo.do_ere_emission) tfields += 4;
  if (runinfo.do_dust_emission && runinfo.do_emission_grain) tfields += 4*(runinfo.n_emission_grain_types);
  vector<string> sttype(tfields);
  vector<string> stform(tfields);
  vector<string> stunit(tfields);

  // column #1
  sttype[0] = "wavelength";
  stunit[0] = "micron";

  // column #2
  sttype[1] = "tau_norm";
  stunit[1] = "ratio";

  // columns for totals
  sttype[2] = "Flux_Input";
  sttype[3] = "Flux";
  sttype[4] = "Flux_Unc";
  sttype[5] = "Flux_rt_d";
  sttype[6] = "Flux_rt_d_unc";
  sttype[7] = "Flux_rt_s";
  sttype[8] = "Flux_rt_s_unc";

  int eg_offset = 9;
  // setup output for ERE if present
  if (runinfo.do_ere_emission) {
    eg_offset += 4;
    sttype[9] = "Flux_ere_d";
    sttype[10] = "Flux_ere_d_unc";
    sttype[11] = "Flux_ere_s";
    sttype[12] = "Flux_ere_s_unc";
  }

  // columns for emission/grain types
  if (runinfo.do_dust_emission && runinfo.do_emission_grain) {
    //    grain_emission_type = GetEmissionLabels(); // from grain model
    // temp stuff
    vector<string> grain_emission_type;
    for (i = 0; i < runinfo.n_emission_grain_types; i++) {
      grain_emission_type.push_back("");
      stringstream ss;
      ss << (i+1);
      ss >> grain_emission_type[i];
    }
      
    // end temp stuff
    for (i = 0; i < runinfo.n_emission_grain_types; i++) {
      sttype[eg_offset+(4*(i))] = "Flux_de_d_" + grain_emission_type[i];
      sttype[eg_offset+1+(4*(i))] = "Flux_de_d_unc_" + grain_emission_type[i];
      sttype[eg_offset+2+(4*(i))] = "Flux_de_s_" + grain_emission_type[i];
      sttype[eg_offset+3+(4*(i))] = "Flux_de_s_unc_" + grain_emission_type[i];
    }
  }

  // now setup the needed char** variables for passing to the cfitsio routine
  //   figured this out from the CCFIT code
  //   for some reason, the strings have to be in a vector (not a single string variable)
  char** ttype = new char*[tfields];
  char** tform = new char*[tfields];
  char** tunit = new char*[tfields];

  for (i = 0; i < tfields; i++) {
//     cout << "i = " << i << " ";
    stform[i] = "E16.6";
    if (i > 1) stunit[i] = "ergs s^-1 um^-1";
    ttype[i] = const_cast<char*>(sttype[i].c_str());
    tform[i] = const_cast<char*>(stform[i].c_str());
    tunit[i] = const_cast<char*>(stunit[i].c_str());
//     cout << sttype[i] << " " << stform[i] << " " << stunit[i] << endl;
  }

#ifdef DEBUG_OGR
  cout << "made it past defining table details" << endl;
#endif

  // now convert the fluxes from ratios to luminosities
  int out_sed_lum_offset = 0;
  vector<double> out_wavelength;
  vector<double> out_sed_lum;
  vector<double> flux_total;
  vector<double> flux_total_unc;
  for (i = 0; i < int(runinfo.wavelength.size()); i++) {
    // output wavelength (do in microns)
    out_wavelength.push_back(runinfo.wavelength[i]*Constant::CM_UM);
    out_sed_lum.push_back(runinfo.sed_lum[i]/Constant::CM_UM);
    // stellar direct
    runinfo.out_sed_lum[0][i] *= out_sed_lum[i];
    runinfo.out_sed_lum_unc[0][i] *= out_sed_lum[i];
    // stellar scattered
    runinfo.out_sed_lum[1][i] *= out_sed_lum[i];
    runinfo.out_sed_lum_unc[1][i] *= out_sed_lum[i];
    // start total
    flux_total.push_back(runinfo.out_sed_lum[0][i] + runinfo.out_sed_lum[1][i]);
    flux_total_unc.push_back(runinfo.out_sed_lum_unc[0][i]*runinfo.out_sed_lum_unc[0][i] +
			     runinfo.out_sed_lum_unc[1][i]*runinfo.out_sed_lum_unc[1][i]);
    out_sed_lum_offset = 2;
    if (runinfo.do_ere_emission) {
      // ere direct
      runinfo.out_sed_lum[out_sed_lum_offset][i] *= runinfo.emitted_ere_lum[0][i]/Constant::CM_UM;
      runinfo.out_sed_lum_unc[out_sed_lum_offset][i] *= runinfo.emitted_ere_lum[0][i]/Constant::CM_UM;
      // ere scattered
      runinfo.out_sed_lum[out_sed_lum_offset+1][i] *= runinfo.emitted_ere_lum[0][i]/Constant::CM_UM;
      runinfo.out_sed_lum_unc[out_sed_lum_offset+1][i] *= runinfo.emitted_ere_lum[0][i]/Constant::CM_UM;

      // add the total ere emitted component to the total
      flux_total[i] += runinfo.out_sed_lum[out_sed_lum_offset][i] + runinfo.out_sed_lum[1+out_sed_lum_offset][i];
      flux_total_unc[i] += runinfo.out_sed_lum_unc[out_sed_lum_offset][i]*runinfo.out_sed_lum_unc[out_sed_lum_offset][i] +
	runinfo.out_sed_lum_unc[1+out_sed_lum_offset][i]*runinfo.out_sed_lum_unc[1+out_sed_lum_offset][i];
      flux_total_unc[i] = pow(flux_total_unc[i],0.5);

      out_sed_lum_offset += 2;
    }
    if (runinfo.do_dust_emission) {
      // offset (adjust if ERE also)
      for (j = 0; j < runinfo.n_emission_grain_types; j++) {
	// dust thermal emitted direct
	runinfo.out_sed_lum[out_sed_lum_offset+(2*j)][i] *= runinfo.emitted_lum[0][i]/Constant::CM_UM;
	runinfo.out_sed_lum_unc[out_sed_lum_offset+(2*j)][i] *= runinfo.emitted_lum[0][i]/Constant::CM_UM;
	// dust thermal emitted scattered
	runinfo.out_sed_lum[1+out_sed_lum_offset+(2*j)][i] *= runinfo.emitted_lum[0][i]/Constant::CM_UM;
	runinfo.out_sed_lum_unc[1+out_sed_lum_offset+(2*j)][i] *= runinfo.emitted_lum[0][i]/Constant::CM_UM;
      }
      // add the total thermal emitted component to the total
      flux_total[i] += runinfo.out_sed_lum[out_sed_lum_offset][i] + runinfo.out_sed_lum[1+out_sed_lum_offset][i];
      flux_total_unc[i] += runinfo.out_sed_lum_unc[out_sed_lum_offset][i]*runinfo.out_sed_lum_unc[out_sed_lum_offset][i] +
	runinfo.out_sed_lum_unc[1+out_sed_lum_offset][i]*runinfo.out_sed_lum_unc[1+out_sed_lum_offset][i];
      flux_total_unc[i] = pow(flux_total_unc[i],0.5);
    }
  }

#ifdef DEBUG_OGR
  cout << "made it past filling in table" << endl;
#endif

  // determine the global energy absorbed and emitted
  double total_stellar_energy = NumUtils::integrate<double>(out_wavelength,out_sed_lum);
  double total_rt_direct_energy = NumUtils::integrate<double>(out_wavelength,runinfo.out_sed_lum[0]);
  double total_rt_scat_energy = NumUtils::integrate<double>(out_wavelength,runinfo.out_sed_lum[1]);
  double total_de_direct_energy = 0.0;
  double total_de_scat_energy = 0.0;
  if (runinfo.do_emission_grain) {
    total_de_direct_energy = NumUtils::integrate<double>(out_wavelength,runinfo.out_sed_lum[out_sed_lum_offset]);
    total_de_scat_energy = NumUtils::integrate<double>(out_wavelength,runinfo.out_sed_lum[out_sed_lum_offset+1]);
  }
  if (runinfo.verbose >= 1) {
    cout << "Emitted(total) = " << total_stellar_energy << endl;
    cout << "Emitted(rt direct) = " << total_rt_direct_energy << endl;
    cout << "Emitted(rt scat) = " << total_rt_scat_energy << endl;
    if (runinfo.do_emission_grain) {
      cout << "Emitted(de direct) = " << total_de_direct_energy << endl;
      cout << "Emitted(de scat) = " << total_de_scat_energy << endl;
      cout << "Absorbed(rt) = " << (total_stellar_energy - total_rt_direct_energy - total_rt_scat_energy) << endl;
      cout << "Emitted(de) = " << (total_de_direct_energy - total_de_scat_energy) << endl;
    }
  }

  // filename of the current output file
  string filename = "!" + output.file_base;
  // add the extra string (designates global luminosities
  filename += "_global_lum.table.fits";

  // output results as a FITS table
  fitsfile *out_ptr;   // pointer to the output fits file
  int status = 0;
  fits_create_file(&out_ptr, filename.c_str(), &status);
  check_fits_io(status, "fits_create_file : output_global_results");
  
  fits_create_tbl(out_ptr, ASCII_TBL, 0, tfields, &ttype[0], &tform[0], &tunit[0], "Global_Outputs", &status);
  check_fits_io(status,"fits_create_tbl : output_global_results");
  
  // put parameter file into header
  fits_params_to_header(runinfo.param_filename, out_ptr);

  // return to regularly scheduled work

  double thmass = geometry.total_h_mass;
  fits_write_key(out_ptr, TDOUBLE, "TOT_HMASS", &thmass, "Total H atom mass", &status);
  check_fits_io(status,"fits_write_key : output_global_results");

  if (runinfo.do_dust_emission) {
    fits_write_key(out_ptr, TDOUBLE, "TOT_ABS", &runinfo.total_absorbed_energy, "total absorbed energy", &status);
    check_fits_io(status,"fits_write_key : output_global_results, abs energy");
    fits_write_key(out_ptr, TDOUBLE, "TOT_EMIT", &runinfo.total_emitted_energy, "total emitted energy", &status);
    check_fits_io(status,"fits_write_key : output_global_results, emit energy");

    fits_write_key(out_ptr, TLONG, "N_CELLS", &geometry.num_cells, "# cells in model", &status);
    check_fits_io(status,"fits_write_key : output_global_results, cells");
    fits_write_key(out_ptr, TLONG, "N_ENOUGH", &runinfo.num_cells_enough, "# cells with enough for dust emission", &status);
    check_fits_io(status,"fits_write_key : output_global_results, cells");
    fits_write_key(out_ptr, TLONG, "N_NOT_EN", &runinfo.num_cells_not_enough, "# cells with not enough for dust emission", &status);
    check_fits_io(status,"fits_write_key : output_global_results, cells");
    fits_write_key(out_ptr, TLONG, "N_ZERO", &runinfo.num_cells_zero, "# cells with zero absorbed energy", &status);
    check_fits_io(status,"fits_write_key : output_global_results, cells");
    fits_write_key(out_ptr, TLONG, "N_FEWWAV", &runinfo.num_cells_too_few_waves, "# cells with too few wavelengths for dust emission", &status);
    check_fits_io(status,"fits_write_key : output_global_results, cells");
  }

  // write the columns
  fits_write_col(out_ptr, TDOUBLE, 1, 1, 0, runinfo.n_waves, &out_wavelength[0], &status);
  check_fits_io(status,"fits_write_col : output global results, wavelength column");

  fits_write_col(out_ptr, TFLOAT, 2, 1, 0, runinfo.n_waves, &runinfo.tau_to_tau_ref[0], &status);
  check_fits_io(status,"fits_write_col : output global results, tau_to_tau_ref");

  fits_write_col(out_ptr, TDOUBLE, 3, 1, 0, runinfo.n_waves, &out_sed_lum[0], &status);
  check_fits_io(status,"fits_write_col : output global results, out_sed_lum");

  fits_write_col(out_ptr, TDOUBLE, 4, 1, 0, runinfo.n_waves, &flux_total[0], &status);
  check_fits_io(status,"fits_write_col : output global results, flux_total");

  fits_write_col(out_ptr, TDOUBLE, 5, 1, 0, runinfo.n_waves, &flux_total_unc[0], &status);
  check_fits_io(status,"fits_write_col : output global results, flux_total_unc");

  int n_col = 1;
  if (runinfo.do_ere_emission) n_col += 1;
  if (runinfo.do_dust_emission) n_col += runinfo.n_emission_grain_types;
  for (i = 0; i < n_col; i++) {
    fits_write_col(out_ptr, TDOUBLE, (6+(4*i)), 1, 0, runinfo.n_waves, &runinfo.out_sed_lum[2*i][0], &status);
    fits_write_col(out_ptr, TDOUBLE, (7+(4*i)), 1, 0, runinfo.n_waves, &runinfo.out_sed_lum_unc[2*i][0], &status);
    fits_write_col(out_ptr, TDOUBLE, (8+(4*i)), 1, 0, runinfo.n_waves, &runinfo.out_sed_lum[(2*i)+1][0], &status);
    fits_write_col(out_ptr, TDOUBLE, (9+(4*i)), 1, 0, runinfo.n_waves, &runinfo.out_sed_lum_unc[(2*i)+1][0], &status);
    check_fits_io(status,"fits_write_col : output global results, grain emission components");
  }
  
  // close FITS File
  fits_close_file(out_ptr, &status);
  check_fits_io(status,"fits_close_file : output global results");
  
}
