// ======================================================================
//   Procedure to output the global results of dirty.
//
// 2008 Jan/KDG - written
// ======================================================================
#include "output_global_results.h"

void output_global_results (runinfo_struct& runinfo,
			    output_struct& output)

{
  int i,j = 0;

  // setup the columns of the output table
  int tfields = 7; // min number of fields
  if (runinfo.do_emission_grain) tfields += 4*(runinfo.n_emission_grain_types);
  vector<string> sttype(tfields);
  vector<string> stform(tfields);
  vector<string> stunit(tfields);

  // column #1
  sttype[0] = "wavelength";
  stunit[0] = "micron";

  // columns for totals
  sttype[1] = "Flux";
  sttype[2] = "Flux_Unc";
  sttype[3] = "Flux_rt_d";
  sttype[4] = "Flux_rt_d_unc";
  sttype[5] = "Flux_rt_s_";
  sttype[6] = "Flux_rt_s_unc";

  // columns for emission/grain types
  if (runinfo.do_emission_grain) {
    int eg_offset = 7;
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
      sttype[eg_offset+(4*(i-1))] = "Flux_de_d_" + grain_emission_type[i];
      sttype[eg_offset+1+(4*(i-1))] = "Flux_de_d_unc_" + grain_emission_type[i];
      sttype[eg_offset+2+(4*(i-1))] = "Flux_de_s_" + grain_emission_type[i];
      sttype[eg_offset+3+(4*(i-1))] = "Flux_de_s_unc_" + grain_emission_type[i];
    }
  }

  // now setup the needed char** variables for passing to the cfitsio routine
  //   figured this out from the CCFIT code
  //   for some reason, the strings have to be in a vector (not a single string variable)
  char** ttype = new char*[tfields];
  char** tform = new char*[tfields];
  char** tunit = new char*[tfields];

  for (i = 0; i < tfields; i++) {
    stform[i] = "E16.6";
    if (i > 0) stunit[i] = "ergs s^-1";
    ttype[i] = const_cast<char*>(sttype[i].c_str());
    tform[i] = const_cast<char*>(stform[i].c_str());
    tunit[i] = const_cast<char*>(stunit[i].c_str());
//     cout << sttype[i] << " " << stform[i] << " " << stunit[i] << endl;
  }

  // now convert the fluxes from ratios to luminosities
  int out_sed_lum_offset = 0;
  vector<double> flux_total;
  vector<double> flux_total_unc;
  for (i = 0; i < int(runinfo.wavelength.size()); i++) {
    // stellar direct
    runinfo.out_sed_lum[0][i] *= runinfo.sed_lum[i];
    runinfo.out_sed_lum_unc[0][i] *= runinfo.sed_lum[i];
    // stellar scattered
    runinfo.out_sed_lum[1][i] *= runinfo.sed_lum[i];
    runinfo.out_sed_lum_unc[1][i] *= runinfo.sed_lum[i];
    // start total
    flux_total.push_back(runinfo.out_sed_lum[0][i] + runinfo.out_sed_lum[1][i]);
    cout << runinfo.wavelength[i] << " ";
    cout << flux_total[i] << " ";
    flux_total_unc.push_back(runinfo.out_sed_lum_unc[0][i]*runinfo.out_sed_lum_unc[0][i] +
			     runinfo.out_sed_lum_unc[1][i]*runinfo.out_sed_lum_unc[1][i]);
    // offset (adjust if ERE also)
    out_sed_lum_offset = 2;
    for (j = 0; j < runinfo.n_emission_grain_types; j++) {
      // dust thermal emitted direct
      runinfo.out_sed_lum[out_sed_lum_offset+(2*j)][i] *= runinfo.emitted_lum[j][i];
      runinfo.out_sed_lum_unc[out_sed_lum_offset+(2*j)][i] *= runinfo.emitted_lum[j][i];
      // dust thermal emitted scattered
      runinfo.out_sed_lum[1+out_sed_lum_offset+(2*j)][i] *= runinfo.emitted_lum[j][i];
      runinfo.out_sed_lum_unc[1+out_sed_lum_offset+(2*j)][i] *= runinfo.emitted_lum[j][i];
    }
    // add the total thermal emitted component to the total
    flux_total[i] += runinfo.out_sed_lum[out_sed_lum_offset][i] + runinfo.out_sed_lum[1+out_sed_lum_offset][i];
    flux_total_unc[i] += runinfo.out_sed_lum_unc[out_sed_lum_offset][i]*runinfo.out_sed_lum_unc[out_sed_lum_offset][i] +
      runinfo.out_sed_lum_unc[1+out_sed_lum_offset][i]*runinfo.out_sed_lum_unc[1+out_sed_lum_offset][i];
    cout << flux_total[i] << endl;
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
  
  // write the columns
  fits_write_col(out_ptr, TFLOAT, 1, 1, 0, runinfo.n_waves, &runinfo.wavelength[0], &status);
  check_fits_io(status,"fits_write_col : output global results, wavelength column");

  fits_write_col(out_ptr, TDOUBLE, 2, 1, 0, runinfo.n_waves, &flux_total[0], &status);
  check_fits_io(status,"fits_write_col : output global results, wavelength column");
  
  // close FITS File
  fits_close_file(out_ptr, &status);
  check_fits_io(status,"fits_close_file : output global results");
  
}
