// ======================================================================
//   Procedure to output the global results of dirty.
//
// 2008 Jan/KDG - written
// ======================================================================
#include "output_global_results.h"

void output_global_results (runinfo_struct& runinfo,
			    output_struct& output)

{
  int i = 0;

  // setup the columns of the output table
  int tfields = 2; // min number of fields
  vector<string> sttype(tfields);
  vector<string> stform(tfields);
  vector<string> stunit(tfields);

  // column #1
  sttype[0] = "wavelength";
  stform[0] = "E16.6";
  stunit[0] = "cm";

  // column #2
  sttype[1] = "Flux";
  stform[1] = "E16.6";
  stunit[1] = "ratio";

  // now setup the needed char** variables for passing to the cfitsio routine
  //   figured this out from the CCFIT code
  //   for some reason, the strings have to be in a vector (not a single string variable)
  char** ttype = new char*[tfields];
  char** tform = new char*[tfields];
  char** tunit = new char*[tfields];

  for (i = 0; i < tfields; i++) {
    ttype[i] = const_cast<char*>(sttype[i].c_str());
    tform[i] = const_cast<char*>(stform[i].c_str());
    tunit[i] = const_cast<char*>(stunit[i].c_str());
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
  
  // close FITS File
  fits_close_file(out_ptr, &status);
  check_fits_io(status,"fits_close_file : output global results");
  
}
