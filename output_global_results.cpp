// ======================================================================
//   Procedure to output the global results of dirty.
//
// 2008 Jan/KDG - written
// ======================================================================
#include "output_results.h"

void output_global_results (runinfo_struct& runinfo)

{
  int i = 0;

  // filename of the current output file
  string filename = "!" + output.file_base;
  // add the extra string (designates global luminosities
  filename += "_global_lum.table.fits";

  // setup the columns of the output table
  int tfields = 6; // min number of fields
  // these should be char *[] - need to figure out how to declare this
  vector<string> ttype;  // column name
  vector<string> tform;  // column form (cfitsio definitions)
  vector<string> tunit;  // units of column
  
  

  // output results as a FITS table
  fitsfile *out_ptr;   // pointer to the output fits file
  int status = 0;
  fits_create_tbl(&out_ptr, ASCII_TBL, 0, ttype, tform, tunit, "Global_Outputs", &status);
  check_fits_io(status,"fits_create_tbl : output_global_results");
  
  // write the columns
  fits_write_col(&out_ptr, TFLOAT, 1, 1, 0, runinfo.n_waves, &runinfo.wavelength[0], &status);
  check_fits_io(status,"fits_write_col : output global results, wavelength column");

  // close FITS File
  fits_close_file(out_ptr, &status);
  check_fits_io(status,"fits_close_file : output global results");

}
