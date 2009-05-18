// ======================================================================
//   Procedure to output the parameter file to a header
//
// 2009 May/KDG - written
// ======================================================================
#include "fits_params_to_header.h"

int fits_params_to_header(string param_filename,
			  fitsfile *out_ptr)

{
  ifstream param_file(param_filename.c_str());

  string line;
  int status = 0;

  while (getline(param_file,line)) {
    
    //    if (!line.length()) continue;

    fits_write_comment(out_ptr, line.c_str(), &status);
    check_fits_io(status,"fits_write_comment : params to header");
  }

  return 0;
}
