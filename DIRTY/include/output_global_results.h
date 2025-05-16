#ifndef _DIRTY_OUTPUT_GLOBAL_RESULTS_
#define _DIRTY_OUTPUT_GLOBAL_RESULTS_

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "Constants.h"
#include "NumUtils.h"
#include "debug.h"
#include "fitsio.h"
#include "geometry_def.h"
#include "output_def.h"
#include "runinfo_def.h"

/* #define OGR_MAX_CHAR_LEN 20 */

//**********************************************************************
// external function definitions

extern int check_fits_io(int status, const char text[100]);

extern int fits_params_to_header(string param_filename, fitsfile *out_ptr);

#endif
