#ifndef _DIRTY_OUTPUT_RESULTS_
#define _DIRTY_OUTPUT_RESULTS_

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "Constants.h"
#include "StringManip.h"
#include "debug.h"
#include "fitsio.h"
#include "geometry_def.h"
#include "output_def.h"
#include "runinfo_def.h"

//**********************************************************************
// external function definitions

extern int check_fits_io(int status, const char text[100]);

extern int fits_params_to_header(string param_filename, fitsfile *out_ptr);

#endif
