#ifndef _DIRTY_OUTPUT_RESULTS_
#define _DIRTY_OUTPUT_RESULTS_

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

#include "geometry_def.h"
#include "output_def.h"
#include "runinfo_def.h"
#include "fitsio.h"
#include "constants.h"
#include "debug.h"

//**********************************************************************
// external function definitions

extern int check_fits_io(int status,
			 char text[100]);

#endif
