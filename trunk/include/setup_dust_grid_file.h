#ifndef _DIRTY_SETUP_DUST_GRID_FILE_
#define _DIRTY_SETUP_DUST_GRID_FILE_

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

#include "ConfigFile.h"
#include "check_input_param.h"
#include "geometry_def.h"
#include "photon_data.h"
#include "fitsio.h"
#include "debug.h"

//**********************************************************************
// external function definitions

//**********************************************************************
// external function definitions

extern int check_fits_io(int status,
			 char text[100]);

extern void setup_dust_grid_subdivide_overdense_cells (geometry_struct& geometry,
						       int spherical_clumps);

#endif
