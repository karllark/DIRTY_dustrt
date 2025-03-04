#ifndef _DIRTY_SETUP_DUST_GRID_FILE_
#define _DIRTY_SETUP_DUST_GRID_FILE_

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "ConfigFile.h"
#include "check_input_param.h"
#include "debug.h"
#include "fitsio.h"
#include "geometry_def.h"
#include "photon_data.h"

//**********************************************************************
// external function definitions

extern int check_fits_io (int status, const char text[100]);

extern void setup_dust_grid_check_grid (geometry_struct &geometry, int cur_grid, int par_grid,
                                        vector<int> &par_idim);

extern void setup_dust_grid_subdivide_overdense_cells (geometry_struct &geometry,
                                                       int spherical_clumps);

#endif
