#ifndef _DIRTY_SETUP_DUST_GRID_
#define _DIRTY_SETUP_DUST_GRID_

#include <iostream>
#include <cmath>
//#include <cstring>

#include "ConfigFile.h"
#include "DataFile.h"
#include "check_input_param.h"
#include "geometry_def.h"
#include "photon_data.h"
#include "random_dirty.h"
#include "constants.h"
#include "debug.h"

//**********************************************************************
// external function definitions

extern void setup_dust_grid_sphere (ConfigFile& param_data,
				    geometry_struct& geometry,
				    random_dirty& random_obj);

extern void setup_dust_grid_shell (ConfigFile& param_data,
				   geometry_struct& geometry,
				   random_dirty& random_obj);

extern void setup_dust_grid_slab (ConfigFile& param_data,
				  geometry_struct& geometry,
				  random_dirty& random_obj);

extern void verify_dust_grid (geometry_struct& geometry);

#endif
