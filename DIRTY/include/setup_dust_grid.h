#ifndef _DIRTY_SETUP_DUST_GRID_
#define _DIRTY_SETUP_DUST_GRID_

#include <cmath>
#include <iostream>
// #include <cstring>

#include "ConfigFile.h"
#include "Constants.h"
#include "DataFile.h"
#include "check_input_param.h"
#include "debug.h"
#include "geometry_def.h"
#include "photon_data.h"
#include "random_dirty.h"

//**********************************************************************
// external function definitions

extern void setup_dust_grid_sphere(ConfigFile &param_data, geometry_struct &geometry, random_dirty &random_obj);

extern void setup_dust_grid_shell(ConfigFile &param_data, geometry_struct &geometry, random_dirty &random_obj);

extern void setup_dust_grid_slab(ConfigFile &param_data, geometry_struct &geometry, random_dirty &random_obj);

extern void setup_dust_grid_dexp_disk(ConfigFile &param_data, geometry_struct &geometry, random_dirty &random_obj);

extern void setup_dust_grid_file(ConfigFile &param_data, geometry_struct &geometry);

extern void verify_dust_grid(geometry_struct &geometry);

#endif
