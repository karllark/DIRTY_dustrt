#ifndef _DIRTY_SETUP_DUST_GRID_SPHERE_
#define _DIRTY_SETUP_DUST_GRID_SPHERE_

#include <cmath>
#include <iostream>

#include "ConfigFile.h"
#include "check_input_param.h"
#include "debug.h"
#include "geometry_def.h"
#include "photon_data.h"
#include "random_dirty.h"

//**********************************************************************
// external function definitions

extern void setup_dust_grid_subdivide_overdense_cells (geometry_struct &geometry, int spherical_clumps);

#endif
