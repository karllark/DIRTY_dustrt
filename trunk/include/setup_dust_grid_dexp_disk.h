#ifndef _DIRTY_SETUP_DUST_GRID_DEXP_DISK_
#define _DIRTY_SETUP_DUST_GRID_DEXP_DISK_

#include <iostream>
#include <cmath>

#include "ConfigFile.h"
#include "check_input_param.h"
#include "geometry_def.h"
#include "photon_data.h"
#include "random_dirty.h"
#include "debug.h"

//**********************************************************************
// external function definitions

extern void setup_dust_grid_subdivide_overdense_cells (geometry_struct& geometry,
						       int spherical_clumps);

#endif
