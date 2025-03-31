#ifndef _DIRTY_SETUP_DUST_GRID_DEXP_DISK_
#define _DIRTY_SETUP_DUST_GRID_DEXP_DISK_

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

extern void setup_dust_grid_subdivide_overdense_cells(geometry_struct &geometry, int spherical_clumps);

// determines the vector of positions and grid numbers given the photon
// position and direction
extern void determine_photon_position_index_initial(geometry_struct &geometry, photon_data &photon);

// determines the photon trajectory (returns the distance and tau traveled)
// extern double calc_photon_trajectory (photon_data& photon,
// 				      geometry_struct& geometry,
// 				      double target_tau,
// 				      double target_dist,
// 				      int& escape,
// 				      double& tau_traveled);
#endif
