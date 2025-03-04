#ifndef _DIRTY_CALC_PHOTON_TRAJECTORY_
#define _DIRTY_CALC_PHOTON_TRAJECTORY_

#include <cmath>
#include <iostream>

#include "debug.h"
#include "geometry_def.h"
#include "photon_data.h"
#include "roundoff_err.h"

// determines the vector of positions and grid numbers given the photon
// position and direction
extern void determine_grid_position_index (geometry_struct &geometry, photon_data &photon);

// function to calculate the distance traveled inside a cell
extern double calc_delta_dist (photon_data &photon, geometry_struct &geometry, double target_tau,
                               double target_dist, int &escape, double &tau_traveled);

extern void determine_photon_position_index (geometry_struct &geometry, photon_data &photon);

#endif
