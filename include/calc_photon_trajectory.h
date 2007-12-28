#ifndef _DIRTY_CALC_PHOTON_TRAJECTORY_
#define _DIRTY_CALC_PHOTON_TRAJECTORY_

#include <iostream>
#include <cmath>

#include "debug.h"
#include "geometry_def.h"
#include "photon_data.h"

// determines the vector of positions and grid numbers given the photon
// position and direction
extern void determine_grid_position_index (geometry_struct& geometry,
					   photon_data& photon);

// function to calculate the distance traveled inside a cell
extern double calc_delta_dist (photon_data& photon,
			       geometry_struct& geometry,
			       double target_tau,
			       int& escape,
			       double& tau_traveled);


#endif
