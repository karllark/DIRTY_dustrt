#ifndef _DIRTY_NEW_PHOTON_DIFFUSE_SOURCE_
#define _DIRTY_NEW_PHOTON_DIFFUSE_SOURCE_

#include <iostream>
#include <cmath>

#include "geometry_def.h"
#include "random_dirty.h"
#include "photon_data.h"
#include "debug.h"

// determines the vector of positions and grid numbers given the photon
// position and direction
extern void determine_photon_position_index_initial (geometry_struct& geometry,
						     photon_data& photon);

// determines the photon trajectory (returns the distance and tau traveled)
extern double calc_photon_trajectory (photon_data& photon,
				      geometry_struct& geometry,
				      double target_tau,
				      int& escape,
				      double& tau_traveled);

// just applies the rotation
extern void rotate_zaxis_for_observer (float transform[3][3],
				       photon_data& photon);

#endif
