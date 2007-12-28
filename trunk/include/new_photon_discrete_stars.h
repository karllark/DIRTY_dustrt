#ifndef _DIRTY_NEW_PHOTON_DISCRETE_STARS__
#define _DIRTY_NEW_PHOTON_DISCRETE_STARS__

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

#endif
