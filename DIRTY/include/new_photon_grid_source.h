#ifndef _DIRTY_NEW_PHOTON_GRID_SOURCE_
#define _DIRTY_NEW_PHOTON_GRID_SOURCE_

#include <iostream>
#include <cmath>

#include "geometry_def.h"
#include "random_dirty.h"
#include "photon_data.h"
#include "runinfo_def.h"
#include "debug.h"

// determines the vector of positions and grid numbers given the photon
// position and direction
extern void determine_photon_position_index_initial (geometry_struct& geometry,
						     photon_data& photon);

#endif
