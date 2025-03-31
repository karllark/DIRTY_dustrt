#ifndef _DIRTY_NEW_PHOTON_DEXP_DISK__
#define _DIRTY_NEW_PHOTON_DEXP_DISK__

#include <cmath>
#include <iostream>

#include "debug.h"
#include "geometry_def.h"
#include "photon_data.h"
#include "random_dirty.h"

// determines the vector of positions and grid numbers given the photon
// position and direction
extern void determine_photon_position_index_initial (geometry_struct &geometry,
                                                     photon_data &photon);

#endif
