#ifndef _DIRTY_CALC_DELTA_DIST__
#define _DIRTY_CALC_DELTA_DIST__

#include <cmath>
#include <iostream>

#include "debug.h"
#include "geometry_def.h"
#include "photon_data.h"
#include "roundoff_err.h"

// determines the photon trajectory (returns the distance and tau traveled)
extern double calc_photon_trajectory(photon_data &photon, geometry_struct &geometry, double target_tau,
                                     double target_dist, int &escape, double &tau_traveled, int repeat_boundary);

extern void determine_photon_position_index(geometry_struct &geometry, photon_data &photon);

#endif
