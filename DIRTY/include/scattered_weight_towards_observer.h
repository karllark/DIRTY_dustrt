#ifndef _DIRTY_SCATTERED_WEIGHT_TOWARDS_OBSERVER__
#define _DIRTY_SCATTERED_WEIGHT_TOWARDS_OBSERVER__

#include <cmath>
#include <iostream>

#include "debug.h"
#include "geometry_def.h"
#include "photon_data.h"

// determines the photon trajectory (returns the distance and tau traveled)

// determines the photon trajectory (returns the distance and tau traveled)
extern double calc_photon_trajectory(photon_data &photon, geometry_struct &geometry, double target_tau,
                                     double target_dist, int &escape, double &tau_traveled, int repeat_boundary);
#endif
