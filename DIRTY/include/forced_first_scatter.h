#ifndef _DIRTY_FORCED_FIRST_SCATTER_
#define _DIRTY_FORCED_FIRST_SCATTER_

#include <cmath>
#include <iostream>

#include "debug.h"
#include "geometry_def.h"
#include "photon_data.h"
#include "random_dirty.h"
#include "roundoff_err.h"

//**********************************************************************
// external function definitions

extern void determine_photon_position_index_initial(geometry_struct &geometry, photon_data &photon);

// determines the photon trajectory (returns the distance and tau traveled)
extern double calc_photon_trajectory(photon_data &photon, geometry_struct &geometry, double target_tau,
                                     double target_dist, int &escape, double &tau_traveled, int repeat_boundary);

#endif
