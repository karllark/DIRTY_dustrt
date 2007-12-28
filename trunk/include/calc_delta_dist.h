#ifndef _DIRTY_CALC_DELTA_DIST__
#define _DIRTY_CALC_DELTA_DIST__

#include <iostream>
#include <cmath>

#include "geometry_def.h"
#include "photon_data.h"
#include "roundoff_err.h"
#include "debug.h"

// determines the photon trajectory (returns the distance and tau traveled)
extern double calc_photon_trajectory (photon_data& photon,
				      geometry_struct& geometry,
				      double target_tau,
				      int& escape,
				      double& tau_traveled);

extern void determine_photon_position_index (geometry_struct& geometry,
					     photon_data& photon);

#endif
