#ifndef _DIRTY_FORCED_FIRST_SCATTER_
#define _DIRTY_FORCED_FIRST_SCATTER_

#include <iostream>
#include <cmath>

#include "geometry_def.h"
#include "photon_data.h"
#include "random_dirty.h"
#include "roundoff_err.h"
#include "debug.h"

//**********************************************************************
// external function definitions

// determines the photon trajectory (returns the distance and tau traveled)
extern double calc_photon_trajectory (photon_data& photon,
				      geometry_struct& geometry,
				      double target_tau,
				      int& escape,
				      double& tau_traveled);

#endif
