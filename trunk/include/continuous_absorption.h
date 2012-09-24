#ifndef _DIRTY_CONTINUOUS_ABSORPTION_
#define _DIRTY_CONTINUOUS_ABSORPTION_

#include <iostream>
#include <cmath>
#include <cassert>

#include "geometry_def.h"
#include "photon_data.h"
#include "roundoff_err.h"
#include "debug.h"

//**********************************************************************
// external function definitions

// function to calculate the distance traveled inside a cell
extern double calc_delta_dist (photon_data& photon,
                               geometry_struct& geometry,
                               double target_tau,
                               int& escape,
                               double& tau_traveled);

#endif
