#ifndef _DIRTY_NEXT_SCATTER_
#define _DIRTY_NEXT_SCATTER_

#include <iostream>
#include <cmath>

#include "geometry_def.h"
#include "photon_data.h"
#include "random_dirty.h"
#include "roundoff_err.h"
#include "debug.h"

#define MAX_NUM_SCAT 50

//**********************************************************************
// external function definitions

// determines the photon trajectory (returns the distance and tau traveled)
extern double calc_photon_trajectory (photon_data& photon,
				      geometry_struct& geometry,
				      double target_tau,
				      int& escape,
				      double& tau_traveled);

extern void deposit_energy (const photon_data& photon,
                     geometry_struct& geometry,
                     const double photon_weight);

extern void move_photon (photon_data& photon,
                  const photon_data& dummy_photon,
                  geometry_struct& geometry,
                  const double target_tau,
                  double& tau_traveled);

#endif
