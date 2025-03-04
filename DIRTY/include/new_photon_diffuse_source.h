#ifndef _DIRTY_NEW_PHOTON_DIFFUSE_SOURCE_
#define _DIRTY_NEW_PHOTON_DIFFUSE_SOURCE_

#include <cmath>
#include <iostream>

#include "debug.h"
#include "geometry_def.h"
#include "photon_data.h"
#include "random_dirty.h"

// determines the vector of positions and grid numbers given the photon
// position and direction
extern void determine_photon_position_index_initial(geometry_struct& geometry,
                                                    photon_data& photon);

// determines the photon trajectory (returns the distance and tau traveled)
extern double calc_photon_trajectory(photon_data& photon,
                                     geometry_struct& geometry,
                                     double target_tau, double target_dist,
                                     int& escape, double& tau_traveled,
                                     int repeat_boundary);

// just applies the rotation
extern void rotate_zaxis_for_observer(float transform[3][3],
                                      photon_data& photon);

#endif
