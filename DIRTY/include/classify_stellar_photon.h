// ======================================================================
//   Header file for classify stellar photon procedure.
// Include files and function definitions.
//
// 2005 May/KDG - written
// ======================================================================

#ifndef _DIRTY_CLASSIFY_STELLAR_PHOTON_
#define _DIRTY_CLASSIFY_STELLAR_PHOTON_

#include <cmath>
#include <iostream>

#include "debug.h"
#include "geometry_def.h"
#include "output_def.h"
#include "photon_data.h"
#include "random_dirty.h"
#include "runinfo_def.h"

//**********************************************************************
// external function definitions

extern void rotate_zaxis_for_observer(float transform[3][3],
                                      photon_data &photon);

extern void determine_photon_position_index_initial(geometry_struct &geometry,
                                                    photon_data &photon);

extern double stellar_weight_towards_observer(photon_data photon,
                                              geometry_struct &geometry,
                                              float observer_position[3]);

extern void compute_observer_trans_matrix(output_struct &output,
                                          geometry_struct &geometry, int i);
#endif
