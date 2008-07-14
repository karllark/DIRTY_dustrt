// ======================================================================
//   Header file for classify scattered photon procedure.  
// Include files and function definitions.
//
// 2006 Apr/KDG - written
// ======================================================================

#ifndef _DIRTY_CLASSIFY_SCATTERED_PHOTON_
#define _DIRTY_CLASSIFY_SCATTERED_PHOTON_

#include <iostream>
#include <cmath>

#include "output_def.h"
#include "geometry_def.h"
#include "photon_data.h"
#include "runinfo_def.h"
#include "debug.h"

//**********************************************************************
// external function definitions

/* extern void determine_photon_position_index_initial (geometry_struct& geometry, */
/*   					             photon_data& photon); */

extern double scattered_weight_towards_observer (photon_data photon,
						 geometry_struct& geometry,
						 float observer_position[3]);

extern void rotate_zaxis_for_observer (float transform[3][3],
				       photon_data& photon);

#endif
