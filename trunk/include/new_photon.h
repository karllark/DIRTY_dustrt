#ifndef _DIRTY_NEW_PHOTON_
#define _DIRTY_NEW_PHOTON_

#include <iostream>
#include <cmath>

#include "geometry_def.h"
#include "random_dirty.h"
#include "photon_data.h"
#include "runinfo_def.h"
#include "debug.h"

// new photon for discrete stars
extern void new_photon_discrete_stars (photon_data& photon,
				       geometry_struct& geometry,
				       random_dirty random_obj);

// new photon for diffuse sources
extern void new_photon_diffuse_source (photon_data& photon,
				       geometry_struct& geometry,
				       random_dirty random_obj);

// new photon for grid sources
extern void new_photon_grid_source (photon_data& photon,
				    geometry_struct& geometry,
				    runinfo_struct& runinfo,
				    random_dirty random_obj);

// new photon for dexp disk
extern void new_photon_dexp_disk (photon_data& photon,
				  geometry_struct& geometry,
				  random_dirty random_obj);

// new photon for power sphere
extern void new_photon_pow_sphere (photon_data& photon,
				     geometry_struct& geometry,
				     random_dirty random_obj);

#endif
