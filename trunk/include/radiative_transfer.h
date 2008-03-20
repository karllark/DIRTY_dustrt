// ======================================================================
//   Header file for radiative transfer procedure.  
// Include files and function definitions.
//
// 2004 Dec/KDG - written
// 2005 May/KDG - added output
// ======================================================================

#ifndef _DIRTY_RADIATIVE_TRANSFER_
#define _DIRTY_RADIATIVE_TRANSFER_

#include <iostream>
#include <cmath>

#include "geometry_def.h"
#include "runinfo_def.h"
#include "output_def.h"
#include "photon_data.h"
#include "random_dirty.h"
#include "debug.h"

//**********************************************************************
// external function definitions

// initialize output structures
extern void initialize_output (output_struct& output,
			       geometry_struct& geometry);

// generates photons for the single star case
extern void new_photon (photon_data& photon,
			geometry_struct& geometry,
			runinfo_struct& runinfo,
			random_dirty random_obj);

// determines the first scattering site (forced)
extern void forced_first_scatter (geometry_struct& geometry,
				  photon_data& photon,
				  random_dirty& random_obj);


// determines the next scattering site (2nd, 3rd, etc.)
extern int next_scatter (geometry_struct& geometry,
			 photon_data& photon,
			 random_dirty& random_obj);

// scatters the photon into a new direction
extern void scatter_photon (geometry_struct& geometry,
			    photon_data& photon,
			    random_dirty& random_obj);

// classifies stellar photon into output
extern void classify_stellar_photon (output_struct& output,
				     photon_data& photon,
				     geometry_struct& geometry,
				     runinfo_struct& runinfo);

// classifies scattered photon into output
extern void classify_scattered_photon (output_struct& output,
				       photon_data& photon,
				       geometry_struct& geometry,
				       runinfo_struct& runinfo);


#endif
