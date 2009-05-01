// ======================================================================
//   Header file for radiative transfer many waves procedure.  
// Include files and function definitions.
//
// 2008 Dec/KDG - written
// ======================================================================

#ifndef _DIRTY_RADIATIVE_TRANSFER_MANY_WAVES_
#define _DIRTY_RADIATIVE_TRANSFER_MANY_WAVES_

#include <iostream>
#include <cmath>

#include "geometry_def.h"
#include "runinfo_def.h"
#include "output_def.h"
#include "photon_data.h"
#include "random_dirty.h"
#include "debug.h"
#include "rt_types.h"

//**********************************************************************
// external function definitions

// get the dust scattering parameters for the current wavelength
extern void get_dust_scat_parameters (int i,
				      runinfo_struct& runinfo,
				      geometry_struct& geometry);

// do the radiative transfer
extern void radiative_transfer (geometry_struct& geometry,
				runinfo_struct& runinfo,
				output_struct& output,
				photon_data& photon,
				random_dirty random_obj);

// output the results of the calculations
extern void output_results (output_struct& output,
			    geometry_struct& geometry,
			    runinfo_struct& runinfo,
			    int index);

// setup absorbed energy grid
extern void setup_absorbed_energy_grid (geometry_struct& geometry,
					runinfo_struct& runinfo,
					int doing_dust_emission);

// store absorbed energy grid (memory or disk)
extern void store_absorbed_energy_grid (geometry_struct& geometry,
					runinfo_struct& runinfo,
					output_struct& output,
					int index,
					int doing_dust_emission);

// check absorbed energy grid
extern void check_absorbed_energy_grid (geometry_struct& geometry,
					runinfo_struct& runinfo);

#endif
