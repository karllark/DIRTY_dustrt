// ======================================================================
//   Header file for DIRTY model.  Include files and function
// definitions.
//
// 2003 Apr/KDG - started
// 2003 Jun/KDG - updated with active grid tracking
// 2004 Dec/KDG - updated with photon scattering
// 2005 May/KDG - updated with output stuff
// 2007 Apr/KDG - added dust grain model
// ======================================================================

#ifndef _DIRTY_
#define _DIRTY_

#include <iostream>
#include <cmath>

#include "debug.h"
#include "geometry_def.h"
#include "output_def.h"
#include "runinfo_def.h"
#include "photon_data.h"
#include "random_dirty.h"
#include "ConfigFile.h"
#include "GrainModel.h"   

//**********************************************************************
// external function definitions

// get the run parameters
extern void get_run_parameters (ConfigFile& param_data,
				output_struct& output,
				geometry_struct& geometry,
				runinfo_struct& runinfo);

// get the dust parameters
extern void get_dust_parameters (ConfigFile& param_data,
				 GrainModel& CurGrainModel,
				 runinfo_struct& runinfo);

// ger the sed parameters
extern void get_sed_parameters (ConfigFile& param_data,
				runinfo_struct& runinfo);

// sets up the dust grid 
extern void setup_dust_grid (ConfigFile& param_data,
			     geometry_struct& geometry,
			     photon_data& photon,
			     random_dirty& random_obj);

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
					int wave_index,
					int doing_dust_emission);

// store absorbed energy grid (memory or disk)
extern void store_absorbed_energy_grid (geometry_struct& geometry,
					runinfo_struct& runinfo,
					output_struct& output,
					int index,
					int doing_dust_emission);

// get the thermal dust emission in each grid cell
extern void get_dust_thermal_emission (geometry_struct& geometry,
				       runinfo_struct& runinfo, 
				       GrainModel& CurGrainModel);

// setup the emitted grid for monte carlo radiative transfer
extern void setup_emitted_grid_for_montecarlo (geometry_struct& geometry,
					       runinfo_struct& runinfo,
					       GrainModel& CurGrainModel);

// output global integrated quantities
extern void output_global_results (runinfo_struct& runinfo);

#endif
