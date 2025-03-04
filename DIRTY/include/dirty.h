// ======================================================================
//   Header file for DIRTY model.  Include files and function
// definitions.
//
// 2003 Apr/KDG - started
// 2003 Jun/KDG - updated with active grid tracking
// 2004 Dec/KDG - updated with photon scattering
// 2005 May/KDG - updated with output stuff
// 2007 Apr/KDG - added dust grain model
// 2008 Mar/KDG - add details for global output
// ======================================================================

#ifndef _DIRTY_
#define _DIRTY_

#include <cmath>
#include <iostream>

#include "ConfigFile.h"
#include "DirtyFailure.h"
#include "GrainModel.h"
#include "debug.h"
#include "geometry_def.h"
#include "output_def.h"
#include "photon_data.h"
#include "random_dirty.h"
#include "rt_types.h"
#include "runinfo_def.h"

//**********************************************************************
// external function definitions

// get the run parameters
extern void get_run_parameters(ConfigFile &param_data, output_struct &output, geometry_struct &geometry,
                               runinfo_struct &runinfo);

// get the dust parameters
extern void get_dust_parameters(ConfigFile &param_data, GrainModel &CurGrainModel, geometry_struct &geometry,
                                runinfo_struct &runinfo);

// ger the sed parameters
extern void get_sed_parameters(ConfigFile &param_data, runinfo_struct &runinfo, GrainModel &CurGrainModel);

// sets up the dust grid
extern void setup_dust_grid(ConfigFile &param_data, geometry_struct &geometry, photon_data &photon,
                            random_dirty &random_obj);

// do the radiative transfer over all the wavelengths
extern void radiative_transfer_many_waves(geometry_struct &geometry, runinfo_struct &runinfo, output_struct &output,
                                          photon_data &photon, random_dirty &random_obj, int rt_type, int iter_num);

// setup the ERE dust emission output
extern void setup_ere_dust_emission_output(output_struct &ere_output, output_struct &output);

// get the ere dust emission in each grid cell
extern void get_dust_ere_emission(geometry_struct &geometry, runinfo_struct &runinfo);

// get the thermal dust emission in each grid cell
extern void get_dust_thermal_emission(geometry_struct &geometry, runinfo_struct &runinfo, GrainModel &CurGrainModel,
                                      DirtyFailure *Failure);

// setup for the thermal dust emission output (used for separate grain/emission
// output)
extern void setup_thermal_dust_emission_output(runinfo_struct &runinfo, output_struct &de_output, output_struct &output,
                                               photon_data &photon);

// setup the emitted grid for monte carlo radiative transfer
extern void setup_emitted_grid_for_montecarlo(geometry_struct &geometry, runinfo_struct &runinfo,
                                              GrainModel &CurGrainModel);

// determine if energy is conserved well enough
extern void check_de_energy_conservation(runinfo_struct &runinfo, int &iter_done);

// output global integrated quantities
extern void output_global_results(runinfo_struct &runinfo, output_struct &output, geometry_struct &geometry);

// output the model grid
extern void output_model_grid(geometry_struct &geometry, output_struct &output, runinfo_struct &runinfo);

#endif
