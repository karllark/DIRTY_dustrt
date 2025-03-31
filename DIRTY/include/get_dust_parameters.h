#ifndef _DIRTY_GET_DUST_PARAMETERS_
#define _DIRTY_GET_DUST_PARAMETERS_

#include <algorithm>
#include <cstring>

#include "ConfigFile.h"
#include "DataFile.h"
#include "GrainModel.h"
#include "check_input_param.h"
#include "geometry_def.h"
#include "runinfo_def.h"

//**********************************************************************
// external function definitions

extern void get_wave_grid (ConfigFile &param_data, runinfo_struct &runinfo);

#endif
