#ifndef _DIRTY_GET_DUST_PARAMETERS_
#define _DIRTY_GET_DUST_PARAMETERS_

#include <cstring>

#include "ConfigFile.h"
#include "DataFile.h"
#include "runinfo_def.h"
#include "check_input_param.h"
#include "GrainModel.h"   

//**********************************************************************
// external function definitions

extern void get_wave_grid (ConfigFile& param_data,
			   runinfo_struct& runinfo);

#endif
