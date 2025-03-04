#ifndef _DIRTY_GET_DUST_THERMAL_EMISSION_
#define _DIRTY_GET_DUST_THERMAL_EMISSION_

#include "Constants.h"
#include "DirtyFailure.h"
#include "GrainModel.h"
#include "NumUtils.h"
#include "debug.h"
#include "geometry_def.h"
#include "runinfo_def.h"

//**********************************************************************
// external function definitions
extern int ComputeDustEmission(vector<float> &J, GrainModel &GrainModel, vector<vector<double>> &EmmittedEnergy,
                               bool &DoStochastic, bool &effective_grain_heating, float &_FailureSz, int &_FailureComp,
                               vector<float> &_transitionSz);

#endif
