#ifndef _DIRTY_GET_DUST_THERMAL_EMISSION_
#define _DIRTY_GET_DUST_THERMAL_EMISSION_

#include "geometry_def.h"
#include "runinfo_def.h"
#include "NumUtils.h"
#include "GrainModel.h"
#include "Constants.h"
#include "debug.h"
#include "DirtyFailure.h"

//**********************************************************************
// external function definitions
extern int ComputeDustEmission (vector <float> & J, GrainModel & GrainModel, 
				vector <vector<double> > & EmmittedEnergy, 
				bool & DoStochastic, float & _FailureSz, int & _FailureComp);

#endif
