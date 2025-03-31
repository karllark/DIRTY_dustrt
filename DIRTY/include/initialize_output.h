// ======================================================================
//   Header file for initialize output procedure.
// Include files and function definitions.
//
// 2005 May/KDG - written
// 2008 Jun/KDG - added computer matrix def
// ======================================================================
#ifndef _DIRTY_INITIALIZE_OUTPUT_
#define _DIRTY_INITIALIZE_OUTPUT_

#include <cmath>
#include <iostream>

#include "debug.h"
#include "geometry_def.h"
#include "output_def.h"

//**********************************************************************
// external function definitions

extern void compute_observer_trans_matrix (output_struct &output, geometry_struct &geometry, int i);

#endif
