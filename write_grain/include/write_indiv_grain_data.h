// ======================================================================
//   Header file for write_indiv_grain_data
//
// 2014 Oct/KDG - started
// ======================================================================

#include <cmath>
#include <fstream>
#include <iostream>

#include "ConfigFile.h"
#include "DataFile.h"
#include "DirtyFlags.h"
#include "FitGrains.h"
#include "ISRF.h"
#include "check_input_param.h"

//**********************************************************************
// external function definitions
extern int ComputeDustEmission_igrain(vector<float> &J, vector<float> &CAbs, vector<float> &Wave,
                                      vector<float> &Enthalpy, vector<float> &CalorimetryTGrid,
                                      vector<double> &EquilibriumEmission, vector<double> &StochasticEmission,
                                      float &TauHeating, float &TauCooling, float &TauScaling,
                                      bool &StochasticallyHeated, string ModelName, string ComponentName,
                                      int ComponentIndex, float Size, int SizeIndex, string RadiationFieldType,
                                      float RadiationFieldScale, float RadiationFieldTemperature);
