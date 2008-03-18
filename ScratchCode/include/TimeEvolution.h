#ifndef _TIMEEVOLUTION_H
#define _TIMEEVOLUTION_H

#include <fstream>
#include <iostream>
#include <vector> 

#include "Constants.h"
#include "HeatUtils.h"
#include "NumUtils.h"
#include "mt_ran.h"

using namespace std;
//using namespace NumUtils; 

void TimeEvolution(vector <float>& wave, vector <float>& J, vector <float>& CAbs, 
		   vector <float>& Temperature, vector <float>& SpecificHeatCapacity, 
		   vector <float>& SpecificEnthalpy, float density, float size,
		   vector <float>& Time, vector <float>& TempHistory);

#endif 
