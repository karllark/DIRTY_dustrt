#ifndef _EQTMP_H
#define _EQTMP_H

#include <iostream>
#include <vector> 

#include "Constants.h"
#include "NumUtils.h"

using namespace std; 
using namespace NumUtils;

double EqTemp(vector <double> wave, vector <double> J, vector <double> Q, 
	      double Tl=1.0, double Th=2500.0);
float EqTemp(vector <float> wave, vector <float> J, vector <float> Q, 
	     float Tl=1.0, float Th=2500.0);

#endif
