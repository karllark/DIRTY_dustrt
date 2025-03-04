/**************************************************************
 * ISRF -- ISRF Class, header file.
 *
 * History:
 *     Written by:   Putrid, Jan 2007, adapted from genISRF()
 *                   Contact: misselt@as.arizona.edu
 *
 **************************************************************/

#ifndef _ISRF_H_
#define _ISRF_H_

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "Constants.h"
#include "NumUtils.h"

using namespace std;

class ISRF
{
private:
  int nWave;
  float XMMP;
  vector<float> theISRF;
  vector<float> wave;

public:
  ISRF (vector<float> in_wave, float in_XMMP);
  ~ISRF () {};

  inline vector<float>
  getISRF (void)
  {
    return theISRF;
  }
  inline vector<float>
  getISRFwave (void)
  {
    return wave;
  }
  inline int
  getN (void)
  {
    return nWave;
  }
  inline float
  getScale (void)
  {
    return XMMP;
  }
};

#endif
