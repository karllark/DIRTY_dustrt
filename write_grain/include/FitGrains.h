#ifndef _FITGRAINS_H_
#define _FITGRAINS_H_

#include <cctype>
#include <iostream>

#include "ConfigFile.h"
#include "Grain.h"
#include "StringManip.h"

using namespace std;

class FitGrains : protected Grain {

public:
  // Empty constructor for vectorization
  FitGrains();
  ~FitGrains() {};

  // Function defines - FitGrains exclusive members
  void MakeFitGrains(ConfigFile &ModelConfigFile, vector<float> &MasterWave);
  vector<float> getWave(void);
  int getNComp(void);

  // These functions are 'pass through' to the relevant Grain objects
  string getModelName(void);
  string ComponentName(int cmp);
  string OpticalFilename(int cmp);
  int nSize(int _cmp);
  float Size(int _cmp, int _szid);
  float Density(int _cmp);
  int nWave(int _cmp);
  vector<float> Wave(int _cmp);
  vector<float> Size(int _cmp);
  vector<float> CAbs(int _cmp, int _szid);
  vector<float> CSca(int _cmp, int _szid);
  vector<float> phFunc(int _cmp, int _szid);
  vector<float> CAbs(int _cmp, float _sz);
  vector<float> CSca(int _cmp, float _sz);
  vector<float> phFunc(int _cmp, float _sz);
  vector<float> wCAbs(int _cmp, int _wid);
  vector<float> CalTemp(int _cmp);
  vector<float> SpecificHeatCapacity(int _cmp);
  vector<float> HeatCapacity(int _szid, int _cmp);
  vector<float> SpecificEnthalpy(int _cmp);
  vector<float> Enthalpy(int _szid, int _cmp);

  // float getTau( float a_wave );
  // vector <float> getTau( void ) { return Tau; }

private:
  // Don't define/store a top level wavelength tabulation. Each component
  // can potentially store it's own wavelength scale.
  int nComp;
  string ModelName;
  vector<Grain> Component;
};

#endif
