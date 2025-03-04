// Jun 3-4 New code to handle failure modes -KAM
#ifndef DIRTYFAILURE_H
#define DIRTYFAILURE_H

#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "DirtyFlags.h"

using namespace std;

class DirtyFailure {

private:
  string sFileName;                  // Output log name
  long nFailure;                     // Number of failures recorded
  vector<long> maxis;                // nth Failures m location
  vector<long> iaxis;                // nth Failures i location
  vector<long> jaxis;                // nth Failures j location
  vector<long> kaxis;                // nth Failures k location
  vector<int> FailureFlag;           // nth Failures integer ID
  vector<string> FailureDescription; // nth Failures string descriptor
  string UnknownFailureDescription;  // ...
  int FailureModeNotFound;           // Unkown failure code
  int nWave;                         // Number of wavelengths
  vector<double> AbsorbedEnergy;
  vector<vector<double>> RadiationField; // Radiation field in failed bin
  vector<double>::iterator iter1, iter2;

  vector<string> GrainModelName;
  vector<int> GrainComponent;
  vector<float> GrainSize;

  // Map return values to string descriptors
  map<int, string> FlagMap;
  map<int, string>::iterator mIter;

public:
  DirtyFailure(string &_FailureLogName, int &_nWave);
  ~DirtyFailure(void) {};

  // Given a failure flag, return a string description of the failure
  inline string getFailureDescription(int _flag) {
    mIter = FlagMap.find(_flag);
    if (mIter != FlagMap.end()) // We've found a defined failure mode
      return (mIter->second);
    else // Don't know what this is...
      return UnknownFailureDescription;
  }
  // Given a string description, return the corresponding flag.
  inline int getFailureFlag(string _description) {
    bool _found = false;
    mIter = FlagMap.begin();
    while (mIter != FlagMap.end() && !_found) {
      if (mIter->second != _description)
        mIter++;
      else
        _found = true;
    }
    if (mIter != FlagMap.end())
      return mIter->first;
    else
      return FailureModeNotFound;
  }

  void AddFailure(int _flag);
  inline void AddCellBook(int &_maxis, int &_iaxis, int &_jaxis, int &_kaxis) {
    maxis[nFailure - 1] = _maxis;
    iaxis[nFailure - 1] = _iaxis;
    jaxis[nFailure - 1] = _jaxis;
    kaxis[nFailure - 1] = _kaxis;
  }
  inline void AddGrainInfo(string _GrainModelName, float &_GrainSize,
                           int &_GrainComponent) {
    GrainModelName[nFailure - 1] = _GrainModelName;
    GrainSize[nFailure - 1] = _GrainSize;
    GrainComponent[nFailure - 1] = _GrainComponent;
  }
  inline void AddEnergyInfo(double &_AbsorbedEnergy,
                            vector<float> &_RadiationField) {
    AbsorbedEnergy[nFailure - 1] = _AbsorbedEnergy;
    copy(_RadiationField.begin(), _RadiationField.end(),
         RadiationField[nFailure - 1].begin());
  }

  // Dump the error log
  void WriteFailureLog(void);

}; // End of DirtyFailure class definition

#endif
