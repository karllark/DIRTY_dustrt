// Jun 3-4 2009: Modified to be DirtyFailure aware.  Dust emission failures will pass failure
// information
//               back up the call stack.
#include "GrainModel.h"
#include "HeatUtils.h"
#include "NumUtils.h"

#include "DirtyFailure.h"

// #include "DirtyFlags.h"
extern int StochasticHeating_bttw(vector<float> &wave, vector<float> &J, vector<float> &cabs, vector<float> &Temp,
                                  vector<float> &Enthalpy, float &EAbs, float &TMin, float &TMax,
                                  vector<double> &Emission, vector<double> &ProbabilityDistribution,
                                  vector<float> &temp);

// int ComputeDustEmission (vector <float> & J, GrainModel & GrainModel,
// 			 vector <vector<double> > & EmmittedEnergy, bool & DoStochastic,
// 			 bool & UseEffective_Grain,
// 			 float & _FailureSz, int & _FailureComp, vector <float> & _transitionSz)

int ComputeDustEmission_igrain(vector<float> &J, vector<float> &CAbs, vector<float> &Wave, vector<float> &Enthalpy,
                               vector<float> &CalorimetryTGrid, vector<double> &EquilibriumEmission,
                               vector<double> &StochasticEmission, float &TauHeating, float &TauCooling,
                               float &TauScaling, bool &StochasticallyHeated, string ModelName, string ComponentName,
                               int ComponentIndex, float Size, int SizeIndex, string RadiationFieldType,
                               float RadiationFieldScale, float RadiationFieldTemperature)

{

    int status;
    uint _nw = Wave.size();
    float EAbs, EquilibriumGrainTemperature;
    float _mintemp, _maxtemp;
    float MeanPhotonEnergy;
    float Tu;
    vector<float> _cJprod(_nw);
    vector<double> ProbabilityDistribution;
    vector<float> temp;
    stringstream ss;
    string pdistFilename;
    ofstream pdistStream;

    _mintemp = 1.0;
    _maxtemp = 3000.0;
    StochasticallyHeated = false;

    if (J.size() != _nw || J.size() != CAbs.size())
    {
        return Flags::FCDE_VECTOR_SIZE_MISMATCH;
    }

    for (uint _wv = 0; _wv < _nw; ++_wv)
    {
        _cJprod[_wv] = CAbs[_wv] * J[_wv];
        StochasticEmission[_wv] = 0.0;
        EquilibriumEmission[_wv] = 0.0;
    }

    status = EqTemp(Wave, J, CAbs, EAbs, EquilibriumGrainTemperature, _mintemp, _maxtemp);

    // Get the heating and cooling timescales
    TauHeating = HeatUtils::getTauAbs(J, CAbs, Wave, MeanPhotonEnergy);
    Tu = HeatUtils::getTmean(CalorimetryTGrid, Enthalpy, MeanPhotonEnergy);
    TauCooling = HeatUtils::getTauRad(CAbs, Wave, MeanPhotonEnergy, Tu);

    if (TauScaling * TauHeating > TauCooling)
        StochasticallyHeated = true;

    EquilibriumEmission = NumUtils::prod_bbodyCGS<double>(Wave, EquilibriumGrainTemperature, CAbs);

    if (StochasticallyHeated)
    {
        // cout << "Should do stochastic" << endl;
        status = StochasticHeating_bttw(Wave, _cJprod, CAbs, CalorimetryTGrid, Enthalpy, EAbs, _mintemp, _maxtemp,
                                        StochasticEmission, ProbabilityDistribution, temp);
        if (status != Flags::FSUCCESS)
            cout << " Huh..." << status << endl;

        ss << ModelName << "_c_" << ComponentIndex + 1 << "_" << ComponentName << "_s_" << SizeIndex + 1 << "_"
           << RadiationFieldType << "_" << RadiationFieldScale << "_PDIST.dat";
        pdistFilename = ss.str();
        pdistStream.open(pdistFilename.c_str());
        pdistStream << "# " << ModelName << endl;
        pdistStream << "# Component, Size Index, Size: " << ComponentName << " " << SizeIndex << " " << Size << endl;
        pdistStream << "# Tau Heat, Cool, Ratio, Termination Ratio: " << TauHeating << " " << TauCooling << " "
                    << TauCooling / TauHeating << " " << TauScaling << endl;
        pdistStream << "# Radition Field: " << RadiationFieldType << " " << RadiationFieldScale << " "
                    << RadiationFieldTemperature << endl;
        pdistStream << "# Temperature  PDist" << endl;
        for (uint i = 0; i < ProbabilityDistribution.size(); ++i)
            pdistStream << temp[i] << " " << ProbabilityDistribution[i] << endl;
    }

    if (status != Flags::FSUCCESS)
        return status; // we've had a heating failure.

    return Flags::FSUCCESS;
}
