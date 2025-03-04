// Jun 3-4 2009.  New code to handle failure modes. -KAM
#include "DirtyFailure.h"

// Constructor.
DirtyFailure::DirtyFailure(string &_FailureLogName, int &_nWave)
{
    nFailure = 0;
    sFileName = _FailureLogName;
    FailureModeNotFound = -99;
    nWave = _nWave;

    // Initialize map using DirtyFlags definitions
    FlagMap[Flags::FSUCCESS] = Flags::SUCCESS;
    FlagMap[Flags::FCDE_VECTOR_SIZE_MISMATCH] = Flags::CDE_VECTOR_SIZE_MISMATCH;
    FlagMap[Flags::FEQ_VECTOR_SIZE_MISMATCH] = Flags::EQ_VECTOR_SIZE_MISMATCH;
    FlagMap[Flags::FEQ_ZEROBOUND] = Flags::EQ_ZEROBOUND;
    FlagMap[Flags::FEQ_ZEROBOUND_LOLIM] = Flags::EQ_ZEROBOUND_LOLIM;
    FlagMap[Flags::FEQ_ZEROBOUND_HILIM] = Flags::EQ_ZEROBOUND_HILIM;
    FlagMap[Flags::FST_EXCEED_MAXBINS] = Flags::ST_EXCEED_MAXBINS;
    FlagMap[Flags::FST_ZERO_PROBABILITY] = Flags::ST_ZERO_PROBABILITY;
    FlagMap[Flags::FST_SMALL_PROBABILITY_HI] = Flags::ST_SMALL_PROBABILITY_HI;
    FlagMap[Flags::FST_SMALL_PROBABILITY_LO] = Flags::ST_SMALL_PROBABILITY_LO;
    FlagMap[Flags::FST_EXCEED_CONVERGENCE_LOOP_COUNT] = Flags::ST_EXCEED_CONVERGENCE_LOOP_COUNT;
    FlagMap[Flags::FST_STM_ZERO_COOLING] = Flags::ST_STM_ZERO_COOLING;
    FlagMap[Flags::FST_CTM_ZERO_ENTH_DEN] = Flags::ST_CTM_ZERO_ENTH_DEN;
    FlagMap[Flags::FST_CTM_TRANSITIONID_EXCEEDS_WAVEID] = Flags::ST_CTM_TRANSITIONID_EXCEEDS_WAVEID;
    FlagMap[Flags::FST_CTM_WAVE_ZERO] = Flags::ST_CTM_WAVE_ZERO;
    FlagMap[Flags::FGDTE_ZERO_ABSORBED_ENERGY] = Flags::GDTE_ZERO_ABSORBED_ENERGY;
    FlagMap[Flags::FGDTE_POOR_ENERGY_CONSERVATION] = Flags::GDTE_POOR_ENERGY_CONSERVATION;
    FlagMap[Flags::FGDTE_NONFINITE_EMITTED_ENERGY_PRE] = Flags::GDTE_NONFINITE_EMITTED_ENERGY_PRE;
    FlagMap[Flags::FGDTE_NONFINITE_EMITTED_ENERGY_PST] = Flags::GDTE_NONFINITE_EMITTED_ENERGY_PST;

    // Make sure our FailureModeNotFound id is not in the map.
    while (FlagMap.find(FailureModeNotFound) != FlagMap.end())
        FailureModeNotFound--;
}

// Add an error
void DirtyFailure::AddFailure(int _flag)
{
    nFailure++;
    FailureFlag.push_back(_flag);
    FailureDescription.push_back(getFailureDescription(_flag));

    GrainModelName.resize(nFailure);
    GrainSize.resize(nFailure);
    GrainComponent.resize(nFailure);

    maxis.resize(nFailure);
    iaxis.resize(nFailure);
    jaxis.resize(nFailure);
    kaxis.resize(nFailure);

    AbsorbedEnergy.resize(nFailure);
    RadiationField.resize(nFailure);
    RadiationField[nFailure - 1].resize(nWave);
}

void DirtyFailure::WriteFailureLog(void)
{
    // Open the output file.
    ofstream failurelog(sFileName.c_str());
    if (failurelog)
    {
        // Write some header info....
        for (long i = 0; i < nFailure; ++i)
        {
            failurelog << "**********************************************************"
                          "************"
                       << endl;
            failurelog << "Begin failure " << i + 1 << " of " << nFailure << endl;
            failurelog << "Failed on cell " << iaxis[i] << "," << jaxis[i] << "," << kaxis[i] << " of " << maxis[i]
                       << endl;
            failurelog << "--> Failure info: " << endl;
            failurelog << "    Failure Flag - " << FailureFlag[i] << endl;
            failurelog << "    Failure Description - " << FailureDescription[i] << endl;
            failurelog << "--> Grain info: " << endl;
            failurelog << "    Grain Model Name - " << GrainModelName[i] << endl;
            failurelog << "    Size - " << GrainSize[i] << " cm  of component " << GrainComponent[i] << endl;
            failurelog << "--> Energy/Radiation Field info: " << endl;
            failurelog << "    Total energy absorbed in bin-  " << AbsorbedEnergy[i] << endl;
            failurelog << "    Radiation field" << endl;
            failurelog << "    ---------------" << endl;
            iter1 = RadiationField[i].begin();
            iter2 = RadiationField[i].end();
            for (; iter1 != iter2; iter1++)
                failurelog << "    " << *iter1 << endl;
            failurelog << "End failure " << i + 1 << " of " << nFailure << endl;
            failurelog << "**********************************************************"
                          "************"
                       << endl;
        }
        failurelog.close();
    } // else... I don't know...
}
