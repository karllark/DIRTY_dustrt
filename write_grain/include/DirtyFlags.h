#ifndef _DIRTYFLAGS_H
#define _DIRTYFLAGS_H

// A common place to hold return value string descriptors. Allows us to
// consistently map strings with integer return values.

#include <string>

namespace Flags
{

// Good return flag, all
const int FSUCCESS = 0;
const std::string SUCCESS = "Success";

// ComputeDustEmission() flags
const int FCDE_VECTOR_SIZE_MISMATCH = 1;
const std::string CDE_VECTOR_SIZE_MISMATCH = "ComputeDustEmission(): Vector size mismatch, J and wave.";

// EqTemp() flags
const int FEQ_VECTOR_SIZE_MISMATCH = 2;
const int FEQ_ZEROBOUND = 3;
const int FEQ_ZEROBOUND_LOLIM = 4;
const int FEQ_ZEROBOUND_HILIM = 5;
const std::string EQ_VECTOR_SIZE_MISMATCH = "EqTemp(): Vector size mismatch, J and Q.";
const std::string EQ_ZEROBOUND = "EqTemp(): Zero bound failure - upper and lower limits equal.";
const std::string EQ_ZEROBOUND_LOLIM = "EqTemp(): Zero bound failure - upper "
                                       "and lower limits equal at lower bound.";
const std::string EQ_ZEROBOUND_HILIM = "EqTemp(): Zero bound failure - upper "
                                       "and lower limits equal at upper bound.";

// StochasticHeating() flags.
const int FST_EXCEED_MAXBINS = 6;
const int FST_ZERO_PROBABILITY = 7;
const int FST_SMALL_PROBABILITY_HI = 8;
const int FST_SMALL_PROBABILITY_LO = 9;
const int FST_EXCEED_CONVERGENCE_LOOP_COUNT = 10;
const std::string ST_EXCEED_MAXBINS = "StochasticHeating(): Exceed maximum bin count.";
const std::string ST_ZERO_PROBABILITY = "StochasticHeating(): Transition probability zero everywhere.";
const std::string ST_SMALL_PROBABILITY_HI = "StochasticHeating(): Transition probablity less than cutoff everywhere, "
                                            "setting upper bound.";
const std::string ST_SMALL_PROBABILITY_LO = "StochasticHeating(): Transition probablity less than cutoff everywhere, "
                                            "setting lower bound.";
const std::string ST_EXCEED_CONVERGENCE_LOOP_COUNT =
    "StochasticHeating(): Exceeded maximum loop count in convergence loop.";
// --> ComputeTransitionMatrix
const int FST_CTM_ZERO_ENTH_DEN = 11;
const int FST_CTM_TRANSITIONID_EXCEEDS_WAVEID = 12;
const int FST_CTM_WAVE_ZERO = 13;
const std::string ST_CTM_ZERO_ENTH_DEN = "StochasticHeating()->ComputeTransitionMatrix(): Enthalpy denominator is "
                                         "zero.";
const std::string ST_CTM_TRANSITIONID_EXCEEDS_WAVEID =
    "StochasticHeating()->ComputeTransitionMatrix(): Found index of transition "
    "wavelength that exceeds maximum wave size.";
const std::string ST_CTM_WAVE_ZERO = "StochasticHeating()->ComputeTransitionMatrix(): wave[id]-wave[id-1] == 0 "
                                     "WTF?";
// --> SolveTransitionMatrix
const int FST_STM_ZERO_COOLING = 14;
const std::string ST_STM_ZERO_COOLING = "StochasticHeating()->SolveTransitionMatrix(): Found a zero cooling "
                                        "element; T-Matrix unsolvable.";

// get_dust_thermal_emission flags
const int FGDTE_ZERO_ABSORBED_ENERGY = 15;
const int FGDTE_POOR_ENERGY_CONSERVATION = 16;
const int FGDTE_NONFINITE_EMITTED_ENERGY_PRE = 17;
const int FGDTE_NONFINITE_EMITTED_ENERGY_PST = 18;
const std::string GDTE_ZERO_ABSORBED_ENERGY = "get_dust_thermal_emission(): No energy absorbed.";
const std::string GDTE_POOR_ENERGY_CONSERVATION =
    "get_dust_thermal_emission(): Energy conservation doesn't meet tolerance.";
const std::string GDTE_NONFINITE_EMITTED_ENERGY_PRE =
    "get_dust_thermal_emission(): Computed emitted energy is not finite, pre "
    "n_H.";
const std::string GDTE_NONFINITE_EMITTED_ENERGY_PST =
    "get_dust_thermal_emission(): Computed emitted energy is not finite, post "
    "n_H.";

} // namespace Flags

#endif
