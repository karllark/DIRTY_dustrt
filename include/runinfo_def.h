// ======================================================================
//   Header file with runinfo definitions.
//
// 2006 Nov/KDG - written
// ======================================================================

#ifndef _DIRTY_RUNINFO_DEF_
#define _DIRTY_RUNINFO_DEF_

#include <vector>
#include <string>
//#include "vect_utils.h"

// structure with multiwavlength information for the 
// cases where empirical scattering parameters are being
// used
//
// lots more will go here as things are fleshed out

struct runinfo_struct {
  int verbose; // verbosity information (0=nothing)
  int do_dust_emission; // do dust thermal emission (0=no)
  int do_emission_grain; // do emission by grain type (0=no)
                         // if set, then the emission is split by grain type 
                         // and equilibrium/non-equilibrium

  int empir_dust;  // set to 1 if empirical dust properties being used
  int model_dust;  // set to 1 if dust grain model properties being used
  int n_waves;
  std::vector<float> wavelength;
  std::vector<float> albedo;
  std::vector<float> g;
  std::vector<float> tau_to_tau_ref;
  std::vector<float> tau_to_h;
  std::vector<float> ave_C_abs;

  // SED info
  std::string sed_type;
  std::vector<double> sed_lum;

  // emitted energy info (by dust grain component)
  std::vector< std::vector<double> > emitted_lum;

  // ERE model info
  float ere_efficiency;  // efficiency of conversion of input to output photons
  float ere_excitation_wavelength;   // wavelength below which photons excite ERE
  float ere_peak_wavelength;  // peak wavelength of ERE emission
  float ere_fwhm; // FWHM of ERE emission

  // output global values
  std::vector< std::vector<float> > out_sed_lum;  // allow for direct/scattered lum
  std::vector< std::vector<float> > out_emitted_lum;  // allow to be by grain/emission type
};

#endif
