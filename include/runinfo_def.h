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
  int dust_thermal_emission;  // set when doing thermal dust emission (used when do_emission_grain=1)
  int n_emission_grain_types;  // number of grain emission types

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

  // output images
  int do_image_output;  // set if output of images at each wavelength are desired

  // output global values
  int do_global_output;  // set if global, multiwavelength output is desired
  int out_sed_lum_offset; // offset for vector to all for ERE/Dust direct/scattered to be saved as well as stellar direct/scattered
  std::vector< std::vector<float> > out_sed_lum;  // allow for direct/scattered lum
  std::vector< std::vector<float> > out_sed_lum_unc;  // allow for direct/scattered lum uncertainties
/*   std::vector< std::vector<float> > out_emitted_lum;  // allow to be by grain/emission type */
};

#endif
