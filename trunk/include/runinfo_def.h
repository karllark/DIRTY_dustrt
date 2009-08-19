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

using namespace std;

struct runinfo_struct {
  string param_filename; // name of file with all the parameters

  int verbose; // verbosity information (0=nothing)
  int do_ere_emission;  // do ERE emission (0=no)
  int do_dust_emission; // do dust thermal emission (0=no)
  int do_stochastic_dust_emission; // do stochastic_dust thermal emission (0=no)
  int do_emission_grain; // do emission by grain type (0=no)
                         // if set, then the emission is split by grain type 
                         // and equilibrium/non-equilibrium
  int dust_thermal_emission;  // set when doing thermal dust emission (used when do_emission_grain=1)
  int n_emission_grain_types;  // number of grain emission types
  float energy_conserve_target;  // energy conservation target

  int rt_check_converged;  // switch to see if we should check for rt convergence
  float rt_converge_target; // target for rt convergance (fractional for total scattered)

  int empir_dust;  // set to 1 if empirical dust properties being used
  int model_dust;  // set to 1 if dust grain model properties being used
  int n_waves;
  bool effective_grain_heating; // true if using effective grain in heating - defaults false.
  std::vector<float> wavelength;
  std::vector<float> albedo;
  std::vector<float> g;
  std::vector<float> tau_to_tau_ref;
  std::vector<float> tau_to_h;
  std::vector<float> ave_C_abs;

  // SED info
  std::string sed_type;
  std::vector<double> sed_lum;

  // emitted ere energy info
  int emitted_ere_energy_grid_initialized;
  std::vector< std::vector<double> > emitted_ere_lum;

  // emitted energy info (by dust grain component)
  std::vector< std::vector<double> > emitted_lum;

  // absorbed energy
  std::vector<double> absorbed_energy;  // aborbed energy by wavelength
  double total_absorbed_energy;   // total absorbed energy (computed in get_thermal_dust_emission.cpp)

  // emitted energy
  double total_emitted_energy;   // total emitted energy

  // cell info
  long num_cells_enough;  // cells with enough energy for dust emission
  long num_cells_not_enough; // cells with no emission
  long num_cells_zero; // cells with zero energy
  long num_cells_too_few_waves; // cells with too few wavelengths for dust emission

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
  std::vector< std::vector<double> > out_sed_lum;  // allow for direct/scattered lum
  std::vector< std::vector<double> > out_sed_lum_unc;  // allow for direct/scattered lum uncertainties
/*   std::vector< std::vector<float> > out_emitted_lum;  // allow to be by grain/emission type */
};

#endif
