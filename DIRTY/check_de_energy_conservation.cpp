// ======================================================================
//   Procedure to check the energy conservation of the dust emission.
// If better than a user input limit, then stop the iteration.  This
// really will check that the total dust emission is no longer changing,
// but it has the same effect as energy conservation.
//
// 2008 Jun/KDG - written
// 2008 Sept 2 - added implicit casts for NumUtils::integrate call -KAM
// ======================================================================
#include "check_de_energy_conservation.h"
#define DEBUG_CDEC

void check_de_energy_conservation(runinfo_struct &runinfo, int &iter_done)

{
  // now convert the fluxes from ratios to luminosities
  int out_sed_lum_offset = 0;
  vector<double> out_wavelength;

  vector<double> de_direct;
  vector<double> de_scattered;
  int i = 0;
  for (i = 0; i < int(runinfo.wavelength.size()); i++) {
    // output wavelength (do in microns)
    out_wavelength.push_back(runinfo.wavelength[i] * Constant::CM_UM);

    // offset (adjust if ERE also)
    out_sed_lum_offset = 2;

    // de direct
    de_direct.push_back(runinfo.out_sed_lum[out_sed_lum_offset][i] *
                        runinfo.emitted_lum[0][i] / Constant::CM_UM);
    // de scattered
    de_scattered.push_back(runinfo.out_sed_lum[1 + out_sed_lum_offset][i] *
                           runinfo.emitted_lum[0][i] / Constant::CM_UM);
  }

  // determine the global energy absorbed and emitted
  double total_de_direct_energy =
      NumUtils::integrate<double>(out_wavelength, de_direct);
  double total_de_scat_energy =
      NumUtils::integrate<double>(out_wavelength, de_scattered);

  double new_total_emit = total_de_direct_energy + total_de_scat_energy;

  // #ifdef debug_CDEC
  cout << "Absorbed(rt) = " << runinfo.total_absorbed_energy << endl;
  cout << "Emitted(de) = " << new_total_emit << endl;
  cout << "Emitted(de) previous = " << runinfo.total_emitted_energy << endl;
  // cout << "Energy Conservation = " <<
  // new_total_emit/runinfo.total_absorbed_energy << endl; #endif

  float emitted_ratio = runinfo.total_emitted_energy / new_total_emit;
  cout << "previous/new emitted(de) = " << emitted_ratio << endl;

  //   float energy_ratio = total_emit/runinfo.total_absorbed_energy;
  //   if (energy_ratio > 1.0) {
  //     cout << "more energy emitted than absorbed...check this out." << endl;
  //     cout << "energy ratio (emit/abs) = " << energy_ratio << endl;
  //     exit(8);
  //   }

  cout << "fabs(1.0 - emitted_ratio) = " << fabs(1.0 - emitted_ratio) << endl;
  cout << "target = " << runinfo.energy_conserve_target << endl;
  if (runinfo.energy_conserve_target >= 0.99)
    iter_done = 1;
  else if (fabs(1.0 - emitted_ratio) <= runinfo.energy_conserve_target)
    iter_done = 1;
  else
    cout << "Another interation required (ratio = " << emitted_ratio << ")"
         << endl;

  runinfo.total_emitted_energy = new_total_emit;

  // Compare dust emitted (escaped) energy to stellar absorbed energy
  if (!iter_done && (runinfo.energy_conserve_target2 > 0)) {
    double energy_conserve_ratio =
        new_total_emit / runinfo.total_absorbed_stellar_energy;
    cout << "energy_conserve_target2 = " << runinfo.energy_conserve_target2
         << endl;
    cout << "energy_conserve_ratio = " << energy_conserve_ratio << endl;
    if (fabs(1.0 - energy_conserve_ratio) <= runinfo.energy_conserve_target2) {
      iter_done = 1;
      cout << "energy_conserve_target2 is met" << endl;
    }
  }
}
