#include "check_absorbed_energy_grid.h"
// #define DEBUG_SAEG
//  2008 sept 2 - added implicit cast for NumUtils::integrate call.
void
check_absorbed_energy_grid (geometry_struct &geometry, runinfo_struct &runinfo)

{
  // loop through the dust density grid and convert to radiation field density
  int i, j, k, m = 0;
  uint x;
  vector<double> abs_energy;
  vector<double> tmp_wave;
  abs_energy.resize (runinfo.wavelength.size (), 0.0);
  tmp_wave.resize (runinfo.wavelength.size (), 0.0);
  for (x = 0; x < runinfo.wavelength.size (); x++)
    {
      // loop over all the defined grids
      double total_abs_energy = 0.0;
      total_abs_energy = 0.0;

      for (m = 0; m < int (geometry.grids.size ()); m++)
        {
          // loop of the cells in this grid
          for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
            for (j = 0; j < geometry.grids[m].index_dim[1]; j++)
              for (i = 0; i < geometry.grids[m].index_dim[0]; i++)
                {
                  // 	    if
                  // (geometry.grids[m].grid(i,j,k).absorbed_energy[x] > 0.) {
                  // 	      cout <<
                  // geometry.grids[m].grid(i,j,k).absorbed_energy[x]
                  // << endl; 	      exit(8);
                  // 	    }
                  //	    abs_energy[x] +=
                  // geometry.grids[m].grid(i,j,k).absorbed_energy[x];
                  if (geometry.grids[m].grid (i, j, k).dust_tau_per_pc > 0.0)
                    {
                      abs_energy[x] += geometry.grids[m].grid (i, j, k).absorbed_energy[x]
                                       * Constant::FPI * runinfo.ave_C_abs[x]
                                       * geometry.grids[m].grid (i, j, k).num_H;
                      if (geometry.grids[m].grid (i, j, k).num_H <= 0.0)
                        {
                          cout << geometry.grids[m].grid (i, j, k).num_H << endl;
                          cout << m << " " << i << " " << j << " " << k << endl;
                          cout << "x = " << x << endl;
                          exit (8);
                        }
                    }
                  // 	    if (x == 1) cout <<
                  // geometry.grids[m].grid(i,j,k).num_H << "
                  // ";
                }
        }
      tmp_wave[x] = runinfo.wavelength[x];
      total_abs_energy += abs_energy[x];
      cout << "caeg; x = " << x << " ";
      cout << runinfo.wavelength[x] * Constant::CM_UM << " ";
      cout << runinfo.sed_lum[x] << " ";
      cout << abs_energy[x] << " ";
      cout << runinfo.out_sed_lum[0][x] << " ";
      cout << runinfo.out_sed_lum[1][x] << " ";
      cout << runinfo.sed_lum[x] * (1.0 - runinfo.out_sed_lum[0][x] - runinfo.out_sed_lum[1][x])
           << " ";
      cout << abs_energy[x]
                  / (runinfo.sed_lum[x]
                     * (1.0 - runinfo.out_sed_lum[0][x] - runinfo.out_sed_lum[1][x]))
           << " ";
      cout << endl;
    }

  double tot_abs_energy = NumUtils::integrate<double> (tmp_wave, abs_energy);
  cout << "total_abs_energy = " << tot_abs_energy << endl;

  exit (8);
}
