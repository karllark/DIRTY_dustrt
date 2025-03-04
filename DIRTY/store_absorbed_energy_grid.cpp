// ======================================================================
//   Procedure to store absorbed energy grid for the current wavelength.
// Should allow either memory (speed) or disk (space) to be used.
//
// 2007 Jun/KDG - written
// 2007 Sep/KDG - added conversion of absorbed energy to radiation field
//                density
// ======================================================================
// notes:
// to convert the absorbed energy to radiation field density
//
// first get the absorbed energy per H atom
//   E_abs = E_model*(input_luminosity/wavelength[in cm])/#H
//   #H = n_H*V = (N_H/L)*V = tau_mod*(V/L)*(tau/N_H)^-1
// where
//   tau/N_H is from dust grain model
//
//   J(lambda) = E_abs(lambda)/[4*pi*ave_C_abs(lambda)]
// where
//   ave_C_abs(lambda) is from the model and is per H atom
//
// compute/keep tau_mod*(V/L)
//   add the tau/N_H division at the final stage
//   this avoids round off error (hopefully)
// ======================================================================
#include "store_absorbed_energy_grid.h"
#include "compat.h"
// #define DEBUG_SAEG

void store_absorbed_energy_grid(geometry_struct &geometry,
                                runinfo_struct &runinfo, output_struct &output,
                                int index, int doing_emission)

{
  // total energy absorbed in photons
  double total_energy_absorbed_photons = 0.0;
  // total energy absorbed at this wavelength
  runinfo.absorbed_energy[geometry.abs_energy_wave_index] = 0.0;

  // constant for num_H calculation
  double num_H_const = 1. / (runinfo.tau_to_h[index] * (Constant::PC_CM));

  // total H mass (so we can get the dust mass assuming a gas-to-dust ratio)
  double total_H_mass = 0.0;

  // loop through the dust density grid and convert to radiation field density
  int i, j, k, m = 0;
  // loop over all the defined grids
  double vol = 1.0;
  double size_x, size_y, size_z;
  for (m = 0; m < int(geometry.grids.size()); m++) {

    //     for (i = 0; i < 3; i++) {
    //       vol *= (Constant::PC_CM)*(geometry.grids[m].positions[i][1] -
    //       geometry.grids[m].positions[i][0]); cout <<
    //       geometry.grids[m].positions[i][0] << " "; cout <<
    //       geometry.grids[m].positions[i][1] << " "; cout << vol << endl;
    //     }
    //     cout << vol << endl;

#ifdef DEBUG_SAEG
    cout << "grid num = " << m << endl;
#endif

    // loop of the cells in this grid
    for (k = 0; k < geometry.grids[m].index_dim[2]; k++) {
      size_z = (Constant::PC_CM) * (geometry.grids[m].positions[2][k + 1] -
                                    geometry.grids[m].positions[2][k]);
      // cout << "z " << size_z << " ";
      // cout << geometry.grids[m].positions[2][k+1] << " ";
      // cout << geometry.grids[m].positions[2][k] << endl;
      for (j = 0; j < geometry.grids[m].index_dim[1]; j++) {
        size_y = (Constant::PC_CM) * (geometry.grids[m].positions[1][j + 1] -
                                      geometry.grids[m].positions[1][j]);
        for (i = 0; i < geometry.grids[m].index_dim[0]; i++) {
          size_x = (Constant::PC_CM) * (geometry.grids[m].positions[0][i + 1] -
                                        geometry.grids[m].positions[0][i]);

          // determine the volume in cm^3
          vol = size_x * size_y * size_z;
          // if (vol <= 0.0) {
          //   cout << size_x << " ";
          //   cout << size_y << " ";
          //   cout << size_z << " ";
          // cout << "vol ; " << vol;
          //   exit(0);
          // }

          // only do the computation the dust optical depth is positive
          // *and* when there is a nonzero amount of energy absorbed
          //   easy way to avoid subgrids and the edges (no dust)
          if (geometry.grids[m].grid(i, j, k).dust_tau_per_pc > 0.0) {
            // # of H atoms in cell
            geometry.grids[m].grid(i, j, k).num_H =
                num_H_const * vol *
                geometry.grids[m].grid(i, j, k).dust_tau_per_pc *
                geometry.tau_to_tau_ref;
            // cout << "; num_h = " << geometry.grids[m].grid(i,j,k).num_H << ";
            // ";
          }

          if ((geometry.grids[m].grid(i, j, k).dust_tau_per_pc > 0.0) &&
              (geometry.grids[m]
                   .grid(i, j, k)
                   .absorbed_energy[geometry.abs_energy_wave_index] > 0.0)) {

#ifdef DEBUG_SAEG
            cout << "dust tau/pc = "
                 << geometry.grids[m].grid(i, j, k).dust_tau_per_pc << endl;
            cout << "# H [total] = " << geometry.grids[m].grid(i, j, k).num_H
                 << endl;
            cout << "Raw E_abs = "
                 << geometry.grids[m]
                        .grid(i, j, k)
                        .absorbed_energy[geometry.abs_energy_wave_index]
                 << endl;
            cout << "sed_lum = " << runinfo.sed_lum[index] << endl;
            cout << "total # photons = " << output.outputs[0].total_num_photons
                 << endl;
            cout << "ave C_abs = " << runinfo.ave_C_abs[index] << endl;
#endif
            total_H_mass +=
                geometry.grids[m].grid(i, j, k).num_H * Constant::HMASS_CGS;

            // convert the absorbed energy to radiation field density
            double j_temp =
                geometry.grids[m]
                    .grid(i, j, k)
                    .absorbed_energy[geometry.abs_energy_wave_index];
            double j_temp_x2 =
                geometry.grids[m]
                    .grid(i, j, k)
                    .absorbed_energy_x2[geometry.abs_energy_wave_index];
            total_energy_absorbed_photons += j_temp;
            double flux_mult_factor = 0.0;
            if (doing_emission)
              if (runinfo.dust_thermal_emission)
                flux_mult_factor = (runinfo.emitted_lum[0][index] /
                                    output.outputs[0].total_num_photons);
              else
                flux_mult_factor = (runinfo.emitted_ere_lum[0][index] /
                                    output.outputs[0].total_num_photons);
            else
              flux_mult_factor = (runinfo.sed_lum[index] /
                                  output.outputs[0].total_num_photons);

            j_temp *= flux_mult_factor;
            j_temp_x2 *= flux_mult_factor * flux_mult_factor;

            flux_mult_factor = (geometry.grids[m].grid(i, j, k).num_H * 4.0 *
                                (Constant::PI)*runinfo.ave_C_abs[index]);
            j_temp /= flux_mult_factor;
            j_temp_x2 /= (flux_mult_factor * flux_mult_factor);

            // if (j_temp > 0) {
            //   cout << "emitted_lum = " << runinfo.sed_lum[index] << "; ";
            //   cout << "num_H = " << geometry.grids[m].grid(i,j,k).num_H << ";
            //   "; cout << "C_abs = " << runinfo.ave_C_abs[index] << "; "; cout
            //   << "j_temp" << " "; cout << j_temp << " "; cout << j_temp_x2 <<
            //   endl;
            // }

            // check for roundoff error before converting to double
            if (j_temp < 1e-38) {
#ifdef DEBUG_SAEG
              cout << "roundoff error warning in store_absorbed_energy_grid."
                   << endl;
              cout << "j_temp = " << j_temp << endl;
              cout << "float(j_temp) = " << float(j_temp) << endl;
#endif
              //	      exit(8);
            } else if (!finite(j_temp)) {
              cout << "j_temp is non finite" << endl;
              cout << "j_temp = " << j_temp << endl;
              exit(8);
            }
            // 	    cout <<
            // geometry.grids[m].grid(i,j,k).absorbed_energy[geometry.abs_energy_wave_index]
            // << " "; 	    cout << j_temp << endl; 	    exit(8);
            geometry.grids[m]
                .grid(i, j, k)
                .absorbed_energy[geometry.abs_energy_wave_index] =
                float(j_temp);
            geometry.grids[m]
                .grid(i, j, k)
                .absorbed_energy_x2[geometry.abs_energy_wave_index] =
                float(j_temp_x2);
            // reset the number of photon to last contribute to this cell
            geometry.grids[m].grid(i, j, k).last_photon_number = -1;
            geometry.grids[m].grid(i, j, k).last_photon_absorbed_energy = 0.0;

#ifdef DEBUG_SAEG
            cout << "J_temp = " << j_temp << endl;
            cout << "float(J_temp) = " << float(j_temp) << endl;
            cout << "J(lambda) = "
                 << geometry.grids[m]
                        .grid(i, j, k)
                        .absorbed_energy[geometry.abs_energy_wave_index]
                 << endl;
            // cout << "prev_J(lambda) = " <<
            // geometry.grids[m].grid(i,j,k).save_radiation_field_density <<
            // endl;
#endif
          }

          // add back the previously saved radiation field density
          if (doing_emission) {
            geometry.grids[m]
                .grid(i, j, k)
                .absorbed_energy[geometry.abs_energy_wave_index] +=
                geometry.grids[m].grid(i, j, k).save_radiation_field_density
                    [geometry.abs_energy_wave_index];
            geometry.grids[m]
                .grid(i, j, k)
                .absorbed_energy_x2[geometry.abs_energy_wave_index] +=
                geometry.grids[m].grid(i, j, k).save_radiation_field_density_x2
                    [geometry.abs_energy_wave_index];
            geometry.grids[m]
                .grid(i, j, k)
                .absorbed_energy_num_photons[geometry.abs_energy_wave_index] +=
                geometry.grids[m]
                    .grid(i, j, k)
                    .save_radiation_field_density_num_photons
                        [geometry.abs_energy_wave_index];
          }

          if (geometry.grids[m].grid(i, j, k).num_H < 0.0)
            exit(8);

          // sum the absorbed energy for this wavelength
          runinfo.absorbed_energy[geometry.abs_energy_wave_index] +=
              geometry.grids[m]
                  .grid(i, j, k)
                  .absorbed_energy[geometry.abs_energy_wave_index] *
              geometry.grids[m].grid(i, j, k).num_H;

#ifdef DEBUG_SAEG
          cout << "final J(lambda) = "
               << geometry.grids[m]
                      .grid(i, j, k)
                      .absorbed_energy[geometry.abs_energy_wave_index]
               << endl;
#endif
          // 	    int ans;
          // 	    cin >> ans;
        }
      }
    }
  }

  //   cout << endl << "index = " << geometry.abs_energy_wave_index << " ";
  //   cout << " " << "index2 = " << index << " ";
  //   cout << "total photons = " << output.outputs[0].total_num_photons;
  //   cout << endl << "total photons absorbed = " <<
  //   total_energy_absorbed_photons << endl; cout << "ratio = " <<
  //   total_energy_absorbed_photons/output.outputs[0].total_num_photons <<
  //   endl; cout << "energy = " <<
  //   (total_energy_absorbed_photons/output.outputs[0].total_num_photons)*runinfo.sed_lum[index]
  //   << endl; exit(8);

  geometry.total_h_mass = total_H_mass;
  //   cout << "total H mass [g] = " << total_H_mass << endl;

  if (geometry.abs_energy_storage_type == 1) {
    cout << "no code written for storage_absorbed_energy_grid for disk storage"
         << endl;
    exit(8);
  }
  // nothing needed for geometry.abs_energy_storage_type == 0 (memory)

  //   cout << "total abs energy = " <<
  //   runinfo.absorbed_energy[geometry.abs_energy_wave_index] << endl;
}
