// ======================================================================
//   Procedure to setup the absorbed energy grid for the current wavelength.
// Should allow either memory (speed) or disk (space) to be used.
//
// 2007 Jun/KDG - written
// ======================================================================
#include "setup_absorbed_energy_grid.h"

void
setup_absorbed_energy_grid (geometry_struct &geometry, runinfo_struct &runinfo, int doing_emission)

{
    if ((!geometry.abs_energy_grid_initialized) && (!doing_emission))
        runinfo.absorbed_energy.resize (runinfo.n_waves, 0.0);

    int i, j, k, m, n = 0;
    // loop over all the defined grids
    for (m = 0; m < int (geometry.grids.size ()); m++)
        {
            // loop of the cells in this grid
            for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
                for (j = 0; j < geometry.grids[m].index_dim[1]; j++)
                    for (i = 0; i < geometry.grids[m].index_dim[0]; i++)
                        {
                            if (doing_emission)
                                {
                                    // save the existing value of the radiation field (needed for
                                    // dust or ere emission) but only on the first iteration
                                    // (i.e., save the stellar RT radiation field)
                                    for (n = 0; n < runinfo.n_waves; n++)
                                        {
                                            if (doing_emission == 1)
                                                {
                                                    geometry.grids[m].grid (i, j, k).save_radiation_field_density[n]
                                                        = geometry.grids[m].grid (i, j, k).absorbed_energy[n];
                                                    geometry.grids[m].grid (i, j, k).save_radiation_field_density_x2[n]
                                                        = geometry.grids[m].grid (i, j, k).absorbed_energy_x2[n];
                                                    geometry.grids[m]
                                                        .grid (i, j, k)
                                                        .save_radiation_field_density_num_photons[n]
                                                        = geometry.grids[m]
                                                              .grid (i, j, k)
                                                              .absorbed_energy_num_photons[n];
                                                }
                                            // 	      if ((k == 5) && (j == 5) && (i == 5)) {
                                            // 		cout << n << " " <<
                                            // geometry.grids[m].grid(i,j,k).save_radiation_field_density[n]
                                            // << endl;
                                            // 	      }
                                            // zero out the current absorbed energy value
                                            geometry.grids[m].grid (i, j, k).absorbed_energy[n] = 0.0;
                                            geometry.grids[m].grid (i, j, k).absorbed_energy_x2[n] = 0.0;
                                            geometry.grids[m].grid (i, j, k).absorbed_energy_num_photons[n] = 0;
                                        }
                                }
                            else
                                {
                                    // if this is the first time or if memory storage requested
                                    // then push zero into the absorbed_energy variable in each
                                    // grid cell
                                    if (geometry.abs_energy_storage_type == 1)
                                        { // disk memory
                                            cout << "need new code for disk storage of absorbed "
                                                    "energy grid"
                                                 << endl;
                                            cout << "ever since limited polychromaticism was "
                                                    "implemented (27 "
                                                    "Mar 2009 - KDG)"
                                                 << endl;
                                            exit (8);
                                            // 	      if
                                            // (!geometry.abs_energy_grid_initialized) {
                                            // 		geometry.grids[m].grid(i,j,k).absorbed_energy.push_back(0.0);
                                            // 		geometry.grids[m].grid(i,j,k).absorbed_energy_num_photons.push_back(0);
                                            // 		geometry.grids[m].grid(i,j,k).save_radiation_field_density.push_back(0.0);
                                            // 		geometry.grids[m].grid(i,j,k).save_radiation_field_density_num_photons.push_back(0);
                                            // 	      } else {
                                            // 		geometry.grids[m].grid(i,j,k).absorbed_energy[0]
                                            // = 0.0;
                                            // 		geometry.grids[m].grid(i,j,k).absorbed_energy_num_photons[0]
                                            // = 0;
                                            // 	      }
                                        }
                                    else if (!geometry.abs_energy_grid_initialized)
                                        {
                                            geometry.grids[m].grid (i, j, k).absorbed_energy.resize (runinfo.n_waves,
                                                                                                     0.0);
                                            geometry.grids[m].grid (i, j, k).absorbed_energy_x2.resize (runinfo.n_waves,
                                                                                                        0.0);
                                            geometry.grids[m].grid (i, j, k).absorbed_energy_num_photons.resize (
                                                runinfo.n_waves, 0);
                                            geometry.grids[m].grid (i, j, k).last_photon_number = -1;
                                            geometry.grids[m].grid (i, j, k).last_photon_absorbed_energy = 0.0;

                                            if (doing_emission == 1)
                                                {
                                                    geometry.grids[m]
                                                        .grid (i, j, k)
                                                        .save_radiation_field_density.resize (runinfo.n_waves, 0.0);
                                                    geometry.grids[m]
                                                        .grid (i, j, k)
                                                        .save_radiation_field_density_x2.resize (runinfo.n_waves, 0.0);
                                                    geometry.grids[m]
                                                        .grid (i, j, k)
                                                        .save_radiation_field_density_num_photons.resize (
                                                            runinfo.n_waves, 0);
                                                }
                                        }
                                    else
                                        {
                                            for (n = 0; n < runinfo.n_waves; n++)
                                                {
                                                    geometry.grids[m].grid (i, j, k).absorbed_energy[n] = 0.0;
                                                    geometry.grids[m].grid (i, j, k).absorbed_energy_x2[n] = 0.0;
                                                    geometry.grids[m].grid (i, j, k).absorbed_energy_num_photons[n] = 0;
                                                }
                                        }
                                }
                        }
        }
    // set to know that the absorbed energy grid has been initialized the first
    // time
    if (!geometry.abs_energy_grid_initialized)
        geometry.abs_energy_grid_initialized = 1;
}
