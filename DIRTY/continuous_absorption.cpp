/*
 * Functions related to continuous absorption. deposit_energy() deposits
 * the energy along the path using a path_tau structures. move_photon()
 * finds the scattering position given a target_tau and the path_tau
 * structures in the dummy photon.
 *
 * Note:
 * This has only been tested with a single-level grid (no subdivision)
 *
 * Reference:
 * http://dirty.as.arizona.edu/twiki/bin/view/DIRTY/ContinuousAbsorption
 *
 * 2012 Apr/KHL - written
 * 2015 Aug/KDG - cube size updated to allow for non-uniform grid spacings
 */
#include "continuous_absorption.h"
#define DEBUG_CONT_ABS

void
deposit_energy (const photon_data &dummy_photon, geometry_struct &geometry, const double photon_weight)
{
    /*
     * Continuous absorption
     * tau_entering/tau_leaving:
     *   total tau traveled when the photon packet enters/leaves the grid cell
     * prob_entering/prob_leaving:
     *   probability that the photon packet enters/leaves the grid cell
     */
    const double abs_weight_init = (1. - geometry.albedo) * photon_weight; // not the reduced scat_weight
    double tau_entering = 0.;
    double prob_entering = 1.;

    for (int i = 0; i < dummy_photon.path_cur_cells; i++)
        {
            cout << i << endl;
            // find the absorbed weight
            double tau_leaving = tau_entering + dummy_photon.path_tau[i];
            double prob_leaving = exp (-tau_leaving);
            double abs_weight = abs_weight_init * (prob_entering - prob_leaving);

            // deposit the energy
            grid_cell &this_cell = geometry.grids[dummy_photon.path_pos_index[0][i]].grid (
                dummy_photon.path_pos_index[1][i], dummy_photon.path_pos_index[2][i],
                dummy_photon.path_pos_index[3][i]);
            this_cell.absorbed_energy[geometry.abs_energy_wave_index] += abs_weight;
            this_cell.absorbed_energy_x2[geometry.abs_energy_wave_index] += abs_weight * abs_weight;
            this_cell.absorbed_energy_num_photons[geometry.abs_energy_wave_index]++;

            // move to the next grid cell
            tau_entering = tau_leaving;
            prob_entering = prob_leaving;
        }
}

void
move_photon (photon_data &photon, const photon_data &dummy_photon, geometry_struct &geometry, const double target_tau,
             double &tau_traveled)
{
#ifdef DEBUG_CONT_ABS
    assert (photon.current_grid_num == 0);
    assert (photon.grid_number[0] == dummy_photon.path_pos_index[0][0]);
    assert (photon.position_index[0][0] == dummy_photon.path_pos_index[1][0]);
    assert (photon.position_index[0][1] == dummy_photon.path_pos_index[2][0]);
    assert (photon.position_index[0][2] == dummy_photon.path_pos_index[3][0]);
#endif

    // ensure that the photon is in the grid cell
    one_grid &current_grid = geometry.grids[photon.grid_number[photon.current_grid_num]];
    for (int j = 0; j < 3; j++)
        {
            // compare photon.position with dummy_photon.path_pos_index[j][0] (1st
            // path segment)
            bool b1 = photon.position[j] < current_grid.positions[j][dummy_photon.path_pos_index[j + 1][0]];
            bool b2 = photon.position[j] > current_grid.positions[j][dummy_photon.path_pos_index[j + 1][0] + 1];
            if (b1 || b2)
                {
                    double correct_pos
                        = current_grid.positions[j][dummy_photon.path_pos_index[j + 1][0] + ((b1) ? 0 : 1)];
                    double cube_size
                        = current_grid.positions[j][dummy_photon.path_pos_index[j + 1][0] + ((b1) ? 0 : 1) + 1]
                          - current_grid.positions[j][dummy_photon.path_pos_index[j + 1][0] + ((b1) ? 0 : 1)];
                    double frac_miss = (photon.position[j] - correct_pos) / cube_size;
                    cout << "move_photon(): photon #" << photon.number << ":" << photon.num_scat
                         << " not in orig cell, frac_miss = " << frac_miss << endl;
                    // quit if frac_miss is too big
                    assert (fabs (frac_miss) < ROUNDOFF_ERR_INDEX);
                    photon.position[j] = correct_pos;
                }
        }

    double tau_left = target_tau;

    // walk through the path until we find the scattering cell
    for (int i = 0; i < dummy_photon.path_cur_cells; i++)
        {
            // if (tau_left - photon.path_tau[i] <= ROUNDOFF_ERR_TRIG) {
            if (tau_left <= dummy_photon.path_tau[i])
                {
                    // photon will scatter in this cell

                    // position index of the scattering cell
                    int k = dummy_photon.current_grid_num;
                    photon.grid_number[k] = dummy_photon.path_pos_index[0][i];
                    photon.position_index[k][0] = dummy_photon.path_pos_index[1][i];
                    photon.position_index[k][1] = dummy_photon.path_pos_index[2][i];
                    photon.position_index[k][2] = dummy_photon.path_pos_index[3][i];

                    // figure out the position when the photon enters the scattering cell
                    // (if i == 0, the current cell is already the scattering cell)
                    if (i != 0)
                        {
                            // find out the distance traveled
                            int enter_index = -1; // on which side the photon enters this grid cell
                            double distance_traveled = -1;

                            // find the side the photon enters this grid cell
                            for (int j = 0; j < 3; j++)
                                {
                                    // check criteria: path_pos_index has changed in that
                                    // direction
                                    if (dummy_photon.path_pos_index[j + 1][i]
                                        != dummy_photon.path_pos_index[j + 1][i - 1])
                                        {
#ifdef DEBUG_CONT_ABS
                                            // make sure that path_pos_index has not changed in other
                                            // directions
                                            assert (i - 1 >= 0);
                                            assert (dummy_photon.path_pos_index[(j + 1) % 3 + 1][i]
                                                    == dummy_photon.path_pos_index[(j + 1) % 3 + 1][i - 1]);
                                            assert (dummy_photon.path_pos_index[(j + 2) % 3 + 1][i]
                                                    == dummy_photon.path_pos_index[(j + 2) % 3 + 1][i - 1]);
#endif

                                            // find the position on the side of the grid cell
                                            double dir_cosine = photon.dir_cosines[j];
                                            double new_pos;
                                            if (dir_cosine > 0.0)
                                                {
                                                    // if the photon is moving forward, it hits the back
                                                    // of the grid cell
                                                    new_pos = geometry.grids[photon.grid_number[k]]
                                                                  .positions[j][photon.position_index[k][j]];
                                                }
                                            else if (dir_cosine < 0.0)
                                                {
                                                    // if the photon is moving backward, it hits the
                                                    // front of the grid cell
                                                    new_pos = geometry.grids[photon.grid_number[k]]
                                                                  .positions[j][photon.position_index[k][j] + 1];
                                                }
                                            else
                                                {
                                                    cout << "we shouldn't reach here (dir_cosine is 0)" << endl;
                                                    exit (2);
                                                }

                                            // move photon in the entering direction
                                            enter_index = j;
                                            distance_traveled = (new_pos - photon.position[j]) / dir_cosine;
                                            photon.position[j] = new_pos;
                                            break;
                                        }
                                }
#ifdef DEBUG_CONT_ABS
                            assert (enter_index >= 0);
                            assert (distance_traveled >= 0);
#endif

                            // move the photon in the other two directions
                            for (int j = 0; j < 3; j++)
                                {
                                    if (j == enter_index)
                                        continue;
                                    photon.position[j] += photon.dir_cosines[j] * distance_traveled;
                                }
                        }

                    // move photon to the actual scattering location
                    /*
                    // using calc_delta_dist(): slow - calc_delta_dist() costs 80% of the
                    time of this function double delta_tau; int escape;
                    photon.path_cur_cells = -1;
                    // set to -1 *not* to save cells traversed calc_delta_dist(photon,
                    geometry, tau_left, escape, delta_tau); tau_traveled += delta_tau;
                    */
                    {
                        // alternatively, directly calculate the scattering position
                        int k = photon.current_grid_num;
                        vector<int> &position_index = photon.position_index[k];
                        one_grid &current_grid = geometry.grids[photon.grid_number[k]];
                        double dust_tau_ref_per_pc
                            = current_grid.grid (position_index[0], position_index[1], position_index[2])
                                  .dust_tau_per_pc;
                        double dust_tau_per_pc = dust_tau_ref_per_pc * geometry.tau_to_tau_ref;

                        for (int j = 0; j < 3; j++)
                            {
                                double distance_j = photon.dir_cosines[j] * tau_left / dust_tau_per_pc;
                                photon.position[j] += distance_j;

                                // ensure that the photon is still in the grid cell
                                bool b1 = photon.position[j] < current_grid.positions[j][position_index[j]];
                                bool b2 = photon.position[j] > current_grid.positions[j][position_index[j] + 1];
                                if (b1 || b2)
                                    {
                                        double correct_pos
                                            = current_grid.positions[j][position_index[j] + ((b1) ? 0 : 1)];
                                        double cube_size
                                            = current_grid.positions[j][position_index[j] + ((b1) ? 0 : 1) + 1]
                                              - current_grid.positions[j][position_index[j] + ((b1) ? 0 : 1)];
                                        double frac_miss = (photon.position[j] - correct_pos) / cube_size;
                                        cout << "move_photon(): photon #" << photon.number << ":" << photon.num_scat
                                             << " not in dest cell, frac_miss = " << frac_miss << endl;
                                        // quit if frac_miss is too big
                                        assert (fabs (frac_miss) < ROUNDOFF_ERR_INDEX);
                                        photon.position[j] = correct_pos;
                                    }
                            }

                        tau_traveled += tau_left;
                    }

                    // done!
                    return;
                }

            // prepare to move on to the next cell
            tau_left -= dummy_photon.path_tau[i];
            tau_traveled += dummy_photon.path_tau[i];
        }

    // if the sum of path_tau[] is smaller than target_tau,
    // the code reaches here and the photon escapes
}
