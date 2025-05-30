// ======================================================================
//   Function to calculate the stellar weight in the direction of the
// observer.
//
// 2007 Feb/KDG - written (adapted from scattered_weight_towards_observer)
// ======================================================================
#include "stellar_weight_towards_observer.h"
// #define DEBUG_STWTO
// #define OUTNUM 1

double stellar_weight_towards_observer(photon_data photon, geometry_struct &geometry, float observer_position[3])

{
    // determine the vector, distance, and direction cosines from the photon
    // birth site to the observer
    double birth_to_obs[3];
    double dist_birth_to_obs = 0.0;
    double dir_cosines_birth_to_obs[3];
    int i;
    for (i = 0; i < 3; i++)
    {
        birth_to_obs[i] = observer_position[i] - photon.position[i];
        dist_birth_to_obs += pow(birth_to_obs[i], 2);
    }
    dist_birth_to_obs = sqrt(dist_birth_to_obs);
    for (i = 0; i < 3; i++)
        dir_cosines_birth_to_obs[i] = birth_to_obs[i] / dist_birth_to_obs;

#ifdef DEBUG_STWTO
    if (photon.number == OUTNUM)
    {
        cout << "birth_to_obs vector = ";
        for (i = 0; i < 3; i++)
            cout << birth_to_obs[i] << " ";
        cout << endl;
        cout << "birth_to_obs dir_cosines = ";
        for (i = 0; i < 3; i++)
            cout << dir_cosines_birth_to_obs[i] << " ";
        cout << endl;
    }
#endif

    // determine the optical depth to the surface from the birth site
    // towards the observer
    for (i = 0; i < 3; i++)
        photon.dir_cosines[i] = dir_cosines_birth_to_obs[i];
    double target_tau = 1e20;
    double target_dist = 1e10 * geometry.radius;
    if (geometry.internal_observer == 1)
        target_dist = dist_birth_to_obs;
    int escape = 0;
    double tau_birth_to_obs = 0.0;

    // check that the direction to the observer will not force an
    // immediate exit from the grid, step back slightly to avoid this
    for (i = 0; i < 3; i++)
    {
        if ((photon.dir_cosines[i] < 0.) &&
            (photon.position[i] == geometry.grids[photon.current_grid_num].positions[i][0]))
            photon.position[i] += 0.01 * (geometry.grids[photon.current_grid_num].positions[i][1] -
                                          geometry.grids[photon.current_grid_num].positions[i][0]);
        else if ((photon.dir_cosines[i] > 0.) &&
                 (photon.position[i] == geometry.grids[photon.current_grid_num]
                                            .positions[i][geometry.grids[photon.current_grid_num].index_dim[i]]))
            photon.position[i] -= 0.01 * (geometry.grids[photon.current_grid_num]
                                              .positions[i][geometry.grids[photon.current_grid_num].index_dim[i]] -
                                          geometry.grids[photon.current_grid_num]
                                              .positions[i][geometry.grids[photon.current_grid_num].index_dim[i] - 1]);
    }

#ifdef DEBUG_STWTO
    if (photon.number == OUTNUM)
    {
        cout << "stellar weight" << endl;
    }
#endif
    // redetermine the location of the photon in the grid
    determine_photon_position_index_initial(geometry, photon);

#ifdef DEBUG_STWTO
    if (photon.number == OUTNUM)
    {
        cout << "initial indexs done" << endl;
        cout << "num_current_grids = " << photon.num_current_grids << endl;
        cout << "current_grid_num = " << photon.current_grid_num << endl;
        for (i = 0; i < 3; i++)
            cout << photon.position_index[photon.current_grid_num][i] << " ";
        cout << endl;
        int k = photon.current_grid_num;
        int grid_val = photon.grid_number[k];
        cout << geometry.grids[grid_val]
                    .grid(photon.position_index[k][0], photon.position_index[k][1], photon.position_index[k][2])
                    .dust_tau_per_pc
             << endl;
        cout << "done" << endl;
    }
#endif

    double distance_traveled = 0.0;
    distance_traveled = calc_photon_trajectory(photon, geometry, target_tau, target_dist, escape, tau_birth_to_obs, 0);

#ifdef DEBUG_STWTO
    if (photon.number == OUTNUM)
    {
        cout << "tau/distance traveled = " << tau_birth_to_obs << " ";
        cout << distance_traveled << endl;
        cout << "Tau from birth site to surface = " << tau_birth_to_obs << endl;
    }
#endif

    double stellar_weight = photon.stellar_weight;

    // calculate the portion of the probability from the tau
    double stellar_weight_tau = exp(-tau_birth_to_obs);

#ifdef DEBUG_STWTO
    if (photon.number == OUTNUM)
    {
        cout << "weight(tau) = " << stellar_weight_tau << endl;
    }
#endif

    // udpate stellar weight
    stellar_weight *= stellar_weight_tau;

    return (stellar_weight);
}
