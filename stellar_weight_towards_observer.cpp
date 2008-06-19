// ======================================================================
//   Function to calculate the stellar weight in the direction of the 
// observer.
//
// 2007 Feb/KDG - written (adapted from scattered_weight_towards_observer)
// ======================================================================
#include "stellar_weight_towards_observer.h"
//#define DEBUG_STWTO

double stellar_weight_towards_observer (photon_data photon,
					geometry_struct& geometry,
					float observer_position[3])
  
{
  // determine the vector, distance, and direction cosines from the photon
  // birth site to the observer
  double birth_to_obs[3];
  double dist_birth_to_obs = 0.0;
  double dir_cosines_birth_to_obs[3];
  int i;
  for (i = 0; i < 3; i++) {
    birth_to_obs[i] = observer_position[i] - photon.position[i];
    dist_birth_to_obs += pow(birth_to_obs[i],2);
  }
  dist_birth_to_obs = sqrt(dist_birth_to_obs);
  for (i = 0; i < 3; i++)
    dir_cosines_birth_to_obs[i] = birth_to_obs[i]/dist_birth_to_obs;

#ifdef DEBUG_STWTO
  cout << "birth_to_obs vector = ";
  for (i = 0; i < 3; i++) 
    cout << birth_to_obs[i] << " ";
  cout << endl;
  cout << "birth_to_obs dir_cosines = ";
  for (i = 0; i < 3; i++) 
    cout << dir_cosines_birth_to_obs[i] << " ";
  cout << endl;
#endif

  // determine the optical depth to the surface from the birth site
  // towards the observer
  for (i = 0; i < 3; i++)
    photon.dir_cosines[i] = dir_cosines_birth_to_obs[i];
  double target_tau = 1e20;
  int escape = 0;
  double tau_birth_to_obs = 0.0;
#ifdef DEBUG_STWTO
  cout << "stellar weight" << endl;
#endif
  // redetermine the location of the photon in the grid
  determine_photon_position_index_initial(geometry, photon);

  double distance_traveled = 0.0;
  distance_traveled = calc_photon_trajectory(photon, geometry, target_tau, escape, tau_birth_to_obs);
#ifdef DEBUG_STWTO
  cout << "tau/distance traveled = " << tau_birth_to_obs << " ";
  cout << distance_traveled << endl;
  cout << "Tau from birth site to surface = " << tau_birth_to_obs << endl;
#endif

  double stellar_weight = photon.stellar_weight;

  // calculate the portion of the probability from the tau
  double stellar_weight_tau = exp(-tau_birth_to_obs);

#ifdef DEBUG_STWTO
  cout << "weight(tau) = " << stellar_weight_tau << endl;
#endif

  // udpate stellar weight
  stellar_weight *= stellar_weight_tau;

  return(stellar_weight);
}
