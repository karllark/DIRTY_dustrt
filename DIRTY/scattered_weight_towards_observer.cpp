// ======================================================================
//   Function to calculate the scattered weight in the direction of the
// observer.
//
// 2006 Apr/KDG - written
// ======================================================================
#include "scattered_weight_towards_observer.h"
// #define OUTNUM 73145
// #define DEBUG_SWTO

double
scattered_weight_towards_observer (photon_data photon, geometry_struct &geometry, float observer_position[3])

{
  // determine the vector, distance, and direction cosines from the scattering
  // site to the observer
  double scat_to_obs[3];
  double dist_scat_to_obs = 0.0;
  double dir_cosines_scat_to_obs[3];
  int i;
  for (i = 0; i < 3; i++)
    {
      scat_to_obs[i] = observer_position[i] - photon.position[i];
      dist_scat_to_obs += pow (scat_to_obs[i], 2);
    }
  dist_scat_to_obs = sqrt (dist_scat_to_obs);
  for (i = 0; i < 3; i++)
    dir_cosines_scat_to_obs[i] = scat_to_obs[i] / dist_scat_to_obs;

#ifdef DEBUG_SWTO
  if (photon.number == OUTNUM)
    {
      cout << "pos vector = ";
      for (i = 0; i < 3; i++)
        cout << photon.position[i] << " ";
      cout << endl;
      cout << "dir cosines vector = ";
      for (i = 0; i < 3; i++)
        cout << photon.dir_cosines[i] << " ";
      cout << endl;
      cout << "scat_to_obs vector = ";
      for (i = 0; i < 3; i++)
        cout << scat_to_obs[i] << " ";
      cout << endl;
      cout << "scat_to_obs dir_cosines = ";
      for (i = 0; i < 3; i++)
        cout << dir_cosines_scat_to_obs[i] << " ";
      cout << endl;
      cout << "cos_alpha contributions = ";
      for (i = 0; i < 3; i++)
        cout << dir_cosines_scat_to_obs[i] * photon.dir_cosines[i] << " ";
      cout << endl;
    }
#endif

  // determine the angle between the photon direction before the scattering and
  // the direction to the observer
  double cos_alpha = 0.0;
  for (i = 0; i < 3; i++)
    cos_alpha += dir_cosines_scat_to_obs[i] * photon.dir_cosines[i];
  // check the resulting value is reasonable
  if (fabs (cos_alpha) > 1.0)
    {
      cout << "ERROR: scattered_weight_towards_observer" << endl;
      cout << "cos_alpha = " << cos_alpha << " is out of bounds [-1,1]" << endl;
      exit (8);
    }

  //   double scat_weight = 0.0;
  // scat_weight += dir_cosines_scat_to_obs[0]*photon.dir_cosines[0];
  // scat_weight += dir_cosines_scat_to_obs[1]*photon.dir_cosines[1];
  //   scat_weight += dir_cosines_scat_to_obs[2]*photon.dir_cosines[2];

  // determine the optical depth to the surface from the scattering site
  // towards the observer
  for (i = 0; i < 3; i++)
    photon.dir_cosines[i] = dir_cosines_scat_to_obs[i];
  double target_tau = 1e20;
  double target_dist = 1e10 * geometry.radius;
  if (geometry.internal_observer == 1)
    target_dist = dist_scat_to_obs;
  int escape = 0;
  double tau_scat_to_obs = 0.0;
  double distance_traveled = 0.0;
  photon.current_grid_num = 0; // set to the base grid to start tarjectory correctly
#ifdef DEBUG_SWTO
  if (photon.number == OUTNUM)
    {
      cout << "starting calc_photon_traj..." << endl;
    }
#endif
  try
    {
      distance_traveled
          = calc_photon_trajectory (photon, geometry, target_tau, target_dist, escape, tau_scat_to_obs, 0);
    }
  catch (std::out_of_range)
    {
      cout << "photon # = " << photon.number << endl;
      cout << "num scat = " << photon.num_scat << endl;
      cout << "scattered_weight_towards_observer: calc_photon_trajectory - out "
              "of range."
           << endl;
      exit (8);
    }
#ifdef DEBUG_SWTO
  if (photon.number == OUTNUM)
    {
      cout << "Tau from scattering site to surface = " << tau_scat_to_obs << endl;
      cout << "cosine of angle = " << cos_alpha << endl;
    }
#endif

  double scat_weight = photon.scat_weight;

  double scat_weight_angle = 0.0;
  // calculate the portion of the probability from the angle
  if (geometry.g == -2)
    { // model phase function
      uint i1, i2, i3;
      i1 = 0;
      i3 = geometry.phi.size ();
      i2 = (i1 + i3) / 2;
      //     cout << i3 << endl;
      while ((i3 - i1) > 1)
        {
          if (cos_alpha < geometry.phi_angle[i2])
            i1 = i2;
          else
            i3 = i2;
          //       cout << cos_alpha << " ";
          //       cout << geometry.phi_angle[i1] << " " <<
          //       geometry.phi_angle[i3]
          //       << " "; cout << i1 << " " << i3 << endl;
          i2 = (i1 + i3) / 2;
        }
      scat_weight_angle = geometry.phi[i1]
                          + ((cos_alpha - geometry.phi_angle[i3]) / (geometry.phi_angle[i1] - geometry.phi_angle[i3]))
                                * (geometry.phi[i1] - geometry.phi[i3]);
      //     cout << ((cos_alpha -
      //     geometry.phi_angle[i3])/(geometry.phi_angle[i1]
      //     - geometry.phi_angle[i3])) <<  endl; cout << cos_alpha << endl;
      //     cout
      //     << geometry.phi_angle[i1] << endl; cout << geometry.phi[i1] << " "
      //     << geometry.phi[i3] << endl; cout << scat_weight_angle << endl;
      //     double tmp_g = 0.384; cout << ((1.0 - pow(tmp_g,2))/
      // 	     (4.0*M_PI*pow(1.0 + pow(tmp_g,2) -
      // 			   2.0*tmp_g*cos_alpha,1.5))) << endl;
      //     exit(8);
    }
  else
    {
      scat_weight_angle = ((1.0 - pow (geometry.g, 2))
                           / (4.0 * M_PI * pow (1.0 + pow (geometry.g, 2) - 2.0 * geometry.g * cos_alpha, 1.5)));
    }
  scat_weight_angle *= geometry.solid_angle;

  // calculate the portion of the probability from the tau
  double scat_weight_tau = exp (-tau_scat_to_obs);

#ifdef DEBUG_SWTO
  if (photon.number == OUTNUM)
    {
      cout << "weight(angle) = " << scat_weight_angle << endl;
      cout << "weight(tau) = " << scat_weight_tau << endl;
      cout << "geometry.solid_angle = " << geometry.solid_angle << endl;
    }
#endif

  // udpate scattered weight
  scat_weight *= scat_weight_angle * scat_weight_tau;
  //   scat_weight = cos_alpha;

  //   if (scat_weight > 1.) {
  //     exit(8);
  //   }

  return (scat_weight);
}
