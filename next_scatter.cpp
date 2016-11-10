// ======================================================================
//   Procedure to scatter a photon into a new direction.  This procedure
// returns a 0 if the photon scatters inside the dust and a 1 if it 
// escapes.
//
// 2004 Dec/KDG - written
// 03 Sep 2015/KDG - added sampling from an exp(tau/100) distribution
//                   1/2 of the time to better sample high optical depth scattering
// ======================================================================
#include "next_scatter.h"
#define DEBUG_NS

int next_scatter (geometry_struct& geometry,
		  photon_data& photon,
		  random_dirty& random_obj)

{
  // **test***
  // float save_scat_bias_fraction = geometry.scat_bias_fraction;
  // geometry.scat_bias_fraction = 0.0;
  
  // find path_tau[]
  photon_data dummy_photon = photon;
  dummy_photon.current_grid_num = 0;  // set to the base grid to start tarjectory correctly
  dummy_photon.path_cur_cells = 0; // set to 0 to save cells traversed
  
  double target_tau = 1e20;
  int escape = 0;
  double tau_path = 0.0;
  calc_photon_trajectory(dummy_photon, geometry, target_tau, escape, tau_path);
  //double bias_norm = 1.0/(1.0 + tau_path);
  // **testing***
  double bias_norm = 1.0/10.0;

  // determine the optical depth to the next scattering
  target_tau = 0.0;

  double ran_num = random_obj.random_num();
  double ran_num2 = random_obj.random_num();
  if (ran_num >= geometry.scat_bias_fraction) { // classical scattering
    target_tau = -log(ran_num2);
  } else { // biased based on tau_path
    target_tau = (-1./bias_norm)*log(ran_num2);
  }

  photon.target_tau = target_tau;

  // check to see if we will start in a subgrid
  if (photon.current_grid_num > 0) {
#ifdef DEBUG_NS
    if (photon.number == OUTNUM) {
      cout << "starting in a subgrid" << endl;
    }
#endif
    photon.current_grid_num = 0;
  }
  
  // determine the site of the next scattering
  double distance_traveled = 0.0;
  double tau_traveled = 0.0;
  escape = 0;
  photon.path_cur_cells = 0;  // set to 0 to save cells tranversed

  distance_traveled = calc_photon_trajectory(photon, geometry, target_tau, escape, tau_traveled);
#ifdef DEBUG_NS
  if (photon.number == OUTNUM) {
    cout << "ns cpt done; ";
    cout << "distance_traveled = " << distance_traveled << endl;
    cout << "target_tau = " << target_tau << endl;
    cout << "photon.scat_weight = " << photon.scat_weight << endl;
  }
#endif

//   int j = 0;
//   for (j = 0; j < 3; j++) 
//     cout << photon.position[j] << " ";
//   cout << endl;

  escape = 0;
  // check if the photon has left the dust
  if ((target_tau - tau_traveled)/geometry.tau > ROUNDOFF_ERR_TRIG)
    escape = 1;

  // check if the photon has scattered enough and there is just no significant weight left
  if (photon.num_scat > geometry.max_num_scat)
    escape = 1;

  // cout << target_tau - tau_traveled << " ";
  // cout << geometry.max_num_scat << " ";
  // cout << geometry.tau << " ";

  // cout << photon.number << " ";
  // cout << photon.num_scat << " ";
  // cout << escape << " ";
  // cout << target_tau << " ";
  // cout << tau_traveled << " ";
  // cout << distance_traveled << endl;

  if (!escape) {
    // update the scattered weight for biasing
    double biased_weight_factor = 0.0;
    biased_weight_factor = (1.0 - (geometry.scat_bias_fraction)) +
      (geometry.scat_bias_fraction*bias_norm*exp((1.-bias_norm)*target_tau));

//     cout << 1.0/biased_weight_factor << " ";

    photon.scat_weight /= biased_weight_factor;

//     cout << photon.number << " ";
//     cout << target_tau << " ";
//     cout << photon.scat_weight << " ";
//     cout << 1.0/biased_weight_factor << " ";
//     cout << endl;
  }

//   cout << endl;

#ifdef DEBUG_NS
  if (photon.number == OUTNUM) {
    cout << "ns escape = " << escape << endl;
  }
#endif
#ifdef DEBUG_NS
  if (photon.number == OUTNUM) {
    if ((target_tau - tau_traveled) < -ROUNDOFF_ERR_TRIG) {
      cout << "*****error*****next_scatter*****" << endl;
      cout << "target_tau = " << target_tau << endl;
      cout << "tau_traveled = " << tau_traveled << endl;
      cout << "diff = " << target_tau - tau_traveled << endl;
    }
  }
#endif

  // **test***
  // geometry.scat_bias_fraction = save_scat_bias_fraction;
  
  // return escape (1 = yes, 0 = no)
  return(escape);
}
