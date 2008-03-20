#ifndef _DIRTY_PHOTON_DATA_
#define _DIRTY_PHOTON_DATA_

#include <vector>

// defines the necessary information on a photon

struct photon_data { 
  double position[3]; // the physical location [parsec]
  double dir_cosines[3];    // the direction cosines describing the direction the photon is moving

  double birth_position[3];  // the birth physical location [parsec]
  // probability of the photon being emitted by different grain/emission types (used for dust thermal emission)
  std::vector<double> birth_photon_type_prob;

  std::vector< std::vector<int> > position_index;  // vector of indexes for each grid photon is in
  std::vector<long> grid_number;      // in which grid is the photon located
  int num_current_grids;       // number of current grids the photon is in
  int current_grid_num;    // number of current grid photon is in

  double scat_weight;     // scattered portion of the weight
  double stellar_weight;  // stellar portion of the weight

  long num_scat;  // number of scatterings

  long number;  // photon number (nth photon)
  double first_tau;  // optical depth to the surface along the birth direction of the photon
};

#endif
