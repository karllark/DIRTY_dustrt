// ======================================================================
//   Procedure to generate a photon and set up the stats on it.
// This procedure is used for a single star located somewhere in the
// dust distribution
//
// 2003 Jun/KDG - written
// 2004 Dec/KDG - various variable initializations
// 2006 Aug/KDG - changed from "1star" to "discrete_stars" to hand
//                1 to n stars with specific positions and luminosities
// ======================================================================
#include "new_photon_discrete_stars.h"
// #define DEBUG_NPDS

void new_photon_discrete_stars(photon_data &photon, geometry_struct &geometry,
                               random_dirty &random_obj)

{
  // setup the weights
  photon.stellar_weight = 1.0;
  photon.scat_weight = 0.0;

  // initialize statistics variables
  photon.num_scat = 0;

  // direction of photon; assuming an isotropic source
  // in direction cosines...
  double phi = M_PI * (2.0 * random_obj.random_num() - 1.0);
  photon.dir_cosines[2] = 2.0 * random_obj.random_num() - 1.0;
  double temp = sqrt(1.0 - pow(photon.dir_cosines[2], 2));
  photon.dir_cosines[0] = cos(phi) * temp;
  photon.dir_cosines[1] = sin(phi) * temp;

  int i = 0;
#ifdef DEBUG_NPDS
  for (i = 0; i < 3; i++)
    cout << photon.dir_cosines[i] << " ";
  cout << "starting dir cosines" << endl;
#endif

  // determine which star the photon will be emitted from
  int pos_index = 0; // default is first star
  if (geometry.num_stars > 1) {
    double ran_num = random_obj.random_num();
    while ((i < geometry.num_stars) &&
           (ran_num > geometry.star_positions[4][i]))
      i++;
    pos_index = i;
    //     cout << "random # = " << ran_num << endl;
    //     cout << "geometry.star_positions = " <<
    //     geometry.star_positions[4][i-1]; cout << " " <<
    //     geometry.star_positions[4][i] << endl; exit(8);
  }

  // set the photon position to the star position
  for (i = 0; i < 3; i++) {
    photon.position[i] = geometry.star_positions[i][pos_index];
    // save the birth photon position
    photon.birth_position[i] = photon.position[i];
  }

  // now determine the position indexes of the photon
  determine_photon_position_index_initial(geometry, photon);
}
