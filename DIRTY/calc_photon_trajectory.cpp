// ======================================================================
//   Function to calculate the trajectory of a photon.  The stopping
// point is either reaching the target optical depth (target_tau) or 
// reaching the edge of the distribution (marked by negative optical
// depth between (0 and, but not including -1: integer negative 
// optical depths denote subgrid indexes.
//
// This function returns the distance traveled.
//
// 2003 Jun/KDG - written
// ======================================================================
#include "calc_photon_trajectory.h"
//#define PHOTON_POS
//#define OUTNUM 0
//#define DEBUG_CPT

double calc_photon_trajectory (photon_data& photon,
			       geometry_struct& geometry,
			       double target_tau,
			       int& escape,
			       double& tau_traveled)

{
  double tau_left = target_tau;  // reduce till zero = done
#ifdef DEBUG_CPT
  if (photon.number == OUTNUM) {
    cout << "oooooooooooooooooooooooooooooooooo" << endl;
    cout << "start cpt; ";
    cout << "tau_left = " << tau_left << endl;
    cout.flush();
  }
#endif
  double distance_traveled = 0.0;
  double delta_tau = 0.0;
  double delta_dist = 0.0;

  // move through the grid until escaping or reaching the target tau
  //   need to also handle the case where the photon starts several subgrids down
  //   multiple calls to calc_photon_trajectory needed
  while ((tau_left > ROUNDOFF_ERR_TRIG) && (!escape)) {

#ifdef DEBUG_CPT
    if (photon.number == OUTNUM) {
      cout << "before delta distance & tau_left & delta_tau & escape = ";
      cout << delta_dist << " ";
      cout << tau_left << " ";
      cout << delta_tau << " ";
      cout << escape << endl;
      cout << "current_grid_num, num_current grids = ";
      cout << photon.current_grid_num << " ";
      cout << photon.num_current_grids << endl;
    }
#endif
    delta_dist = calc_delta_dist(photon, geometry, tau_left, escape, delta_tau);
#ifdef PHOTON_POS
    int i = 0;
    for (i = 0; i < 3; i++)
      cout << photon.position[i] << " ";
    cout << endl;
#endif    
    tau_traveled += delta_tau;
    tau_left -= delta_tau;
    distance_traveled += delta_dist;
#ifdef DEBUG_CPT
    if (photon.number == OUTNUM) {
      cout << "cpt: photon.path_cur_cells = " << photon.path_cur_cells << endl;
      cout << "delta distance & tau_left & delta_tau & escape = ";
      cout << delta_dist << " ";
      cout << tau_left << " ";
      cout << delta_tau << " ";
      cout << escape << endl;
      cout << "current_grid_num, num_current grids = ";
      cout << photon.current_grid_num << " ";
      cout << photon.num_current_grids << endl;
    }
#endif
  }

  // check if the photon is in the right grid cell
  // this should be a temporary check...
  // make sure the photon is in this grid
//   int k = photon.current_grid_num;
//   int cur_grid_num = photon.grid_number[k];
//   int i = 0;
//   for (i = 0; i < 3; i++) {
//     if ((photon.position[i] < geometry.grids[cur_grid_num].positions[i][0]) ||
// 	(photon.position[i] > geometry.grids[cur_grid_num].positions[i][geometry.grids[cur_grid_num].index_dim[i]])) {
//       cout << "outside of bounds of current grid (in calc_photon_trajectory.cpp)" << endl;
//       cout << "This should not happen." << endl;
//       cout << "photon # = " << photon.number << endl;
//       cout << "i = " << i << endl;
//       cout << "photon.position[i] = " << photon.position[i] << endl;
//       cout << "photon.birth_position[i] = " << photon.birth_position[i] << endl;
//       cout << "photon.dir_consines[i] = " << photon.dir_cosines[i] << endl;
//       cout << "photon.num_scat = " << photon.num_scat << endl;
//       cout << "cur_grid_num = " << cur_grid_num << endl;
//       cout << "min/max of grid = " << geometry.grids[cur_grid_num].positions[i][0] << " ";
//       cout << geometry.grids[cur_grid_num].positions[i][geometry.grids[cur_grid_num].index_dim[i]] << endl;
//       cout << "size of cube [pc] = " << geometry.grids[cur_grid_num].phys_cube_size[i] << endl;
//       cout << "proposed index = " << int((photon.position[i] - geometry.grids[cur_grid_num].positions[i][0])/
// 					 geometry.grids[cur_grid_num].phys_cube_size[i]) << endl;
//       cout << "real value of pindex = " << (photon.position[i] - geometry.grids[cur_grid_num].positions[i][0])/
// 	geometry.grids[cur_grid_num].phys_cube_size[i] << endl;
//       cout << "out of a possible " << geometry.grids[cur_grid_num].index_dim[i] << endl;
//       cout << "out of a possible2 " << geometry.grids[cur_grid_num].positions[i].size() << endl;
//       int parent_grid_num = geometry.grids[cur_grid_num].parent_grid_num;
//       cout << "parent grid cell min = " << geometry.grids[parent_grid_num].positions[i][photon.position_index[k-1][i]] << endl;
//       cout << "parent grid cell max = " << geometry.grids[parent_grid_num].positions[i][photon.position_index[k-1][i]+1] << endl;
//       exit(8);
//     }
//   }

  return(distance_traveled);
}
