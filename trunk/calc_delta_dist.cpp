// ======================================================================
//   Function to calculate the distance the photon travels in a cell.
//
// 2003 Jun/KDG - written
// 2007 Feb/KDG - fixed case where dir cosines are zero -> nan distance traveled
//                now zero distance traveled
// ======================================================================
#include "calc_delta_dist.h"
//#define OUTNUM 108473
//#define DEBUG_CDD

double calc_delta_dist (photon_data& photon,
			geometry_struct& geometry,
			double target_tau,
			int& escape,
			double& tau_traveled)
  
{
#ifdef DEBUG_CDD
  if (photon.number == OUTNUM) {
    cout << "starting cdd.." << endl;
  }
#endif
  tau_traveled = 0.0;
  double distance_traveled = 0.0;  // resulting distance traveled
  int exit_cell = 1;  // tells if the cell has been exited
  int min_index = 0;

  // save the current grid number, the grid number, and the dust_tau_per_pc
  int k = photon.current_grid_num;
  long grid_val = photon.grid_number[k];
#ifdef DEBUG_CDD
  if (photon.number == OUTNUM) {
    cout << "k = " << k << endl;
    cout << "grid_val = " << grid_val << endl;
    cout << "max # grids = " << geometry.grids.size() << endl;
    cout << photon.position_index[0].size() << endl;
    cout << geometry.max_grid_depth << endl;
    cout << photon.position_index[k][0] << endl;
    cout << "end testing" << endl;
  }
#endif
  float dust_tau_ref_per_pc = geometry.grids[grid_val].grid(photon.position_index[k][0],photon.position_index[k][1],photon.position_index[k][2]).dust_tau_per_pc;
  float dust_tau_per_pc = dust_tau_ref_per_pc*geometry.tau_to_tau_ref;

  // print the photon statistics
#ifdef DEBUG_CDD
  int di = 0;
  if (photon.number == OUTNUM) {
    cout << endl;
    cout << "cdd photon in" << endl;
    cout << "number = " << photon.number << endl;
    cout << "pos = ";
    for (di = 0; di < 3; di++)
      cout << photon.position[di] << " ";
    cout << endl;
    cout << "dir_cos = ";
    for (di = 0; di < 3; di++)
      cout << photon.dir_cosines[di] << " ";
    cout << endl;
    cout << "pos idx = ";
    for (di = 0; di < 3; di++)
      cout << photon.position_index[k][di] << " ";
    cout << endl;
    cout << "k = " << k << " < " << photon.num_current_grids-1 << endl;
    cout << "dust_tau_per_pc = " << dust_tau_per_pc << endl;
    cout << "number of grids = " << geometry.grids.size() << endl;
    cout << "photon.current_grid_num = " << photon.current_grid_num << endl;
    cout << "photon.grid_number[photon.current_grid_num] = " << photon.grid_number[photon.current_grid_num] << endl;
  }
#endif

  if ((k < (photon.num_current_grids-1)) || (dust_tau_ref_per_pc <= -1.0)) {
#ifdef DEBUG_CDD
    if (photon.number == OUTNUM) {
      cout << endl << "=============" << endl;
      cout << "k, num_current_grids, dust_tau_per_pc, dust_tau_ref_per_pc = " << k << " " << photon.num_current_grids;
      cout << " " << dust_tau_per_pc << " " << dust_tau_ref_per_pc << endl;
    }
#endif
    photon.current_grid_num++;
    // determine the position indexes for this subgrid (any more nested subgrids)
    if (dust_tau_ref_per_pc <= -1.0) {
#ifdef DEBUG_CDD
      if (photon.number == OUTNUM) {
	cout << "new "; cout.flush();
      }
      
      // add in handling of subgrids
      // check if photon has already been in this grid
      // basically, are we restarting the photon tracking in a subgrid?

#endif
      photon.grid_number[photon.current_grid_num] = -int(dust_tau_ref_per_pc);

      // check that the grid number is an actual integer
      // can have problems if tau is accidently set negative
      if ((photon.grid_number[photon.current_grid_num] + dust_tau_ref_per_pc) != 0.0) {
	cout << "problem with subgrid designation of in parent grid." << endl;
	cout << "dust_tau_ref_per_pc = " << dust_tau_ref_per_pc << endl;
	cout << "photon.grid_number[photon.current_grid_num] = " << photon.grid_number[photon.current_grid_num] << endl;
	exit(8);
      }

      if (photon.grid_number[photon.current_grid_num] >= int(geometry.grids.size())) {
	cout << "grid desired does not exist" << endl;
	cout << "grid desired = " << photon.grid_number[photon.current_grid_num] << endl;
	cout << "# of grids = " << geometry.grids.size() << endl;
	cout << "dust_tau_ref_per_pc = " << dust_tau_ref_per_pc << endl;
	exit(8);
      }

      determine_photon_position_index(geometry, photon);
    }
#ifdef DEBUG_CDD
    if (photon.number == OUTNUM) {
      cout << "subgrid found" << endl;
      cout << "target_tau = " << target_tau << endl;
      cout << "tau_traveled = " << tau_traveled << endl;
    }
#endif
    
    distance_traveled = calc_photon_trajectory(photon, geometry, target_tau, escape, tau_traveled);
    
    int escape_grid = 0;
    if (escape) { // determine which of x,y,z exited first
      photon.current_grid_num--;
      photon.num_current_grids--;
      int i = 0;
      escape = 0;
      escape_grid = 1;
      exit_cell = 1;
      for (i = 0; i < 3; i++) {
	if (((photon.dir_cosines[i] > 0.0) && 
	     (photon.position[i] == geometry.grids[photon.grid_number[photon.current_grid_num]].positions[i][photon.position_index[k][i]+1])) ||
	    ((photon.dir_cosines[i] < 0.0) && 
	     (photon.position[i] == geometry.grids[photon.grid_number[photon.current_grid_num]].positions[i][photon.position_index[k][i]])))
	  min_index = i;
#ifdef DEBUG_CDD
	if (photon.number == OUTNUM) {
	  cout << "position, cell wall+/- = ";
	  cout << photon.position[i] << " ";
	  cout << geometry.grids[photon.grid_number[photon.current_grid_num]].positions[i][photon.position_index[k][i]] << " ";
	  cout << geometry.grids[photon.grid_number[photon.current_grid_num]].positions[i][photon.position_index[k][i]+1] << endl;
	}
#endif
      }
    }
    // now determine if the photon has traveled far enough
    if (fabs(target_tau - tau_traveled) < ROUNDOFF_ERR_TRIG) {
      if (!escape_grid) exit_cell = 0;
      escape = 1;
    }

#ifdef DEBUG_CDD
    if (photon.number == OUTNUM) {
      cout << "target_tau = " << target_tau << endl;
      cout << "tau_traveled = " << tau_traveled << endl;
      cout << "diff = " << target_tau - tau_traveled << endl;
      cout << "current_grid_num, num_current grids = ";
      cout << photon.current_grid_num << " ";
      cout << photon.num_current_grids << endl;
      cout << "min_index = " << min_index << endl;
      cout << "=============" << endl;
    }
#endif

  } else { // otherwise determine the distance in this cell

    double delta_position[3];  // x,y,z distance to the edge of the cell (parallel to the x,y,z axes)
    double delta_distance[3];  // distance the photon would need to travelel to exit the cell in the x,y,z directions

  // determine the x,y,z distances to the edge of the cell 
  //   this is parallel to the x,y,z axes
    int i = 0;
#ifdef DEBUG_CDD
    if (photon.number == OUTNUM) {
      cout << "exit_cell = " << exit_cell << endl;
    }
#endif
    for (i = 0; i < 3; i++) {
      if (photon.dir_cosines[i] > 0.0) {
	delta_position[i] = geometry.grids[photon.grid_number[photon.current_grid_num]].positions[i][photon.position_index[k][i]+1] -
	  photon.position[i];					 
      } else if (photon.dir_cosines[i] < 0.0) {
	delta_position[i] = photon.position[i] - 
	  geometry.grids[photon.grid_number[photon.current_grid_num]].positions[i][photon.position_index[k][i]];
      } else {  // photon.dir_cosines[i] == 0.0
	delta_position[i] = 0.0;
      }
      
      // calculate the distance the photon would  have to travel in each x,y,z direction
      // to exit the cell
      if (photon.dir_cosines[i] != 0.0) 
        delta_distance[i] = fabs(delta_position[i]/photon.dir_cosines[i]);
      else
        delta_distance[i] = 0.0;
#ifdef DEBUG_CDD
      if (photon.number == OUTNUM) {
	cout << "i = " << i << " ";
	cout << "delt pos dir = ";
	cout << delta_position[i] << " ";
	cout << delta_distance[i] << " ";
	cout << photon.dir_cosines[i] << " ";
	cout << geometry.grids[photon.grid_number[photon.current_grid_num]].positions[i][photon.position_index[k][i]+1] << " ";
	cout << photon.position[i] << " ";
	cout << endl;
      }
#endif
    }
  
    // determine the first axis the photon exits along
    //   minumum distance traveled
    double min_dist = 1e20;
    min_index = -1;
    for (i = 0; i < 3; i++) 
      if ((delta_distance[i] < min_dist) && (delta_distance[i] >= 0.) &&
	  (!((photon.dir_cosines[i] == 0.0) && (delta_distance[i] == 0.0)))) {
	min_dist = delta_distance[i];
	min_index = i;
      }

    if (min_index == -1) {
      // find the first non-zero dir cosine to send the photon along
      min_index = 0;
      while (photon.dir_cosines[min_index] == 0.0)
         min_index++;
      if ((min_index == 2) && (photon.dir_cosines[min_index] == 0.0)) {
         cout << "no nonzero distances traveled in this cell." << endl;
         exit(8);
      }
    }

#ifdef DEBUG_CDD
    if (photon.number == OUTNUM) {
      cout << "index transversed = " << min_index << endl;
    }
#endif    
    // save the current cell dust coeff (tau/pc)
    // and distance traveled
    distance_traveled = delta_distance[min_index];
    
    // determine the optical depth transversed before the photon exits the cell
    tau_traveled = dust_tau_per_pc*distance_traveled;
    
    // if this tau_traveled is larger than the amount of tau_left, then decrease the distance traveled
    if (tau_traveled > target_tau) {
      distance_traveled *= target_tau/tau_traveled;
      tau_traveled = target_tau;
      exit_cell = 0;
    }

#ifdef DEBUG_CDD
    if (photon.number == OUTNUM) {
      cout << "***single cell***" << endl;
      cout << "target_tau = " << target_tau << endl;
      cout << "tau_traveled = " << tau_traveled << endl;
      cout << "diff = " << target_tau - tau_traveled << endl;
      cout << "current_grid_num, num_current grids = ";
      cout << photon.current_grid_num << " ";
      cout << photon.num_current_grids << endl;
      cout << "min_index = " << min_index << endl;
    }
#endif
  // update the photon position for the distance_traveled
    for (i = 0; i < 3; i++)
      photon.position[i] += distance_traveled*photon.dir_cosines[i];

  }

  // end single cell, rest should be cell or subgrid compatable

  // update the position_index if the cell is exited
  // and check if the photon has escaped the grid
  // make sure the position is exactly the edge of the cell
  if (exit_cell) {
    if (photon.dir_cosines[min_index] > 0.0) {
      photon.position_index[k][min_index]++;
      photon.position[min_index] = geometry.grids[photon.grid_number[photon.current_grid_num]].positions[min_index][photon.position_index[k][min_index]];
    } else {
      photon.position[min_index] = geometry.grids[photon.grid_number[photon.current_grid_num]].positions[min_index][photon.position_index[k][min_index]];
      photon.position_index[k][min_index]--;
    }

    // now determine if the photon has left the grid
#ifdef DEBUG_CDD
    if (photon.number == OUTNUM) {
      cout << "escape (pre check1/2) = " << escape << endl;
    }
#endif
    if (!escape)
      if ((photon.position_index[k][min_index] >= geometry.grids[photon.grid_number[photon.current_grid_num]].index_dim[min_index]) ||
	  (photon.position_index[k][min_index] < 0))
	escape = 1;

#ifdef DEBUG_CDD
    if (photon.number == OUTNUM) {
      cout << "escape (post check1) = " << escape << endl;
    }
#endif
    
    // determine if the photon has left the model while still inside the grid (-1 > dust_tau_per_pc > 0)
    if (!escape) {
//       if ((photon.position_index[k][min_index] < geometry.grids[photon.grid_number[photon.current_grid_num]].index_dim[min_index]) &&
// 	  (photon.position_index[k][min_index] >= 0)) {
	dust_tau_ref_per_pc = geometry.grids[photon.grid_number[photon.current_grid_num]].
	  grid(photon.position_index[k][0],photon.position_index[k][1],photon.position_index[k][2]).dust_tau_per_pc;
	if ((dust_tau_ref_per_pc > -1.0) && (dust_tau_ref_per_pc < 0.0))
	  escape = 1;
//       }
    }
#ifdef DEBUG_CDD
    if (photon.number == OUTNUM) {
      cout << "escape (post check2) = " << escape << endl;
    }
#endif

  }

#ifdef DEBUG_CDD
  if (photon.number == OUTNUM) {
    cout << "dist & tau traveled = ";
    cout << distance_traveled << " ";
    cout << tau_traveled << endl;
    cout << "tau target = " << target_tau << endl;
    cout << "pos = ";
    for (di = 0; di < 3; di++)
      cout << photon.position[di] << " ";
    cout << endl;
    cout << "pos idx = ";
    for (di = 0; di < 3; di++)
      cout << photon.position_index[k][di] << " ";
    cout << endl;
    cout << "exit_cell = " << exit_cell << endl;
    cout << "escape = " << escape << endl;
  }
#endif

  return(distance_traveled);
}
