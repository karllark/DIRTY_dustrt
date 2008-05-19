// ======================================================================
//   Procedure to generate a photon and set up the stats on it.          
// This procedure is used for photons emitted from the grid (needed for 
// dust emission cases).
//
// 2007 Nov/KDG - written
// 2008 Mar/KDG - added setting of birth position (forgot it)
// ======================================================================
#include "new_photon_grid_source.h"
//#define DEBUG_NPGS

void new_photon_grid_source (photon_data& photon,
			     geometry_struct& geometry,
			     runinfo_struct& runinfo,
			     random_dirty random_obj)
  

{
  // setup the weights
  photon.stellar_weight = 1.0;
  photon.scat_weight = 0.0;

  // initialize statistics variables
  photon.num_scat = 0;

  // direction of photon; assuming an isotropic source
  // in direction cosines...
  double phi = M_PI*(2.0*random_obj.random_num() - 1.0);
  photon.dir_cosines[2] = 2.0*random_obj.random_num() - 1.0;
  double temp = sqrt(1.0 - pow(photon.dir_cosines[2],2));
  photon.dir_cosines[0] = cos(phi)*temp;
  photon.dir_cosines[1] = sin(phi)*temp;

  int i = 0;
#ifdef DEBUG_NPGS
  for (i = 0; i < 3; i++)
    cout << photon.dir_cosines[i] << " ";
  cout << "starting dir cosines" << endl;
#endif

  // determine where the photon is emitted from
  // loop through the grids and determine which grid using a random number
  double ran_val = random_obj.random_num();
  int grid_num = -1;
  uint m = 0;
  int x = geometry.wave_index;

  // check if the we are at the beginning of the grid
  if (ran_val <= geometry.grids[m].grid(geometry.grids[m].index_dim[0]-1,
				       geometry.grids[m].index_dim[1]-1,
				       geometry.grids[m].index_dim[2]-1).emitted_energy[0][x])
    grid_num = 0;
 
#ifdef DEBUG_NPGS
 cout << "grid_num (1st check) = " << grid_num << endl;
#endif
  while ((grid_num == -1) && (m < geometry.grids.size())) {
    m++;
#ifdef DEBUG_NPGS
    cout << "(2nd check) m = " << m << endl;
#endif
    if ((ran_val > geometry.grids[m-1].grid(geometry.grids[m-1].index_dim[0]-1,geometry.grids[m-1].index_dim[1]-1,
					    geometry.grids[m-1].index_dim[2]-1).emitted_energy[0][x]) &&
	ran_val <= geometry.grids[m].grid(geometry.grids[m].index_dim[0]-1,geometry.grids[m].index_dim[1]-1,
					    geometry.grids[m].index_dim[2]-1).emitted_energy[0][x])
      grid_num = m ;
  } 
  
#ifdef DEBUG_NPGS
  cout << "grid_num = " << grid_num << endl;
#endif
  if (grid_num == -1) {
    cout << "Can't emitted photon from grid, random number of " << ran_val << " larger than anything in grid." << endl;
    exit(8);
  }

  // then go through the grid and find the cell
  // then randomly distribute the photon in the cell

  // find z plane
  int z_val = -1;
  int k = 0;

  // check if the we are in the first z plane
  if (ran_val <= geometry.grids[grid_num].grid(geometry.grids[grid_num].index_dim[0]-1,
					       geometry.grids[grid_num].index_dim[1]-1,
					       k).emitted_energy[0][x])

    z_val = 0;
  
  while ((z_val == -1) && (k < geometry.grids[grid_num].index_dim[2]-1)) {
    k++;
    if ((ran_val > geometry.grids[grid_num].grid(geometry.grids[grid_num].index_dim[0]-1,
						 geometry.grids[grid_num].index_dim[1]-1,k-1).emitted_energy[0][x]) &&
	(ran_val <= geometry.grids[grid_num].grid(geometry.grids[grid_num].index_dim[0]-1,
						  geometry.grids[grid_num].index_dim[1]-1,k).emitted_energy[0][x]))
      z_val = k;
  }

  if (z_val == -1) {
    cout << "Can't emitted photon from grid, random number of " << ran_val << " larger than anything z dim of designated grid." << endl;
    exit(8);
  }

  // find y row
  int y_val = -1;
  int j = 0;

  // check if the we are in the first y column
  if (ran_val <= geometry.grids[grid_num].grid(geometry.grids[grid_num].index_dim[0]-1,
					       j,z_val).emitted_energy[0][x])

    y_val = 0;

  while ((y_val == -1) && (j < geometry.grids[grid_num].index_dim[1]-1)) {
    j++;
    if ((ran_val > geometry.grids[grid_num].grid(geometry.grids[grid_num].index_dim[0]-1,
						 j-1,z_val).emitted_energy[0][x]) &&
	(ran_val <= geometry.grids[grid_num].grid(geometry.grids[grid_num].index_dim[0]-1,
						  j,z_val).emitted_energy[0][x]))
      y_val = j;
  }
  if (y_val == -1) {
    cout << "Can't emitted photon from grid, random number of " << ran_val << " larger than anything y dim of designated grid." << endl;
    exit(8);
  }

  // find x pixel
  i = 0;
  int x_val = -1;

  // check if the we are in the first y column
  if (ran_val <= geometry.grids[grid_num].grid(i,y_val,z_val).emitted_energy[0][x])
    x_val = 0;

  while ((x_val == -1) && (i < geometry.grids[grid_num].index_dim[0]-1)) {
    i++;
    if ((ran_val > geometry.grids[grid_num].grid(i-1,y_val,z_val).emitted_energy[0][x]) &&
	(ran_val <= geometry.grids[grid_num].grid(i,y_val,z_val).emitted_energy[0][x]))
      x_val = i;
  }
  if (x_val == -1) {
    cout << "Can't emitted photon from grid, random number of " << ran_val << " larger than anything x dim of designated grid." << endl;
    exit(8);
  }

#ifdef DEBUG_NPGS
  cout << "ran val & grid_cell vals = ";
  cout << x_val << " ";
  cout << y_val << " ";
  cout << z_val << " ";
  cout << grid_num << " ";
  cout << ran_val << " ";
  cout << geometry.grids[grid_num].grid(x_val,y_val,z_val).emitted_energy[0][x] - ran_val << " ";
  if (x_val > 0) 
    cout << geometry.grids[grid_num].grid(x_val-1,y_val,z_val).emitted_energy[0][x] << " ";
  else
    cout << 0.0 << endl;
  cout << geometry.grids[grid_num].grid(x_val,y_val,z_val).emitted_energy[0][x] << endl;
#endif

  // determine where in the cell to emit the photon (random)
  photon.position[0] = geometry.grids[grid_num].positions[0][x_val] +
    random_obj.random_num()*(geometry.grids[grid_num].positions[0][x_val+1] - geometry.grids[grid_num].positions[0][x_val]);
  photon.position[1] = geometry.grids[grid_num].positions[1][y_val] +
    random_obj.random_num()*(geometry.grids[grid_num].positions[1][y_val+1] - geometry.grids[grid_num].positions[1][y_val]);
  photon.position[2] = geometry.grids[grid_num].positions[2][z_val] +
    random_obj.random_num()*(geometry.grids[grid_num].positions[2][z_val+1] - geometry.grids[grid_num].positions[2][z_val]);

  for (i = 0; i < 3; i++)
    photon.birth_position[i] = photon.position[i];

  // save the birth grain/emission type probabilities if needed for dust thermal emission part of dirty
  if (runinfo.dust_thermal_emission && runinfo.do_emission_grain)
    for (i = 0; i < runinfo.n_emission_grain_types; i++) {
      photon.birth_photon_type_prob[i] = geometry.grids[grid_num].grid(x_val,y_val,z_val).emitted_energy[i][geometry.wave_index];
//       cout << "i = " << i << "; " << photon.birth_photon_type_prob[i] << endl;
    }


  // now determine the position indexes of the photon
  determine_photon_position_index_initial(geometry, photon);
  
#ifdef DEBUG_NPGS
  cout << "done with npgs." << endl;
#endif

}
