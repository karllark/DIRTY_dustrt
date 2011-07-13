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

// function object for searching grid layer according to the emitted_energy
struct grid_comparator: public std::binary_function<one_grid, double, bool>
{
  int x;
  grid_comparator(int wave_index): x(wave_index) {}
  inline bool operator()(one_grid &grid, double value) const
  {
    return grid.grid.back().emitted_energy[0][x] < value;
  }
};

// function object for searching grid cell according to the emitted_energy
struct cell_comparator: public std::binary_function<grid_cell, double, bool>
{
  int x;
  cell_comparator(int wave_index): x(wave_index) {}
  inline bool operator()(grid_cell &cell, double value) const
  {
    return cell.emitted_energy[0][x] < value;
  }
};

void new_photon_grid_source (photon_data& photon,
			     geometry_struct& geometry,
			     runinfo_struct& runinfo,
			     random_dirty& random_obj)
  

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
  int x = geometry.wave_index;

  // check which grid are we in
  vector<one_grid>::iterator selected_grid;
  if (ran_val <= geometry.grids[0].grid.back().emitted_energy[0][x]) {
    // dedicated test for the main grid because very often we end up here
    selected_grid = geometry.grids.begin();
    grid_num = 0;
  } else {
    // binary search, O(N) => O(lnN)
    selected_grid = lower_bound(geometry.grids.begin() + 1, geometry.grids.end(), ran_val, grid_comparator(x));
    grid_num = selected_grid - geometry.grids.begin();
  }

  // use binary search to check which cell are we in
  NumUtils::Cube<grid_cell>::iterator selected_cell = lower_bound(selected_grid->grid.begin(), selected_grid->grid.end(), ran_val, cell_comparator(x));
  int cell_num = selected_cell - selected_grid->grid.begin();
  int x_val, y_val, z_val;
  selected_grid->grid.get_xyz(cell_num, x_val, y_val, z_val);

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
