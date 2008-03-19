// ======================================================================
//   Header file for geometry_def.h (grid cell definition)  
// Include files and function definitions.
//
// 2007 Mar/KDG - changed absorbed energy from float to double
// 2008 Mar/KDG - added num_H to save the number of H atoms per cell
// ======================================================================
#ifndef _DIRTY_GRID_CELL_
#define _DIRTY_GRID_CELL_

// defines the contents of a grid cell
struct grid_cell {
  float dust_tau_per_pc;      // dust optical depth per parsec tau/pc
  vector<double> absorbed_energy;  // energy absorbed in this cell at this wavelength
                                   // can handle multiwavelengths if desired
  double save_radiation_field_density;    // save the existing radiation field density (needed for the self-absorption calculation)
  vector< vector<double> > emitted_energy;  // energy emitted in this cell as a function of wavelength
                                            // with ability to have emitted energy split into components
  double num_H;  // number of HI atoms in this cell
};

#endif
