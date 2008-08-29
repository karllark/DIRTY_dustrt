// ======================================================================
//   Header file for classify stellar photon procedure.  
// Include files and function definitions.
//
// 2003 Jun/KDG - written
// 2004 Dec/KDG - added photon tracking
// 2005 May/KDG - added angular radius
// 2008 Aug/KDG - added double exp disk
// ======================================================================

#ifndef _DIRTY_GEOMETRY_DEF_
#define _DIRTY_GEOMETRY_DEF_

#include "NumUtils.h"
//#include "vect_utils.h"
#include "grid_cell.h"

#define MAX_OBSERVERS 100
#define MAX_MULTIPLE_STARS 100

// types of photon sources (for new_photon.cc)
#define NEW_PHOTON_DISCRETE_STARS 1
#define NEW_PHOTON_DIFFUSE_ISOTROPIC 2
#define NEW_PHOTON_DIFFUSE_FILE 3
#define NEW_PHOTON_GRID 4
#define NEW_PHOTON_DEXP_DISK 5

// info defining a single grid (can be nested inside of other single grids)
struct one_grid {
  int parent_grid_num;
  long index_dim[3];
  double phys_grid_size[3];
  double phys_cube_size[3];
  vector<vector<double> > positions;
  NumUtils::Cube<grid_cell> grid;
};

// info fully defining a geometry
struct geometry_struct {
  vector<one_grid> grids;  // vector of grids
  int max_grid_depth;   // maximum number of nesting grids

  long num_cells; // number of cells

  // absorbed energy storage info
  int abs_energy_storage_type;  // how to store the absorbed energy
                                // 0 = memory, 1 = disk
  int abs_energy_grid_initialized;
  int abs_energy_wave_index;    // wavelength index of absorbed energy grid
                                // 0 if disk is used, otherwise equal to wavelength index

  int emitted_energy_grid_initialized;

  // global dust properties
  //   at a later date, may want to make these cell dependent
  double albedo;   // dust scattering albedo
  double g;        // dust scattering phase function asymmetry (g = <cos(theta)>)

  double total_h_mass; // total hydrogen mass (where dust emission was done)

  int wave_index;  // wavelength index (needed for new_photon_grid_source)

  float distance;  // distance of model from observer
  int num_observers;  // number of different observer directions
  float observer_angles[2][MAX_OBSERVERS]; // theta,phi for each observer

  string source_type;   // type of source (stars, diffuse, etc.)
  double total_source_luminosity;  // in ergs s^-1 A^-1
  int new_photon_source_type; // integer to control how to emit photons

  int num_stars;   // number of stars
  double star_positions[5][MAX_MULTIPLE_STARS];  // position of stars in physical units
          // luminosity is the 4th element, 5th element the running sum of the luminosity
          // between 0 and 1 for the determination of which star emits the current photon

  // dexp stellar variables
  double stellar_scaleheight;
  double stellar_scalelength;
  double stellar_emit_constant_z;
  double stellar_emit_constant_xy;

  int randomize_observer;  // randomize observer position to integrate over 4pi

  // vectors for the diffuse source (ISRF)
  vector<double> diffuse_source_theta;  // in radians
  vector<double> diffuse_source_phi;    // cos(phi)
  vector<double> diffuse_source_intensity;  // assumed to be in ergs s^-1 sr^-1 A^-1 cm^-2
  vector<double> diffuse_source_sum_intensity; // between 0 and 1

  string type;   // type of geometry [controls which routine is called to setup dust grid]
  float radius;   // radius of model
  double angular_radius; // radius of model in angular units
  float tau;      // radial optical depth of model at reference wavelength
  float tau_wave;  // reference wavelength for tau
  float tau_to_tau_ref;  // ratio of current tau to reference tau
  float wavelength;   // wavelenght of current tau
  float max_tau_per_cell;  // the maximum tau per cell (used in subdividing cells)
  float filling_factor;  // filling factor of high density clumps;
  float density_ratio;   // density ratio of low/high clumps
  float clump_densities[2];  // densities of clumps in tau/pc
  float solid_angle;  // solid angle for scattering probabilities

  long n_photons;  // number of photons to run
};

#endif
