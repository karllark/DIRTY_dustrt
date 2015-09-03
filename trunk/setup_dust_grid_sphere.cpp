// ======================================================================
// setup a spherical clumpy dust grid
//
// KDG 30 Jun 2006 - split spherical part from main setup_dust_grid code
// KDG 13 Apr 2007 - removed geometry.solid_angle (see setup_dust_grid.cc)
// KDG 16 Jun 2008 - fixed max_grid_depth calculation
// ======================================================================
#include "setup_dust_grid_sphere.h"

void setup_dust_grid_sphere (ConfigFile& param_data,
			     geometry_struct& geometry,
			     random_dirty& random_obj)

{
  // get the parameters needed for a spherical geometry from the ConfigFile

  // model radius
  geometry.radius = param_data.FValue("Geometry","radius");
  check_input_param("radius",geometry.radius,0.0,1e38);

  // angular radius needs to be large enough to allow for any rotation and still
  // have all the photons encompassed in the final image
  geometry.angular_radius = atan(1.45*geometry.radius/geometry.distance);
  
  // radial optical depth
  geometry.tau = param_data.FValue("Geometry","tau");
  check_input_param("tau",geometry.tau,0.0,1000.);
#ifdef DEBUG_SDG
  cout << "input tau = " << geometry.tau << endl;
#endif
  // maximum optical depth per cell (controls when a cell is subdivided)
  geometry.max_tau_per_cell = param_data.FValue("Geometry","max_tau_per_cell");
  check_input_param("max_tau_per_cell",geometry.max_tau_per_cell,0.0,1e10);

  // filling factor of high density material
  geometry.filling_factor = param_data.FValue("Geometry","filling_factor");
  check_input_param("filling_factor",geometry.filling_factor,0.0,1.0);

  // ratio of low/high density clumps
  geometry.density_ratio = param_data.FValue("Geometry","density_ratio");
  check_input_param("density_ratio",geometry.density_ratio,0.0,1.0);

  // spherical or cubical clumps
  string clump_type = param_data.SValue("Geometry","clump_type");
  int spherical_clumps = 0;
  if (clump_type == "sphere") {
    spherical_clumps = 1;
    // adjust filling factor to account of a sphere inscribed in a cube
    geometry.filling_factor /= (4./3.)*M_PI*pow(0.5,3.0);
//     cout << "new filling factor = " << geometry.filling_factor << endl;
  }

  // size of grid on one side (all sides equal)
  int grid_size = param_data.IValue("Geometry","grid_size");
  check_input_param("grid_size",grid_size,0,1000);

  // set the maximum grid depth
  geometry.max_grid_depth = 1;

  // declare main grid
  one_grid main_grid;

  // setup size of main grid
  main_grid.index_dim[0] = grid_size;
  main_grid.index_dim[1] = grid_size;
  main_grid.index_dim[2] = grid_size;

  // fill position arrays with
  vector<double> x_pos(main_grid.index_dim[0]+1);
  vector<double> y_pos(main_grid.index_dim[1]+1);
  vector<double> z_pos(main_grid.index_dim[2]+1);
  int i;
  double tmp_val;
  for (i = 0; i <= main_grid.index_dim[0]; i++) {
    tmp_val = double(i)*(2.0*geometry.radius)/double(main_grid.index_dim[0]) - geometry.radius;
    x_pos[i] = tmp_val;
    y_pos[i] = tmp_val;
    z_pos[i] = tmp_val;
#ifdef DEBUG_SDG
    cout << "xyz grid position = ";
    cout << x_pos[i] << " ";
    cout << y_pos[i] << " ";
    cout << z_pos[i] << " ";
    cout << "; i = " << i << endl;
#endif
  }

  // add position arrays to main grid
  main_grid.positions.push_back(x_pos);
  main_grid.positions.push_back(y_pos);
  main_grid.positions.push_back(z_pos);

  // give the size of the grid
  main_grid.phys_grid_size[0] = 2.0*geometry.radius;
  main_grid.phys_grid_size[1] = main_grid.phys_grid_size[0];
  main_grid.phys_grid_size[2] = main_grid.phys_grid_size[0];

  // give the size of the cubes
  // main_grid.phys_cube_size[0] = main_grid.phys_grid_size[0]/main_grid.index_dim[0];
  // main_grid.phys_cube_size[1] = main_grid.phys_grid_size[1]/main_grid.index_dim[1];
  // main_grid.phys_cube_size[2] = main_grid.phys_grid_size[2]/main_grid.index_dim[2];

  // allocate main grid
  main_grid.grid.CSize(main_grid.index_dim[0],main_grid.index_dim[1],main_grid.index_dim[2]);

  // fill main grid with dust density
  geometry.clump_densities[0] = geometry.tau/
    (geometry.radius*(geometry.filling_factor + geometry.density_ratio*(1.0 - geometry.filling_factor)));
  geometry.clump_densities[1] = geometry.density_ratio*geometry.clump_densities[0];

  int j,k;
  float radius = 0.0;
  float x_val = 0.0;
  float y_val = 0.0;
  float z_val = 0.0;
  for (k = 0; k < main_grid.index_dim[2]; k++) {
    z_val = (main_grid.positions[2][k] + main_grid.positions[2][k+1])/2.0;
    for (j = 0; j < main_grid.index_dim[1]; j++) {
      y_val = (main_grid.positions[1][j] + main_grid.positions[1][j+1])/2.0;
      for (i = 0; i < main_grid.index_dim[0]; i++) {
	x_val = (main_grid.positions[0][i] + main_grid.positions[0][i+1])/2.0;
	radius = sqrt(x_val*x_val + y_val*y_val + z_val*z_val);
	if (radius <= geometry.radius)
	  if (random_obj.random_num() <= geometry.filling_factor)
	    main_grid.grid(i,j,k).dust_tau_per_pc = geometry.clump_densities[0];
	  else
	    main_grid.grid(i,j,k).dust_tau_per_pc = geometry.clump_densities[1];
	else
	  main_grid.grid(i,j,k).dust_tau_per_pc = -0.5;  // this means the edge of the dust

	// temp
// 	if (k == 3) main_grid.grid(i,j,k).dust_tau_per_pc = 100.0;

//	main_grid.grid(i,j,k).absorbed_energy = 0.0;
      }
    }
  }

  // identify this grid as main grid
  main_grid.parent_grid_num = -1;

  // add main grid to grids vector
  geometry.grids.push_back(main_grid);

  // subdivide all overdense cells
  setup_dust_grid_subdivide_overdense_cells(geometry, spherical_clumps);
}
