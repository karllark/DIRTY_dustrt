// ======================================================================
// setup a spherical clumpy dust grid
//
// KDG 30 Jun 2006 - split spherical part from main setup_dust_grid code
// KDG 13 Apr 2007 - removed geometry.solid_angle (see setup_dust_grid.cc)
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
  geometry.max_grid_depth = 2;

  // declare main grid
  one_grid main_grid;

  // setup size of main grid
  main_grid.index_dim[0] = grid_size;
  main_grid.index_dim[1] = grid_size;
  main_grid.index_dim[2] = grid_size;

  // fill position arrays with
  vector<float> x_pos(main_grid.index_dim[0]+1);
  vector<float> y_pos(main_grid.index_dim[1]+1);
  vector<float> z_pos(main_grid.index_dim[2]+1);
  int i;
  float tmp_val;
  for (i = 0; i <= main_grid.index_dim[0]; i++) {
    tmp_val = float(i)*(2.0*geometry.radius)/float(main_grid.index_dim[0]) - geometry.radius;
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
  main_grid.phys_cube_size[0] = main_grid.phys_grid_size[0]/main_grid.index_dim[0];
  main_grid.phys_cube_size[1] = main_grid.phys_grid_size[1]/main_grid.index_dim[1];
  main_grid.phys_cube_size[2] = main_grid.phys_grid_size[2]/main_grid.index_dim[2];

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

  // loop through the main_grid (grids[0]) and subgrid any overdense cells
  float x_tau = 0.0;
  int cur_subgrid_num = 1;
  int subdivide = 0;
  for (k = 0; k < main_grid.index_dim[2]; k++)
    for (j = 0; j < main_grid.index_dim[1]; j++)
      for (i = 0; i < main_grid.index_dim[0]; i++) {
	x_tau = (geometry.grids[0].positions[0][i+1] - geometry.grids[0].positions[0][i])*
	  geometry.grids[0].grid(i,j,k).dust_tau_per_pc;

	subdivide = 0;
	if (x_tau > geometry.max_tau_per_cell) subdivide = 1;
	if ((spherical_clumps) && (geometry.grids[0].grid(i,j,k).dust_tau_per_pc == geometry.clump_densities[0])) subdivide = 1;

	if (subdivide) {
// 	  cout << i << ",";
// 	  cout << j << ",";
// 	  cout << k << " needs a subgrid ";
// 	  cout << x_tau << " ";

	  one_grid subgrid;
	  if (x_tau > geometry.max_tau_per_cell) 
	    subgrid.index_dim[0] = int(x_tau/geometry.max_tau_per_cell) + 1;
	  else if (spherical_clumps) {
	    subgrid.index_dim[0] = 10;  // make a sphere
	  }

	  subgrid.index_dim[1] = subgrid.index_dim[0];
	  subgrid.index_dim[2] = subgrid.index_dim[0];
// 	  cout << subgrid.index_dim[0] << endl;

	  vector<float> x_subpos(subgrid.index_dim[0]+1);
	  vector<float> y_subpos(subgrid.index_dim[1]+1);
	  vector<float> z_subpos(subgrid.index_dim[2]+1);

	  subgrid.phys_grid_size[0] = geometry.grids[0].positions[0][i+1] - geometry.grids[0].positions[0][i];
	  subgrid.phys_grid_size[1] = geometry.grids[0].positions[1][j+1] - geometry.grids[0].positions[1][j];
	  subgrid.phys_grid_size[2] = geometry.grids[0].positions[2][k+1] - geometry.grids[0].positions[2][k];

	  subgrid.phys_cube_size[0] = subgrid.phys_grid_size[0]/subgrid.index_dim[0];
	  subgrid.phys_cube_size[1] = subgrid.phys_grid_size[1]/subgrid.index_dim[1];
	  subgrid.phys_cube_size[2] = subgrid.phys_grid_size[2]/subgrid.index_dim[2];

	  int m;
	  for (m = 0; m <= subgrid.index_dim[0]; m++) {
	    x_subpos[m] = geometry.grids[0].positions[0][i] + (float(m)/subgrid.index_dim[0])*subgrid.phys_grid_size[0];
	    y_subpos[m] = geometry.grids[0].positions[1][j] + (float(m)/subgrid.index_dim[1])*subgrid.phys_grid_size[1];
	    z_subpos[m] = geometry.grids[0].positions[2][k] + (float(m)/subgrid.index_dim[2])*subgrid.phys_grid_size[2];
// 	    cout << "size = " << x_subpos[m] << " ";
// 	    cout << y_subpos[m] << " ";
// 	    cout << z_subpos[m] << endl;
	  }
	  subgrid.positions.push_back(x_subpos);
	  subgrid.positions.push_back(y_subpos);
	  subgrid.positions.push_back(z_subpos);

	  subgrid.grid.CSize(subgrid.index_dim[0],subgrid.index_dim[1],subgrid.index_dim[2]);

	  float dust_tau_per_pc = geometry.grids[0].grid(i,j,k).dust_tau_per_pc/subgrid.index_dim[0];
	  float index_radius;  // radius of subgrid position in index values
	  int n,o;
	  for (o = 0; o < subgrid.index_dim[2]; o++) 
	    for (n = 0; n < subgrid.index_dim[1]; n++) 
	      for (m = 0; m < subgrid.index_dim[0]; m++) {
		index_radius = sqrt(pow(float(o) - subgrid.index_dim[2]/2.,2.0) + pow(float(n) - subgrid.index_dim[1]/2.,2.0) +
				    pow(float(m) - subgrid.index_dim[0]/2.,2.0));
		if (index_radius < subgrid.index_dim[0]/2.)
		  subgrid.grid(m,n,o).dust_tau_per_pc = dust_tau_per_pc;
		else
		  subgrid.grid(m,n,o).dust_tau_per_pc = dust_tau_per_pc*geometry.density_ratio;

// 		subgrid.grid(m,n,o).absorbed_energy = 0.0;
	      }

	  // setup ties to main grid
	  subgrid.parent_grid_num = 0;
	  geometry.grids[0].grid(i,j,k).dust_tau_per_pc = -cur_subgrid_num;
	  cur_subgrid_num++;

	  geometry.grids.push_back(subgrid);
	}
      }

}
