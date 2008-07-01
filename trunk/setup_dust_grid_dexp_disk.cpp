// ======================================================================
// setup a double exponential clumpy dust grid
//   (think spiral galaxy)
//
// KDG Jun 2008 - adapted from setup_dust_shell
// ======================================================================
#include "setup_dust_grid_dexp_disk.h"
#define DEBUG_SDGDD

void setup_dust_grid_dexp_disk (ConfigFile& param_data,
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

  // dust scalelength
  double dust_scalelength = param_data.FValue("Geometry","dust_scalelength");
  check_input_param("dust_scalelength",dust_scalelength,0.0,geometry.radius);

  // dust scaleheight
  double dust_scaleheight = param_data.FValue("Geometry","dust_scaleheight");
  check_input_param("dust_scaleheight",dust_scaleheight,0.0,geometry.radius);
  
  // dust vertical truncation
  double dust_vertical_trunc = param_data.FValue("Geometry","dust_vertical_trunc");
  check_input_param("dust_vertical_trunc",dust_vertical_trunc,0.0,geometry.radius);
  
  // radial optical depth
  geometry.tau = param_data.FValue("Geometry","tau");
  check_input_param("tau",geometry.tau,0.0,1000.);
#ifdef DEBUG_SDGDD
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

  // min size of grid on one side (assume xy plane)
  float min_grid_size = param_data.FValue("Geometry","min_grid_size");
  check_input_param("min_grid_size",min_grid_size,0,geometry.radius);

  // starting size of grid on one side (assume xy plane)
  float start_grid_size = param_data.FValue("Geometry","start_grid_size");
  check_input_param("start_grid_size",start_grid_size,0,geometry.radius);

  // set the maximum grid depth
  geometry.max_grid_depth = 1;

  // declare main grid
  one_grid main_grid;

  // determine the main grid size
  //  equal in xy, smaller in z (unless vertical truncation = radius of model)
  int xy_grid_size = 2*int(geometry.radius/start_grid_size);
  int z_grid_size = 2*int(dust_vertical_trunc/start_grid_size);
  // make sure this defines a cube
  if (((0.5*xy_grid_size*start_grid_size) != geometry.radius) ||
      ((0.5*z_grid_size*start_grid_size) != dust_vertical_trunc)) {
    cout << "non cubical cells required -> not allowed (maybe check if this could be allowed?)" << endl;
    cout << "start_grid_size [pc] = " << start_grid_size << endl;
    cout << "dust_vertical_trunc [pc] = " << dust_vertical_trunc << endl;
    cout << "radius [pc] = " << geometry.radius << endl;
    cout << "xy_grid_size = " << xy_grid_size << endl;
    cout << "xy_grid_size (real) = " << 2.*geometry.radius/start_grid_size << endl;
    cout << "z_grid_size = " << z_grid_size << endl;
    cout << "z_grid_size (real) = " << 2.*dust_vertical_trunc/start_grid_size << endl;
    exit(8);
  }
#ifdef DEBUG_SDGDD
  cout << "grid size = " << xy_grid_size << " ";
  cout << xy_grid_size << " ";
  cout << z_grid_size << endl;
#endif
  
  // setup size of main grid
  main_grid.index_dim[0] = xy_grid_size;
  main_grid.index_dim[1] = xy_grid_size;
  main_grid.index_dim[2] = z_grid_size;

  // fill position arrays with
  vector<double> x_pos(main_grid.index_dim[0]+1);
  vector<double> y_pos(main_grid.index_dim[1]+1);
  vector<double> z_pos(main_grid.index_dim[2]+1);
  int i;
  float tmp_val;
  // do z
  for (i = 0; i <= main_grid.index_dim[2]; i++) {
    tmp_val = double(i)*(2.0*dust_vertical_trunc)/double(main_grid.index_dim[2]) - dust_vertical_trunc;
    z_pos[i] = tmp_val;
  }
  // do xy
  for (i = 0; i <= main_grid.index_dim[0]; i++) {
    tmp_val = double(i)*(2.0*geometry.radius)/double(main_grid.index_dim[0]) - geometry.radius;
    x_pos[i] = tmp_val;
    y_pos[i] = tmp_val;
#ifdef DEBUG_SDGDD
    cout << "xyz grid position = ";
    cout << x_pos[i] << " ";
    cout << y_pos[i] << " ";
    if (i <= main_grid.index_dim[2])
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
  main_grid.phys_grid_size[2] = 2.0*dust_vertical_trunc;

  // give the size of the cubes
  main_grid.phys_cube_size[0] = main_grid.phys_grid_size[0]/main_grid.index_dim[0];
  main_grid.phys_cube_size[1] = main_grid.phys_grid_size[1]/main_grid.index_dim[1];
  main_grid.phys_cube_size[2] = main_grid.phys_grid_size[2]/main_grid.index_dim[2];

  // allocate main grid
  main_grid.grid.CSize(main_grid.index_dim[0],main_grid.index_dim[1],main_grid.index_dim[2]);

  // fill main grid with dust density 
  // (now these are just the factors and the density is imposed later)
  geometry.clump_densities[0] = (geometry.filling_factor + geometry.density_ratio*(1.0 - geometry.filling_factor));
  geometry.clump_densities[1] = geometry.density_ratio*geometry.clump_densities[0];

  int j,k;

  double tmp_density = 0.0;
//   double tmp_den_constant = (geometry.tau*(radial_density_poly+1.))/
//     (pow(outer_radius,radial_density_poly+1.) - pow(inner_radius,radial_density_poly+1.));
  double radius = 0.0;
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
// 	  tmp_density = tmp_den_constant*pow(radius,radial_density_poly);
	if (radius < geometry.radius) {
	  if (random_obj.random_num() <= geometry.filling_factor)
	    main_grid.grid(i,j,k).dust_tau_per_pc = tmp_density*geometry.clump_densities[0];
	  else
	    main_grid.grid(i,j,k).dust_tau_per_pc = tmp_density*geometry.clump_densities[1];
	  
// #ifdef DEBUG_SDGDD
// 	  cout << main_grid.grid(i,j,k).dust_tau_per_pc << " ";
// #endif
	} else
	  main_grid.grid(i,j,k).dust_tau_per_pc = -0.5;  // this means the edge of the dust
      }
    }
  }
  
  // identify this grid as main grid
  main_grid.parent_grid_num = -1;
  
  // add main grid to grids vector
  geometry.grids.push_back(main_grid);
  
//   // now subdivide the center cubes if the inner radius is inside the cube dimension
//   if (subdivide_radius > 0.0) {
//     if (subdivide_radius < pow(3.,0.5)*(main_grid.phys_cube_size[0]/2.)) subdivide_radius = main_grid.phys_cube_size[0];
//     int subdivide = 0;
//     int subdivide_any = 0;
//     int m = 0;  // only main_grid for now
//     int cur_subgrid_num = int(geometry.grids.size());
//     for (k = 0; k < geometry.grids[m].index_dim[2]; k++) {
//       z_val = (geometry.grids[m].positions[2][k] + geometry.grids[m].positions[2][k+1])/2.0;
//       for (j = 0; j < geometry.grids[m].index_dim[1]; j++) {
// 	y_val = (geometry.grids[m].positions[1][j] + geometry.grids[m].positions[1][j+1])/2.0;
// 	for (i = 0; i < geometry.grids[m].index_dim[0]; i++) {
// 	  x_val = (geometry.grids[m].positions[0][i] + geometry.grids[m].positions[0][i+1])/2.0;
// 	  radius = sqrt(x_val*x_val + y_val*y_val + z_val*z_val);
// 	  if (radius <= subdivide_radius) {
// 	    cout << "inner radius = " << inner_radius;
// 	    cout << "; radius = " << radius << "; ";
// 	    cout << i << " ";
// 	    cout << j << " ";
// 	    cout << k << endl;
// 	    subdivide = 1;
// 	    subdivide_any = 1;
// 	  } else subdivide = 0;

// 	  if (subdivide) {
// 	    one_grid subgrid;
// 	    subgrid.index_dim[0] = int(2.*geometry.grids[0].phys_cube_size[0]/inner_radius) + 1;
// #ifdef DEBUG_SDG
// 	    cout << "subgird size = " << subgrid.index_dim[0] << endl;
// #endif
	    
// 	    subgrid.index_dim[1] = subgrid.index_dim[0];
// 	    subgrid.index_dim[2] = subgrid.index_dim[0];
	    
// 	    vector<double> x_subpos(subgrid.index_dim[0]+1);
// 	    vector<double> y_subpos(subgrid.index_dim[1]+1);
// 	    vector<double> z_subpos(subgrid.index_dim[2]+1);
	    
// 	    subgrid.phys_grid_size[0] = geometry.grids[m].positions[0][i+1] - geometry.grids[m].positions[0][i];
// 	    subgrid.phys_grid_size[1] = geometry.grids[m].positions[1][j+1] - geometry.grids[m].positions[1][j];
// 	    subgrid.phys_grid_size[2] = geometry.grids[m].positions[2][k+1] - geometry.grids[m].positions[2][k];
	    
// 	    subgrid.phys_cube_size[0] = subgrid.phys_grid_size[0]/subgrid.index_dim[0];
// 	    subgrid.phys_cube_size[1] = subgrid.phys_grid_size[1]/subgrid.index_dim[1];
// 	    subgrid.phys_cube_size[2] = subgrid.phys_grid_size[2]/subgrid.index_dim[2];

// 	    int l;
// 	    for (l = 0; l <= subgrid.index_dim[0]; l++) {
// 	      x_subpos[l] = geometry.grids[m].positions[0][i] + (double(l)/subgrid.index_dim[0])*subgrid.phys_grid_size[0];
// 	      y_subpos[l] = geometry.grids[m].positions[1][j] + (double(l)/subgrid.index_dim[1])*subgrid.phys_grid_size[1];
// 	      z_subpos[l] = geometry.grids[m].positions[2][k] + (double(l)/subgrid.index_dim[2])*subgrid.phys_grid_size[2];
// 	    }
// 	    subgrid.positions.push_back(x_subpos);
// 	    subgrid.positions.push_back(y_subpos);
// 	    subgrid.positions.push_back(z_subpos);

// 	    subgrid.grid.CSize(subgrid.index_dim[0],subgrid.index_dim[1],subgrid.index_dim[2]);

// 	    double tradius;  // radius of subgrid position in index values
// 	    float tz_val,ty_val,tx_val = 0.;
// 	    int n,o;
// 	    for (o = 0; o < subgrid.index_dim[2]; o++) {
// 	      tz_val = (subgrid.positions[2][o] + subgrid.positions[2][o+1])/2.0;
// 	      for (n = 0; n < subgrid.index_dim[1]; n++) {
// 		ty_val = (subgrid.positions[2][n] + subgrid.positions[2][n+1])/2.0;
// 		for (l = 0; l < subgrid.index_dim[0]; l++) {
// 		  tx_val = (subgrid.positions[2][l] + subgrid.positions[2][l+1])/2.0;
// 		  tradius = sqrt(tx_val*tx_val + ty_val*ty_val + tz_val*tz_val);
// 		  tmp_density = tmp_den_constant*pow(tradius,radial_density_poly);
// 		  if (tradius < inner_radius)
// 		    subgrid.grid(l,n,o).dust_tau_per_pc = geometry.grids[m].grid(i,j,k).dust_tau_per_pc*0.;
// 		  else
// 		    subgrid.grid(l,n,o).dust_tau_per_pc = tmp_density*geometry.clump_densities[0];
// 		}
// 	      }
// 	    }

// 	    // setup ties to parent grid
// 	    subgrid.parent_grid_num = m;
// 	    geometry.grids[m].grid(i,j,k).dust_tau_per_pc = -cur_subgrid_num;
// 	    cur_subgrid_num++;

// 	    geometry.grids.push_back(subgrid);
// 	  }
// 	}
//       }
//     }
//     if (subdivide_any)
//       geometry.max_grid_depth++;
//   }

  // subdivide all overdense cells
  setup_dust_grid_subdivide_overdense_cells(geometry, spherical_clumps);

  exit(8);
}
