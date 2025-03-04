// ======================================================================
// setup a shell clumpy dust grid
//
// KDG 6 Jul 2006 - adapted from setup_dust_sphere
// KDG 13 Apr 2007 - removed geometry.solid_angle (see setup_dust_grid.cc)
// KDG 15 May 2008 - moved the subdivision of grid cells to a separate routine
// KDG 13 Jun 2008 - added polynomial density profile (except n=-1)
// KDG 16 Jun 2008 - fixed max_grid_depth calculation
// ======================================================================
#include "setup_dust_grid_shell.h"
// #define DEBUG_SDG

void
setup_dust_grid_shell (ConfigFile &param_data, geometry_struct &geometry, random_dirty &random_obj)

{
  // get the parameters needed for a spherical geometry from the ConfigFile

  // model radius
  geometry.radius = param_data.FValue ("Geometry", "radius");
  check_input_param ("radius", geometry.radius, 0.0, 1e38);

  // angular radius needs to be large enough to allow for any rotation and
  // still have all the photons encompassed in the final image
  geometry.angular_radius = atan (1.45 * geometry.radius / geometry.distance);

  // shell inner radius
  double very_inner_radius = param_data.FValue ("Geometry", "very_inner_radius");
  check_input_param ("very_inner_radius", very_inner_radius, 0.0, geometry.radius);

  // shell inner radius
  double inner_radius = param_data.FValue ("Geometry", "inner_radius");
  check_input_param ("inner_radius", inner_radius, 0.0, geometry.radius);

  // shell outer radius
  double outer_radius = param_data.FValue ("Geometry", "outer_radius");
  check_input_param ("outer_radius", outer_radius, 0.0, geometry.radius);

  // shell subdivide radius
  float subdivide_radius = param_data.FValue ("Geometry", "subdivide_radius");
  check_input_param ("subdivide_radius", subdivide_radius, 0.0, geometry.radius);

  // check the outer_radius is greater than the inner radius
  if (outer_radius <= inner_radius)
    {
      cout << "outer shell radius is less than or equal to the inner shell "
              "radius"
           << endl;
      cout << "inner radius = " << inner_radius << endl;
      cout << "outer radius = " << outer_radius << endl;
      exit (8);
    }

  // radial optical depth
  geometry.tau = param_data.FValue ("Geometry", "tau");
  check_input_param ("tau", geometry.tau, 0.0, 1000.);
#ifdef DEBUG_SDG
  cout << "input tau = " << geometry.tau << endl;
#endif

  // polynomial for the shell radial density profile
  double radial_density_poly = param_data.FValue ("Geometry", "shell_density_poly");
  check_input_param ("shell density polynomial", radial_density_poly, -100., 100);
#ifdef DEBUG_SDG
  cout << "poly = " << radial_density_poly << endl;
#endif
  if (radial_density_poly == -1.)
    {
      cout << "shell density polynomial is equal to -1 - new code needed!" << endl;
      exit (8);
    }

  // maximum optical depth per cell (controls when a cell is subdivided)
  geometry.max_tau_per_cell = param_data.FValue ("Geometry", "max_tau_per_cell");
  check_input_param ("max_tau_per_cell", geometry.max_tau_per_cell, 0.0, 1e10);
  geometry.max_tau_per_cell_x = geometry.max_tau_per_cell;
  geometry.max_tau_per_cell_y = geometry.max_tau_per_cell;
  geometry.max_tau_per_cell_z = geometry.max_tau_per_cell;

  // filling factor of high density material
  geometry.filling_factor = param_data.FValue ("Geometry", "filling_factor");
  check_input_param ("filling_factor", geometry.filling_factor, 0.0, 1.0);

  // ratio of low/high density clumps
  geometry.density_ratio = param_data.FValue ("Geometry", "density_ratio");
  check_input_param ("density_ratio", geometry.density_ratio, 0.0, 1.0);

  // spherical or cubical clumps
  string clump_type = param_data.SValue ("Geometry", "clump_type");
  int spherical_clumps = 0;
  if (clump_type == "sphere")
    {
      spherical_clumps = 1;
      // adjust filling factor to account of a sphere inscribed in a cube
      geometry.filling_factor /= (4. / 3.) * M_PI * pow (0.5, 3.0);
      //     cout << "new filling factor = " << geometry.filling_factor <<
      //     endl;
    }

  // size of grid on one side (all sides equal)
  int grid_size = param_data.IValue ("Geometry", "grid_size");
  check_input_param ("grid_size", grid_size, 0, 1000);

  // set the maximum grid depth
  geometry.max_grid_depth = 1;

  // declare main grid
  one_grid main_grid;

  // setup size of main grid
  main_grid.index_dim[0] = grid_size;
  main_grid.index_dim[1] = grid_size;
  main_grid.index_dim[2] = grid_size;

  // fill position arrays with
  vector<double> x_pos (main_grid.index_dim[0] + 1);
  vector<double> y_pos (main_grid.index_dim[1] + 1);
  vector<double> z_pos (main_grid.index_dim[2] + 1);
  int i;
  float tmp_val;
  for (i = 0; i <= main_grid.index_dim[0]; i++)
    {
      tmp_val = double (i) * (2.0 * geometry.radius) / double (main_grid.index_dim[0])
                - geometry.radius;
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
  main_grid.positions.push_back (x_pos);
  main_grid.positions.push_back (y_pos);
  main_grid.positions.push_back (z_pos);

  // give the size of the grid
  main_grid.phys_grid_size[0] = 2.0 * geometry.radius;
  main_grid.phys_grid_size[1] = main_grid.phys_grid_size[0];
  main_grid.phys_grid_size[2] = main_grid.phys_grid_size[0];

  // give the size of the cubes
  // main_grid.phys_cube_size[0] =
  // main_grid.phys_grid_size[0]/main_grid.index_dim[0];
  // main_grid.phys_cube_size[1] =
  // main_grid.phys_grid_size[1]/main_grid.index_dim[1];
  // main_grid.phys_cube_size[2] =
  // main_grid.phys_grid_size[2]/main_grid.index_dim[2];

  // allocate main grid
  main_grid.grid.CSize (main_grid.index_dim[0], main_grid.index_dim[1], main_grid.index_dim[2]);

  // fill main grid with dust density
  // (now these are just the factors and the density is imposed later)
  geometry.clump_densities[0]
      = 1.0 / (geometry.filling_factor + geometry.density_ratio * (1.0 - geometry.filling_factor));
  geometry.clump_densities[1] = geometry.density_ratio * geometry.clump_densities[0];

  int j, k;

  // these computations need to be done with radii less than 1
  // I don't understand why, but I checked this in IDL and get the same results
  // 20 Aug 2008 - KDG (should check this more, but don't know what the next
  // step is)

  // normalized all the radii to the outer radius (see above comment)
  very_inner_radius /= outer_radius;
  inner_radius /= outer_radius;
  //  extended_outer_radius /= outer_radius;
  double save_outer_radius = outer_radius;
  outer_radius = 1.0;

  double tmp_density = 0.0;
  double tmp_den_constant1 = 0.0;
  double tmp_den_constant2 = 0.0;
  if (very_inner_radius == inner_radius)
    {
      tmp_den_constant2 = (geometry.tau * (radial_density_poly + 1.))
                          / (pow (outer_radius, radial_density_poly + 1.)
                             - pow (inner_radius, radial_density_poly + 1.));
    }
  else
    {
      // 1st term
      tmp_den_constant2 = pow (inner_radius, radial_density_poly) * 0.5
                          * (inner_radius * inner_radius - very_inner_radius * very_inner_radius)
                          / (inner_radius - very_inner_radius);
      // 2nd term
      tmp_den_constant2 -= inner_radius * inner_radius * very_inner_radius;
      // 3rd term
      tmp_den_constant2 += (pow (outer_radius, radial_density_poly + 1.)
                            - pow (inner_radius, radial_density_poly + 1.))
                           / (radial_density_poly + 1.);
      tmp_den_constant2 = geometry.tau / tmp_den_constant2;
      tmp_den_constant1 = tmp_den_constant2 * pow (inner_radius, radial_density_poly)
                          / (inner_radius - very_inner_radius);
    }

  //   cout << tmp_den_constant1 << endl;
  //   cout << tmp_den_constant2 << endl;
  //   exit(8);

  double radius = 0.0;
  float x_val = 0.0;
  float y_val = 0.0;
  float z_val = 0.0;
  double min_good_radius = 100000.;
  double max_good_radius = 0.;
  for (k = 0; k < main_grid.index_dim[2]; k++)
    {
      z_val = (main_grid.positions[2][k] + main_grid.positions[2][k + 1]) / 2.0;
      for (j = 0; j < main_grid.index_dim[1]; j++)
        {
          y_val = (main_grid.positions[1][j] + main_grid.positions[1][j + 1]) / 2.0;
          for (i = 0; i < main_grid.index_dim[0]; i++)
            {
              x_val = (main_grid.positions[0][i] + main_grid.positions[0][i + 1]) / 2.0;
              radius = sqrt (x_val * x_val + y_val * y_val + z_val * z_val);
              radius /= save_outer_radius;
              if (radius > geometry.radius / save_outer_radius)
                main_grid.grid (i, j, k).dust_tau_per_pc = -0.5; // this means the edge of the dust
              else if ((radius > outer_radius) && (radius < geometry.radius / save_outer_radius))
                main_grid.grid (i, j, k).dust_tau_per_pc = 1e-20; // negligible amount of dust
              else
                {
                  if (radius < very_inner_radius)
                    tmp_density = 0.0;
                  else if (radius < inner_radius)
                    tmp_density = tmp_den_constant1 * (radius - very_inner_radius);
                  else if (radius <= outer_radius)
                    tmp_density = tmp_den_constant2 * pow (radius, radial_density_poly);
                  tmp_density /= save_outer_radius;

                  if (random_obj.random_num () <= geometry.filling_factor)
                    main_grid.grid (i, j, k).dust_tau_per_pc
                        = tmp_density * geometry.clump_densities[0];
                  else
                    {
                      if (radius < min_good_radius)
                        min_good_radius = radius;
                      if (radius > max_good_radius)
                        max_good_radius = radius;
                      main_grid.grid (i, j, k).dust_tau_per_pc
                          = tmp_density * geometry.clump_densities[1];
                    }
                  // 	  cout << main_grid.grid(i,j,k).dust_tau_per_pc << " ";
                  // 	  cout << radius << " ";
                  // 	  cout << pow(radius,radial_density_poly) << " ";
                  // 	  cout << endl;
                }

#ifdef DEBUG_SDG
              cout << main_grid.grid (i, j, k).dust_tau_per_pc << " ";
#endif
            }
        }
    }
    //   exit(8);

#ifdef DEBUG_SDG
  cout << min_good_radius << " " << inner_radius << endl;
  cout << max_good_radius << " " << outer_radius << endl;
  //  cout << tmp_den_constant2 << " " << tmp_den_constant << endl;
#endif

  //   double tmp_den_constant2 = (geometry.tau*(radial_density_poly+1.))/
  //     (pow(max_good_radius,radial_density_poly+1.) -
  //     pow(min_good_radius,radial_density_poly+1.));

  // check that the density constant is not significantly different
  //   if ((fabs(tmp_den_constant2 - tmp_den_constant)/tmp_den_constant) >
  //   0.01)
  //   {
  //     // if it is different, adjust the numbers in the main grid
  //     double tmp_ratio = tmp_den_constant2/tmp_den_constant;
  //     for (k = 0; k < main_grid.index_dim[2]; k++)
  //       for (j = 0; j < main_grid.index_dim[1]; j++)
  // 	for (i = 0; i < main_grid.index_dim[0]; i++)
  // 	  if (main_grid.grid(i,j,k).dust_tau_per_pc > 0.0)
  // 	    main_grid.grid(i,j,k).dust_tau_per_pc *= tmp_ratio;
  //   }

  // identify this grid as main grid
  main_grid.parent_grid_num = -1;

#ifdef DEBUG_SDG
  cout << "made it past main push" << endl;
  cout << "subdivide_radious = " << subdivide_radius << endl;
#endif

  // add main grid to grids vector
  geometry.grids.push_back (main_grid);

  // useful cube size
  double phys_cube_size = main_grid.phys_grid_size[0] / main_grid.index_dim[0];

  // now subdivide the center cubes if the inner radius is inside the cube
  // dimension
  if (subdivide_radius > 0.0)
    {
      cout << "input subdivide radius = " << subdivide_radius << endl;
      if (subdivide_radius < pow (3., 0.5) * (phys_cube_size / 2.))
        subdivide_radius = phys_cube_size;
      cout << "subdivide radius used (depends on main grid size) = " << subdivide_radius << endl;
      int subdivide = 0;
      int subdivide_any = 0;
      int m = 0; // only main_grid for now
      int cur_subgrid_num = int (geometry.grids.size ());
      int poss_index = 0;
      for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
        {
          z_val = (geometry.grids[m].positions[2][k] + geometry.grids[m].positions[2][k + 1]) / 2.0;
          for (j = 0; j < geometry.grids[m].index_dim[1]; j++)
            {
              y_val = (geometry.grids[m].positions[1][j] + geometry.grids[m].positions[1][j + 1])
                      / 2.0;
              for (i = 0; i < geometry.grids[m].index_dim[0]; i++)
                {
                  x_val
                      = (geometry.grids[m].positions[0][i] + geometry.grids[m].positions[0][i + 1])
                        / 2.0;
                  radius = sqrt (x_val * x_val + y_val * y_val + z_val * z_val);
                  if (radius <= subdivide_radius)
                    {
                      cout << 10. * phys_cube_size
                                  / ((inner_radius - very_inner_radius) * save_outer_radius)
                           << endl;
                      cout << phys_cube_size << endl;
                      cout << (inner_radius - very_inner_radius) * save_outer_radius << endl;
                      poss_index = int (5. * phys_cube_size
                                        / ((inner_radius - very_inner_radius) * save_outer_radius))
                                   + 1;
                      if (poss_index > 1)
                        {
                          subdivide = 1;
                          subdivide_any = 1;
                        }
                      else
                        subdivide = 0;
                    }
                  else
                    subdivide = 0;

                  if (subdivide)
                    {
                      one_grid subgrid;
                      subgrid.index_dim[0] = poss_index;
#ifdef DEBUG_SDG
#endif
                      cout << "very inner radius = " << very_inner_radius * save_outer_radius
                           << " ";
                      cout << "inner radius = " << inner_radius * save_outer_radius;
                      cout << "; radius = " << radius << "; ";
                      cout << i << " ";
                      cout << j << " ";
                      cout << k << endl;
                      // cout << "geometry.grids[0].phys.cube_size[0] = " <<
                      // geometry.grids[0].phys_cube_size[0] << endl;
                      cout << "subgird size = " << subgrid.index_dim[0] << endl;

                      subgrid.index_dim[1] = subgrid.index_dim[0];
                      subgrid.index_dim[2] = subgrid.index_dim[0];

                      vector<double> x_subpos (subgrid.index_dim[0] + 1);
                      vector<double> y_subpos (subgrid.index_dim[1] + 1);
                      vector<double> z_subpos (subgrid.index_dim[2] + 1);

                      subgrid.phys_grid_size[0] = geometry.grids[m].positions[0][i + 1]
                                                  - geometry.grids[m].positions[0][i];
                      subgrid.phys_grid_size[1] = geometry.grids[m].positions[1][j + 1]
                                                  - geometry.grids[m].positions[1][j];
                      subgrid.phys_grid_size[2] = geometry.grids[m].positions[2][k + 1]
                                                  - geometry.grids[m].positions[2][k];

                      // subgrid.phys_cube_size[0] =
                      // subgrid.phys_grid_size[0]/subgrid.index_dim[0];
                      // subgrid.phys_cube_size[1] =
                      // subgrid.phys_grid_size[1]/subgrid.index_dim[1];
                      // subgrid.phys_cube_size[2] =
                      // subgrid.phys_grid_size[2]/subgrid.index_dim[2];

                      // 	    cout << subgrid.phys_cube_size[0] << endl;
                      // 	    cout << very_inner_radius << endl;

                      int l;
                      for (l = 0; l <= subgrid.index_dim[0]; l++)
                        {
                          x_subpos[l]
                              = geometry.grids[m].positions[0][i]
                                + (double (l) / subgrid.index_dim[0]) * subgrid.phys_grid_size[0];
                          y_subpos[l]
                              = geometry.grids[m].positions[1][j]
                                + (double (l) / subgrid.index_dim[1]) * subgrid.phys_grid_size[1];
                          z_subpos[l]
                              = geometry.grids[m].positions[2][k]
                                + (double (l) / subgrid.index_dim[2]) * subgrid.phys_grid_size[2];
                        }
                      subgrid.positions.push_back (x_subpos);
                      subgrid.positions.push_back (y_subpos);
                      subgrid.positions.push_back (z_subpos);

                      subgrid.grid.CSize (subgrid.index_dim[0], subgrid.index_dim[1],
                                          subgrid.index_dim[2]);

                      double tradius; // radius of subgrid position in index
                                      // values
                      float tz_val, ty_val, tx_val = 0.;
                      int n, o;
                      for (o = 0; o < subgrid.index_dim[2]; o++)
                        {
                          tz_val = (subgrid.positions[2][o] + subgrid.positions[2][o + 1]) / 2.0;
                          // 	      cout << "o = " << o << endl;
                          for (n = 0; n < subgrid.index_dim[1]; n++)
                            {
                              ty_val
                                  = (subgrid.positions[1][n] + subgrid.positions[1][n + 1]) / 2.0;
                              for (l = 0; l < subgrid.index_dim[0]; l++)
                                {
                                  tx_val = (subgrid.positions[0][l] + subgrid.positions[0][l + 1])
                                           / 2.0;
                                  tradius
                                      = sqrt (tx_val * tx_val + ty_val * ty_val + tz_val * tz_val);
                                  tradius /= save_outer_radius;
                                  if (tradius < very_inner_radius)
                                    {
                                      tmp_density = 0.0;
                                      subgrid.grid (l, n, o).dust_tau_per_pc = 0.0;
                                    }
                                  else if (tradius < inner_radius)
                                    {
                                      tmp_density
                                          = tmp_den_constant1 * (tradius - very_inner_radius);
                                      //  		    subgrid.grid(l,n,o).dust_tau_per_pc
                                      //  =
                                      //  geometry.grids[m].grid(i,j,k).dust_tau_per_pc;
                                    }
                                  else if (tradius <= outer_radius)
                                    {
                                      tmp_density
                                          = tmp_den_constant2 * pow (tradius, radial_density_poly);
                                      //  		    subgrid.grid(l,n,o).dust_tau_per_pc
                                      //  =
                                      //  geometry.grids[m].grid(i,j,k).dust_tau_per_pc;
                                    }
                                  else if (tradius > outer_radius)
                                    {
                                      tmp_density = 1.0;
                                      subgrid.grid (l, n, o).dust_tau_per_pc = -0.5;
                                    }
                                  tmp_density /= save_outer_radius;

                                  if (subgrid.grid (l, n, o).dust_tau_per_pc != -0.5)
                                    subgrid.grid (l, n, o).dust_tau_per_pc
                                        = tmp_density * geometry.clump_densities[0];
                                  // 		  if
                                  // ((subgrid.grid(l,n,o).dust_tau_per_pc ==
                                  // 0.0)) 		    cout << "0 ";
                                  // else cout <<
                                  // subgrid.grid(l,n,o).dust_tau_per_pc << " /
                                  // ";
                                  //  		    cout << "1 ";
                                }
                              // 		cout << endl;
                            }
                        }

                      // setup ties to parent grid
                      subgrid.parent_grid_num = m;
                      geometry.grids[m].grid (i, j, k).dust_tau_per_pc = -cur_subgrid_num;
                      cur_subgrid_num++;

                      geometry.grids.push_back (subgrid);

                      // 	    exit(8);
                    }
                }
            }
        }
      if (subdivide_any)
        geometry.max_grid_depth++;
    }

#ifdef DEBUG_SDG
  cout << "starting setup_dust_grid_shell" << endl;
#endif

  // subdivide all overdense cells
  setup_dust_grid_subdivide_overdense_cells (geometry, spherical_clumps);

#ifdef DEBUG_SDG
  cout << "finished setup_dust_grid_shell" << endl;
#endif
}
