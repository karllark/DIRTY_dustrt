// ======================================================================
// setup a double exponential clumpy dust grid
//   (think spiral galaxy)
//
// KDG Jun 2008 - adapted from setup_dust_shell
// ======================================================================
#include "setup_dust_grid_dexp_disk.h"
// #define DEBUG_SDGDD

void
setup_dust_grid_dexp_disk (ConfigFile &param_data, geometry_struct &geometry,
                           random_dirty &random_obj)

{
  // get the parameters needed for a spherical geometry from the ConfigFile

  // model radius
  geometry.radius = param_data.FValue ("Geometry", "radius");
  check_input_param ("radius", geometry.radius, 0.0, 1e38);

  // angular radius needs to be large enough to allow for any rotation and still
  // have all the photons encompassed in the final image
  geometry.angular_radius = atan (1.1 * geometry.radius / geometry.distance);

  // dust scalelength
  double dust_scalelength = param_data.FValue ("Geometry", "dust_scalelength");
  check_input_param ("dust_scalelength", dust_scalelength, 0.0, geometry.radius);

  // dust scaleheight
  double dust_scaleheight = param_data.FValue ("Geometry", "dust_scaleheight");
  check_input_param ("dust_scaleheight", dust_scaleheight, 0.0, geometry.radius);

  // dust vertical truncation
  double dust_vertical_trunc = param_data.FValue ("Geometry", "dust_vertical_trunc");
  check_input_param ("dust_vertical_trunc", dust_vertical_trunc, 0.0, geometry.radius);

  // stellar scalelength
  geometry.stellar_scalelength = param_data.FValue ("Geometry", "stellar_scalelength");
  check_input_param ("stellar_scalelength", geometry.stellar_scalelength, 0.0, geometry.radius);

  // stellar scaleheight
  geometry.stellar_scaleheight = param_data.FValue ("Geometry", "stellar_scaleheight");
  check_input_param ("stellar_scaleheight", geometry.stellar_scaleheight, 0.0, geometry.radius);

  // constant needed for z calculation
  geometry.stellar_emit_constant_z
      = 1.0 / (1.0 - exp (-1.0 * dust_vertical_trunc / geometry.stellar_scaleheight));

  // radial optical depth
  geometry.tau = param_data.FValue ("Geometry", "tau");
  check_input_param ("tau", geometry.tau, 0.0, 1000.);
#ifdef DEBUG_SDGDD
  cout << "input tau = " << geometry.tau << endl;
#endif

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
      //     cout << "new filling factor = " << geometry.filling_factor << endl;
    }

  // min size of grid on one side (assume xy plane)
  float min_grid_size = param_data.FValue ("Geometry", "min_grid_size");
  check_input_param ("min_grid_size", min_grid_size, 0, geometry.radius);

  // starting size of grid on one side (assume xy plane)
  float start_grid_size = param_data.FValue ("Geometry", "start_grid_size");
  check_input_param ("start_grid_size", start_grid_size, 0, geometry.radius);

  // constant needed for xy calculation (not used as went to vector method)
  //   geometry.stellar_emit_constant_xy = 1.0/(1.0 -
  //   exp(-1.0*geometry.radius/geometry.stellar_scalelength)*
  // 					   ((geometry.radius/geometry.stellar_scalelength)
  // + 1.0)); vector of values for xy calculation
  geometry.stellar_emit_n_xy = 100 * int (geometry.radius / min_grid_size) + 1;
  geometry.stellar_emit_xy_vals.resize (geometry.stellar_emit_n_xy);
  int ii = 0;
  double tmp_min_rad = 0.0;
  double tmp_rad = 0.0;
  double tmp_sum = 0.0;
  geometry.stellar_emit_xy_vals[0] = 0.0;
  for (ii = 1; ii < geometry.stellar_emit_n_xy; ii++)
    {
      tmp_min_rad
          = ((double (ii) - 1.0) / double (geometry.stellar_emit_n_xy - 1)) * geometry.radius;
      tmp_rad = ((double (ii)) / double (geometry.stellar_emit_n_xy - 1)) * geometry.radius;
      geometry.stellar_emit_xy_vals[ii] = exp (-1.0 * tmp_rad / geometry.stellar_scalelength);
      geometry.stellar_emit_xy_vals[ii] *= (tmp_rad * tmp_rad - tmp_min_rad * tmp_min_rad);
      tmp_sum += geometry.stellar_emit_xy_vals[ii];
      geometry.stellar_emit_xy_vals[ii] = tmp_sum;
    }
  for (ii = 0; ii < geometry.stellar_emit_n_xy; ii++)
    geometry.stellar_emit_xy_vals[ii] /= tmp_sum;

  // set the maximum grid depth
  geometry.max_grid_depth = 1;

  // declare main grid
  one_grid main_grid;

  // determine the main grid size
  //  equal in xy, smaller in z (unless vertical truncation = radius of model)
  int xy_grid_size = 2 * int (geometry.radius / start_grid_size);
  int z_grid_size = 2 * int (dust_vertical_trunc / start_grid_size);
  // make sure this defines a cube
  //   if (((0.5*xy_grid_size*start_grid_size) != geometry.radius) ||
  //       ((0.5*z_grid_size*start_grid_size) != dust_vertical_trunc)) {
  //     cout << "non cubical cells required -> not allowed (maybe check if this
  //     could be allowed?)" << endl; cout << "start_grid_size [pc] = " <<
  //     start_grid_size << endl; cout << "dust_vertical_trunc [pc] = " <<
  //     dust_vertical_trunc << endl; cout << "radius [pc] = " <<
  //     geometry.radius << endl; cout << "xy_grid_size = " << xy_grid_size <<
  //     endl; cout << "xy_grid_size (real) = "
  //     << 2.*geometry.radius/start_grid_size << endl; cout << "z_grid_size = "
  //     << z_grid_size << endl; cout << "z_grid_size (real) = "
  //     << 2.*dust_vertical_trunc/start_grid_size << endl; exit(8);
  //   }
  // arbitrarily make the z a finer grid
  float subdiv_z = param_data.FValue ("Geometry", "subdivide_z_factor");
  if (isnan (subdiv_z))
    subdiv_z = 1;
  check_input_param ("subdivide_z_factor", subdiv_z, 1, 10);
  z_grid_size *= subdiv_z;
  //   cout << subdiv_z << endl;
  //   exit(8);

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
  vector<double> x_pos (main_grid.index_dim[0] + 1);
  vector<double> y_pos (main_grid.index_dim[1] + 1);
  vector<double> z_pos (main_grid.index_dim[2] + 1);
  int i;
  float tmp_val;
  // do z
  for (i = 0; i <= main_grid.index_dim[2]; i++)
    {
      tmp_val = double (i) * (2.0 * dust_vertical_trunc) / double (main_grid.index_dim[2])
                - dust_vertical_trunc;
      z_pos[i] = tmp_val;
    }
  // do xy
  for (i = 0; i <= main_grid.index_dim[0]; i++)
    {
      tmp_val = double (i) * (2.0 * geometry.radius) / double (main_grid.index_dim[0])
                - geometry.radius;
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
  main_grid.positions.push_back (x_pos);
  main_grid.positions.push_back (y_pos);
  main_grid.positions.push_back (z_pos);

  // give the size of the grid
  main_grid.phys_grid_size[0] = 2.0 * geometry.radius;
  main_grid.phys_grid_size[1] = main_grid.phys_grid_size[0];
  main_grid.phys_grid_size[2] = 2.0 * dust_vertical_trunc;

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
      = (geometry.filling_factor + geometry.density_ratio * (1.0 - geometry.filling_factor));
  geometry.clump_densities[1] = geometry.density_ratio * geometry.clump_densities[0];

  int j, k;

  double tmp_density = 0.0;
  double tmp_den_constant
      = geometry.tau
        / (dust_scaleheight * (1.0 - exp (-1. * dust_vertical_trunc / dust_scaleheight)));
  float x_val = 0.0;
  float y_val = 0.0;
  float z_val = 0.0;
  for (k = 0; k < main_grid.index_dim[2]; k++)
    {
      z_val = (main_grid.positions[2][k] + main_grid.positions[2][k + 1]) / 2.0;
      for (j = 0; j < main_grid.index_dim[1]; j++)
        {
          y_val = (main_grid.positions[1][j] + main_grid.positions[1][j + 1]) / 2.0;
          for (i = 0; i < main_grid.index_dim[0]; i++)
            {
              x_val = (main_grid.positions[0][i] + main_grid.positions[0][i + 1]) / 2.0;
              float xy_val = sqrt (x_val * x_val + y_val * y_val);
              tmp_density = tmp_den_constant * exp (-1.0 * abs (z_val) / dust_scaleheight)
                            * exp (-1.0 * xy_val / dust_scalelength);
              if ((xy_val > geometry.radius) || (fabs (z_val) > dust_vertical_trunc))
                {
                  main_grid.grid (i, j, k).dust_tau_per_pc = 1e-20;
                }
              else if (random_obj.random_num () <= geometry.filling_factor)
                main_grid.grid (i, j, k).dust_tau_per_pc
                    = tmp_density * geometry.clump_densities[0];
              else
                main_grid.grid (i, j, k).dust_tau_per_pc
                    = tmp_density * geometry.clump_densities[1];
              // 	if ((i == main_grid.index_dim[0]/2) && (k ==
              // main_grid.index_dim[2]/2)) { 	  cout << i << " " << j << " " << k << " :
              // "; 	  cout << main_grid.grid(i,j,k).dust_tau_per_pc << " "; 	  cout <<
              // x_val << " " << y_val << " " << xy_val << endl;
              // 	}
            }
        }
    }

  // identify this grid as main grid
  main_grid.parent_grid_num = -1;

  // add main grid to grids vector
  geometry.grids.push_back (main_grid);

  // subdivide if the min_grid_size < start_grid_size
  if (min_grid_size < start_grid_size)
    {

      int subdivide = 0;
      int subdivide_any = 0;
      int m = 0; // only main_grid for now
      int cur_subgrid_num = int (geometry.grids.size ());
      float start_min_frac = start_grid_size / min_grid_size;
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
                  float xy_val = sqrt (x_val * x_val + y_val * y_val);

                  float z_frac = fabs (z_val) / dust_scaleheight;
                  float xy_frac = fabs (xy_val) / dust_scalelength;

                  int poss_index = 0;
                  float frac_sub = z_frac;
                  if ((z_frac < start_min_frac) && (xy_frac < start_min_frac)
                      && (geometry.grids[m].grid (i, j, k).dust_tau_per_pc > 0.0))
                    {
                      if (xy_frac < z_frac)
                        frac_sub = xy_frac;
                      poss_index = int (start_min_frac / frac_sub);

                      if (poss_index > 1)
                        {
#ifdef DEBUG_SDGDD
                          cout << "subdivide ";
                          cout << z_frac << " " << xy_frac << endl;
                          cout << i << " " << j << " " << k << endl;
#endif
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
                      if (poss_index > int (start_min_frac))
                        poss_index = int (start_min_frac);
                      subgrid.index_dim[0] = poss_index;
#ifdef DEBUG_SDGDD
                      cout << "subgrid x cube size = "
                           << geometry.grids[m].phys_grid_size[0] / float (poss_index) << endl;
                      cout << "subgird x cubes (1 dim) = " << subgrid.index_dim[0] << endl;
#endif

                      subgrid.index_dim[1] = subgrid.index_dim[0];
                      subgrid.index_dim[2] = int (float (subgrid.index_dim[0]) * subdiv_z);
                      subgrid.index_dim[0] = 2;
                      subgrid.index_dim[1] = 2;

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

                      int l;
#ifdef DEBUG_SDGDD
                      cout << cur_subgrid_num << " ";
                      for (l = 0; l < 3; l++)
                        cout << subgrid.phys_grid_size[l] << " ";
                      cout << endl;
#endif

                      // ensure that the first position exactly equals the edge of the
                      // parent cell
                      //   no roundoff error problems (hopefully)
                      x_subpos[0] = geometry.grids[m].positions[0][i];
                      y_subpos[0] = geometry.grids[m].positions[1][j];
                      z_subpos[0] = geometry.grids[m].positions[2][k];

                      for (l = 1; l < subgrid.index_dim[0]; l++)
                        x_subpos[l]
                            = geometry.grids[m].positions[0][i]
                              + (double (l) / subgrid.index_dim[0]) * subgrid.phys_grid_size[0];

                      for (l = 1; l < subgrid.index_dim[1]; l++)
                        y_subpos[l]
                            = geometry.grids[m].positions[1][j]
                              + (double (l) / subgrid.index_dim[1]) * subgrid.phys_grid_size[1];

                      for (l = 1; l < subgrid.index_dim[2]; l++)
                        z_subpos[l]
                            = geometry.grids[m].positions[2][k]
                              + (double (l) / subgrid.index_dim[2]) * subgrid.phys_grid_size[2];

                      // ensure that the last position exactly equals the edge of the
                      // parent cell
                      //   no roundoff error problems (hopefully)
                      x_subpos[subgrid.index_dim[0]] = geometry.grids[m].positions[0][i + 1];
                      y_subpos[subgrid.index_dim[1]] = geometry.grids[m].positions[1][j + 1];
                      z_subpos[subgrid.index_dim[2]] = geometry.grids[m].positions[2][k + 1];

                      subgrid.positions.push_back (x_subpos);
                      subgrid.positions.push_back (y_subpos);
                      subgrid.positions.push_back (z_subpos);

                      subgrid.grid.CSize (subgrid.index_dim[0], subgrid.index_dim[1],
                                          subgrid.index_dim[2]);

                      float tz_val, ty_val, tx_val = 0.;
                      int n, o;
                      for (o = 0; o < subgrid.index_dim[2]; o++)
                        {
                          tz_val = (subgrid.positions[2][o] + subgrid.positions[2][o + 1]) / 2.0;
                          for (n = 0; n < subgrid.index_dim[1]; n++)
                            {
                              ty_val
                                  = (subgrid.positions[1][n] + subgrid.positions[1][n + 1]) / 2.0;
                              for (l = 0; l < subgrid.index_dim[0]; l++)
                                {
                                  tx_val = (subgrid.positions[0][l] + subgrid.positions[0][l + 1])
                                           / 2.0;
                                  float txy_val = sqrt (tx_val * tx_val + ty_val * ty_val);
                                  tmp_density = tmp_den_constant
                                                * exp (-1.0 * abs (tz_val) / dust_scaleheight)
                                                * exp (-1.0 * txy_val / dust_scalelength);
                                  if ((txy_val > geometry.radius)
                                      || (fabs (tz_val) > dust_vertical_trunc))
                                    {
                                      subgrid.grid (l, n, o).dust_tau_per_pc = 1e-20;
                                    }
                                  else
                                    {
                                      if (random_obj.random_num () <= geometry.filling_factor)
                                        subgrid.grid (l, n, o).dust_tau_per_pc
                                            = tmp_density * geometry.clump_densities[0];
                                      else
                                        subgrid.grid (l, n, o).dust_tau_per_pc
                                            = tmp_density * geometry.clump_densities[1];
                                    }
                                }
                            }
                        }

                      // setup ties to parent grid
                      subgrid.parent_grid_num = m;
                      geometry.grids[m].grid (i, j, k).dust_tau_per_pc = -cur_subgrid_num;
                      cur_subgrid_num++;

                      geometry.grids.push_back (subgrid);
                    }
                }
            }
        }
      if (subdivide_any)
        geometry.max_grid_depth++;
    }

  // subdivide all overdense cells
  setup_dust_grid_subdivide_overdense_cells (geometry, spherical_clumps);

  //******************************************************************
  // test that the geometry has the right vertical optical depth
  // ***does not work*** argh!

  //   photon_data tmp_photon;

  //   // setup the weights
  //   tmp_photon.number = 0;
  //   tmp_photon.stellar_weight = 1.0;
  //   tmp_photon.scat_weight = 0.0;

  //   // initialize statistics variables
  //   tmp_photon.num_scat = 0;

  //   // direction of tmp_photon; assuming an isotropic source
  //   // in direction cosines...
  //   tmp_photon.dir_cosines[0] = 0.0;
  //   tmp_photon.dir_cosines[1] = 0.0;
  //   tmp_photon.dir_cosines[2] = 1.0;

  //   // direction of photon; assuming an isotropic source
  //   // in direction cosines...
  //   double phi = M_PI*(2.0*random_obj.random_num() - 1.0);
  //   tmp_photon.dir_cosines[2] = 2.0*random_obj.random_num() - 1.0;
  //   double temp = sqrt(1.0 - pow(tmp_photon.dir_cosines[2],2));
  //   tmp_photon.dir_cosines[0] = cos(phi)*temp;
  //   tmp_photon.dir_cosines[1] = sin(phi)*temp;

  //   // set the tmp_photon position to the star position
  //   for (i = 0; i < 3; i++) {
  //     tmp_photon.position[i] = 1.0;
  //     tmp_photon.birth_position[i] = tmp_photon.position[i];
  //   }

  //   // setup the photon structure with positions for the maximum grid depth
  //   vector<int> one_index(3);
  //   for (i = 0; i < geometry.max_grid_depth; i++) {
  //     tmp_photon.position_index.push_back(one_index);
  //     tmp_photon.grid_number.push_back(long(0));
  //   }

  //   // add in the initialization for the photon path variables
  //   tmp_photon.path_max_cells = 10;
  //   vector<int> tmp_index(tmp_photon.path_max_cells);
  //   for (i = 0; i < 4; i++)
  //     tmp_photon.path_pos_index.push_back(tmp_index);
  //   tmp_photon.path_tau.resize(tmp_photon.path_max_cells);

  //   // now determine the position indexes of the tmp_photon
  //   determine_photon_position_index_initial(geometry, tmp_photon);

  //   // now move the tmp_photon to the edge
  //   double target_tau = 1e20;
  //   int escape = 0;
  //   double distance_traveled = 0.0;
  //   double tau_to_surface = 0.0;
  //   tmp_photon.path_cur_cells = 0;  // set to -1 *not* to save cells
  //   tranversed

  //   distance_traveled = calc_photon_trajectory(tmp_photon, geometry,
  //   target_tau, escape, tau_to_surface);

  //   for (i = 0; i < tmp_photon.path_cur_cells; i++) {
  //     cout << tmp_photon.path_tau[i] << " ";
  //   }
  //   cout << endl;

  //   cout << "optical depth to surface = " << tau_to_surface << endl;
  //   cout << "escape = " << escape << endl;
  //  exit(8);
  //******************************************************************
}
