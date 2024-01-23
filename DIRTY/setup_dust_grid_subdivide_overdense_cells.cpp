// ======================================================================
// Subdivide overdense cells in the grid
//
// KDG 15 May 2008 - Written (taken from setup_dust_grid_shell)
// KDG 16 Jun 2008 - fixed max_grid_depth calculation & fixed assignment of grid
// number to parent KDG 14 Jul 2008 - fixed error of interating over recently
// created subgrids KDG 17 Jun 2015 - modified code to all for different
// subdivision for different axes
// ======================================================================
#include "setup_dust_grid_subdivide_overdense_cells.h"
// #define DEBUG_SDGSOC

void setup_dust_grid_subdivide_overdense_cells(geometry_struct& geometry,
                                               int spherical_clumps)

{
  int i, j, k, m = 0;

  // loop over all existing grids and subgrid any overdense cells
  float x_tau, y_tau, z_tau = 0.0;
  int cur_subgrid_num = int(geometry.grids.size());
  int subdivide = 0;
  int subdivide_any = 0;

  geometry.num_cells = 0;

  long num_cells_orig = 0;
  long num_cells_subdivide = 0;

  int n_grids = int(geometry.grids.size());
  for (m = 0; m < n_grids; m++) {
    for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
      for (j = 0; j < geometry.grids[m].index_dim[1]; j++)
        for (i = 0; i < geometry.grids[m].index_dim[0]; i++) {
          num_cells_orig++;

          x_tau = (geometry.grids[m].positions[0][i + 1] -
                   geometry.grids[m].positions[0][i]) *
                  geometry.grids[m].grid(i, j, k).dust_tau_per_pc;
          y_tau = (geometry.grids[m].positions[1][j + 1] -
                   geometry.grids[m].positions[1][j]) *
                  geometry.grids[m].grid(i, j, k).dust_tau_per_pc;
          z_tau = (geometry.grids[m].positions[2][k + 1] -
                   geometry.grids[m].positions[2][k]) *
                  geometry.grids[m].grid(i, j, k).dust_tau_per_pc;

          subdivide = 0;
          if (x_tau > geometry.max_tau_per_cell_x) subdivide = 1;
          if (y_tau > geometry.max_tau_per_cell_y) subdivide = 1;
          if (z_tau > geometry.max_tau_per_cell_z) subdivide = 1;
          if ((spherical_clumps) &&
              (geometry.grids[0].grid(i, j, k).dust_tau_per_pc ==
               geometry.clump_densities[0]))
            subdivide = 1;

          if (subdivide) {
            num_cells_subdivide++;

#ifdef DEBUG_SDGSOC
            cout << m << ",";
            cout << i << ",";
            cout << j << ",";
            cout << k << " needs a subgrid; ";
            cout << "max_tau_x = " << geometry.max_tau_per_cell_x << endl;
            cout << "cell tau_x = " << x_tau << endl;
#endif

            // note that we've subdivided at least one grid cell
            subdivide_any = 1;

            one_grid subgrid;
            if (spherical_clumps) {
              subgrid.index_dim[0] = 10;  // make a sphere
              subgrid.index_dim[1] = subgrid.index_dim[0];
              subgrid.index_dim[2] = subgrid.index_dim[0];
            } else {
              if (x_tau > geometry.max_tau_per_cell_x)
                subgrid.index_dim[0] =
                    int(x_tau / float(geometry.max_tau_per_cell_x)) + 1;
              else
                subgrid.index_dim[0] = 1;
              if (y_tau > geometry.max_tau_per_cell_y)
                subgrid.index_dim[1] =
                    int(y_tau / float(geometry.max_tau_per_cell_y)) + 1;
              else
                subgrid.index_dim[1] = 1;
              if (z_tau > geometry.max_tau_per_cell_z)
                subgrid.index_dim[2] =
                    int(z_tau / float(geometry.max_tau_per_cell_z)) + 1;
              else
                subgrid.index_dim[2] = 1;

              //  	      cout << x_tau/float(geometry.max_tau_per_cell_x)
              //  << " ";
              // 	      cout << y_tau/float(geometry.max_tau_per_cell_y)
              // << " "; 	      cout << z_tau/float(geometry.max_tau_per_cell_z) <<
              // endl;

              //  	      cout <<
              //  int(x_tau/float(geometry.max_tau_per_cell_x)) << " ";
              // 	      cout <<
              // int(y_tau/float(geometry.max_tau_per_cell_y)) << " "; 	      cout <<
              // int(z_tau/float(geometry.max_tau_per_cell_z)) << endl;

              // 	      cout << subgrid.index_dim[0] << " ";
              // 	      cout << subgrid.index_dim[1] << " ";
              // 	      cout << subgrid.index_dim[2] << endl;

              // 	      cout << x_tau << " ";
              // 	      cout << y_tau << " ";
              // 	      cout << z_tau << endl;

              // 	      cout << geometry.max_tau_per_cell_x << " ";
              // 	      cout << geometry.max_tau_per_cell_y << " ";
              // 	      cout << geometry.max_tau_per_cell_z << endl;
              // 	      exit(8);
            }

            vector<double> x_subpos(subgrid.index_dim[0] + 1);
            vector<double> y_subpos(subgrid.index_dim[1] + 1);
            vector<double> z_subpos(subgrid.index_dim[2] + 1);

            subgrid.phys_grid_size[0] = geometry.grids[m].positions[0][i + 1] -
                                        geometry.grids[m].positions[0][i];
            subgrid.phys_grid_size[1] = geometry.grids[m].positions[1][j + 1] -
                                        geometry.grids[m].positions[1][j];
            subgrid.phys_grid_size[2] = geometry.grids[m].positions[2][k + 1] -
                                        geometry.grids[m].positions[2][k];

            // subgrid.phys_cube_size[0] =
            // subgrid.phys_grid_size[0]/subgrid.index_dim[0];
            // subgrid.phys_cube_size[1] =
            // subgrid.phys_grid_size[1]/subgrid.index_dim[1];
            // subgrid.phys_cube_size[2] =
            // subgrid.phys_grid_size[2]/subgrid.index_dim[2];

            int l;
            // ensure that the first position exactly equals the edge of the
            // parent cell
            //   no roundoff error problems (hopefully)
            x_subpos[0] = geometry.grids[m].positions[0][i];
            y_subpos[0] = geometry.grids[m].positions[1][j];
            z_subpos[0] = geometry.grids[m].positions[2][k];

            for (l = 1; l < subgrid.index_dim[0]; l++)
              x_subpos[l] = geometry.grids[m].positions[0][i] +
                            (double(l) / subgrid.index_dim[0]) *
                                subgrid.phys_grid_size[0];

            for (l = 1; l < subgrid.index_dim[1]; l++)
              y_subpos[l] = geometry.grids[m].positions[1][j] +
                            (double(l) / subgrid.index_dim[1]) *
                                subgrid.phys_grid_size[1];

            for (l = 1; l < subgrid.index_dim[2]; l++)
              z_subpos[l] = geometry.grids[m].positions[2][k] +
                            (double(l) / subgrid.index_dim[2]) *
                                subgrid.phys_grid_size[2];

            // ensure that the last position exactly equals the edge of the
            // parent cell
            //   no roundoff error problems (hopefully)
            x_subpos[subgrid.index_dim[0]] =
                geometry.grids[m].positions[0][i + 1];
            y_subpos[subgrid.index_dim[1]] =
                geometry.grids[m].positions[1][j + 1];
            z_subpos[subgrid.index_dim[2]] =
                geometry.grids[m].positions[2][k + 1];

            subgrid.positions.push_back(x_subpos);
            subgrid.positions.push_back(y_subpos);
            subgrid.positions.push_back(z_subpos);

            subgrid.grid.CSize(subgrid.index_dim[0], subgrid.index_dim[1],
                               subgrid.index_dim[2]);

            float dust_tau_per_pc =
                geometry.grids[m].grid(i, j, k).dust_tau_per_pc;
            float index_radius;  // radius of subgrid position in index values
            int n, o;
            for (o = 0; o < subgrid.index_dim[2]; o++)
              for (n = 0; n < subgrid.index_dim[1]; n++)
                for (l = 0; l < subgrid.index_dim[0]; l++) {
                  geometry.num_cells++;
                  if (spherical_clumps) {
                    index_radius =
                        sqrt(pow(float(o) - subgrid.index_dim[2] / 2., 2.0) +
                             pow(float(n) - subgrid.index_dim[1] / 2., 2.0) +
                             pow(float(l) - subgrid.index_dim[0] / 2., 2.0));
                    if (index_radius < subgrid.index_dim[0] / 2.)
                      subgrid.grid(l, n, o).dust_tau_per_pc = dust_tau_per_pc;
                    else
                      subgrid.grid(l, n, o).dust_tau_per_pc =
                          dust_tau_per_pc * geometry.density_ratio;
                  } else {
                    subgrid.grid(l, n, o).dust_tau_per_pc = dust_tau_per_pc;
                  }
                }

            // setup ties to parent grid
            subgrid.parent_grid_num = m;
            geometry.grids[m].grid(i, j, k).dust_tau_per_pc = -cur_subgrid_num;
            cur_subgrid_num++;

            geometry.grids.push_back(subgrid);
          } else {
            geometry.num_cells++;
            // calculate number of H atoms
            // 	    if (geometry.grids[m].grid(i,j,k).dust_tau_per_pc > 0.0)
            // 	      geometry.grids[m].grid(i,j,k).num_H =
            // num_H_const*vol*geometry.grids[m].grid(i,j,k).dust_tau_per_pc*geometry.tau_to_tau_ref;
            // 	    else
            geometry.grids[m].grid(i, j, k).num_H = 0.0;
          }
        }
  }

  if (subdivide_any) geometry.max_grid_depth++;

  cout << "Result of optional subdivide overdense cells" << endl;
  cout << "total number of cells = " << geometry.num_cells << endl;
  cout << "number of orig cells = " << num_cells_orig << endl;
  cout << "number of orig subdiv cells = " << n_grids << endl;
  cout << "number of new subdiv cells = " << num_cells_subdivide << endl;
}
