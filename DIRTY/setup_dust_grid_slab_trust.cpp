// ======================================================================
// setup a slab of dust - specific for TRUST with custom grid
//
// KDG 28 Feb 2007 - adapted from setup_dust_shell
// KDG 13 Apr 2007 - removed geometry.solid_angle (see setup_dust_grid.cc)
// KDG  8 May 2008 - fixed error in setting up the dust grid physical dimensions
//                   for the case where the xy and z sizes are unequal
// KDG 16 Jun 2008 - fixed max_grid_depth calculation
// KDG 16 Jun 2015 - added reading of different values for x,y,z max cell tau
// KDG 28 Jul 2015 - adapted from standard slab setup for TRUST benchmark
//                   mainly to allow for a log spacing of the z axis
// ======================================================================
#include "setup_dust_grid_slab.h"
// #define DEBUG_SDG

void setup_dust_grid_slab_trust(ConfigFile &param_data, geometry_struct &geometry, random_dirty &random_obj)

{
    // get the parameters needed for a spherical geometry from the ConfigFile

    // model x,y,z size
    float size_xy = param_data.FValue("Geometry", "size_xy");
    check_input_param("size_xy", size_xy, 0.0, 1e38);
    float size_z = param_data.FValue("Geometry", "size_z");
    check_input_param("size_z", size_z, 0.0, 1e38);

    // radius of geometry (needed for output)
    geometry.radius = size_xy / 2.;
    if (size_z / 2. > geometry.radius)
        geometry.radius = size_z / 2.;

    // angular radius needs to be large enough to allow for any rotation and still
    // have all the photons encompassed in the final image
    geometry.angular_radius = atan(1.5 * geometry.radius / (geometry.distance - geometry.radius));
    //  geometry.angular_radius = atan(1.5*geometry.radius/(geometry.distance -
    //  geometry.radius));

    // slab start/end
    float slab_z1 = param_data.FValue("Geometry", "slab_z1");
    check_input_param("slab_z1", slab_z1, -1. * size_z / 2., size_z / 2.);
    float slab_z2 = param_data.FValue("Geometry", "slab_z2");
    check_input_param("slab_z2", slab_z2, -1. * size_z / 2., size_z / 2.);

    // check that slab_z1 < slab_z2
    if (slab_z1 > slab_z2)
    {
        cout << "slab start is at a greater distance than the slab end" << endl;
        cout << "slab z1 = " << slab_z1 << endl;
        cout << "slab z2 = " << slab_z2 << endl;
        exit(8);
    }

    // slab optical depth
    geometry.tau = param_data.FValue("Geometry", "slab_tau");
    check_input_param("slab tau", geometry.tau, 0.0, 1000.);
#ifdef DEBUG_SDG
    cout << "input slab tau = " << geometry.tau << endl;
#endif

    // nonslab density ratio
    float nonslab_density_ratio = param_data.FValue("Geometry", "nonslab_density_ratio");
    //  check_input_param("nonslab density
    //  ratio",nonslab_density_ratio,0.99e-7,1.);
    check_input_param("nonslab density ratio", nonslab_density_ratio, 0.0, 1.);

    // get the number of bins for each axis
    int nbins_x = param_data.IValue("Geometry", "slab_nbins_x");
    check_input_param("slab_nbins_x", nbins_x, 1, 1000);
    int nbins_y = param_data.IValue("Geometry", "slab_nbins_y");
    check_input_param("slab_nbins_y", nbins_y, 1, 1000);
    int nbins_z = param_data.IValue("Geometry", "slab_nbins_z");
    check_input_param("slab_nbins_z", nbins_z, 1, 1000);

    // set the maximum grid depth
    geometry.max_grid_depth = 1;

    // declare main grid
    one_grid main_grid;

    // setup size of main grid
    main_grid.index_dim[0] = nbins_x;
    main_grid.index_dim[1] = nbins_y;
    main_grid.index_dim[2] = nbins_z + 1; // add one for the non-slab region

    // fill position arrays with the physical dimensions
    vector<double> x_pos(main_grid.index_dim[0] + 1);
    vector<double> y_pos(main_grid.index_dim[1] + 1);
    vector<double> z_pos(main_grid.index_dim[2] + 1);
    int i;

    // x values
    for (i = 0; i <= main_grid.index_dim[0]; i++)
        x_pos[i] = double(i) * (size_xy) / double(main_grid.index_dim[0]) - size_xy / 2.;

    // y values
    for (i = 0; i <= main_grid.index_dim[1]; i++)
        y_pos[i] = double(i) * (size_xy) / double(main_grid.index_dim[1]) - size_xy / 2.;

    // z values (linear)
    //   float deltaz = slab_z2 - slab_z1;
    //   for (i = 0; i <= main_grid.index_dim[2]-1; i++)
    //     z_pos[i] = slab_z1 +
    //     double(i)*(deltaz)/double(main_grid.index_dim[2]-1);
    //   z_pos[main_grid.index_dim[2]] = size_z/2.;

    // z values (do as a log10 spacing)
    // vector<double> z_pos_log(main_grid.index_dim[2]+1);
    float log_slab_z1 = log10(fabs(slab_z1));
    float log_slab_z2 = log10(fabs(slab_z2));
    float log_deltaz = (log_slab_z2 - log_slab_z1) / double(main_grid.index_dim[2] - 1);
    for (i = 0; i <= (main_grid.index_dim[2] - 1); i++)
        z_pos[i] = -1.0 * pow(10.0, log_slab_z1 + double(i) * log_deltaz);
    z_pos[main_grid.index_dim[2]] = size_z / 2.;

    // for (i = 0; i <= (main_grid.index_dim[2]); i++)
    //   cout << i << " " << z_pos[i]  << endl;
    // exit(0);

    // add position arrays to main grid
    main_grid.positions.push_back(x_pos);
    main_grid.positions.push_back(y_pos);
    main_grid.positions.push_back(z_pos);

    // give the size of the grid
    main_grid.phys_grid_size[0] = size_xy;
    main_grid.phys_grid_size[1] = size_xy;
    main_grid.phys_grid_size[2] = size_z;

    // give the size of the cubes
    // main_grid.phys_cube_size[0] =
    // main_grid.phys_grid_size[0]/main_grid.index_dim[0];
    // main_grid.phys_cube_size[1] =
    // main_grid.phys_grid_size[1]/main_grid.index_dim[1];
    // main_grid.phys_cube_size[2] =
    // main_grid.phys_grid_size[2]/main_grid.index_dim[2];

    // allocate main grid
    main_grid.grid.CSize(main_grid.index_dim[0], main_grid.index_dim[1], main_grid.index_dim[2]);

    // densities (in tau/pc)
    float slab_density = geometry.tau / (slab_z2 - slab_z1);
    float nonslab_density = nonslab_density_ratio * slab_density;

    int j, k;
    float z_val = 0.0;
    for (k = 0; k < main_grid.index_dim[2]; k++)
    {
        z_val = (main_grid.positions[2][k] + main_grid.positions[2][k + 1]) / 2.0;
        for (j = 0; j < main_grid.index_dim[1]; j++)
        {
            // y_val = (main_grid.positions[1][j] + main_grid.positions[1][j+1])/2.0;
            for (i = 0; i < main_grid.index_dim[0]; i++)
            {
                // x_val = (main_grid.positions[0][i] +
                // main_grid.positions[0][i+1])/2.0;
                //  slab position should be a negative z if behind stars (observer at
                //  positive z)
                if ((z_val >= slab_z1) && (z_val <= slab_z2))
                {
                    main_grid.grid(i, j, k).dust_tau_per_pc = slab_density;
                }
                else
                    main_grid.grid(i, j, k).dust_tau_per_pc = nonslab_density;
            }
        }
    }

    // identify this grid as main grid
    main_grid.parent_grid_num = -1;

    // add main grid to grids vector
    geometry.grids.push_back(main_grid);

    geometry.max_tau_per_cell_x = 1e6;
    geometry.max_tau_per_cell_y = 1e6;
    geometry.max_tau_per_cell_z = 1e6;
    geometry.filling_factor = 0.0;
    geometry.density_ratio = 1.0;
    int spherical_clumps = 0;

    // subdivide all overdense cells
    setup_dust_grid_subdivide_overdense_cells(geometry, spherical_clumps);
}
