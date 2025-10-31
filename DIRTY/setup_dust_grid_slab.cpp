// ======================================================================
// setup a shell clumpy dust grid
//
// KDG 28 Feb 2007 - adapted from setup_dust_shell
// KDG 13 Apr 2007 - removed geometry.solid_angle (see setup_dust_grid.cc)
// KDG  8 May 2008 - fixed error in setting up the dust grid physical dimensions
//                   for the case where the xy and z sizes are unequal
// KDG 16 Jun 2008 - fixed max_grid_depth calculation
// KDG 16 Jun 2015 - added reading of different values for x,y,z max cell tau
// ======================================================================
#include "setup_dust_grid_slab.h"
#include "compat.h"
// #define DEBUG_SDG

void setup_dust_grid_slab(ConfigFile &param_data, geometry_struct &geometry, random_dirty &random_obj)

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

    // maximum optical depth per cell (controls when a cell is subdivided)
    geometry.max_tau_per_cell = param_data.FValue("Geometry", "max_tau_per_cell");
    // get if the max_tau_per_cell is going to be input for each dimension if it
    // is the same for all 3
    if (finite(geometry.max_tau_per_cell))
    {
        check_input_param("max_tau_per_cell", geometry.max_tau_per_cell, 0.0, 1e10);
        geometry.max_tau_per_cell_x = geometry.max_tau_per_cell;
        geometry.max_tau_per_cell_y = geometry.max_tau_per_cell;
        geometry.max_tau_per_cell_z = geometry.max_tau_per_cell;
    }
    else
    {
        geometry.max_tau_per_cell_x = param_data.FValue("Geometry", "max_tau_per_cell_x");
        check_input_param("max_tau_per_cell_x", geometry.max_tau_per_cell_x, 0.0, 1e10);

        geometry.max_tau_per_cell_y = param_data.FValue("Geometry", "max_tau_per_cell_y");
        check_input_param("max_tau_per_cell_y", geometry.max_tau_per_cell_y, 0.0, 1e10);

        geometry.max_tau_per_cell_z = param_data.FValue("Geometry", "max_tau_per_cell_z");
        check_input_param("max_tau_per_cell_z", geometry.max_tau_per_cell_z, 0.0, 1e10);
    }

    // cout << geometry.max_tau_per_cell_x << " ";
    // cout << geometry.max_tau_per_cell_y << " ";
    // cout << geometry.max_tau_per_cell_z << endl;

    // filling factor of high density material
    geometry.filling_factor = param_data.FValue("Geometry", "filling_factor");
    check_input_param("filling_factor", geometry.filling_factor, 0.0, 1.0);

    // ratio of low/high density clumps
    geometry.density_ratio = param_data.FValue("Geometry", "density_ratio");
    check_input_param("density_ratio", geometry.density_ratio, 0.0, 1.0);

    // spherical or cubical clumps
    string clump_type = param_data.SValue("Geometry", "clump_type");
    int spherical_clumps = 0;
    if (clump_type == "sphere")
    {
        spherical_clumps = 1;
        // adjust filling factor to account of a sphere inscribed in a cube
        geometry.filling_factor /= (4. / 3.) * M_PI * pow(0.5, 3.0);
        //     cout << "new filling factor = " << geometry.filling_factor << endl;
    }

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

    // fill position arrays with
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

    // choose between linear and log spacing (log spacing used for TRUST slab benchmark)
    string slab_zspacing = param_data.SValue("Geometry", "slab_zspacing");
    if (slab_zspacing == "log")
    {
        float log_slab_z1 = log10(fabs(slab_z1));
        float log_slab_z2 = log10(fabs(slab_z2));
        float log_deltaz = (log_slab_z2 - log_slab_z1) / double(main_grid.index_dim[2] - 1);
        for (i = 0; i <= (main_grid.index_dim[2] - 1); i++)
            z_pos[i] = -1.0 * pow(10.0, log_slab_z1 + double(i) * log_deltaz);
        z_pos[main_grid.index_dim[2]] = size_z / 2.;
    }
    else
    {
        float deltaz = slab_z2 - slab_z1;
        for (i = 0; i <= main_grid.index_dim[2]-1; i++)
            z_pos[i] = slab_z1 +
            double(i)*(deltaz)/double(main_grid.index_dim[2]-1);
         z_pos[main_grid.index_dim[2]] = size_z/2.;
    }

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

    // fill main grid with dust density
    geometry.clump_densities[0] =
        geometry.tau /
        ((slab_z2 - slab_z1) * (geometry.filling_factor + geometry.density_ratio * (1.0 - geometry.filling_factor)));
    geometry.clump_densities[1] = geometry.density_ratio * geometry.clump_densities[0];

    // determine the nonslap density (compute here to keep consistant between a
    // homogenous slab and a clumpy slab)
    double nonslab_density = nonslab_density_ratio * geometry.tau / (slab_z2 - slab_z1);

    int j, k;
    // float radius = 0.0;
    double x_val = 0.0;
    double y_val = 0.0;
    double z_val = 0.0;
    for (k = 0; k < main_grid.index_dim[2]; k++)
    {
        z_val = (main_grid.positions[2][k] + main_grid.positions[2][k + 1]) / 2.0;
        for (j = 0; j < main_grid.index_dim[1]; j++)
        {
            y_val = (main_grid.positions[1][j] + main_grid.positions[1][j + 1]) / 2.0;
            for (i = 0; i < main_grid.index_dim[0]; i++)
            {
                x_val = (main_grid.positions[0][i] + main_grid.positions[0][i + 1]) / 2.0;
                // slab position should be a negative z if behind stars (observer at
                // positive z)
                if ((z_val >= slab_z1) && (z_val <= slab_z2))
                {
                    if (random_obj.random_num() <= geometry.filling_factor)
                        main_grid.grid(i, j, k).dust_tau_per_pc = geometry.clump_densities[0];
                    else
                        main_grid.grid(i, j, k).dust_tau_per_pc = geometry.clump_densities[1];
                }
                else
                    main_grid.grid(i, j, k).dust_tau_per_pc = nonslab_density;
                //  	cout << main_grid.grid(i,j,k).dust_tau_per_pc << " ";

                //  	main_grid.grid(i,j,k).absorbed_energy = 0.0;
            }
            //        cout << endl;
        }
    }

    // identify this grid as main grid
    main_grid.parent_grid_num = -1;

    // add main grid to grids vector
    geometry.grids.push_back(main_grid);

    // subdivide all overdense cells
    setup_dust_grid_subdivide_overdense_cells(geometry, spherical_clumps);
}
