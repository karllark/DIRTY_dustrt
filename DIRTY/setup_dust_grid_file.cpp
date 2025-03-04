// ======================================================================
// setup a dust grid contained in a file (arbitrary grid)
//
// KDG 11 Jan 2009 - started
// KDG 6-7 Oct 2009 - major work
// ======================================================================
#include "setup_dust_grid_file.h"
// #define DEBUG_SDGF

void
setup_dust_grid_file (ConfigFile &param_data, geometry_struct &geometry)

{
  // maximum optical depth per cell (controls when a cell is subdivided)
  geometry.max_tau_per_cell = param_data.FValue ("Geometry", "max_tau_per_cell");
  check_input_param ("max_tau_per_cell", geometry.max_tau_per_cell, 0.0, 1e10);
  geometry.max_tau_per_cell_x = geometry.max_tau_per_cell;
  geometry.max_tau_per_cell_y = geometry.max_tau_per_cell;
  geometry.max_tau_per_cell_z = geometry.max_tau_per_cell;

  // get the filename of the position file
  string pos_filename = param_data.SValue ("Geometry", "type_file_pos");
  // check that the file exists
  ifstream pos_file (pos_filename.c_str ());
  if (pos_file.fail ())
    {
      cout << "Geometry position file (" << pos_filename << ") does not exist." << endl;
      exit (8);
    }
  pos_file.close ();

  // get the filename of the tau/pc file
  string tau_pc_filename = param_data.SValue ("Geometry", "type_file_tau_pc");
  // check that the file exists
  ifstream tau_pc_file (tau_pc_filename.c_str ());
  if (tau_pc_file.fail ())
    {
      cout << "Geometry tau/pc file (" << tau_pc_filename << ") does not exist." << endl;
      exit (8);
    }
  tau_pc_file.close ();

  // open the position and tau/pc files for reading
  int status = 0;
  int hdutype = 0;           // type of header
  fitsfile *pos_file_ptr;    // pointer to pos FITS file
  fitsfile *tau_pc_file_ptr; // pointer to tau/pc FITS file

  fits_open_file (&pos_file_ptr, pos_filename.c_str (), READONLY,
                  &status); // open the file
  check_fits_io (status, "fits_open_file : pos_file");

  fits_open_file (&tau_pc_file_ptr, tau_pc_filename.c_str (), READONLY,
                  &status); // open the file
  check_fits_io (status, "fits_open_file : tau_pc_file");

  // determine the number of grids in the files (1 = main_grid, rest are
  // subgrids)
  int num_grids = 0;
  fits_get_num_hdus (pos_file_ptr, &num_grids, &status);

  int num_grids_tau = 0;
  fits_get_num_hdus (pos_file_ptr, &num_grids_tau, &status);

  // make sure the tau_pc file has the same number of grids
  if (num_grids != num_grids_tau)
    {
      cout << "Number of grids in position and tau/pc files does not match" << endl;
      cout << "# grids positon, tau = " << num_grids << ", " << num_grids_tau << endl;
      exit (8);
    }

  // useful variables
  long pos_naxes[2];
  char comment[72];

  double nulval = 0.0;
  int anynul = 0;
  long fpixel[2];
  long cfpixel[3];

  // need
  // geometry.radius
  // geometry.angular_radius
  // geometry.tau
  // geometry.max_tau_per_cell
  // geometry.max_grid_depth

  // loop over the grids reading them into the necessary internal variables
  int i, k, l, m = 0;
#ifdef DEBUG_SDGF
  cout << "num_grids = " << num_grids << endl;
#endif
  for (i = 0; i < num_grids; i++)
    {
#ifdef DEBUG_SDGF
      cout << endl << "grid # = " << (i + 1) << endl;
#endif

      // declare grid
      one_grid subgrid;

      // get the size of the grid
      fits_read_key (pos_file_ptr, TLONG, "NAXIS1", &pos_naxes[0], comment, &status);
      fits_read_key (pos_file_ptr, TLONG, "NAXIS2", &pos_naxes[1], comment, &status);
      check_fits_io (status, "fits_read_key : pos_file naxes1&2");

      // get the detailed size of each dimension
      fits_read_key (tau_pc_file_ptr, TLONG, "NAXIS1", &subgrid.index_dim[0], comment, &status);
      fits_read_key (tau_pc_file_ptr, TLONG, "NAXIS2", &subgrid.index_dim[1], comment, &status);
      fits_read_key (tau_pc_file_ptr, TLONG, "NAXIS3", &subgrid.index_dim[2], comment, &status);
      check_fits_io (status, "fits_read_key : tau_pc_file x,y,zsize");

      // check that the cube size is less than the positions
      if (!((subgrid.index_dim[0] < pos_naxes[0]) && (subgrid.index_dim[1] < pos_naxes[0])
            && (subgrid.index_dim[2] < pos_naxes[0])))
        {
          cout << "At least one of the cube dimensions exceeds the position "
                  "dimensions"
               << endl;
          cout << "position dimension = " << pos_naxes[0] << endl;
          exit (8);
        }

      // create the position arrays
      vector<double> x_pos (subgrid.index_dim[0] + 1);
      vector<double> y_pos (subgrid.index_dim[1] + 1);
      vector<double> z_pos (subgrid.index_dim[2] + 1);

      // read in the positions
      nulval = 0.0;
      anynul = 0;
      fpixel[0] = 1;
      fpixel[1] = 1;
      fits_read_pix (pos_file_ptr, TDOUBLE, &fpixel[0], subgrid.index_dim[0] + 1, &nulval, &x_pos[0], &anynul, &status);
      check_fits_io (status, "setup_dust_grid_file, fits_read_pix: pos[0]");
#ifdef DEBUG_SDGF
      for (k = 0; k < (subgrid.index_dim[1] + 1); k++)
        cout << x_pos[k] << " ";
      cout << endl;
#endif

      fpixel[0] = 1;
      fpixel[1] = 2;
      fits_read_pix (pos_file_ptr, TDOUBLE, &fpixel[0], subgrid.index_dim[1] + 1, &nulval, &y_pos[0], &anynul, &status);
      check_fits_io (status, "setup_dust_grid_file, fits_read_pix: pos[1]");
#ifdef DEBUG_SDGF
      for (k = 0; k < (subgrid.index_dim[1] + 1); k++)
        cout << y_pos[k] << " ";
      cout << endl;
#endif

      fpixel[0] = 1;
      fpixel[1] = 3;
      fits_read_pix (pos_file_ptr, TDOUBLE, &fpixel[0], subgrid.index_dim[2] + 1, &nulval, &z_pos[0], &anynul, &status);
      check_fits_io (status, "setup_dust_grid_file, fits_read_pix: pos[3]");
#ifdef DEBUG_SDGF
      for (k = 0; k < (subgrid.index_dim[2] + 1); k++)
        cout << z_pos[k] << " ";
      cout << endl;
#endif

      // add position arrays to main grid
      subgrid.positions.push_back (x_pos);
      subgrid.positions.push_back (y_pos);
      subgrid.positions.push_back (z_pos);

      // determine some needed quanties
#ifdef DEBUG_SDGF
      cout << "dim, size, cube size" << endl;
#endif
      for (k = 0; k < 3; k++)
        {
          subgrid.phys_grid_size[k] = subgrid.positions[k][subgrid.index_dim[k]] - subgrid.positions[k][0];
          // subgrid.phys_cube_size[k] =
          // subgrid.phys_grid_size[k]/subgrid.index_dim[k];
#ifdef DEBUG_SDGF
          cout << k << " " << subgrid.phys_grid_size[k] << endl;
#endif
        }

      // determine useful quantities for geometry record
      if (i == 0)
        {
          geometry.radius = subgrid.phys_grid_size[0] / 2.;
          for (k = 1; k < 3; k++)
            if (subgrid.phys_grid_size[k] / 2. > geometry.radius)
              geometry.radius = subgrid.phys_grid_size[k] / 2.;

          geometry.angular_radius = atan (1.5 * geometry.radius / (geometry.distance - geometry.radius));

#ifdef DEBUG_SDGF
          cout << "radius [pc] = " << geometry.radius << endl;
          cout << "angular radius [rad] = " << geometry.angular_radius << endl;
#endif

          geometry.tau = 0.0;
          fits_read_key (tau_pc_file_ptr, TFLOAT, "RAD_TAU", &geometry.tau, comment, &status);
          if (status != 0)
            {
              geometry.tau = 0.0;
              status = 0;
            }
          check_input_param ("arbitrary file: radial optical depth", geometry.tau, 0.0, 1000.);

          // variables that are not used, but are written to an output FITS
          // file need to be set to avoid bad float to string conversions in
          // creating the FITS file
          geometry.density_ratio = 1.0;
          geometry.clump_densities[0] = 0.0;
          geometry.clump_densities[1] = 0.0;

          // max grid depth
          fits_read_key (tau_pc_file_ptr, TLONG, "GRDDEPTH", &geometry.max_grid_depth, comment, &status);
          if (status != 0)
            {
              if (num_grids == 1)
                geometry.max_grid_depth = 1;
              else
                geometry.max_grid_depth = 2;
              status = 0;
            }
          check_input_param ("arbitrary file: max_grid_depth", geometry.tau, 0.0, 1e6);
        }

      // now get the cube of tau/pc
      // create a 3d matrix to read the tau/pc cube into
      NumUtils::Cube<float> tmp_tau;
      tmp_tau.CSize (subgrid.index_dim[0], subgrid.index_dim[1], subgrid.index_dim[2]);

      cfpixel[0] = 1;
      cfpixel[1] = 1;
      cfpixel[2] = 1;
      fits_read_pix (tau_pc_file_ptr, TFLOAT, &cfpixel[0],
                     subgrid.index_dim[0] * subgrid.index_dim[1] * subgrid.index_dim[2], &nulval, &tmp_tau[0], &anynul,
                     &status);
      check_fits_io (status, "setup_dust_grid_file, fits_read_pix: tau");
#ifdef DEBUG_SDGF
      cout << "tau/pc along z-axis" << endl;
      for (k = 0; k < subgrid.index_dim[2]; k++)
        cout << tmp_tau (int (subgrid.index_dim[2] / 2), int (subgrid.index_dim[2] / 2), k) << " ";
      cout << endl;
#endif
      // allocate the grid cells
      subgrid.grid.CSize (subgrid.index_dim[0], subgrid.index_dim[1], subgrid.index_dim[2]);

      // now loop over the grid and put the tau/pc into the grid cell
      for (k = 0; k < subgrid.index_dim[2]; k++)
        for (l = 0; l < subgrid.index_dim[1]; l++)
          for (m = 0; m < subgrid.index_dim[0]; m++)
            subgrid.grid (m, l, k).dust_tau_per_pc = tmp_tau (m, l, k);

      // identify the parent grid
      if (i == 0)
        subgrid.parent_grid_num = -1;
      else
        {
          fits_read_key (tau_pc_file_ptr, TINT, "PAR_GRID", &subgrid.parent_grid_num, comment, &status);
          check_fits_io (status, "fits_read_key : tau_pc_file parent_grid number (PAR_GRID)");
          check_input_param ("arbitrary file input: parent_grid number", subgrid.parent_grid_num, 0, (num_grids - 1));
        }

      geometry.grids.push_back (subgrid);

      // move to the next subgrid, unless we are at the end
      if (i < (num_grids - 1))
        {
          fits_movrel_hdu (pos_file_ptr, 1, &hdutype, &status);
          fits_movrel_hdu (tau_pc_file_ptr, 1, &hdutype, &status);
          check_fits_io (status, "fits_movrel_hdu: pos & tau_pc files");
        }
    }

  // close both files
  fits_close_file (pos_file_ptr, &status);
  check_fits_io (status, "fits_close_file : pos_file");
  fits_close_file (tau_pc_file_ptr, &status);
  check_fits_io (status, "fits_close_file : tau_pc_file");

  // check to make sure that all subgrids designated exist
  if (num_grids > 1)
    {
      vector<int> par_idim;
      par_idim.reserve (3);
      for (i = 0; i < 3; i++)
        par_idim[i] = 0;
      setup_dust_grid_check_grid (geometry, 0, -1, par_idim);
    }

  int spherical_clumps = 0;
  // subdivide all overdense cells
  setup_dust_grid_subdivide_overdense_cells (geometry, spherical_clumps);
}
