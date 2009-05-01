// ======================================================================
// setup a dust grid contained in a file (arbitrary grid)
//
// KDG 11 Jan 2009 - written
// ======================================================================
#include "setup_dust_grid_file.h"
// #define DEBUG_SDG

void setup_dust_grid_file (ConfigFile& param_data,
			   geometry_struct& geometry)

{
  // get the filename of the position file
  string pos_filename = param_data.SValue("Geometry","type_file_pos");
  // check that the file exists
  ifstream pos_file(pos_filename.c_str());
  if (pos_file.fail()) {
    cout << "Geometry position file (" << pos_filename << ") does not exist." << endl;
    exit(8);
  }
  pos_file.close();

  // get the filename of the tau/pc file
  string tau_pc_filename = param_data.SValue("Geometry","type_file_tau_pc");
  // check that the file exists
  ifstream tau_pc_file(tau_pc_filename.c_str());
  if (tau_pc_file.fail()) {
    cout << "Geometry tau/pc file (" << tau_pc_filename << ") does not exist." << endl;
    exit(8);
  }
  tau_pc_file.close();
  
  // open the position and tau/pc files for reading
  int status = 0;
  int hdutype = 0;   // type of header
  fitsfile *pos_file_ptr;   // pointer to pos FITS file
  fitsfile *tau_pc_file_ptr;  // pointer to tau/pc FITS file
  
  fits_open_file(&pos_file_ptr, pos_filename.c_str(), READONLY, &status);   // open the file
  check_fits_io(status, "fits_open_file : pos_file");

  fits_open_file(&tau_pc_file_ptr, tau_pc_filename.c_str(), READONLY, &status);   // open the file
  check_fits_io(status, "fits_open_file : tau_pc_file");

  // determine the number of grids in the file (1 = main_grid, rest are subgrids)
  int num_grids = 0;
  fits_get_num_hdus(pos_file_ptr, &num_grids, &status);

  // make sure the tau_pc file has the same number of grids
  fits_get_num_hdus(pos_file_ptr, &num_grids, &status);

  // move to the 1st extension (unless only primary exists) 
  // primary is for common information (all grid info in extensions)
  if (num_grids > 1) {
    fits_movabs_hdu(pos_file_ptr, 2, &hdutype, &status);
    fits_movabs_hdu(tau_pc_file_ptr, 2, &hdutype, &status);
    check_fits_io(status, "fits_movabs_hdu: pos & tau_pc files");
    num_grids--;
  }

  // useful variables
  long pos_naxes[2];
  char comment[72];

  // need
  // geometry.radius
  // geometry.angular_radius
  // geometry.tau
  // geometry.max_tau_per_cell
  // geometry.max_grid_depth

  // loop over the grids reading them into the necessary internal variables
  int i = 0;
  cout << "num_grids = " << num_grids << endl;
  for (i = 0; i < num_grids; i++) {
    // declare grid
    one_grid subgrid;

    // get the size of the grid
    fits_read_key(pos_file_ptr, TLONG, "NAXIS1", &pos_naxes[0], comment, &status);
    fits_read_key(pos_file_ptr, TLONG, "NAXIS2", &pos_naxes[1], comment, &status);
    check_fits_io(status, "fits_read_key : pos_file naxes1&2");

    // get the detailed size of each dimension
    fits_read_key(pos_file_ptr, TLONG, "XSIZE", &subgrid.index_dim[0], comment, &status);
    fits_read_key(pos_file_ptr, TLONG, "YSIZE", &subgrid.index_dim[1], comment, &status);
    fits_read_key(pos_file_ptr, TLONG, "ZSIZE", &subgrid.index_dim[2], comment, &status);
    check_fits_io(status, "fits_read_key : pos_file x,y,zsize");

    // create a 3d matrix to read the position info into
    NumUtils::Matrix<float> tmp_pos;
//     tmp_pos.MSize(3,max_index+1);

//     fits_read_img(pos_file_ptr, TFLOAT, 1, 

    
  }

  // close both files
  fits_close_file(pos_file_ptr, &status);
  check_fits_io(status, "fits_close_file : pos_file");
  fits_close_file(tau_pc_file_ptr, &status);
  check_fits_io(status, "fits_close_file : tau_pc_file");

  exit(0);
  
  int spherical_clumps = 0;
  // subdivide all overdense cells
  setup_dust_grid_subdivide_overdense_cells(geometry, spherical_clumps);

}
