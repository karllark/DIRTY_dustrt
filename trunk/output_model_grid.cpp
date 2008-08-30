// ======================================================================
//   Procedure to output the model grid of dirty.
//
// 2008 Aug/KDG - written
// ======================================================================
#include "output_model_grid.h"

void output_model_grid (geometry_struct& geometry,
			output_struct& output)

{
  // filename of the current output file
  string filename = "!" + output.file_base;
  filename += "_rad_field.fits";

  string filename_tau = "!" + output.file_base;
  filename_tau += "_tau_ref_per_pc.fits";

  string filename_pos = "!" + output.file_base;
  filename_pos += "_pos.fits";
      
#ifdef DEBUG_OMG
  cout << "filename for output (rad_field) = " << filename << endl;
  cout << "filename for output (tau_ref_per_pc) = " << filename_tau << endl;
#endif
      
  // create a FITS file with extensions to fill with the output of the model
  fitsfile *out_ptr;   // pointer to the output fits file
  int status = 0;
  fits_create_file(&out_ptr,filename.c_str(), &status);
  check_fits_io(status, "fits_create_file : output_model_grid");

  // create a FITS file with extensions to fill with the output of the model
  fitsfile *out_tau_ptr;   // pointer to the output fits file
  fits_create_file(&out_tau_ptr,filename_tau.c_str(), &status);
  check_fits_io(status, "fits_create_file : output_model_grid (tau)");

  // create a FITS file with extensions to fill with the output of the model
  fitsfile *out_pos_ptr;   // pointer to the output fits file
  fits_create_file(&out_pos_ptr,filename_pos.c_str(), &status);
  check_fits_io(status, "fits_create_file : output_model_grid (pos)");
      
  // the primary header is a blank image (very small)
  fits_create_img(out_ptr, 8, 0, 0, &status);
  check_fits_io(status,"fits_create_img : output_model_grid");

  // the primary header is a blank image (very small)
  fits_create_img(out_tau_ptr, 8, 0, 0, &status);
  check_fits_io(status,"fits_create_img : output_model_grid (tau)");

  // the primary header is a blank image (very small)
  fits_create_img(out_pos_ptr, 8, 0, 0, &status);
  check_fits_io(status,"fits_create_img : output_model_grid (pos)");
      
  // populate the primary header with the details of the run
      
  // final stuff for primary header
  fits_write_comment(out_ptr, "**---------------------------------**",&status);
  fits_write_comment(out_ptr, "Output of the DIRTY model",&status);
  fits_write_comment(out_ptr, "Karl D. Gordon & Karl A. Misselt", &status);
  fits_write_comment(out_ptr, "version v2.0prealpha (Aug 2008)", &status);
  fits_write_comment(out_ptr, "**---------------------------------**",&status);
  check_fits_io(status,"fits_write_comment : output_model_grid");

  // final stuff for primary header
  fits_write_comment(out_tau_ptr, "**---------------------------------**",&status);
  fits_write_comment(out_tau_ptr, "Output of the DIRTY model",&status);
  fits_write_comment(out_tau_ptr, "Karl D. Gordon & Karl A. Misselt", &status);
  fits_write_comment(out_tau_ptr, "version v2.0prealpha (Aug 2008)", &status);
  fits_write_comment(out_tau_ptr, "**---------------------------------**",&status);
  check_fits_io(status,"fits_write_comment : output_model_grid (tau)");

  // final stuff for primary header
  fits_write_comment(out_pos_ptr, "**---------------------------------**",&status);
  fits_write_comment(out_pos_ptr, "Output of the DIRTY model",&status);
  fits_write_comment(out_pos_ptr, "Karl D. Gordon & Karl A. Misselt", &status);
  fits_write_comment(out_pos_ptr, "version v2.0prealpha (Aug 2008)", &status);
  fits_write_comment(out_pos_ptr, "**---------------------------------**",&status);
  check_fits_io(status,"fits_write_comment : output_model_grid (pos)");
      
  int i,j,k = 0;
  uint m = 0;
  for (m = 0; m < geometry.grids.size(); m++) {
    // create a 4d matrix to copy the grid info into for output
//     NumUtils::FourCube<float> tmp_rad_field;
//     tmp_rad_field.FCSize(geometry.grids[m].index_dim[0],geometry.grids[m].index_dim[1],geometry.grids[m].index_dim[2]);

    // create a 3d matrix to copy the grid info into for output
    NumUtils::Cube<float> tmp_tau;
    tmp_tau.CSize(geometry.grids[m].index_dim[0],geometry.grids[m].index_dim[1],geometry.grids[m].index_dim[2]);

    // create a 3d matrix to copy the grid info into for output
    // assumes this is a cube - will fail for the dexp geometry (can be fixed, but need to code correctly)
    NumUtils::Matrix<float> tmp_pos;
    tmp_pos.MSize(3,geometry.grids[m].index_dim[0]+1);

    // get pos data
    for (i = 0; i <= geometry.grids[m].index_dim[0]; i++)
      for (k = 0; k < 3; k++)
	tmp_pos(k,i) = geometry.grids[m].positions[k][i];
    
    // loop of the cells in this grid
    for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
      for (j = 0; j < geometry.grids[m].index_dim[1]; j++)
	for (i = 0; i < geometry.grids[m].index_dim[0]; i++) {
	  tmp_tau(i,j,k) = geometry.grids[m].grid(i,j,k).dust_tau_per_pc;
// 	  tmp_rad_field(i,j,k,n) = geometry.grids[m].grid(i,j,k).absorbed_energy[n];
	}

    // create and output each grid
    fits_create_img(out_tau_ptr, -32, 3, geometry.grids[m].index_dim, &status);
    check_fits_io(status,"fits_create_image : output_model_grid (tau)");
      
    fits_write_img(out_tau_ptr, TFLOAT, 1, geometry.grids[m].index_dim[0]*geometry.grids[m].index_dim[1]*geometry.grids[m].index_dim[2], 
		   &tmp_tau[0], &status);

    // create and output each grid
    long tmp_pos_index[2];
    tmp_pos_index[0] = 3;
    tmp_pos_index[1] = geometry.grids[m].index_dim[0] + 1;
    fits_create_img(out_pos_ptr, -32, 2, tmp_pos_index, &status);
    check_fits_io(status,"fits_create_image : output_model_grid (pos)");
      
    fits_write_img(out_pos_ptr, TFLOAT, 1, tmp_pos_index[0]*tmp_pos_index[1], 
		   &tmp_pos[0], &status);

  }

  // close FITS File
  fits_close_file(out_ptr, &status);
  check_fits_io(status,"fits_close_file : output_model_grid");

  // close FITS File
  fits_close_file(out_tau_ptr, &status);
  check_fits_io(status,"fits_close_file : output_model_grid (tau)");

  // close FITS File
  fits_close_file(out_pos_ptr, &status);
  check_fits_io(status,"fits_close_file : output_model_grid (pos)");

}
