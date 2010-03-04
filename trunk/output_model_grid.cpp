// ======================================================================
//   Procedure to output the model grid of dirty.
//
// 2008 Aug/KDG - written
// ======================================================================
#include "output_model_grid.h"

void output_model_grid (geometry_struct& geometry,
			output_struct& output,
			runinfo_struct& runinfo)

{
  // filename of the current output file
  string filename = "!" + output.file_base;
  filename += "_rad_field.fits";

  string filename_unc = "!" + output.file_base;
  filename_unc += "_rad_field_unc.fits";

  // filename of the wavelength grid output file
  string filename_wave = "!" + output.file_base;
  filename_wave += "_wave_grid.fits";

  string filename_tau = "!" + output.file_base;
  filename_tau += "_tau_ref_per_pc.fits";

  string filename_pos = "!" + output.file_base;
  filename_pos += "_pos.fits";
      
#ifdef DEBUG_OMG
  cout << "filename for output (rad_field) = " << filename << endl;
  cout << "filename for output (tau_ref_per_pc) = " << filename_tau << endl;
#endif
      
  // create a FITS file with extensions to fill with the output of the model
  int status = 0;
  fitsfile *out_ptr;   // pointer to the output fits file
  fits_create_file(&out_ptr,filename.c_str(), &status);
  check_fits_io(status, "fits_create_file : output_model_grid (rad_field)");

  fitsfile *out_unc_ptr;   // pointer to the output fits file
  fits_create_file(&out_ptr,filename_unc.c_str(), &status);
  check_fits_io(status, "fits_create_file : output_model_grid (rad_field unc)");

  // create a FITS file with extensions to fill with the output of the model
  fitsfile *out_tau_ptr;   // pointer to the output fits file
  fits_create_file(&out_tau_ptr,filename_tau.c_str(), &status);
  check_fits_io(status, "fits_create_file : output_model_grid (tau)");

  // create a FITS file with extensions to fill with the wavelength grid
  fitsfile *out_wave_ptr;   // pointer to the output fits file
  fits_create_file(&out_wave_ptr,filename_wave.c_str(), &status);
  check_fits_io(status, "fits_create_file : output_model_grid (wave)");

  // create a FITS file with extensions to fill with the output of the model
  fitsfile *out_pos_ptr;   // pointer to the output fits file
  fits_create_file(&out_pos_ptr,filename_pos.c_str(), &status);
  check_fits_io(status, "fits_create_file : output_model_grid (pos)");
      
  int i,j,k,n = 0;
  uint m = 0;
  int n_waves = geometry.grids[0].grid(0,0,0).absorbed_energy.size();
  for (m = 0; m < geometry.grids.size(); m++) {
    // create a 4d matrix to copy the grid info into for output
    NumUtils::FourVector<float> tmp_rad_field;
    tmp_rad_field.FVSize(geometry.grids[m].index_dim[0],geometry.grids[m].index_dim[1],geometry.grids[m].index_dim[2],n_waves);

    // create a 4d matrix to copy the grid info into for output
    NumUtils::FourVector<float> tmp_rad_field_unc;
    tmp_rad_field_unc.FVSize(geometry.grids[m].index_dim[0],geometry.grids[m].index_dim[1],geometry.grids[m].index_dim[2],n_waves);

    // create a 3d matrix to copy the grid info into for output
    NumUtils::Cube<float> tmp_tau;
    tmp_tau.CSize(geometry.grids[m].index_dim[0],geometry.grids[m].index_dim[1],geometry.grids[m].index_dim[2]);

    // get the max index size
    int max_index = geometry.grids[m].index_dim[0];
    for (i = 1; i < 3; i++) 
      if (geometry.grids[m].index_dim[i] > max_index) max_index = geometry.grids[m].index_dim[i];

    // create a 3d matrix to copy the grid info into for output
    NumUtils::Matrix<double> tmp_pos;
    tmp_pos.MSize(max_index+1,3);

    // get pos data
    // fill each dimension separately to allow for non-cubes
    for (i = 0; i <= geometry.grids[m].index_dim[0]; i++)
      tmp_pos(i,0) = geometry.grids[m].positions[0][i];
    for (i = 0; i <= geometry.grids[m].index_dim[1]; i++)
      tmp_pos(i,1) = geometry.grids[m].positions[1][i];
    for (i = 0; i <= geometry.grids[m].index_dim[2]; i++)
      tmp_pos(i,2) = geometry.grids[m].positions[2][i];
    
    // loop of the cells in this grid
    float rad_unc = 0.0;
    for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
      for (j = 0; j < geometry.grids[m].index_dim[1]; j++)
	for (i = 0; i < geometry.grids[m].index_dim[0]; i++) {
	  tmp_tau(i,j,k) = geometry.grids[m].grid(i,j,k).dust_tau_per_pc;
	  for (n = 0; n < n_waves; n++) {
	    tmp_rad_field(i,j,k,n) = geometry.grids[m].grid(i,j,k).absorbed_energy[n];
	    if (geometry.grids[m].grid(i,j,k).absorbed_energy_num_photons[n] >= 5) {
	      rad_unc = geometry.grids[m].grid(i,j,k).absorbed_energy_x2[n]/geometry.grids[m].grid(i,j,k).absorbed_energy_num_photons[n] -
		pow(geometry.grids[m].grid(i,j,k).absorbed_energy[n]/geometry.grids[m].grid(i,j,k).absorbed_energy_num_photons[n],2.0);
	      if (rad_unc > 0.0)
		rad_unc = sqrt(rad_unc);
	      else
		rad_unc = 0.0;
	      tmp_rad_field_unc(i,j,k,n) = rad_unc;
	    }
	  }
	}
  
    // create and output each grid (rad_field)
    long four_vector_size[4];
    for (i = 0; i < 3; i++) four_vector_size[i] = geometry.grids[m].index_dim[i];
    four_vector_size[3] = n_waves;
    fits_create_img(out_ptr, -32, 4, four_vector_size, &status);
    check_fits_io(status,"fits_create_image : output_model_grid (rad_field)");
      
    fits_write_img(out_ptr, TFLOAT, 1, geometry.grids[m].index_dim[0]*geometry.grids[m].index_dim[1]*geometry.grids[m].index_dim[2]*n_waves, 
		   &tmp_rad_field[0], &status);

    // create and output each grid (rad_field unc)
    fits_create_img(out_unc_ptr, -32, 4, four_vector_size, &status);
    check_fits_io(status,"fits_create_image : output_model_grid (rad_field unc)");
      
    fits_write_img(out_unc_ptr, TFLOAT, 1, geometry.grids[m].index_dim[0]*geometry.grids[m].index_dim[1]*geometry.grids[m].index_dim[2]*n_waves, 
		   &tmp_rad_field_unc[0], &status);

    // create and output each grid (tau)
    fits_create_img(out_tau_ptr, -32, 3, geometry.grids[m].index_dim, &status);
    check_fits_io(status,"fits_create_image : output_model_grid (tau)");
      
    fits_write_img(out_tau_ptr, TFLOAT, 1, geometry.grids[m].index_dim[0]*geometry.grids[m].index_dim[1]*geometry.grids[m].index_dim[2], 
		   &tmp_tau[0], &status);

    // create and output each grid (positions)
    long tmp_pos_index[2];
    tmp_pos_index[0] = max_index+1;
    tmp_pos_index[1] = 3;
    fits_create_img(out_pos_ptr, -64, 2, tmp_pos_index, &status);
    check_fits_io(status,"fits_create_image : output_model_grid (pos)");
      
    fits_write_img(out_pos_ptr, TDOUBLE, 1, tmp_pos_index[0]*tmp_pos_index[1], 
		   &tmp_pos[0], &status);

//     // output extra info needed to fully define geometry
//     fits_write_key(out_pos_ptr, TLONG, "XSIZE", &geometry.grids[m].index_dim[0], "index dimension of x vals", &status);
//     fits_write_key(out_pos_ptr, TLONG, "YSIZE", &geometry.grids[m].index_dim[1], "index dimension of y vals", &status);
//     fits_write_key(out_pos_ptr, TLONG, "ZSIZE", &geometry.grids[m].index_dim[2], "index dimension of z vals", &status);
//     check_fits_io(status,"fits_write_key : output_model_grid (pos)");

    if (m != 0) {
      fits_write_key(out_tau_ptr, TINT, "PAR_GRID", &geometry.grids[m].parent_grid_num, "grid number of partent", &status);
    } else {
      // populate the primary header with the details of the run
      
      // final stuff for primary header
      fits_write_comment(out_ptr, "**---------------------------------**",&status);
      fits_write_comment(out_ptr, "Output of the DIRTY model",&status);
      fits_write_comment(out_ptr, "Karl D. Gordon & Karl A. Misselt", &status);
      fits_write_comment(out_ptr, "version v2.0prealpha (Oct 2009)", &status);
      fits_write_comment(out_ptr, "**---------------------------------**",&status);
      check_fits_io(status,"fits_write_comment : output_model_grid");

      // final stuff for primary header
      fits_write_comment(out_unc_ptr, "**---------------------------------**",&status);
      fits_write_comment(out_unc_ptr, "Output of the DIRTY model",&status);
      fits_write_comment(out_unc_ptr, "Karl D. Gordon & Karl A. Misselt", &status);
      fits_write_comment(out_unc_ptr, "version v2.0prealpha (Oct 2009)", &status);
      fits_write_comment(out_unc_ptr, "**---------------------------------**",&status);
      check_fits_io(status,"fits_write_comment : output_model_grid (unc)");
      
      // extra information needed for setting up the geometry
      fits_write_key(out_tau_ptr, TFLOAT, "RAD_TAU", &geometry.tau, "radial optical depth", &status);
      fits_write_key(out_tau_ptr, TLONG, "GRDDEPTH", &geometry.max_grid_depth, "maximum depth of the grid", &status);

      // final stuff for primary header
      fits_write_comment(out_tau_ptr, "**---------------------------------**",&status);
      fits_write_comment(out_tau_ptr, "Output of the DIRTY model",&status);
      fits_write_comment(out_tau_ptr, "Karl D. Gordon & Karl A. Misselt", &status);
      fits_write_comment(out_tau_ptr, "version v2.0prealpha (Oct 2009)", &status);
      fits_write_comment(out_tau_ptr, "**---------------------------------**",&status);
      check_fits_io(status,"fits_write_comment : output_model_grid (tau)");
      
      // final stuff for primary header
      fits_write_comment(out_pos_ptr, "**---------------------------------**",&status);
      fits_write_comment(out_pos_ptr, "Output of the DIRTY model",&status);
      fits_write_comment(out_pos_ptr, "Karl D. Gordon & Karl A. Misselt", &status);
      fits_write_comment(out_pos_ptr, "version v2.0prealpha (Oct 2009)", &status);
      fits_write_comment(out_pos_ptr, "**---------------------------------**",&status);
      check_fits_io(status,"fits_write_comment : output_model_grid (pos)");
    }

  }

  // create and output the wavelength
  long wavelength_size[1];
  wavelength_size[0] = n_waves;
  fits_create_img(out_wave_ptr, -32, 1, wavelength_size, &status);
  check_fits_io(status,"fits_create_image : output_model_grid (wavelength grid)");
  
  fits_write_img(out_wave_ptr, TFLOAT, 1, n_waves, &runinfo.wavelength[0], &status);

  // close FITS File
  fits_close_file(out_ptr, &status);
  check_fits_io(status,"fits_close_file : output_model_grid");

  // close FITS File
  fits_close_file(out_unc_ptr, &status);
  check_fits_io(status,"fits_close_file : output_model_grid (unc)");

  // close FITS File
  fits_close_file(out_wave_ptr, &status);
  check_fits_io(status,"fits_close_file : output_model_grid (wave)");

  // close FITS File
  fits_close_file(out_tau_ptr, &status);
  check_fits_io(status,"fits_close_file : output_model_grid (tau)");

  // close FITS File
  fits_close_file(out_pos_ptr, &status);
  check_fits_io(status,"fits_close_file : output_model_grid (pos)");

}
