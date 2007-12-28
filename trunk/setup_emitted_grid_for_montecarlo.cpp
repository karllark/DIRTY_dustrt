// ======================================================================
//   Procedure to setup the emitted energy grid for Monte Carlo 
// radiative transfer.  This means getting the total emitted energy
// per wavelength and normalizing the full grid to run from 0 to 1 
// through the full grid.  The goal is to be able to pick a grid cell
// where a photon is emitted randomly allowing for an arbitrary number
// of photons to be emitted.  This is mainly for the dust emission,
// but also may be 
//
// 2007 Nov/KDG - written
// ======================================================================
#include "setup_emitted_grid_for_montecarlo.h"

void setup_emitted_grid_for_montecarlo (geometry_struct& geometry,
					runinfo_struct& runinfo,
					GrainModel& CurGrainModel)

{
  int i,j,k,m,z = 0;
  uint x = 0;

  int n_emit_components = 1;
  if (runinfo.do_emission_grain) n_emit_components += 2*CurGrainModel.getNComp();

  // set for grid emission
  geometry.new_photon_source_type = NEW_PHOTON_GRID;

  vector< vector<double> > sum_emitted_lum;
  sum_emitted_lum.resize(n_emit_components);
  for (z = 0; z < n_emit_components; z++)
    sum_emitted_lum[z].resize(runinfo.wavelength.size(),0.0);

  // loop over all the defined grids
  for (m = 0; m < int(geometry.grids.size()); m++) {
    // loop of the cells in this grid
    for (k = 0; k < geometry.grids[m].index_dim[2]; k++)
      for (j = 0; j < geometry.grids[m].index_dim[1]; j++)
	for (i = 0; i < geometry.grids[m].index_dim[0]; i++) {
	  
	    for (x = 0; x < runinfo.wavelength.size(); x++) {

	      // normalize if there we need the emission by grain/emission type
	      if (runinfo.do_emission_grain) 
		for (z = 1; z < n_emit_components; z++)
		  geometry.grids[m].grid(i,j,k).emitted_energy[z][x] /= geometry.grids[m].grid(i,j,k).emitted_energy[0][x];

	      sum_emitted_lum[0][x] += geometry.grids[m].grid(i,j,k).emitted_energy[0][x];
	      geometry.grids[m].grid(i,j,k).emitted_energy[0][x] = sum_emitted_lum[0][x]/runinfo.emitted_lum[0][x];

	    }
    
	}
  }

}
