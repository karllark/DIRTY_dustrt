// ======================================================================
//   Procedure which calls the right procedure to generate a new photon.
//
// 2006 Aug/KDG - written
// 21 Nov 2006 KDG - added code to choose between discrete stars and
//                   a diffuse radiation field
//  8 Mar 2007 KDG - fixed and speed up the code for deciding how to emit a photon
// ======================================================================
#include "new_photon.h"

void new_photon (photon_data& photon,
		 geometry_struct& geometry,
		 runinfo_struct& runinfo,
		 random_dirty random_obj)

{
  switch(geometry.new_photon_source_type)
    {
    case NEW_PHOTON_DISCRETE_STARS: 
      new_photon_discrete_stars(photon, geometry, random_obj);
      break;
    case NEW_PHOTON_DIFFUSE_ISOTROPIC:
    case NEW_PHOTON_DIFFUSE_FILE:
      new_photon_diffuse_source(photon, geometry, random_obj);
      break;
    case NEW_PHOTON_GRID:
      new_photon_grid_source(photon, geometry, runinfo, random_obj);
      break;
    case NEW_PHOTON_DEXP_DISK:
      new_photon_dexp_disk(photon, geometry, random_obj);
      break;
    default: 
      cout << "new_photon for input source type (" << geometry.source_type << ") not found [NEW CODE NEEDED]." << endl;
      exit(8);
    }

}
