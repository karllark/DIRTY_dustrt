// ======================================================================
//   Procedure to do Monte Carlo radiative transfer for photons through
// dust.  This procedure assumes the dust grid is already constructed.
//
// 2004 Dec/KDG - written
// 2005 May/KDG - added output stuff
// ======================================================================
#include "radiative_transfer.h"
//#define DEBUG_RT
//#define OUTNUM 4628428

void radiative_transfer (geometry_struct& geometry,
			 runinfo_struct& runinfo,
			 output_struct& output,
			 photon_data& photon,
			 random_dirty random_obj)

{
  int i;

  // initialize (and allocate if needed) output structure
  initialize_output(output, geometry);

  // TBD: add code to dynamically determine if enough photons have been
  //      run depending on the output quantitity desired
  for (i = 0; i < geometry.n_photons; i++) {
    // print a status statement if asked
    if (runinfo.verbose >= 1) {
      if ((geometry.n_photons > 100) && ((i % (geometry.n_photons/20)) == 0)) cout << "current # = " << i << endl;
    }
    // start the photon at the star position
    new_photon(photon, geometry, runinfo, random_obj);
    photon.number = i;
#ifdef DEBUG_RT
    if (photon.number > OUTNUM) cout << photon.number << " ";
    if (photon.number > OUTNUM) cout << "np1* done; "; cout.flush();
#endif

#ifdef SAVE_TRAJ
    cout << photon.number << " ";
    cout << photon.num_scat << " ";
    int k;
    for (k = 0; k < 3; k++) 
      cout << photon.position[k] << " ";
    cout << " save_traj begin [pnum, nscat, pos(1-3)]";
    cout << endl;
#endif

    // find the first scattering site (forced)
    forced_first_scatter(geometry, photon, random_obj);
#ifdef DEBUG_RT
    if (photon.number > OUTNUM) cout << "ffs done; "; cout.flush();
#endif

#ifdef SAVE_TRAJ
    cout << photon.number << " ";
    cout << photon.num_scat << " ";
    for (k = 0; k < 3; k++) 
      cout << photon.position[k] << " ";
    cout << " save_traj [pnum, nscat, pos(1-3)]";
    cout << endl;
#endif

    // classify stellar photon(s)
    classify_stellar_photon(output, photon, geometry, runinfo);
#ifdef DEBUG_RT
    if (photon.number > OUTNUM) cout << "cstp done; "; cout.flush();
#endif

    // loop until photon escapes
    int escape = 0;
    while (!escape) {

      // classify the scattered photon
      classify_scattered_photon(output, photon, geometry, runinfo);
#ifdef DEBUG_RT
      if (photon.number > OUTNUM) cout << "cscp done; "; cout.flush();
#endif

      // scatter the photon into a new direction
      scatter_photon(geometry, photon, random_obj);
#ifdef DEBUG_RT
      if (photon.number > OUTNUM) cout << "sp done; "; cout.flush();
#endif

      // determine the next scattering site or
      // the position where the photon escapes the dust
      escape = next_scatter(geometry, photon, random_obj);
#ifdef DEBUG_RT
      if (photon.number > OUTNUM) cout << "ns done; "; cout.flush();
#endif

#ifdef SAVE_TRAJ
      cout << photon.number << " ";
      cout << photon.num_scat << " ";
      double rad = 0.0;
      for (k = 0; k < 3; k++) {
	cout << photon.position[k] << " ";
	rad += pow(photon.position[k],2);
      }
      cout << escape << " ";
      cout << pow(rad,0.5) << " ";
      cout << " save_traj [pnum, nscat, pos(1-3), escape, radius]";
      cout << endl;
#endif

#ifdef DEBUG_RT
      if (photon.number > OUTNUM) cout << endl;
#endif
    }


  }

}
