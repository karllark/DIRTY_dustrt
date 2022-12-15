// ======================================================================
//   Procedure to do run the radiative transfer for all the wavelengths.
//
// 2008 Dec/KDG - written
// ======================================================================
#include "radiative_transfer_many_waves.h"
//#define DEBUG_MRT
//#define OUTNUM -1
//#define OUTNUM 1559

void radiative_transfer_many_waves (geometry_struct& geometry,
				    runinfo_struct& runinfo,
				    output_struct& output,
				    photon_data& photon,
				    random_dirty& random_obj,
				    int rt_type,
				    int iter_num)

{
#ifdef DEBUG_MRT
      cout << "mrt: mrt begin; ";
      cout.flush();
#endif
  // setup/(re)initialize absorbed energy grid
  setup_absorbed_energy_grid(geometry, runinfo, iter_num);

  // loop over all wavelengths
  int i;
  for (i = 0; i < runinfo.n_waves; i++) {

    int do_rt = 1;
    if (rt_type > REG_RT) {
      switch (rt_type)
	{
	case ERE_RT:
	  if (runinfo.emitted_ere_lum[0][i] < 1e-3*runinfo.out_sed_lum[0][i]*runinfo.sed_lum[i]) do_rt = 0;
	  break;
	case DE_RT:
	  if (runinfo.emitted_lum[0][i] < 1e-20*runinfo.out_sed_lum[0][i]*runinfo.sed_lum[i]) do_rt = 0;
	  break;
	}
    }

    if (geometry.abs_energy_storage_type == 0) {
      geometry.abs_energy_wave_index = i;
    } else {
      geometry.abs_energy_wave_index = 0;
    }

    if (do_rt) {
      if (runinfo.verbose >= 1) {
	cout << "working on wavelength [micron] = " << runinfo.wavelength[i]*Constant::CM_UM << endl;
	cout << "tau = " << runinfo.tau_to_tau_ref[i]*geometry.tau << " ";
	cout << "albedo = " << runinfo.albedo[i] << " ";
	cout << "g = " << runinfo.g[i] << " ";
	cout << endl;
      }

      //     // check the absorbed energy grid (temp needed as energy not conserved)
      //     // KDG - 23 Mar 2008
//       check_absorbed_energy_grid(geometry, runinfo);

      // setup dust grains for this wavelength
      get_dust_scat_parameters(i, runinfo, geometry);

      // do RT part
#ifdef DEBUG_MRT
      cout << "rt: rt begin; ";
      cout.flush();
#endif
      radiative_transfer(geometry, runinfo, output, photon, random_obj);
#ifdef DEBUG_MRT
      cout << "rt: rt end; ";
      cout.flush();
#endif

      // output RT results
      output_results(output, geometry, runinfo, i);

    }
      //     // check the absorbed energy grid (temp needed as energy not conserved)
      //     // KDG - 23 Mar 2008
      //     cout << "**test2**" << endl;
    //    check_absorbed_energy_grid(geometry, runinfo);

      // store the result (either in memory or on disk)
      // remember to zero out the absorbed energy grid
//     if ((runinfo.do_dust_emission) || (runinfo.do_ere_emission)) {
#ifdef DEBUG_MRT
      cout << endl;
      cout << "wave [cm] = " << geometry.wavelength << endl;
      cout << "te: saeg start; ";
      cout.flush();
#endif
			if (runinfo.do_dust_emission == 1) {
        store_absorbed_energy_grid(geometry, runinfo, output, i, iter_num);
			}
#ifdef DEBUG_MRT
      cout << "te: saeg end; ";
      cout.flush();
#endif

      //     // check the absorbed energy grid (temp needed as energy not conserved)
      //     // KDG - 23 Mar 2008
      //check_absorbed_energy_grid(geometry, runinfo);
      //     if (i == 1) exit(8);
      //  	exit(8);
//     }

    // check the absorbed energy grid (temp needed as energy not conserved)
    // KDG - 23 Mar 2008
//     check_absorbed_energy_grid(geometry, runinfo);
  }
}
