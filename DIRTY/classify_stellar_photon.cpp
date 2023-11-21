// ======================================================================
//   Procedure to classify stellar photon(s) into the output images
// and global totals.
//
// 2005 May/KDG - written
// 2006 Apr/KDG - updated to include rotate_zaxis_for_observe
// 2007 Feb/KDG - added calcuations of stellar weight (correct for observer
// position)
// 2007 Dec/KDG - changed denominator of atan from d-z to d
// 2008 Mar/KDG - added output for emission/grain types 2008 May/KDG - changed
// denominator of atan from d to d-z (see fix in setup_dust_grid_slab.cpp)
// ======================================================================
#include "classify_stellar_photon.h"
// #define DEBUG_CSP
// #define OUTNUM 0

void classify_stellar_photon(output_struct& output, photon_data& photon,
                             geometry_struct& geometry, runinfo_struct& runinfo,
                             random_dirty& random_obj)

{
  int i;
  photon_data tmp_photon;
  double save_stellar_weight = 0.0;
  int image_indxs[2];

#ifdef DEBUG_CSP
  if (photon.number == OUTNUM) {
    cout << "in cstp" << endl;
    cout << "photon.path_cur_cells = " << photon.path_cur_cells << endl;
  }
#endif

  // loop over the line-of-sights or albedos or whatever
  for (i = 0; i < output.num_outputs; i++) {
    // setup so the hard stuff is not recomputed for the case where there are
    // multiple outputs, but only 1 observer position
    if ((i == 0) || (geometry.num_observers > 1)) {
      // copy input photon into temporary photon to ensure no change
      tmp_photon = photon;

      // copy the photon birth position to the position location
      int k;
      for (k = 0; k < 3; k++) {
        tmp_photon.position[k] = tmp_photon.birth_position[k];
#ifdef DEBUG_CSP
        if (photon.number == OUTNUM) {
          cout << "pos; k = " << k << " " << tmp_photon.position[k] << endl;
        }
#endif
      }

      // now determine the position indexes of the photon
      determine_photon_position_index_initial(geometry, tmp_photon);

#ifdef DEBUG_CSP
      if (photon.number == OUTNUM) {
        cout << "photon number = " << photon.number << endl;
        cout << "classify_stellar_photon: in stellar_weight = "
             << tmp_photon.stellar_weight << endl;
      }
#endif

      if (geometry.internal_observer == 0) {
        // if averging over the entire 4pi steradians is desired, randomize the
        // observer position
        if (geometry.randomize_observer) {
          geometry.observer_angles[0][i] =
              acos(2.0 * random_obj.random_num() - 1.0);
          geometry.observer_angles[1][i] =
              M_PI * (2.0 * random_obj.random_num() - 1.0);

          compute_observer_trans_matrix(output, geometry, i);
        }
      } else {  // internal observer
        for (int k = 0; k < 3; k++)
          output.outputs[i].observer_position[k] =
              geometry.observer_position[k];
      }
      // determine stellar weight towards observer
      tmp_photon.stellar_weight = stellar_weight_towards_observer(
          tmp_photon, geometry, output.outputs[i].observer_position);

      // save the stellar weight for use when different outputs are for
      // emission/grain types
      save_stellar_weight = tmp_photon.stellar_weight;

#ifdef DEBUG_CSP
      if (photon.number == OUTNUM) {
        cout << "photon.first.tau = " << tmp_photon.first_tau << endl;
        cout << "classify_stellar_photon: stellar_weight = "
             << tmp_photon.stellar_weight << endl;
      }
#endif

      if (geometry.internal_observer == 0) {
        // external observer case
        // transform photon positions so that the line-of-sight is along
        // the z-axis
        rotate_zaxis_for_observer(output.outputs[i].rotate_transform,
                                  tmp_photon);

        // compute x,y angles and image indexs
        double angle;
        for (k = 0; k < 2; k++) {
          angle = atan(tmp_photon.position[k] /
                       (geometry.distance - tmp_photon.position[2]));
          image_indxs[k] = int((1.0 + (angle / geometry.angular_radius)) *
                               output.image_size[k] * 0.5);
          // 	image_indxs[k] = int((1.0 +
          // (tmp_photon.birth_position[k]/(1.1*geometry.radius)))*output.image_size[k]*0.5);
          // check the index is on the image
          if ((image_indxs[k] < 0) ||
              (image_indxs[k] > (output.image_size[k] - 1))) {
            cout << "classify stellar photon: image_indxs[" << k
                 << "] = " << image_indxs[k] << " (beyond image bounds)"
                 << endl;
            cout << "real version of index = "
                 << (1.0 + (angle / geometry.angular_radius)) *
                        output.image_size[k] * 0.5
                 << endl;
            exit(8);
          }
        }
      } else {
        // internal observer case
        double obs_to_birth[3];
        double dist_obs_to_birth = 0.0;
        double dir_obs_to_birth[3];
        int k;
        for (k = 0; k < 3; k++) {
          obs_to_birth[k] =
              output.outputs[i].observer_position[k] - tmp_photon.position[k];
          dist_obs_to_birth += pow(obs_to_birth[k], 2);
        }
        dist_obs_to_birth = sqrt(dist_obs_to_birth);
        for (k = 0; k < 3; k++)
          dir_obs_to_birth[k] = obs_to_birth[k] / dist_obs_to_birth;

        double theta = acos(dir_obs_to_birth[2]);
        double phi = acos(dir_obs_to_birth[0] / sin(theta));
        if ((dir_obs_to_birth[1] / sin(theta)) < 0.) phi *= -1.0;

        image_indxs[0] = int(output.image_size[0] * ((phi + M_PI) / (2.0 * M_PI)));
        image_indxs[1] = int(output.image_size[1] * (1. - dir_obs_to_birth[2]) / 2.0);
      }
#ifdef DEBUG_CSP
      if (photon.number == OUTNUM) {
        cout << "geometry.distance = " << geometry.distance << endl;
        cout << "stellar x,y positions = (" << tmp_photon.position[0] << ","
             << tmp_photon.position[1] << ")" << endl;
        cout << "model radius = " << geometry.radius << endl;
        cout << "stellar photon image indexs = (" << image_indxs[0] << ","
             << image_indxs[1] << ")" << endl;
      }
#endif

    } else {
      // reset the stellar weight for emission/grain output
      tmp_photon.stellar_weight = save_stellar_weight;
    }

    // modify weight by probability the emission was due to a specific
    // grain/emission type only used for dust emission part of dirty
    if (runinfo.dust_thermal_emission && runinfo.do_emission_grain && (i > 0))
      tmp_photon.stellar_weight *= tmp_photon.birth_photon_type_prob[i];

#ifdef DEBUG_CSP
    if (photon.number == OUTNUM) {
      cout << "tmp_photon.stellar_weight = " << tmp_photon.stellar_weight
           << endl;
      cout << "stellar_weight (before) = "
           << output.outputs[i].total_stellar_weight << endl;
      cout << "i = " << i << endl;
    }
#endif
    // update global values
    output.outputs[i].total_num_photons += 1.0;
    output.outputs[i].total_stellar_weight += tmp_photon.stellar_weight;
    output.outputs[i].total_stellar_weight_x2 +=
        pow(tmp_photon.stellar_weight, 2.);

    output.outputs[i].ave_first_tau += tmp_photon.first_tau;
    output.outputs[i].ave_first_tau_x2 += pow(tmp_photon.first_tau, 2.);

    // update the image values
    output.outputs[i].num_stellar_photons_xy(image_indxs[0], image_indxs[1]) +=
        1.0;
    output.outputs[i].stellar_weight_xy(image_indxs[0], image_indxs[1]) +=
        tmp_photon.stellar_weight;
    output.outputs[i].stellar_weight_xy_x2(image_indxs[0], image_indxs[1]) +=
        pow(tmp_photon.stellar_weight, 2.0);

#ifdef DEBUG_CSP
    if (photon.number == OUTNUM) {
      cout << "stellar_weight (after) = "
           << output.outputs[i].total_stellar_weight << endl;
    }
#endif
  }
}
