// ======================================================================
//   Procedure to classify scattered photon(s) into the output images
// and global totals.
//
// 2006 Apr/KDG - written
// 2007 Dec/KDG - changed denominator of atan from d-z to d 
// 2008 Mar/KDG - added output for emission/grain types
// 2008 May/KDG - changed denominator of atan from d to d-z (see fix in setup_dust_grid_slab.cpp)
// ======================================================================
#include "classify_scattered_photon.h"
//#define DEBUG_CSCP

void classify_scattered_photon (output_struct& output,
				photon_data& photon,
				geometry_struct& geometry,
				runinfo_struct& runinfo)

{
#ifdef DEBUG_CSCP
  cout << "starting csp.." << endl;
#endif

  int i;
  photon_data tmp_photon;
  double save_scat_weight = 0.0;
  int image_indxs[2];

  // loop over the line-of-sights or albedos or whatever
  for (i = 0; i < output.num_outputs; i++) {

#ifdef DEBUG_CSCP
    cout << "output # = " << i << endl;
#endif
    // setup so the hard stuff is not recomputed for the case where there are
    // multiple outputs, but only 1 observer position
    if ((i == 0) || (geometry.num_observers > 1)) {

      // copy input photon into temporary photon to ensure no change
      tmp_photon = photon;

#ifdef DEBUG_CSCP
      cout << "starting scat_weight..." << endl;
#endif
      // determine the probability the photon would have scattered to the observer
      tmp_photon.scat_weight = scattered_weight_towards_observer(tmp_photon, geometry, 
								 output.outputs[i].observer_position);

#ifdef DEBUG_CSCP
      cout << "starting rotate..." << endl;
#endif
      // transform photon positions so that the line-of-sight is along the z-axis 
      rotate_zaxis_for_observer(output.outputs[i].rotate_transform,tmp_photon);

      // take into account the albedo for this scattering
      tmp_photon.scat_weight *= geometry.albedo;

      // save the scattered weight for use when different outputs are for emission/grain types
      save_scat_weight = tmp_photon.scat_weight;

      // compute x,y angles and image indexs
      double angle;
      int k;
      for (k = 0; k < 2; k++) {
	angle = atan(tmp_photon.position[k]/(geometry.distance - tmp_photon.position[2]));
	// don't know why having d-z doesn't work and having d does work, but that seems the case
	// this might mean there is some problem with the tranformation, but I can't figure it out.
	// all this done with Chris in testing the slab model for the SPINR Orion data - KDG 18 Dec 2007
	//angle = atan(tmp_photon.position[k]/(geometry.distance));
	image_indxs[k] = int((1.0 + (angle/geometry.angular_radius))*output.image_size[k]*0.5);
	// check the index is on the image
	if ((image_indxs[k] < 0) || (image_indxs[k] > (output.image_size[k]-1))) {
	  cout << "classify_scattered_photon: image_indxs[" << k << "] = " << image_indxs[k] << 
	    " (beyond image bounds)" << endl;
	  cout << "position = " << tmp_photon.position[k] << endl;
	  cout << "angle = " << angle << endl;
	  cout << "model angular radius = " << geometry.angular_radius << endl;
	  cout << "k = " << k << endl;
	  cout << "image_size[i] = " << output.image_size[k] << endl;
	  cout << "photon # = " << photon.number << endl;
	  int m = 0;
	  cout << "tmp_photon.position[] = ";
	  for (m = 0; m < 3; m++)
	    cout << tmp_photon.position[m] << " ";
	  cout << endl;
	  cout << "geometry.distance = " << geometry.distance << endl;
	  //cout << "updated model angular radius = " << updated_angular_radius << endl;
	  exit(8);
	}
      }

#ifdef DEBUG_CSCP
      cout << "scattered x,y positions = (" << tmp_photon.position[0] << "," << tmp_photon.position[1] << ")" << endl;
      cout << "model radius = " << geometry.radius << endl;
      cout << "scattered photon image indexs = (" << image_indxs[0] << "," << image_indxs[1] << ")" << endl;
#endif
    } else {
      // reset the scattered weight for emission/grain output
      tmp_photon.scat_weight = save_scat_weight;
    }

    // modify weight by probability the emission was due to a specific grain/emission type
    // only used for dust emission part of dirty
    if (runinfo.dust_thermal_emission && runinfo.do_emission_grain && (i > 0))
      tmp_photon.scat_weight *= tmp_photon.birth_photon_type_prob[i];

    // update global values
    output.outputs[i].total_num_scattered_photons += 1.0;
    output.outputs[i].total_scattered_weight += tmp_photon.scat_weight;
    output.outputs[i].total_scattered_weight_x2 += pow(tmp_photon.scat_weight,2.);

    // update the image values
    output.outputs[i].num_photons_xy(image_indxs[0],image_indxs[1]) += 1.0;
    output.outputs[i].scattered_weight_xy(image_indxs[0],image_indxs[1]) += tmp_photon.scat_weight;
    output.outputs[i].scattered_weight_xy_x2(image_indxs[0],image_indxs[1]) += pow(tmp_photon.scat_weight,2.0);

  }

}
