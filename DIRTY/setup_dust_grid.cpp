// ======================================================================
// set up the general properties of the dust grid and call the right
// subroutine depending on the type of geometry
//
// 30 Jun 2006 KDG - modified to put spherical portion in seperate function
// 21 Nov 2006 KDG - added diffuse source option
// 27 Feb 2007 KDG - added slab geometry
// 13 Apr 2007 KDG - moved geometry.solid_angle definition to where the
//                   source emission is defined
// KDG 16 Jun 2008 - added subdivide grid cells & fixed max_grid_depth
// calculation
// ======================================================================
#include "setup_dust_grid.h"

void
setup_dust_grid (ConfigFile &param_data, geometry_struct &geometry, photon_data &photon, random_dirty &random_obj)

{
  // model distance
  geometry.distance = param_data.FValue ("Geometry", "distance");
  check_input_param ("distance", geometry.distance, 0.0, 1e38);

  // setup the dust grid proper depending on type
  geometry.type = param_data.SValue ("Geometry", "type");
  //  if (strcmp(geometry.type.c_str(),"sphere")) {
  if (geometry.type == "sphere")
    {
      setup_dust_grid_sphere (param_data, geometry, random_obj);
    }
  else if (geometry.type == "shell")
    {
      setup_dust_grid_shell (param_data, geometry, random_obj);
    }
  else if (geometry.type == "slab")
    {
      setup_dust_grid_slab (param_data, geometry, random_obj);
    }
  else if (geometry.type == "slab_trust")
    {
      setup_dust_grid_slab_trust (param_data, geometry, random_obj);
    }
  else if (geometry.type == "dexp_disk")
    {
      setup_dust_grid_dexp_disk (param_data, geometry, random_obj);
    }
  else if (geometry.type == "file")
    {
      setup_dust_grid_file (param_data, geometry);
    }
  else
    {
      cout << "Setup for input geometry type (" << geometry.type << ") not found [NEW CODE NEEDED]." << endl;
      exit (8);
    }

  // verify that the dust grid has been constructed correctly
  verify_dust_grid (geometry);

  // get the source type
  geometry.source_type = param_data.SValue ("Geometry", "source_type");

  if (geometry.source_type == "stars")
    {
      // set the solid angle (stars emit isotropically)
      geometry.solid_angle = 4. * M_PI;

      geometry.new_photon_source_type = NEW_PHOTON_DISCRETE_STARS;

      // number of stars
      geometry.num_stars = param_data.IValue ("Geometry", "n_stars");
      check_input_param ("n_stars", geometry.num_stars, 1, MAX_MULTIPLE_STARS);

      if (geometry.num_stars == 1)
        {
          // read in star position
          geometry.star_positions[0][0] = param_data.FValue ("Geometry", "star_pos_x");
          check_input_param ("star_pos_x", geometry.star_positions[0][0], -1 * geometry.radius, geometry.radius);
          geometry.star_positions[1][0] = param_data.FValue ("Geometry", "star_pos_y");
          check_input_param ("star_pos_y", geometry.star_positions[1][0], -1 * geometry.radius, geometry.radius);
          geometry.star_positions[2][0] = param_data.FValue ("Geometry", "star_pos_z");
          check_input_param ("star_pos_z", geometry.star_positions[2][0], -1 * geometry.radius, geometry.radius);
          // current set this to a default value (maybe use input for later)
          geometry.total_source_luminosity = 1.;
        }
      else
        {
          // get the filename
          string multiple_stars_filename = param_data.SValue ("Geometry", "star_file");
          // check that the file exists
          ifstream multiple_stars_file (multiple_stars_filename.c_str ());
          if (multiple_stars_file.fail ())
            {
              cout << "Multiple star file w/ positions (" << multiple_stars_filename << ") does not exist." << endl;
              exit (8);
            }
          multiple_stars_file.close ();

          // read in the positions
          vector<double> pos_x, pos_y, pos_z, lum;
          DataFile (multiple_stars_filename, pos_x, pos_y, pos_z, lum);

          // check we got the number of positions we expected
          if (int (pos_x.size ()) != geometry.num_stars)
            {
              cout << "The number of positions in multiple star file (" << pos_x.size () << ")" << endl;
              cout << "does not match the number expected (" << geometry.num_stars << ")" << endl;
              exit (8);
            }

          // test the results and put them into the appropriate locations
          int i;
          for (i = 0; i < geometry.num_stars; i++)
            {
              check_input_param ("multiple star pos_x", pos_x[i], -1 * geometry.grids[0].phys_grid_size[0] / 2.,
                                 geometry.grids[0].phys_grid_size[0] / 2.);
              check_input_param ("multiple star pos_y", pos_y[i], -1 * geometry.grids[0].phys_grid_size[1] / 2.,
                                 geometry.grids[0].phys_grid_size[1] / 2.);
              check_input_param ("multiple star pos_z", pos_z[i], -1 * geometry.grids[0].phys_grid_size[2] / 2.,
                                 geometry.grids[0].phys_grid_size[2] / 2.);
              check_input_param ("multiple star luminosity", lum[i], 0.0, 1e40);

              geometry.star_positions[0][i] = pos_x[i];
              geometry.star_positions[1][i] = pos_y[i];
              geometry.star_positions[2][i] = pos_z[i];
              geometry.star_positions[3][i] = lum[i];
              geometry.total_source_luminosity += lum[i];
              geometry.star_positions[4][i] = geometry.total_source_luminosity;
            }
          // now divide the 5th position by the total luminosity to get a
          // number between 0 and 1 for the easy random generation of photons
          for (i = 0; i < geometry.num_stars; i++)
            {
              geometry.star_positions[4][i] /= geometry.total_source_luminosity;
              //       cout << i << " " << geometry.star_positions[4][i] <<
              //       endl;
            }
        }
    }
  else if (geometry.source_type == "diffuse")
    {
      // set the solid angle
      // assuming the diffuse field convers 4*pi sr
      // if not, then this should be set to the angular area of the emitted
      // diffuse field this would be the case for a diffuse field illuminating
      // the source from just one direction
      //   not trivial given that the solid angle would have to be fairly small
      //   for this to work found this during testing of the isotropic field
      //   (scattered/input luminosity > 1)
      geometry.solid_angle = 4. * M_PI;

      // get the filename
      string source_filename = param_data.SValue ("Geometry", "source_file");

      if (source_filename == "isotropic")
        {
          geometry.new_photon_source_type = NEW_PHOTON_DIFFUSE_ISOTROPIC;
          // current set this to a default value (maybe use input for later)
          geometry.total_source_luminosity = 4.0 * M_PI;
        }
      else
        {
          geometry.new_photon_source_type = NEW_PHOTON_DIFFUSE_FILE;
          // check that the file exists
          ifstream source_file (source_filename.c_str ());
          if (source_file.fail ())
            {
              cout << "Diffuse source file (" << source_filename << ") does not exist." << endl;
              exit (8);
            }
          source_file.close ();

          // read in the radiation field
          vector<double> theta, phi, intensity, sr_area;
          DataFile (source_filename, theta, phi, intensity, sr_area);

          // test the results and put them into the appropriate locations
          int i;
          // double sum_lum = 0.0;  // running sum of the luminosity
          for (i = 0; i < int (theta.size ()); i++)
            {
              check_input_param ("diffuse source theta", theta[i], -1 * M_PI / 2., M_PI / 2.);
              check_input_param ("diffuse source phi", phi[i], 0.0, 2. * M_PI);
              check_input_param ("diffuse source intensity", intensity[i], 0, 1e40);
              check_input_param ("diffuse source sr area", sr_area[i], 0., 4. * M_PI);

              geometry.diffuse_source_theta.push_back (theta[i]);
              geometry.diffuse_source_phi.push_back (phi[i]);
              geometry.diffuse_source_intensity.push_back (intensity[i]);
              geometry.total_source_luminosity += intensity[i] * sr_area[i];
              geometry.diffuse_source_sum_intensity.push_back (geometry.total_source_luminosity);
            }
          // now divide the sum_intensity by the total luminosity to get a
          // number between 0 and 1 for the easy random generation of photons
          for (i = 0; i < int (theta.size ()); i++)
            {
              geometry.diffuse_source_sum_intensity[i] /= geometry.total_source_luminosity;
            }

          // input source intensity in ergs s^-1 A^-1 cm^-2 sr^-1
          // need to get ride of cm^-2 and sr^-1 for total luminosity

          // multiply by the size of each cell in sr (assumed to be equal in
          // size) no longer needed as this the size of each cell is explicitly
          // in calculation geometry.total_source_luminosity
          // *= 4.*M_PI/theta.size();
          //       cout << geometry.total_source_luminosity << endl;
        }
      if (geometry.type == "sphere")
        {
          // muliple by the size of the nebula in cm^2
          geometry.total_source_luminosity *= M_PI * pow (0.95 * geometry.radius * Constant::PC_CM, 2.0);
          //      cout << geometry.total_source_luminosity << endl;
        }
      else
        {
          cout << "diffuse source from file not setup for anything but sphere" << endl;
          cout << "geometry.type = " << geometry.type << endl;
          cout << "switch to sphere or write some more code" << endl;
          exit (8);
        }
    }
  else if (geometry.source_type == "dexp_disk")
    {
      geometry.solid_angle = 4. * M_PI;
      geometry.new_photon_source_type = NEW_PHOTON_DEXP_DISK;
    }
  else if (geometry.source_type == "pow_sphere")
    {
      geometry.solid_angle = 4. * M_PI;
      geometry.new_photon_source_type = NEW_PHOTON_POW_SPHERE;

      double alpha = param_data.FValue ("Geometry", "pow_sphere_exponent");
      // if (isnan(alpha)) alpha = 0.;
      check_input_param ("pow_sphere_exponent", alpha, -100., 100.);
      geometry.pow_sphere_exponent = alpha;

      double rmax = param_data.FValue ("Geometry", "radius");
      double r1 = param_data.FValue ("Geometry", "pow_sphere_inner_radius");
      // if (isnan(r1)) r1 = 0.;
      check_input_param ("pow_sphere_inner_radius", r1, 0., rmax);
      geometry.pow_sphere_inner_radius = r1;

      double r2 = param_data.FValue ("Geometry", "pow_sphere_outer_radius");
      // if (isnan(r2)) r2 = rmax;
      check_input_param ("pow_sphere_outer_radius", r2, r1, rmax);
      geometry.pow_sphere_outer_radius = r2;

      if (alpha == -3.)
        {
          cout << "pow_sphere: exponent == -3 case not implemented (need "
                  "logarithm)"
               << endl;
          exit (8);
        }
      else if ((alpha < -3.) && (r1 == 0.))
        {
          cout << "pow_sphere: exponent < -3 and outer_radius == 0 (infinity "
                  "at "
                  "origin)"
               << endl;
          exit (8);
        }

      geometry.pow_sphere_constant1 = pow (r2, 3.0 + alpha) - pow (r1, 3.0 + alpha);
      geometry.pow_sphere_constant2 = pow (r1, 3.0 + alpha);
      geometry.pow_sphere_constant3 = 1.0 / (3.0 + alpha);
    }
  else
    {
      cout << "Setup for input source type (" << geometry.source_type << ") not found [NEW CODE NEEDED]." << endl;
      exit (8);
    }

  // set the observer details

  // number of observers (internal or external)
  //  poorly named input parameter name (now that internal observers are
  //  possible)
  geometry.num_observers = param_data.IValue ("Geometry", "n_obs_angles");
  check_input_param ("n_obs_angles", geometry.num_observers, 1, 100);

  // check if an internal observer position is set
  geometry.internal_observer = param_data.IValue ("Geometry", "internal_obs");
  if (geometry.internal_observer == -99)
    geometry.internal_observer = 0;
  check_input_param ("internal_obs", geometry.internal_observer, 0, 1);

  if (geometry.num_observers == 1)
    {
      // external observer: read in observer angles and convert to radians
      if (geometry.internal_observer == 0)
        {
          geometry.observer_angles[0][0] = param_data.FValue ("Geometry", "obs_theta");
          check_input_param ("obs_theta", geometry.observer_angles[0][0], 0., 180.);
          geometry.observer_angles[0][0] *= M_PI / 180.;

          geometry.observer_angles[1][0] = param_data.FValue ("Geometry", "obs_phi");
          check_input_param ("obs_theta", geometry.observer_angles[1][0], 0., 360.);
          geometry.observer_angles[1][0] *= M_PI / 180.;
        }
      else
        { // internal observer
          geometry.observer_position[0] = param_data.FValue ("Geometry", "obs_x");
          check_input_param ("obs_x", geometry.observer_position[0], -1.0 * geometry.radius, geometry.radius);
          geometry.observer_position[1] = param_data.FValue ("Geometry", "obs_y");
          check_input_param ("obs_y", geometry.observer_position[0], -1.0 * geometry.radius, geometry.radius);
          geometry.observer_position[2] = param_data.FValue ("Geometry", "obs_z");
          check_input_param ("obs_z", geometry.observer_position[0], -1.0 * geometry.radius, geometry.radius);
        }
    }
  else
    {
      // multiple internal observers not implemented (yet)
      if (geometry.internal_observer == 1)
        {
          cout << "multiple internal observers not implemented (yet)" << endl;
          exit (8);
        }

      // get the filename
      string multiple_obs_filename = param_data.SValue ("Geometry", "obs_file");
      // check that the file exists
      ifstream multiple_obs_file (multiple_obs_filename.c_str ());
      if (multiple_obs_file.fail ())
        {
          cout << "Multiple obs file w/ positions (" << multiple_obs_filename << ") does not exist." << endl;
          exit (8);
        }
      multiple_obs_file.close ();

      // read in the positions
      vector<double> obs_theta, obs_phi;
      DataFile (multiple_obs_filename, obs_theta, obs_phi);

      // check we got the number of positions we expected
      if (int (obs_theta.size ()) != geometry.num_observers)
        {
          cout << "The number of directions in multiple observers file (" << obs_theta.size () << ")" << endl;
          cout << "does not match the number expected (" << geometry.num_observers << ")" << endl;
          exit (8);
        }

      // test the results and put them into the appropriate locations
      int i;
      for (i = 0; i < geometry.num_observers; i++)
        {
          check_input_param ("multiple obs obs_theta", obs_theta[i], 0.0, 180.0);
          check_input_param ("multiple obs obs_theta", obs_phi[i], 0.0, 360.0);

          geometry.observer_angles[0][i] = obs_theta[i] * M_PI / 180.;
          geometry.observer_angles[1][i] = obs_phi[i] * M_PI / 180.;
        }
    }

  // randomize observer position
  geometry.randomize_observer = param_data.IValue ("Geometry", "randomize_observer");
  if (geometry.randomize_observer == -99)
    geometry.randomize_observer = 0;
  check_input_param ("ranomize_observer", geometry.randomize_observer, 0, 1);

  // setup the photon structure with positions for the maximum grid depth
  vector<int> one_index (3);
  int i = 0;
  for (i = 0; i < geometry.max_grid_depth; i++)
    {
      photon.position_index.push_back (one_index);
      photon.grid_number.push_back (long (0));
    }

  // add in the initialization for the photon path variables
  photon.path_max_cells = 10;
  vector<int> tmp_index (photon.path_max_cells);
  for (i = 0; i < 4; i++)
    photon.path_pos_index.push_back (tmp_index);
  photon.path_tau.resize (photon.path_max_cells);
}
