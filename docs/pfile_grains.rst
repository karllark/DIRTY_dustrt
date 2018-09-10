###########
Dust Grains
###########

This section is started with the single line:

**[Dust Grains]**

Model Type
==========

**type=dust_model [dust_model, single_wavelength]**

If dust_model, then use a wavelength grid (see below)

if single_wavelength, then specify the needed details

.. note::
  More options possible, need to add

**wavelength=0.55 [??, ??]**

Wavelength in microns

**albedo=0.5 [0.0, 1.0]**

Single scattering albedo.

**g=0.7 [-1.0, 1.0]**

Henyey-Greenstein phase function asymmetry.


Wavelength Grid
===============

**wave_type=res [res,file]**

The type of wavelength grid can be set to be either res (by resolution) or file (by an ASCII file).

**wave_type=res**
**wave_min=0.0912**
**wave_max=1000.**
**wave_resolution=10.**

If the wave_type is res, then the wavelength grid is defined by a min/max wavelength (in microns) and desired resolution.

**wave_file=wave_grid_pah_opt.dat  [string]**

If the wave_type is file, then the wavelength grid is set by a a file where
each line in the file is a wavelength point given in microns.

As an example, a wavelength grid optimized for normal ISM dust is given by
wave_grid_pah_opt.dat. This is a grid at a resolution of 10 from
0.0912-1000 microns with an embedded higher resolution grid (resolution = 50)
from (roughly) 3 to 20 microns.
