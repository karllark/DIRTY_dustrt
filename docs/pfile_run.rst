###
Run
###

This section is started with the single line:

**[Run]**

Run Parameters
==============

**verbose=1 [0,1,2]**

The verbosity level (how much text is printed to the screen)
with the higher numbers representing more info.

**num_photons=1e5 [1,1e10]**

The number of photons to run at each wavelength point.

**do_dust_emission=1 [0,1]**

To have the code compute the dust thermal emission (equilibrium and non-equilibrium), set to 1.

**energy_conserve_target=0.05 [0.,1.]**

The conservation target for the dust emission where 1 is no conservation and
0 is exact conservation. This controls the number of iterations needed to
account of dust self-absorption. There is a maximum number of iterations
that is currently hard coded to 10.

Output Details
==============

**do_global_output=1 [0,1]**

**do_image_output=0 [0,1]**

The amount of output can be a FITS table of global luminosities and/or images
at each wavelength. At least one of these two outputs must be set to 1.

An image is created for each wavelength in your wavelength grid; the number
preceding 'um' in the filename indicates the wavelength in microns of each
image. Furthermore, images are split into dust emission images and stellar
emission images. To get the final total image, you need to combine the stellar
radiative transfer with the dust emission. The image naming scheme is as follows:

<filename>_de_ge?_w*_*um.fits

These are the images for the dust emission. The ge1 files are the total dust
emission, the subsequent (ge2, ge3 files) are it split up between grain types
and emission (equilibrium/non-equilibrium).

<filename>_w*_*um.fits

These are the images from the stellar radiative transfer before the dust emission.

The size of the output images (at each wavelength) can be given two ways.

**output_image_size=201 [1,5e4]**

For square images.  Usually used for external observers where the output image is 
a tangent projection.

**output_image_x_size=201 [1,5e4]**

**output_image_y_size=201 [1,5e4]**

For rectangular images.  This option may be particularly useful for internal observers.
For internal observers, the image is given in (phi, theta) angles.  Specifically the 
x-axis=phi and y-axis=cos(theta).

**output_filebase=standard_sphere [any string]**

The size base filename of the images.

**output_type=ratio [ratio,sbrightness]**

The units of the output images are either a ratio to the source luminosity
(basically unitless) or as a surface brightness (cm^-2 sr^-1)

**output_model_grid=1 [0,1]**

Parts of the model grid can be output as a multi-extension FITS files.
 Currently, the files output give the tau/pc, the positions of the cell edges,
 the radiation field, the number of H atoms in each cell, and the wavelength grid.

**do_emission_type=1 [0,1]**

The output of the dust emission can be separated into emission from grains in
equilibrium and non-equilibrium.

**do_emission_grain=1 [0,1]**

The output of the dust emission can be separated into emission from grains
of different compositions (graphite, silicate, PAH, etc.).
