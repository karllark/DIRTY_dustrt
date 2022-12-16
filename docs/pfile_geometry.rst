########
Geometry
########

.. note::
  All this needs to be checked.

This section is started with the single line:

**[Geometry]**

In general, length units are in pc.

Global Details
==============

**distance=100.0  [0,1e38]**

The distance to the model is specified in parsecs. It is important that the distance is large
enough to make sure that the entire model grid is beyond the observer.

**n_obs_angles=1 [1,100]**

The number of observer angles is usually 1, but could be up to 100.
If the number of observer angles is 1, then the single observer (theta,phi) position is input using:

**obs_theta=0. [0.,180.]**

**obs_phi=0. [0.,360.]**

where the angles are input in degrees.
If the number of observer angles is larger than 1, then a file is used.

**obs_file=obs_pos.dat [any string]**

The obs_file has n_obs_angles number of lines with each line containing the observer (theta,phi) pairs separated by a space.
The number of observer (theta,phi) pairs in the file must be the same as the n_obs_angles value.
If the properties of an average line-of-sight are desired, then set:

**randomize_observer=1**

This has the effect of picking a new random observer theta/phi for each photon, averaging over all lines-of-sight.

Source Properties
=================

The source type is set by:

**source_type=stars [stars,diffuse,dexp_disk,pow_sphere]**

stars
-----

If the source_type is "stars", then the number of stars is specified by:

**n_stars=1 [1,100]**

If the number of stars is 1, then the star's position (in parsecs) in the model grid is specified by:

**star_pos_x=0.  [min_x,max_x]**

**star_pos_y=0.  [min_y,max_y]**

**star_pos_z=0.  [min_z,max_z]**

If the number of stars is larger than 1, then the stars' positions are given in a file.

**star_file=star_pos.dat**

The star_file has n_stars number of lines with each line containing 4 numbers:
the stellar positions (x,y,z) in the model and the luminosity of the star.
The luminosity units are not important if you are running a single wavelength only,
just that the relative luminosities between different stars are correct.
The number of star positions and luminosities in the file must be the same as the n_stars value.
If you are running multiple wavelengths and dust emission, then things don't work currently.

diffuse
-------

Info needed.

dexp_disk
---------

A double exponential disk of stars is designated and usually used with the
same type of dust disk. The stellar density falls off as an exponential with
a scale length (xy) and a scale height (z).

**stellar_scalelength=3000.  [0,radius]**

**stellar_scaleheight=300.  [0,radius]**

pow_sphere
----------

An isotropic stellar distribution following the power law r-α. It's a sphere if
the inner_radius is zero; otherwise it's a shell.

**pow_sphere_exponent=0. [-100.,100.]**

**pow_sphere_inner_radius=0. [0.,radius]**

**pow_sphere_outer_radius=1e3 [pow_sphere_inner_radius,radius]**

Dust Distribution
=================

All the dust geometries can have locally varying dust distributions.
The simplest description of this is a two phase dust distribution with the
distribution determined by two parameters, the filling factor of the high
density dust and the density ratio between the low and high density phases.
The clumps can be made spherical by setting the clump type to 'sphere'.
This will create a subgrid for each clump that will roughly resolve the
spherical clump (the rest of the subgrid is filled with the low density dust
- the filling factor is adjusted to account for the lower filling factor of
spherical clumps than cubical clumps in the original gird). For a clump
type that is 'cube', the entire grid cell is filled with high density dust.

The parameters with reasonable values are:

**filling_factor=0.15 [0.,1.]**

**density ratio=0.01 [-0.0,1.0]**

**clump_type=cube [cube,sphere]**

Global Geometries
-----------------

The type of global geometry is picked by the type parameter.

**type=sphere [sphere,shell,slab,dexp_disk,file]**

sphere
~~~~~~

Sphere of dust.  Info needed.

shell
~~~~~

The shell geometry is setup to have a shell that has an evacuated
dust-free inner region. It is sometimes desirable to have a non-instantaneous
ramp up of the dust density at the inner boundary (e.g. circumstellar winds,
etc.). In this case, the very_inner_radius should be set to a value between 0
and the inner_radius. Then, the dust density will ramp up from 0 to the inner
radius value linearly with radius. The radial density profile in the shell is
set by the shell_density_poly where the profile is r^poly. A value of 0 will
provide a uniform density shell. The subdivide_radius parameter provides a
way to provide higher resolution inside of this radius to help resolve quickly
changing shell density profiles as well as a better resolution of the inner
shell boundary. This can be important where the majority of dust is near the
inner boundary (e.g., AGB stars) or where the temperature structure is changing
quickly and large cubic cells are too coarse to resolve the action.
Set subdivide_inner_radius to 0 to disable this option.

**radius=1000. [0.,1e38]**

**very_inner_radius=0. [0.,radius]**

**inner_radius=300. [0.,radius]**

**outer_radius=1000. [0.,radius]**

**subdivide_radius=0. [0.,radius]**

**shell_density_poly=0. [-100.,100.]**

slab
~~~~

.. note::
  This section needs updating for the TRUST slab work.

The slab geometry is setup to have the slab exist in the xy plane with a
set depth in the z plane. The (theta,phi)=(0,0) position of the observer is
on the z-axis. The nonslab density ratio needs to be greater than 0 (but
can be very small) to avoid computational issues when external sources emit
photons which donot intersect the slab. Having some dust in the nonslab
region allows for the scattering code to work without allowing for the
special case where there is no dust along the photon's line-of-sight.

**size_xy=10.  [0.,1e38]**

**size_z=10. [0.,1e38]**

**slab_z1=3. [-size_z/2,size_z/2]**

**slab_z2=4. [-size_x/2,size_z/2]**

size_xy and size_z are the physical sizes of the grid in pc, in the respective directions.
slab_z1/2 are the start/end positions of the slab, in z direction

**nonslab_density_ratio=1e-5 [1e-5,1.]**

The region of the model not part of the slab should have some non-zero
density.

**grid_size=10 [0,1000]**

grid_size is the number of grid cells in the x and y direction.
The number of grid cells in the z direction is calculated as grid_size*int(size_z/size_xy)

dexp_disk
~~~~~~~~~

Double exponential disk for nominally for disk (spiral) galaxy modeling.
The face-on optical depth from the center to infinity along the z-axis is
set by tau which is at the wavelength set by tau_wave. The maximum size of
the disk is set by the radius. The density falls off as an exponential
with a scale length (xy) and a scale height (z). In the z direction
(scale height) the disk is truncated at the vertical height.

**tau=0.5  [0,1000]**

**tau_wave=1.0  [micron, default if not set 0.55]**

**radius=12000.  [0,1e38]**

**dust_scalelength=3000. [0,radius]**

**dust_scaleheight=150. [0,radius]**

**dust_vertical_trunc=1800. [0,radius]**

arbitrary
~~~~~~~~~

Completely arbitrary distributions of dust can be input using the "file" option.
Two files are required to specify the tau/pc in each cell in the model
and the x,y,z coordinates of the cell boundaries.

**type_file_pos=filebase_pos.fits**

**type_file_tau_pc=filebase_tau_ref_per_pc.fits**

The filenames can be any string, the form given above just makes is clear
that the two are associated with each other.

Both files are FITS format files where the main grid is in the 1st
hdu (header data unit, 1st is also known as the primary hdu) and
any subgrids are in subsequent hdus. The 1st hdu should have the
LONG (datatype) keyword "GRDDEPTH" set to the maximum depth of the
grid (e.g., for a single main grid, "GRDEPTH" = 1, for a main grid
with cells subdivided would have "GRDDEPTH" = 2, for a main grid with
subdivided cells and with some of the subdivided cells having subdivided
cells would be "GRDDEPTH" = 3, etc.).

For subgrids, the cell in the main grid that the subgrid is subdividing
should be filled with the negative of the subgrid number (e.g., the 1st
subgrid is numbered 2 and the cell in the main grid should have a value
of -2, the 2nd subgrid is numbered 3 and the cell in the main grid has
a value of -3, etc.). In addition, the header of each subgrid should
include the INT (datatype) keyword "PAR_GRID" with the index of the of
the grid where of the cell it subdivides (e.g., subgrids of cells in
the main grid will have "PAR_GRID" = 1).

A grid tau_ref_per_pc value of -0.5 designates that grid cells outside
of this one are also filled with -0.5 and the model effectively stops
(no more photon propagation).

The tau_ref_per_pc file contains a cube in each hdu where the cube is o
rdered x,y,z. The values in each cell should be FLOAT (datatype).

The pos file contains a (n+1)x3 image in each hdu where the image is
ordered (x,y,z positions)x3. The positions give the edges of each cell
in units of parsecs. The positions should be given as DOUBLE (datatype).
Since it is possible that one of the 3 dimensions has a different length
than the other two, the n is the largest of the dimension of the cube.
For the dimensions that are smaller, the rest of the values should be
filled with zeros.
