.. _paramfile:

##############
Parameter File
##############

Basic Properties
----------------

The DIRTY model is run mainly from a single parameter file. This parameter
file is made up of sections that are identified by a line with:

[Section]

Each parameter is specified by a keyword, the "=" sign, and the parameter value
on its own line in the file. The format is:

keyword=value

It is important that there is no space between the end of the keyword, the "=" sign, and value.

Sections
--------

The allowed values are given in [] in the following sections as:

keyword=value [min,max]

.. toctree::
  :maxdepth: 2

  Geometry: Source/Dust Geometry definitions <pfile_geometry.rst>
  Grains: Type of dust grains <pfile_grains.rst>
  Model BookKeeping: Dust properties <pfile_dust.rst>
  Sed: Source spectral energy distribution <pfile_sed.rst>
  Run: Details of the run (I/O, # photons, etc.) <pfile_run.rst>
