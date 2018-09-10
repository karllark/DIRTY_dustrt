###
SED
###

This sections starts with

**[SED]**

If this section does not exist or has no content, then no SED is used.

**type=bb_file [bb_file, ssp_file]**

The type of file blackbody (bb) or singe stellar population (ssp_file).

**sed_file=filename [any valid filename]**

The SED is usually contained in an ASCII file. The type and name of the file

A bb_file is a file containing the SED of a black body. This is an ASCII file
containing wavelength [microns] and luminosities [ergs/s/Hz]. Such a file can
be created using the IDL program generate_blackbody_input.pro located in the
pro subdirectory.

A ssp_file is a file containing the SED of a single stellar population.
This is an ASCII file containing wavelength [microns] and
luminosities [log(ergs/s/Hz)]. These files are usually generated from
population synthesis models. The scaling of this SED is given by:

**sfr_or_mass=1.0 [0,1e20]**

where a SFR is given for a burst population and a mass is given for constant
star formation. Note that the PEGASE models that have been used with DIRTY have
a burst SFRs normalized to 5e-10 solar masses. The constant SFR PEGASE models
are normalized to 1 solar mass/year. [need to check if this is right].
