// ======================================================================
//   Procedure to get the dust scattering parameters for the wavelength
// currently being computed.
//
// 2006 Nov/KDG - written
// ======================================================================
#include "get_dust_scat_parameters.h"

void
get_dust_scat_parameters (int i, runinfo_struct &runinfo, geometry_struct &geometry)
{
  if ((runinfo.empir_dust == 1) || (runinfo.model_dust == 1))
    {
      geometry.albedo = runinfo.albedo[i];
      geometry.g = runinfo.g[i];
      geometry.tau_to_tau_ref = runinfo.tau_to_tau_ref[i];
      geometry.wavelength = runinfo.wavelength[i];
      geometry.wave_index = i;
    }
  else
    {
      cout << "Code not written get_dust_scat_parameters for non-empir_dust and "
              "non-model_dust"
           << endl;
      exit (8);
    }
}
