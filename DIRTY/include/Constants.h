#ifndef _CONSTANTS_H
#define _CONSTANTS_H


namespace Constant {
 
  // "Physical" constants
  //   updated for TRUST values 12 Jun 2015
  const double LIGHT=2.99792458e10;            // Speed of light cm sec^-1
  const double PI=3.14159265358979323846;      // Pi
  const double BOLTZMAN=1.3806488e-16;         // Boltzman Constant, erg K^-1
  const double PLANCK=6.62606957e-27;          // Planck Constant erg sec
  const double SOL_RAD=6.955e+10;              // Solar radius in cm.
  const double SOL_MASS=1.9891e+33;            // Solar mass in gm.
  const double HMASS_AMU=1.0078250321;         // Hydrogen mass in atomic mass units.
  const double HMASS_CGS=1.67353263e-24;       // Hydrogen mass in gm. 

  // Conversion factors and common constant combinations.
  const double PLANCKLIGHT=(PLANCK)*(LIGHT);   // hc in erg cm 
  const double IPLANCKLIGHT=1.0/(PLANCKLIGHT); // 1/hc
  const double FPI=4.0*(PI);                   // 4pi
  const double FPI_HC=FPI*IPLANCKLIGHT;        // 4pi/hc
  const double FPI_HC2=FPI_HC*IPLANCKLIGHT;    // 4pi/(hc)^2
  const double CGS_AMU=6.02214151e+23;         // gm --> Atomic mass units
  const double AMU_CGS=1.66053886e-24;         // Atomic mass units --> gm
  const double ANG_UM=1.0e-4;                  // Angstrom --> um
  const double ANG_CM=1.0e-8;                  // Angstrom --> cm
  const double UM_CM=1.0e-4;                   // um --> cm
  const double CM_UM=1.0e4;                    // cm --> um
  const double AU_CM=1.49597871e+13;            // AU --> cm
  const double PC_CM=3.08567758e+18;            // Parsec --> cm
  const double JY_CGS=1.0e-23;                 // Constants without wave conversion.
  const double MKS_CGS=10.0; 
  const double U_J=(FPI)/(LIGHT);              // Convert Energy Density to Field
  const double VFAC=(FPI)/3.0;                 // Volume factor (4/3)*pi 
  const double ERG_EV=6.24150974e11;           // Ergs -> electron volts


}

#endif
