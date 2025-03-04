/**************************************************************
 * Grain -- Grain Class, header file.
 *
 * History:
 *     Written by:   Putrid, August 2004
 *                   Contact: misselt@as.arizona.edu
 *                   Re-written completely, Sept 2006.
 *
 **************************************************************/

#ifndef _GRAIN_H_
#define _GRAIN_H_

#include <cmath>
#include <fstream>
#include <iostream>

#include "Constants.h"
#include "NumUtils.h"
#include "StringManip.h"

using namespace std;

class Grain
{ // Begin Grain Class definition

  private:
  protected:
    int nsize; // Number of sizes in file.
    int nwave; // Number of wavelengths in file
    int nTemp; // Number of temperatures in calorimetry grid

    string ComponentName;   // name of the grain component
    string OpticalFile;     // File that Cs where read from.
    string CalorimetryFile; // File that calorimetry was read from
    int nelements;          // Number of discrete elements comprising the grain.
    float density;          // Bulk density of grain material
    float SubT;             // Sublimation Temperature of grain material
    vector<float> amu;      // Mass of ith element in amu
    vector<float> stoich;   // Stoichiometric contribution of ith element
    vector<string> atom;    // Name of ith element (eg. Si, C, O, etc...)

    vector<float> size; // Size grid.
    vector<float> mass; // Mass grid.
    vector<float> wave; // Wave grid.

    vector<vector<float> > CAbs;   // Absorption cross sections
    vector<vector<float> > CSca;   // Scattering cross sections
    vector<vector<float> > phFunc; // Scattering phase funtions

    vector<float> Temperature;           // Temperature grid
    vector<float> SpecificHeatCapacity;  // Specific heat capacity (erg/K/gm)
    vector<vector<float> > HeatCapacity; // Heat capacity (erg/K)
    vector<float> SpecificEnthalpy;      // Specific enthalpy (erg/gm)
    vector<vector<float> > Enthalpy;     // Enthalpy (erg)

  public:
    // Constructors/Destructors - empty so we can vectorize with MakeGrain; sigh.
    Grain ();
    ~Grain () {};

    friend class GrainModel;
    // Member Functions.
    // Really the constructor - make a member function so we can easily vectorize
    void MakeGrain (string const &inComponentName, string const &fOpticalConstants, string const &fCalorimetry,
                    vector<float> &MasterWave, vector<float> &MasterSize, string const &Path, float a_min = -1,
                    float a_max = -1);

    // Return the name of the dust component
    inline string
    getComponentName (void)
    {
        return ComponentName;
    }

    // Return the filename of the optical constants
    inline string
    getOpticalFilename (void)
    {
        return OpticalFile;
    }

    // Return the bulk density of the grain material
    inline float
    getDensity (void)
    {
        return density;
    }

    // Return sublimation temperature of the grain
    inline float
    getSubTemp (void)
    {
        return SubT;
    }

    // Return number of grain sizes in object
    inline int
    getNSize (void)
    {
        return nsize;
    }

    // Return number of wavelengths in object
    inline int
    getNWave (void)
    {
        return nwave;
    }

    // Bounds checking needed....

    // Return the vector of sizes
    inline vector<float>
    getSize (void)
    {
        return size;
    }
    inline float
    getSize (int size_id)
    {
        return size[size_id];
    }

    // Return the vector of masses
    inline vector<float>
    getMass (void)
    {
        return mass;
    }

    // Return the vector of wavelengths
    inline vector<float>
    getWave (void)
    {
        return wave;
    }

    // Return a vector with CAbs tabulated at a specific size.
    inline vector<float>
    getCAbs (int size_id)
    {
        return CAbs[size_id];
    }

    // Return a vector with CSca tabulated at a specific size.
    inline vector<float>
    getCSca (int size_id)
    {
        return CSca[size_id];
    }

    // Return a vector with Scattering Phase Function tabulated at a specific
    // size.
    inline vector<float>
    getphFunc (int size_id)
    {
        return phFunc[size_id];
    }

    // Retrun the temperature grid
    inline vector<float>
    getTemperature (void)
    {
        return Temperature;
    }

    // Return the specific heat capacity of the grain material
    inline vector<float>
    getSpecificHeatCapacity (void)
    {
        return SpecificHeatCapacity;
    }
    inline vector<float>
    getHeatCapacity (int _szid)
    {
        return HeatCapacity[_szid];
    }

    // Return the specific enthalpy of the grain material
    inline vector<float>
    getSpecificEnthalpy (void)
    {
        return SpecificEnthalpy;
    }
    inline vector<float>
    getEnthalpy (int _szid)
    {
        return Enthalpy[_szid];
    }

    // Return a vector with CAbs tabulated at a specific size; here, the size
    // can be arbitrary and the function will return an interpolated CAbs.
    vector<float> getCAbs (float a_size);

    // Return a vector with CSca tabulated at a specific size; here, the size
    // can be arbitrary and the function will return an interpolated CSca.
    vector<float> getCSca (float a_size);

    // Return a vector with phFunc tabulated at a specific size; here the size can
    // be arbitrary and the function will return an interpolated phFunc
    vector<float> getphFunc (float a_size);

    // Return a vector with CAbs tabulated at a specific wavelength.
    inline vector<float>
    w_getCAbs (int wave_id)
    {
        vector<float> retvect (nsize);
        for (int i = 0; i < nsize; i++)
            retvect[i] = CAbs[i][wave_id];
        return retvect;
    }

    // Return a vector with CSca tabulated at a specific wavelength.
    inline vector<float>
    w_getCSca (int wave_id)
    {
        vector<float> retvect (nsize);
        for (int i = 0; i < nsize; i++)
            retvect[i] = CSca[i][wave_id];
        return retvect;
    }
    // Return a vector with phFunc tabulated at a specific wavelength.
    inline vector<float>
    w_getphFunc (int wave_id)
    {
        vector<float> retvect (nsize);
        for (int i = 0; i < nsize; i++)
            retvect[i] = phFunc[i][wave_id];
        return retvect;
    }

}; // End Grain class.

#endif
