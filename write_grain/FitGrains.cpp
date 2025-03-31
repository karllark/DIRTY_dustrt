// *******************************************************************************
// FitGrains -- FitGrains Class
//
// Requires:
//     ConfigFile
//     Grain
//     StringManip
//
// History:
//     Written by:
// *******************************************************************************

#include "FitGrains.h"

FitGrains::FitGrains () {}

void
FitGrains::MakeFitGrains (ConfigFile &_mbkcf, vector<float> &MasterWave)
{

  string _thisSection;
  string _thisComponentName;
  string _thisCrossSectionFile;
  string _thisCalorimetryFile;
  string _TopLevelPath;

  string _CrossSectionsPath;
  string _CalorimetryPath;
  string _ModelDefinitionPath;
  string _FullModelDefinitionFile;

  float _thisA_min;
  float _thisA_max;

  vector<float> _MasterSize;

  // typedef vector <float>::iterator iter;
  // iter iCAbsEff,iCScaEff,iphFuncEff,iTau,iAlbedo;
  // iter ithisCAbs,ithisCSca,ithisphFunc,ithisSize,ithisSizeDist;

  // Make the wavelength vector for the GrainModel the input MasterWave.
  // Generate an energy scale at the same time.
  // copy(MasterWave.begin(),MasterWave.end(),back_inserter(Wave));
  // nWave = Wave.size();
  // vector <float>::iterator _itb,_ite,_it;
  //_itb=Wave.begin();
  //_ite=Wave.end();
  // for (_it=_itb;_it!=_ite;++_it) EScale.push_back(Constant::PLANCKLIGHT/(*_it));

  // Take input size tabulation as gospel for each grain material.

  _MasterSize.push_back (-1);

  // Model Book Keeping.  Paths, Model definition configuration file.
  // Top Level:
  _TopLevelPath = _mbkcf.SValue ("Model Book Keeping", "Path to Dust Properties");
  // CrossSections:
  _CrossSectionsPath = _mbkcf.SValue ("Model Book Keeping", "Cross Section SubDir");
  // Calorimetry:
  _CalorimetryPath = _mbkcf.SValue ("Model Book Keeping", "Calorimetry SubDir");
  // Model Definiiton:
  _ModelDefinitionPath = _mbkcf.SValue ("Model Book Keeping", "Model SubDir");
  // Model Name:
  ModelName = _mbkcf.SValue ("Model Book Keeping", "Model Name");

  // YES! Check and make sure model definition file exists.
  _FullModelDefinitionFile = _TopLevelPath + _ModelDefinitionPath + ModelName;
  if (!StringManip::FileExists (_FullModelDefinitionFile))
    { // NO!
      cout << "****************************************************" << endl;
      cout << "Fit Grains Definition File " << _FullModelDefinitionFile << " not found." << endl;
      cout << "****************************************************" << endl;
      exit (8);
    }
  // YES! Proceed with model definition.

  // Create a ConfigFile Object for the actual model definition.
  ConfigFile _mcf (_FullModelDefinitionFile);

  // Get number of components and resize member vectors as necessary.
  nComp = _mcf.IValue ("Model", "Number of Components");
  Component.resize (nComp);

  // Get each component.
  //  - Wavelength grid will be set by MasterWave.
  //  - Size grid will be set by the Cross Section file for each component.

  // Construct the Size distribution for each component.
  for (int cmp = 0; cmp < nComp; cmp++)
    { // Loop over all components

      _thisA_min = -1;
      _thisA_max = -1;

      // Construct section key for this component.
      // Assume key is of the form "Comp i" where i=1,2,...,nComp
      _thisSection = "Component " + StringManip::vtos (cmp + 1);

      // Get the name of this component
      _thisComponentName = _mcf.SValue (_thisSection, "Component Name");

      // Get grain property definitions.
      _thisCrossSectionFile = _mcf.SValue (_thisSection, "Cross Sections");
      _thisCalorimetryFile = _mcf.SValue (_thisSection, "Calorimetry");

      // Populate Component[i] with the appropriate Grain object.
      Component[cmp].MakeGrain (_thisComponentName, _CrossSectionsPath + _thisCrossSectionFile,
                                _CalorimetryPath + _thisCalorimetryFile, MasterWave, _MasterSize,
                                _TopLevelPath, _thisA_min, _thisA_max);

    } // End Component loop.
}

vector<float>
FitGrains::Wave (int cmp)
{
  return Component[cmp].getWave ();
}

int
FitGrains::nWave (int cmp)
{
  return Component[cmp].getNWave ();
}

int
FitGrains::getNComp (void)
{
  return nComp;
}

string
FitGrains::getModelName (void)
{
  return ModelName;
}

// These functions are 'pass through' to the relevant Grain objects
string
FitGrains::ComponentName (int cmp)
{
  return Component[cmp].getComponentName ();
}

string
FitGrains::OpticalFilename (int cmp)
{
  return Component[cmp].getOpticalFilename ();
}

int
FitGrains::nSize (int _cmp)
{
  return Component[_cmp].getNSize ();
}

vector<float>
FitGrains::Size (int _cmp)
{
  return Component[_cmp].getSize ();
}

float
FitGrains::Size (int _cmp, int _szid)
{
  return Component[_cmp].getSize (_szid);
}

vector<float>
FitGrains::CAbs (int _cmp, int _szid)
{
  return Component[_cmp].getCAbs (_szid);
}

vector<float>
FitGrains::CSca (int _cmp, int _szid)
{
  return Component[_cmp].getCSca (_szid);
}

vector<float>
FitGrains::phFunc (int _cmp, int _szid)
{
  return Component[_cmp].getphFunc (_szid);
}

float
FitGrains::Density (int _cmp)
{
  return Component[_cmp].getDensity ();
}

vector<float>
FitGrains::CAbs (int _cmp, float _sz)
{
  return Component[_cmp].getCAbs (_sz);
}

vector<float>
FitGrains::CSca (int _cmp, float _sz)
{
  return Component[_cmp].getCSca (_sz);
}

vector<float>
FitGrains::phFunc (int _cmp, float _sz)
{
  return Component[_cmp].getphFunc (_sz);
}

vector<float>
FitGrains::wCAbs (int _cmp, int _wid)
{
  return Component[_cmp].w_getCAbs (_wid);
}

vector<float>
FitGrains::CalTemp (int _cmp)
{
  return Component[_cmp].getTemperature ();
}

vector<float>
FitGrains::SpecificHeatCapacity (int _cmp)
{
  return Component[_cmp].getSpecificHeatCapacity ();
}

vector<float>
FitGrains::HeatCapacity (int _szid, int _cmp)
{
  return Component[_cmp].getHeatCapacity (_szid);
}

vector<float>
FitGrains::SpecificEnthalpy (int _cmp)
{
  return Component[_cmp].getSpecificEnthalpy ();
}

vector<float>
FitGrains::Enthalpy (int _szid, int _cmp)
{
  return Component[_cmp].getEnthalpy (_szid);
}

// float FitGrains::getTau( float a_wave );
// vector <float> FitGrains::getTau( void ) { return Tau; }
