#ifndef _GRAINMODEL_H_ 
#define _GRAINMODEL_H_

//*******************************************************************************
  // Grain -- Grain Class, header file.  
  // 
  // VARIABLES: 
  //    TotalNormalization - Grains/Hatom, ie. Sum_i( Int[da_i f_i(a)] ) 
  //    Normalization[] - Normalization of ith Component, ie. Int[da_i f_i(a)]
  // 
  // FUNCTIONS: 
  //    getNormalization - Return TotalNormalization [H^-1]. 
  //    getNormalization(i) - Return Normalization for component i [H^-1].
  //    getCAbsEff     - Get size integrated, component summed absorption cross
  //                     section as a function of wavelength [cm^2].
  //    getCAbsEffNorm - Get size integrated, component summed absorption cross 
  //                     section normalized to number of H atoms [cm^2 H^-1]. 
  //    getCScaEff     - Get size integrated, component summed scattering cross
  //                     section as a function of wavelength [cm^2].
  //    getAsymEff     - Get size integrated, component summed asymmetry  
  //                     parameter as a function of wavelength [Dimensionless].
  //    getAlbedo      - Get albedo (== CScaEff/(CScaEff+CAbsEff) as a function
  //                     of wavelength. 
  //    getTau         - Get optical depth per H column as f(wave) [H^-1].
  //    getTau(W)      - Get optical depth per H column at wavelength W. W must 
  //                     be in [cm]. Returned as [^-1]. 
  //    getMDust       - Get the mass of dust.  [gm H^-1].
  //    getMDust(i)    - Get the mass of dust for component i. [gm H^-1]
  //
  // History: 
  //     Written by:   Putrid, Fall/Winter 2006
  //                   Contact: misselt@as.arizona.edu
  //                   
  //
//*******************************************************************************

#include <iostream>

//#include "Constants.h"
#include "ConfigFile.h"
#include "Grain.h" 
#include "StringManip.h"
#include "EqTemp.h"

using namespace std; 

class GrainModel: protected Grain { 

public:

  // Empty constructor for vectorization in the case that we want to allow 
  // segregation of grains within the model space. 
  GrainModel(); 
  ~GrainModel() {};

  // Create a specific GrainModel object
  void MakeGrainModel ( ConfigFile & ModelConfigFile , vector <float> & MasterWave); 

  inline vector <float> getSizeDistribution ( int _cmp ) { return SizeDistribution[_cmp]; }
  inline vector <float> getCAbsEff( void ) { return CAbsEff; } 
  inline vector <float> getCAbsEffNorm( void ) { 
    vector <float> CAbsEffNorm; 
    transform(CAbsEff.begin(),CAbsEff.end(),back_inserter(CAbsEffNorm),
	      bind2nd(multiplies<float>(),TotalNormalization)); 
    return CAbsEffNorm;
  }
  inline vector <float> getCScaEff( void ) { return CScaEff; } 
  inline vector <float> getAlbedo( void ) { return Albedo; }
  inline vector <float> getphFuncEff( void ) { return phFuncEff; }
  inline vector <float> getWave( void ) { return Wave; } 
  inline vector <float> getEScale( void ) { return EScale; }

  float getTau( float a_wave );  
  inline vector <float> getTau( void ) { return Tau; }
  inline float getMDust( int cmp ) { return DustMass[cmp]; }
  inline float getMDust( void ) { return TotalDustMass; }
  
  inline int getNComp( void ) { return nComp; }
  inline float getNormalization( void ) { return TotalNormalization; }
  inline float getNormalization( int cmp ) { return Normalization[cmp]; }
  //void ComputeEmission ( vector <float> rField );
  // These functions interface with the Grain model to retrieve properties of 
  // individual components of the model. 
  //  - Need bound checking
  inline vector <float> Size ( int _cmp ) { return Component[_cmp].getSize(); }
  inline float Size ( int _cmp, int _szid ) { return Component[_cmp].getSize(_szid); }
  inline vector <float> CAbs ( int _cmp, int _szid ) { return Component[_cmp].getCAbs(_szid); }
  inline int nSize ( int _cmp ) { return Component[_cmp].getNSize(); }
  inline float Density ( int _cmp ) { return Component[_cmp].getDensity(); }
  inline vector <float> CAbs ( int _cmp, float _sz ) { return Component[_cmp].getCAbs(_sz); }
  inline vector <float> wCAbs ( int _cmp, int _wid ) { return Component[_cmp].w_getCAbs(_wid); }

  inline vector <float> CalTemp ( int _cmp ) { return Component[_cmp].getTemperature(); }
  inline vector <float> SpecificHeatCapacity ( int _cmp ) 
    { return Component[_cmp].getSpecificHeatCapacity(); }
  inline vector <float> HeatCapacity ( int _szid, int _cmp ) 
    { return Component[_cmp].getHeatCapacity(_szid); } 
  inline vector <float> SpecificEnthalpy ( int _cmp ) 
    { return Component[_cmp].getSpecificEnthalpy(); }
  inline vector <float> Enthalpy ( int _szid, int _cmp ) 
    { return Component[_cmp].getEnthalpy(_szid); } 

private: 

  int nComp; 
  int nWave; 
  vector <float> Wave; 
  vector <float> EScale; 
  vector <Grain> Component; 
  string ModelName; 
  
  vector <vector<float> > SizeDistribution; 
  vector <float> Temperature;      // Hold a single component temperature at a time...

  //vector <vector<float> > Temperature;      // Hold a single component temperature at a time...
 
  vector <float> Normalization; 
  float TotalNormalization;
  float TotalDustMass; 
  vector <float> DustMass; 

  //vector <float> ComponentTemperature; 

  // Effective quantities; f(wave), integrated and summed over size 
  // distribution for this grain model. 
  vector <float> CAbsEff;   
  vector <float> CScaEff;
  vector <float> Albedo; 
  vector <float> Tau; 
  vector <float> phFuncEff; 

  vector <float> getZDA_sdist(vector <float> coeff, int cmp); 
  
  // Map definitions.  
  // ModelID - Pre-defined allowed model types. 
  // SizeDist - Pre-defined allowed size distribution types
  // iMap will be used to iterate through the maps. 
  map<string,int>::iterator iMap; 
  map<string,int> ModelID; 
  void ModelMapping() { 
    if ( !ModelID.empty()) return; 
    ModelID["MRN"]=0; 
    ModelID["BARE-GR-S"]=1; 
    ModelID["KG-BARE-S"]=2; 
    ModelID["KG-BARE-S1"]=3; 
    ModelID["KG-BARE-S2"]=4; 
    ModelID["KG-BARE-S3"]=5; 
    ModelID["KG-BARE-S4"]=6; 
  }
  map<string,int> SizeDistID; 
  void SizeDistMapping() { 
    if (!SizeDistID.empty()) return; 
    SizeDistID["ZDA"] = 0;
  }

}; 

#endif









