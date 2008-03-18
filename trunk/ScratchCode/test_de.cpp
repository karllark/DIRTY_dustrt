#include "ISRF.h"
#include "StringManip.h"
#include "Grain.h"
#include "EqTemp.h"         // For equilibrium heating
#include "TimeEvolution.h"  // Time evolution
#include "NumUtils.h"       // Numerical stuff - namespaced so can use with older defs. 
#include "mt_ran.h"         //  
#include "GrainModel.h"   

using namespace std; 

int main (int argc, char * argv[]) { 

  // ****************************************************************************
  // Test making a grain. 
  // These determine whether things will be interpolated. 
  // one element, -1 vectors mean take the grids defined in the cross section 
  // file.  Pass vectors tabulated where you want data as MW, MS [IN CGS, IE CM 
  // FOR BOTH!] to interpolate the cross section file onto a new grid. 
  vector <float> MW(1); 
  vector <float> MS(1); 
  MW[0]=-1; 
  MS[0]=-1; 

  // Number of components - only going to fill in one here, but this illustrates
  // the vector nature of the Grain object.
  int ncomp=2; 
  // Declare a vector of grain object. 
  vector <Grain> MyGrain(ncomp); 
  // Instantiate a Grain object and stick it into the first element of the vector. 
  MyGrain[0].MakeGrain("CrossSections/PAHneu_30_mod.srt",
		       "Calorimetry/Graphitic_Calorimetry.dat",
		       MW,MS,"/home/misselt/Science/Dust/OpticalProperties/");
  // Get the wavelenght grid MyGrain[0] is defined on. 
  vector <float> thisWave = MyGrain[0].getWave(); 
  int nWave = MyGrain[0].getNWave(); 
  // Get the size grid MyGrain[0] is defined on. 
  vector <float> thisSize = MyGrain[0].getSize();
  int nSize = MyGrain[0].getNSize();
  // Get an absorption cross section for this grain - 0 is for the first size
  // in the grid.  
  vector <float> thisCAbs(nWave);
  thisCAbs = MyGrain[0].getCAbs(0);
  // Overloaded getCAbs() - with a float size, will either find an exact match,
  // or interpolate from the Grain object size grid. 
  // vector <float> thatCAbs = MyGrain[0].getCAbs((float)0.15e-4); 
  // For the above case, Grain will return an error since PAHs are not defined 
  // to that large a size...
  // See Grain.h/cpp for details on other member functions. 
  // ****************************************************************************

  // ****************************************************************************
  //  How to use the GrainModel class. 
  // First, we'll use a master wavelength grid from the previously defined 
  // Grain object. 
  // MasterWave == thisWave == MyGrain[0].getWave();
 
  // Declare a GrainModel object (single, but this is also vectorizable should 
  // we eventually want to have different grains models in different grids. 
  GrainModel thisGrainModel;
  // Hard coded GrainModel configuration. 
  string aConfigFile="Craptastic.cfg";
  ConfigFile thisConfig(aConfigFile);
  // Instiate a GrainModel - requires a configuration file which will tell the 
  // class where to get the model info and the wavelength grid we want things 
  // defined on. 
  // All the work goes on behind this line. 
  thisGrainModel.MakeGrainModel(thisConfig,thisWave); 
  // Get the effective (size integrated, summed) grain quantities. 
  vector <float> thiseffcabs=thisGrainModel.getCAbsEff(); 
  vector <float> thiseffcsca=thisGrainModel.getCScaEff();
  vector <float> thiseffphfunc=thisGrainModel.getphFuncEff();

  // Define an albedo vector and compute the albedo 
  // Albedo maybe something we should store in the GrainModel object...
  vector <float> albedo; 
  for (int i=0;i<nWave;i++) { 
    albedo.push_back(thiseffcsca[i]/(thiseffcsca[i]+thiseffcabs[i])); 
    cout << thisWave[i] << " " << albedo[i] << " " << thiseffphfunc[i] << endl; 
  }
  // ****************************************************************************

  // ****************************************************************************
  // Example of computing an ISRF
  // We will compute it on thisWave and scale it by 1 MMP (~1.3G)
  // Instantiate the object
  // Again, everything is CGS!
  ISRF MyISRF(thisWave,1.0); 
  // Get the field from the object. 
  vector <float> thisISRF = MyISRF.getISRF(); 
  // ****************************************************************************

  // ****************************************************************************
  // Example of computing the Equilibrium temperature. 
  vector <float> thisTemp; 

  float Tlo,Thi;

  thisCAbs = MyGrain[0].getCAbs(0); 
  thisTemp.push_back(EqTemp(thisWave,thisISRF,thisCAbs));

  for (int sz=1;sz<nSize;sz++) { 
   
    thisCAbs = MyGrain[0].getCAbs(sz); 

    // Should make Tlo,Thi a tunable parameter somehow...
    // Need to experiment on specifying this for stochastic portion. 
    Tlo = (0.3*thisTemp[sz-1]); // < 1.0) ? 1.0 : (thisTemp[sz-1]-10.0); 
    Thi = (4.0*thisTemp[sz-1]); // > 2500.0) ? 2500.0 : (thisTemp[sz-1]+10.0);
    
    thisTemp.push_back(EqTemp(thisWave,thisISRF,thisCAbs,((Tlo<0)?1.0:Tlo),((Thi>2500.0)?2500.0:Thi)));

  }
    
  float mysize=0.0040*1.0e-4;
  int sizeid=NumUtils::index(mysize,thisSize); 

  cout << "EQUILIBRIUM TEMPERATURE AT SIZE " << mysize << "(" 
       << thisSize[sizeid] << ") IS " << thisTemp[sizeid] << endl; 
  // ****************************************************************************

//   // ****************************************************************************
//   // Example of computing the time evolution of a small grain
  
//   // Storage for heat capacity, enthalpy and the temperature they are defined on. 
//   vector <float> SHeatCapacity = MyGrain[0].getSpecificHeatCapacity(); 
//   vector <float> SEnthalpy = MyGrain[0].getSpecificEnthalpy();
//   vector <float> Temperature = MyGrain[0].getTemperature(); 
//   float Density=MyGrain[0].getDensity(); 

//   // Time and Temperature History.
//   vector <float> Time; 
//   vector <float> TH; 

//   // Get Cabs as mysize. 
//   thisCAbs = MyGrain[0].getCAbs(mysize);

//   vector <float> integrand1(thisWave.size());
//   vector <float> integrand2(thisWave.size());

//   TimeEvolution(thisWave,thisISRF,thisCAbs,Temperature,SHeatCapacity,
// 		SEnthalpy,Density,mysize,Time,TH); 

//   for (int i=0;i<Time.size();i++) cout << Time[i] << " " << TH[i] << endl; 

//   // ****************************************************************************

}
