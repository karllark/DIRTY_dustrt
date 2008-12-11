#include "ISRF.h"
#include "StringManip.h"
#include "Grain.h"
#include "EqTemp.h"         // For equilibrium heating
#include "TimeEvolution.h"  // Time evolution
#include "NumUtils.h"       // Numerical stuff - namespaced so can use with older defs. 
#include "mt_ran.h"         //  
#include "GrainModel.h"   
#include <algorithm>
#include <functional>
#include <numeric>

#include "DataFile.h"

extern void ComputeDustEmission (vector <float> & J, 
				 GrainModel & GrainModel,
				 vector <vector<double> > & EmmittedEnergy, 
				 bool & DoStochastic);

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

  // Just to define a wavelength grid. 
  Grain MyGrain; 
  // Instantiate a Grain object. 
  MyGrain.MakeGrain("CrossSections/PAHneu_30_mod.srt",
		    "Calorimetry/Graphitic_Calorimetry_1000.dat",
		    MW,MS,"/home/misselt/Science/Dust/OpticalProperties/");

  // Get the wavelength grid MyGrain[0] is defined on.
  MW.erase(MW.begin(),MW.end());
  MW = MyGrain.getWave();
  int nWave=MW.size(); 
  // Generate an ISRF
  ISRF MyISRF(MW,1.0); 
  // Get the field from the object. 
  vector <float> thisISRF = MyISRF.getISRF();

  // Generate the Grain model
  string aConfigFile="Craptastic.cfg";
  ConfigFile thisConfig(aConfigFile);
  GrainModel thisGrainModel; 
  thisGrainModel.MakeGrainModel(thisConfig,MW);
  
  int ncomp=thisGrainModel.getNComp();
  //   float mysize=0.0040*1.0e-4;
//   int sizeid=NumUtils::index(mysize,thisSize); 
  vector <vector<double> > thisEmission; 
  thisEmission.resize(2*ncomp+1); 
  for (int i=0;i<(2*ncomp+1);++i) thisEmission[i].resize(nWave);
  bool DoStochastic=true; 
  //   cout << "Calling compute emission " <<endl;
  ComputeDustEmission(thisISRF, thisGrainModel, thisEmission, DoStochastic); 

  for (int i=0;i<nWave;i++) { 
    //cout << thisEmission[i].size() << endl; 
    cout << MW[i] << " " << thisEmission[0][i] << " ";
    for (int c=0;c<ncomp;++c) { 
      
      cout << thisEmission[2*c+1][i] << " " << thisEmission[2*c+2][i] << " "; 
    }
    cout << endl; 
  }

//   vector <float> thissize=thisGrainModel.Size(0);
//   for (int i=0;i<thissize.size();++i) cout << i << " " << thissize[i] << endl; 


}
