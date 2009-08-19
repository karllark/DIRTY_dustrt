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
#include <fstream>
#include <iostream>
#include "Constants.h"
#include "DataFile.h"
#include "DirtyFailure.h" 
#include "DirtyFlags.h"

extern int  ComputeDustEmission (vector <float> & J, 
				 GrainModel & GrainModel,
				 vector <vector<double> > & EmmittedEnergy, 
				 bool & DoStochastic, bool & UseEffective_Grain,
	 		         float & _FailureSz, int & _FailureComp);

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
  MyGrain.MakeGrain("CrossSections/PAH_28_1201_neu",
		    "Calorimetry/Graphitic_Calorimetry_1000.dat",
		    MW,MS,"/home/misselt/Dust/OpticalProperties/");

  // Get the wavelength grid MyGrain[0] is defined on.
  MW.erase(MW.begin(),MW.end());
  MW = MyGrain.getWave();
  int nWave=MW.size(); 
  // Generate an ISRF
  ISRF MyISRF(MW,1.0); 
  // Get the field from the object. 
  vector <float> thisISRF = MyISRF.getISRF();
  //for (int i=0;i<MW.size();++i) cout << MW[i] << " " << thisISRF[i] << endl; 
  //exit(8);
  cout << "Done in test...." << endl; 
  // Generate the Grain model
  string aConfigFile="/home/misselt/dirty/ScratchCode/Craptastic.cfg";
  ConfigFile thisConfig(aConfigFile);
  GrainModel thisGrainModel; 

  cout << thisConfig.SValue("Model Book Keeping","Model Name") << endl;
  thisGrainModel.MakeGrainModel(thisConfig,MW);
  
  cout << "Total mass (gm/H): " << thisGrainModel.getMDust() << endl; 
  cout << "Total numer (/H): " << thisGrainModel.getNumber() << endl; 

  bool UseEffective_Grain=false;
  
  int ncomp=thisGrainModel.getNComp();

  int ncomp_effective; 
  if (UseEffective_Grain) ncomp_effective=1; else ncomp_effective=ncomp;

  cout << "Mass, number per H: " << endl; 
  for (int i=0;i<ncomp;i++) cout << thisGrainModel.getMDust(i) << " , " << thisGrainModel.getNumber(i) << endl;  
  
  vector <vector<double> > thisEmission; 
  thisEmission.resize(2*ncomp_effective+1); 
  for (int i=0;i<(2*ncomp_effective+1);++i) thisEmission[i].resize(nWave);
  bool DoStochastic=false; 
  if (UseEffective_Grain) DoStochastic=false; 
  
  string sOutBase="SED_ISRF";
  string sOutput; 

  // Dust masses in gr/H
  float TotalDustMass=thisGrainModel.getMDust(); 
  vector <float> DustMass(ncomp); 
  for (int c=0;c<ncomp_effective;++c) DustMass[c] = thisGrainModel.getMDust(c); 
  int status; 
  thisISRF=MyISRF.getISRF(); 
  float thise; int thisid;
  status = ComputeDustEmission(thisISRF, thisGrainModel, thisEmission, DoStochastic,UseEffective_Grain,thise,thisid); 
  for (int i=0;i<nWave;i++) { 
    
    cout << MW[i] << " " << thisISRF[i] << " " << thisEmission[0][i] << " ";
    
    for (int c=0;c<ncomp_effective;++c) { 
      
      cout << thisEmission[2*c+1][i] << " " << thisEmission[2*c+2][i] << " "; 
      
    }
    cout << endl; 
  }


}
