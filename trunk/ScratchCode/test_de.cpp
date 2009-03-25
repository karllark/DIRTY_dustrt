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
		    MW,MS,"/work/Science/Dust/OpticalProperties/");

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
  string aConfigFile="/work/Science/dirty/ScratchCode/Craptastic.cfg";
  ConfigFile thisConfig(aConfigFile);
  GrainModel thisGrainModel; 
  thisGrainModel.MakeGrainModel(thisConfig,MW);
  
  cout << "Total mass (gm/H): " << thisGrainModel.getMDust() << endl; 
  cout << "Total numer (/H): " << thisGrainModel.getNumber() << endl; 
  
  int ncomp=thisGrainModel.getNComp();

  cout << "Mass, number per H: " << endl; 
  for (int i=0;i<ncomp;i++) cout << thisGrainModel.getMDust(i) << " , " << thisGrainModel.getNumber(i) << endl;  

  exit(8); 

  vector <vector<double> > thisEmission; 
  thisEmission.resize(2*ncomp+1); 
  for (int i=0;i<(2*ncomp+1);++i) thisEmission[i].resize(nWave);
  bool DoStochastic=true; 

  // a range of ISRF scalings. 
  int nisrf=101; 
  vector <float> ISRF_Scaling(nisrf); 
  float isrf_min=0.1; 
  float isrf_max=1.0e6;
  float isrf_delta=(1.0/static_cast<float>(nisrf-1))*log10(isrf_max/isrf_min); 
  ISRF_Scaling[0]=isrf_min; 
  for (int i=1;i<nisrf;++i) ISRF_Scaling[i]=ISRF_Scaling[i-1]*pow(10,isrf_delta);
  ISRF_Scaling[14]=1.0;
  ISRF_Scaling[29]=10.0; 
  ISRF_Scaling[43]=100.0;
  ISRF_Scaling[57]=1000.0;
  ISRF_Scaling[71]=10000.0; 
  ISRF_Scaling[86]=100000.0; 

  string sOutBase="SED_ISRF";
  string sOutput; 

  // Dust masses in gr/H
  float TotalDustMass=thisGrainModel.getMDust(); 
  vector <float> DustMass(ncomp); 
  for (int c=0;c<ncomp;++c) DustMass[c] = thisGrainModel.getMDust(c); 
  
  for (int is=0;is<nisrf;++is) {
    
    cout << "Working on scaling " << is << "(" << StringManip::vtos(ISRF_Scaling[is]) << ")" << endl; 
    sOutput=sOutBase+StringManip::vtos(ISRF_Scaling[is])+"_"+thisGrainModel.getModelName()+".dat"; 
    
    ISRF MyISRF(MW,ISRF_Scaling[is]); 
    thisISRF=MyISRF.getISRF(); 

    //   cout << "Calling compute emission " <<endl;
    ComputeDustEmission(thisISRF, thisGrainModel, thisEmission, DoStochastic); 

    ofstream fOutput(sOutput.c_str()); 
    
    fOutput << "# ISRF " << ISRF_Scaling[is] << " DustModel " << thisGrainModel.getModelName() << endl; 
    fOutput << "# This is the emission (luminosity) per unit dust mass" << endl; 

    // we output the emmision per unit dust mass. 
    for (int i=0;i<nWave;i++) { 
      
      fOutput << MW[i] << " " << thisISRF[i] << " " << thisEmission[0][i]/TotalDustMass << " ";

      for (int c=0;c<ncomp;++c) { 
	
	fOutput << thisEmission[2*c+1][i]/DustMass[c] << " " << thisEmission[2*c+2][i]/DustMass[c] << " "; 

      }
      fOutput << endl; 
    }
    
    fOutput.close(); 
  }


}
