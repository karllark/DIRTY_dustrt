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
  
//   float dtg=_mcf.FValue("Model","Dust to Gas Mass Ratio"); 
//   cout << dtg << endl; 
//   float amin=_mcf.FValue("Component 1","a_min"); 
//   cout << amin << endl; 


  // ****************************************************************************
  // Test the Gaussian model
  int nWave=1200; 
  vector <float> MW(nWave);  
  MW[0]=0.09;         
  MW[nWave-1]=1000.0; 
  float delta=(1.0/static_cast<float>(nWave-1))*log10(MW[nWave-1]/MW[0]); 
  for (int wv=1;wv<nWave-1;++wv) MW[wv]=MW[wv-1]*pow(10,delta); 
  // Put wave length scale on cm. 
  transform(MW.begin(),MW.end(),MW.begin(),bind2nd(multiplies<float>(),Constant::UM_CM));

  string cf="GaussianTest.cfg"; 
  ConfigFile Gconf(cf); 
  
  GrainModel* GaussianModel = new GrainModel(); 
  
  GaussianModel->MakeGrainModel(Gconf,MW); 
  
  cout << "Computed quantities from the normal GrainModel() member functions" << endl; 
  cout << "Total dust mass in gm/H: " << GaussianModel->getMDust() << endl; 
  cout << "Dust to gas mass ratio: " << GaussianModel->getMDust()/(Constant::AMU_CGS*GaussianModel->getMeanMolecularWeight()) << " (Compared to input value of : " << GaussianModel->getDustToGasMassRatio() << endl; 
  cout << "Mass fractions (component mass, fraction)" << endl; 
  for (int c=0;c<GaussianModel->getNComp();++c) cout << GaussianModel->getMDust(c) << " " << GaussianModel->getMDust(c)/GaussianModel->getMDust() << endl; 

  // Get the "extinction curve" (ie tau/H)
  vector <float> tau_gauss=GaussianModel->getTau();

  vector <float> sd;
  vector <float> sz;
  for (int i=0;i<GaussianModel->getNComp();i++) {
    sz=GaussianModel->Size(i);
    sd=GaussianModel->getSizeDistribution(i); 
    for (int j=0;j<sz.size();j++) cout << sz[j] << " " << sd[j] << endl; 
    cout << endl; 
  }

  cout << "*****************  BARE-GR-S **********************" << endl; 
  // For fun, lets get the BARE-GR-S output as well. 
  string aConfigFile="/work/Science/dirty/ScratchCode/Craptastic.cfg";
  ConfigFile thisConfig(aConfigFile);
  GrainModel thisGrainModel; 
  thisGrainModel.MakeGrainModel(thisConfig,MW);
  int nComp=thisGrainModel.getNComp(); 

  cout << "Total mass (gm/H): " << thisGrainModel.getMDust() << endl; 
  cout << "Mass per H, fraction: " << endl; 
  for (int i=0;i<nComp;i++) cout << thisGrainModel.getMDust(i) << " , " << thisGrainModel.getMDust(i)/thisGrainModel.getMDust() << endl; 
  vector <float> tau_bgs=thisGrainModel.getTau();
  cout << "*****************  BARE-GR-S **********************" << endl; 
  cout << endl; 
  cout << "*****************  RV31_BC6-WD01 **********************" << endl; 
  // How about Rv 3.1 MW from WD01? 
  string newConfigFile="/work/Science/dirty/ScratchCode/WGD01.cfg";
  ConfigFile newConfig(newConfigFile);
  GrainModel newGrainModel; 
  newGrainModel.MakeGrainModel(newConfig,MW);
  nComp=newGrainModel.getNComp(); 

  cout << "Total mass (gm/H): " << newGrainModel.getMDust() << endl; 
  cout << "Mass per H, fraction: " << endl; 
  for (int i=0;i<nComp;i++) cout << newGrainModel.getMDust(i) << " , " << newGrainModel.getMDust(i)/newGrainModel.getMDust() << endl; 
  vector <float> tau_wd01=newGrainModel.getTau();
  
  for (int wv=0;wv<nWave;wv++) cout << MW[wv] <<  " " << tau_gauss[wv] << " " << tau_wd01[wv] << endl;
 
  //vector <float> sd;
  //vector <float> sz;
  for (int i=0;i<nComp;i++) {
    sz=newGrainModel.Size(i);
    sd=newGrainModel.getSizeDistribution(i); 
    for (int j=0;j<sz.size();j++) cout << sz[j] << " " << sd[j] << endl; 
    cout << endl; 
  }
  // ****************************************************************************



  // ****************************************************************************
  // Test making a grain. 
  // These determine whether things will be interpolated. 
  // one element, -1 vectors mean take the grids defined in the cross section 
  // file.  Pass vectors tabulated where you want data as MW, MS [IN CGS, IE CM 
  // FOR BOTH!] to interpolate the cross section file onto a new grid. 
  //vector <float> MW(1); 
//   vector <float> MS(1);
//   MW[0]=-1;
//   MS[0]=-1;

//   // Just to define a wavelength grid. 
//   Grain MyGrain; 
//   // Instantiate a Grain object. 
//   MyGrain.MakeGrain("CrossSections/PAHneu_30_mod.srt",
// 		    "Calorimetry/Graphitic_Calorimetry_1000.dat",
// 		    MW,MS,"/work/Science/Dust/OpticalProperties/");

//   // Get the wavelength grid MyGrain[0] is defined on.
//   MW.erase(MW.begin(),MW.end());
//   MW = MyGrain.getWave();
//   nWave=MW.size(); 
//   // Generate an ISRF
//   ISRF MyISRF(MW,1.0); 
//   // Get the field from the object. 
//   vector <float> thisISRF = MyISRF.getISRF();
//   //for (int i=0;i<MW.size();++i) cout << MW[i] << " " << thisISRF[i] << endl; 
//   // energy density: 
//   float ed=(Constant::FPI/Constant::LIGHT)*NumUtils::integrate<float>(MW,thisISRF);  
//    cout << "Energy density in erg/cm^3: " << ed << endl; 
//    //exit(8);
//   cout << "Done in test...." << endl; 
//   // Generate the Grain model
//   string aConfigFile="/work/Science/dirty/ScratchCode/Craptastic.cfg";
//   ConfigFile thisConfig(aConfigFile);
//   GrainModel thisGrainModel; 

//   cout << thisConfig.SValue("Model Book Keeping","Model Name") << endl;
//   thisGrainModel.MakeGrainModel(thisConfig,MW);
  
//   cout << "Total mass (gm/H): " << thisGrainModel.getMDust() << endl; 
//   cout << "Total numer (/H): " << thisGrainModel.getNumber() << endl; 

//   bool UseEffective_Grain=false;
  
//   int ncomp=thisGrainModel.getNComp();

//   int ncomp_effective; 
//   if (UseEffective_Grain) ncomp_effective=1; else ncomp_effective=ncomp;

//   cout << "Mass, number per H: " << endl; 
//   for (int i=0;i<ncomp;i++) cout << thisGrainModel.getMDust(i) << " , " << thisGrainModel.getNumber(i) << endl;  
  
//   vector <vector<double> > thisEmission; 
//   thisEmission.resize(2*ncomp_effective+1); 
//   for (int i=0;i<(2*ncomp_effective+1);++i) thisEmission[i].resize(nWave);
//   bool DoStochastic=false; 
//   if (UseEffective_Grain) DoStochastic=false; 
  
//   string sOutBase="SED_ISRF";
//   string sOutput; 

//   // Dust masses in gr/H
//   float TotalDustMass=thisGrainModel.getMDust(); 
//   vector <float> DustMass(ncomp); 
//   for (int c=0;c<ncomp_effective;++c) DustMass[c] = thisGrainModel.getMDust(c); 
//   int status; 
//   thisISRF=MyISRF.getISRF(); 
//   float thise; int thisid;
//   status = ComputeDustEmission(thisISRF, thisGrainModel, thisEmission, DoStochastic,UseEffective_Grain,thise,thisid); 
//   for (int i=0;i<nWave;i++) { 
    
//     cout << MW[i] << " " << thisISRF[i] << " " << thisEmission[0][i] << " ";
    
//     for (int c=0;c<ncomp_effective;++c) { 
      
//       cout << thisEmission[2*c+1][i] << " " << thisEmission[2*c+2][i] << " "; 
      
//     }
//     cout << endl; 
//   }


}
