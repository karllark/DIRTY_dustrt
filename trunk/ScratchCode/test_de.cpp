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

  // Test dust heating for RT validation
  int nWave; 
  string wavefile="wavegrid_codeval.dat"; 
  vector <double> dMW; 
  DataFile(wavefile, dMW);
  nWave = dMW.size(); 
  transform(dMW.begin(),dMW.end(),dMW.begin(),bind2nd(multiplies<double>(),Constant::UM_CM)); 

  vector <float> MW(dMW.begin(), dMW.end());

  string cf="EffectiveGrainTest.cfg"; 
  ConfigFile GConf(cf); 
  GrainModel EffectiveGrain; 
  EffectiveGrain.MakeGrainModel(GConf, MW); 

  //cout << EffectiveGrain.getNormalization() << endl; 
  //vector <float> avgCAbs=EffectiveGrain.getCAbsEff(); 
  
  int nRF = 5; 
  float array[]={0.01,0.1,1.0,10.0,100.0}; 
  vector <float> RFScale; 
  RFScale.assign(array,array+5);  

  bool DoStochastic=false; 
  bool UseEffectiveGrain=true; 
  
  int status,fcmp;
  float fsz; 
  status=0; fcmp=-1; fsz=-1.0; 
  
  int nEmissionComp=3; 
  vector <vector<double> > EmittedEnergy(nEmissionComp); 
  for (int i=0;i<nEmissionComp;++i) EmittedEnergy[i].resize(nWave,0.0);

  vector <float> thisISRF; 
  
  float thisNorm = Constant::UM_CM/(EffectiveGrain.getNormalization()); 

  for (int r=0;r<nRF;++r) { 
    cout << "*** ISRF SCALE: " << RFScale[r] << endl; 
    ISRF* myISRF = new ISRF(MW,RFScale[r]); 
    thisISRF = myISRF->getISRF(); 
    status=ComputeDustEmission(thisISRF, EffectiveGrain, EmittedEnergy, 
			       DoStochastic, UseEffectiveGrain, fsz, fcmp); 
    for (int wv=0;wv<nWave;++wv) cout << MW[wv]*Constant::CM_UM << " " << EmittedEnergy[0][wv]*thisNorm << " " << EmittedEnergy[1][wv]*thisNorm << " " << EmittedEnergy[2][wv]*thisNorm << endl; 
    cout << " " << endl; 
  }

  // ****************************************************************************
  // Test the Gaussian model
//   int nWave=1200; 
//   vector <float> MW(nWave);  
//   MW[0]=0.09;         
//   MW[nWave-1]=1000.0; 
//   float delta=(1.0/static_cast<float>(nWave-1))*log10(MW[nWave-1]/MW[0]); 
//   for (int wv=1;wv<nWave-1;++wv) MW[wv]=MW[wv-1]*pow(10,delta); 
//   // Put wave length scale on cm. 
//   transform(MW.begin(),MW.end(),MW.begin(),bind2nd(multiplies<float>(),Constant::UM_CM));

//   string cf="GaussianTest.cfg"; 
//   ConfigFile Gconf(cf); 
  
//   GrainModel* GaussianModel = new GrainModel(); 
  
//   GaussianModel->MakeGrainModel(Gconf,MW); 
  
//   cout << "Computed quantities from the normal GrainModel() member functions" << endl; 
//   cout << "Total dust mass in gm/H: " << GaussianModel->getMDust() << endl; 
//   cout << "Dust to gas mass ratio: " << GaussianModel->getMDust()/(Constant::AMU_CGS*GaussianModel->getMeanMolecularWeight()) << " (Compared to input value of : " << GaussianModel->getDustToGasMassRatio() << endl; 
//   cout << "Mass fractions (component mass, fraction)" << endl; 
//   for (int c=0;c<GaussianModel->getNComp();++c) cout << GaussianModel->getMDust(c) << " " << GaussianModel->getMDust(c)/GaussianModel->getMDust() << endl; 

//   // Get the "extinction curve" (ie tau/H)
//   vector <float> tau_gauss=GaussianModel->getTau();

//   vector <float> sd;
//   vector <float> sz;
//   for (int i=0;i<GaussianModel->getNComp();i++) {
//     sz=GaussianModel->Size(i);
//     sd=GaussianModel->getSizeDistribution(i); 
//     for (int j=0;j<sz.size();j++) cout << sz[j] << " " << sd[j] << endl; 
//     cout << endl; 
//   }

//   cout << "*****************  BARE-GR-S **********************" << endl; 
//   // For fun, lets get the BARE-GR-S output as well. 
//   string aConfigFile="/work/Science/dirty/ScratchCode/Craptastic.cfg";
//   ConfigFile thisConfig(aConfigFile);
//   GrainModel thisGrainModel; 
//   thisGrainModel.MakeGrainModel(thisConfig,MW);
//   int nComp=thisGrainModel.getNComp(); 

//   cout << "Total mass (gm/H): " << thisGrainModel.getMDust() << endl; 
//   cout << "Mass per H, fraction: " << endl; 
//   for (int i=0;i<nComp;i++) cout << thisGrainModel.getMDust(i) << " , " << thisGrainModel.getMDust(i)/thisGrainModel.getMDust() << endl; 
//   vector <float> tau_bgs=thisGrainModel.getTau();
//   cout << "*****************  BARE-GR-S **********************" << endl; 
//   cout << endl; 
//   cout << "*****************  RV31_BC6-WD01 **********************" << endl; 
//   // How about Rv 3.1 MW from WD01? 
//   string newConfigFile="/work/Science/dirty/ScratchCode/WGD01.cfg";
//   ConfigFile newConfig(newConfigFile);
//   GrainModel newGrainModel; 
//   newGrainModel.MakeGrainModel(newConfig,MW);
//   nComp=newGrainModel.getNComp(); 

//   cout << "Total mass (gm/H): " << newGrainModel.getMDust() << endl; 
//   cout << "Mass per H, fraction: " << endl; 
//   for (int i=0;i<nComp;i++) cout << newGrainModel.getMDust(i) << " , " << newGrainModel.getMDust(i)/newGrainModel.getMDust() << endl; 
//   vector <float> tau_wd01=newGrainModel.getTau();
  
//   for (int wv=0;wv<nWave;wv++) cout << MW[wv] <<  " " << tau_gauss[wv] << " " << tau_wd01[wv] << endl;
 
//   //vector <float> sd;
//   //vector <float> sz;
//   for (int i=0;i<nComp;i++) {
//     sz=newGrainModel.Size(i);
//     sd=newGrainModel.getSizeDistribution(i); 
//     for (int j=0;j<sz.size();j++) cout << sz[j] << " " << sd[j] << endl; 
//     cout << endl; 
//   }
  // ****************************************************************************




}
