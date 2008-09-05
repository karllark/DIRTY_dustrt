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

extern void ComputeDustEmission (vector <float> & J, 
				 GrainModel & GrainModel,
				 vector <vector<double> > & EmmittedEnergy, 
				 bool DoStochastic);

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


//   vector <float> v1,v2; 
//   for (int i=1;i<=10;++i) v1.push_back((float)i);
//   for (int i=1;i<=20;++i) v2.push_back((float)i); 
//   cout << "Size before assign: " <<  v1.size() << " " << v2.size() << endl; 
//   v2.assign(v1.begin(),v1.end()); 
//   cout << "Size after assign: " << v1.size() << " " << v2.size() << endl; 
//   cout << "****" << endl; 
//   for (int i=0;i<v2.size();++i) cout << i+1 << " " << v1[i] << " " << v2[i] << endl;
//   cout << "****" << endl; 
  
//   vector <bool> vb; 
//   for (int i=0;i<10;++i) vb.push_back(false); 
//   for (int i=0;i<10;++i) cout << i << " " << vb[i] << endl; 
//   transform(vb.begin()+3,vb.end(),vb.begin()+3,bind2nd(multiplies<bool>(),false)); 
//   for (int i=0;i<10;++i) cout << i << " " << vb[i] << endl; 
  //exit(8);

//   vector <float> atest,ainvtest; 
//   for (int i=1;i<=10;++i) atest.push_back((float)i);
//   int idx0=NumUtils::index((float)-1.0,atest); 
//   int idx1=NumUtils::index((float)3.4,atest);
//   int idx2=NumUtils::index((float)12.0,atest);
//   for (int i=0;i<atest.size();i++) cout << i << " " << atest[i] << endl; 
//   cout << endl; 
//   cout << idx0 << " " << idx1 << " " << idx2 << endl << endl; 
//   //exit(8); 
//   ainvtest.assign(atest.begin()+idx1,atest.begin()+idx2); 
//   cout << ainvtest.size() << endl; 
//   for (int i=0;i<ainvtest.size();i++) cout << i << " " << ainvtest[i] << endl; 
//   cout << endl; 
//   idx1=NumUtils::index((float)1.5,atest);
//   idx2=NumUtils::index((float)9.5,atest);
//   cout << idx1 << " " << idx2 << endl << endl; 
//   ainvtest.assign(atest.begin()+idx1,atest.begin()+idx2); 
//   cout << ainvtest.size() << endl; 
//   for (int i=0;i<ainvtest.size();i++) cout << i << " " << ainvtest[i] << endl; 
//   cout << endl; 
//   idx1=NumUtils::index((float)4.8,atest);
//   idx2=NumUtils::index((float)6.7,atest);
//   cout << idx1 << " " << idx2 << endl << endl; 
//   ainvtest.assign(atest.begin()+idx1,atest.begin()+idx2); 
//   cout << ainvtest.size() << endl; 
//   for (int i=0;i<ainvtest.size();i++) cout << i << " " << ainvtest[i] << endl; 
//   cout << endl; 

//   for (int i=0;i<atest.size();++i) atest[i] *= 2.5; 
//   for (int i=0;i<atest.size();++i) cout << i <<" " << atest[i] << endl; 
//   idx1 = NumUtils::index((float)5.6,atest);
//   cout << idx1 << endl; 
//   vector <float>::iterator it; 
//   it=atest.begin()+idx1; 
//   cout << idx1 << " " << *it << " " << *(it-1) << endl; 

//   vector <float>::iterator it1,it2,it3,it4; 
//   it1=atest.begin(); 
//   it2=atest.end(); 
//   for (it3=it1;it3!=it2;++it3) { 
//     ainvtest.push_back(1.0/(*it3));
//   }
//   it3=ainvtest.begin();
//   for (it4=it1;it4!=it2;++it4,++it3) cout << *it4 << " " << *it3 << endl; 
  //exit(8);

  // Number of components - only going to fill in one here, but this illustrates
  // the vector nature of the Grain object.
//   int ncomp=2; 
//   // Declare a vector of grain object. 
  // Just to define a wavelength grid. 
  Grain MyGrain; 
//   // Instantiate a Grain object and stick it into the first element of the vector. 
  MyGrain.MakeGrain("CrossSections/AmC_121_large_longwave",
		    "Calorimetry/Graphitic_Calorimetry_1000.dat",
		    MW,MS,"/work/Science/Dust/OpticalProperties/");
  // Get the wavelenght grid MyGrain[0] is defined on. 
  vector <float> thisWave = MyGrain.getWave(); 
  int nWave = MyGrain.getNWave(); 
//   // Get the size grid MyGrain[0] is defined on. 
//   vector <float> thisSize = MyGrain[0].getSize();
//   int nSize = MyGrain[0].getNSize();
//   // Get an absorption cross section for this grain - 0 is for the first size
//   // in the grid.  
//   vector <float> thisCAbs(nWave);
//   thisCAbs = MyGrain[0].getCAbs(0);
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

//   // Hard coded GrainModel configuration. 
//   string aConfigFile="Craptastic.cfg";
//   ConfigFile thisConfig(aConfigFile);
//   // Instiate a GrainModel - requires a configuration file which will tell the 
//   // class where to get the model info and the wavelength grid we want things 
//   // defined on. 
//   // All the work goes on behind this line. 

//   // Get the effective (size integrated, summed) grain quantities. 
//   vector <float> thiseffcabs=thisGrainModel.getCAbsEff(); 
//   vector <float> thiseffcsca=thisGrainModel.getCScaEff();
//   vector <float> thiseffphfunc=thisGrainModel.getphFuncEff();
//   vector <float> thisEnthalpy,thisSpecificEnthalpy,thisCalTemp; 
//   ncomp = thisGrainModel.getNComp();
//   float trala;
//   for (int c=0;c<1;c++) { 
//     thisSize=thisGrainModel.Size(c); 
//     thisSpecificEnthalpy=thisGrainModel.SpecificEnthalpy(c);
//     thisCalTemp=thisGrainModel.CalTemp(c); 
//     trala = thisGrainModel.Density(c); 
//     for (int sz=0;sz<1;++sz) { 
//       thisEnthalpy=thisGrainModel.Enthalpy(sz,c); 
//       cout << sz << " " << thisSize[sz] << " " << trala << endl; 
//       for (int t=0;t<thisCalTemp.size();++t) { 
// 	cout << thisCalTemp[t] << " " << thisEnthalpy[t] << "  " << Constant::VFAC*trala*thisSize[sz]*thisSize[sz]*thisSize[sz]*thisSpecificEnthalpy[t] << " " << thisSpecificEnthalpy[t] << endl;
//       }
//     }
//   }
//   exit(8);
  // Define an albedo vector and compute the albedo 
  // Albedo maybe something we should store in the GrainModel object...
//   vector <float> albedo; 
//   for (int i=0;i<nWave;i++) { 
//     albedo.push_back(thiseffcsca[i]/(thiseffcsca[i]+thiseffcabs[i])); 
//     //cout << thisWave[i] << " " << albedo[i] << " " << thiseffphfunc[i] << endl; 
//   }
//   // ****************************************************************************

  // ****************************************************************************
  // Example of computing an ISRF
  // We will compute it on thisWave and scale it by 1 MMP (~1.3G)
  // Instantiate the object
  // Again, everything is CGS!
  ISRF MyISRF(thisWave,1.0); 
  // Get the field from the object. 
  vector <float> thisISRF = MyISRF.getISRF(); 

  GrainModel thisGrainModel;
  string aConfigFile="Craptastic.cfg";
  ConfigFile thisConfig(aConfigFile);
  thisGrainModel.MakeGrainModel(thisConfig,thisWave); 
  int ncomp = thisGrainModel.getNComp();

//   for (int i=0;i<thisISRF.size();++i) cout << thisWave[i]*1.0e4 << " " << thisISRF[i] << endl; 
// 	exit(8);
  // ****************************************************************************


  // ****************************************************************************
  // Example of computing the Equilibrium temperature. 
//   vector <float> thisTemp; 

//   float Tlo,Thi;
  
//   for (int c=0;c<3;c++) { 
    
//     thisSize = thisGrainModel.Size(c); 
//     nSize = thisSize.size(); 
//     thisCAbs = thisGrainModel.CAbs(c,0); 
//     thisCAbs = MyGrain[0].getCAbs(0);
//     thisTemp.resize(nSize);
//     thisTemp[0] = EqTemp(thisWave,thisISRF,thisCAbs);
    
//     for (int sz=1;sz<nSize;sz++) { 
      
//       thisCAbs = thisGrainModel.CAbs(c,sz); 
//       //thisCAbs = MyGrain[0].getCAbs(sz); 
//       // Should make Tlo,Thi a tunable parameter somehow...
//       // Need to experiment on specifying this for stochastic portion. 
//       Tlo = (0.3*thisTemp[sz-1]); // < 1.0) ? 1.0 : (thisTemp[sz-1]-10.0); 
//       Thi = (4.0*thisTemp[sz-1]); // > 2500.0) ? 2500.0 : (thisTemp[sz-1]+10.0);
      
//       thisTemp[sz] = EqTemp(thisWave,thisISRF,thisCAbs,((Tlo<0)?1.0:Tlo),((Thi>2500.0)?2500.0:Thi));
//       //      thisTemp.push_back(EqTemp(thisWave,thisISRF,thisCAbs,((Tlo<0)?1.0:Tlo),((Thi>2500.0)?2500.0:Thi)));
      
//     }
//   }

//   float mysize=0.0040*1.0e-4;
//   int sizeid=NumUtils::index(mysize,thisSize); 
  vector <vector<double> > thisEmission; 
  thisEmission.resize(2*ncomp+1); 
  for (int i=0;i<(2*ncomp+1);++i) thisEmission[i].resize(nWave);
  bool DoStochastic=true; 
  cout << "Calling compute emission " <<endl;
  ComputeDustEmission(thisISRF, thisGrainModel, thisEmission, DoStochastic); 

  for (int i=0;i<nWave;i++) { 
    //cout << thisEmission[i].size() << endl; 
    cout << thisWave[i] << " " << thisEmission[0][i] << " ";
    for (int c=0;c<ncomp;++c) { 
      
      cout << thisEmission[2*c+1][i] << " " << thisEmission[2*c+2][i] << " "; 
    }
    cout << endl; 
  }

  thisEmission.resize(2*ncomp+1);
  for (int i=0;i<(2*ncomp+1);++i) thisEmission[i].resize(nWave);
  cout << "Calling compute emission " <<endl;
  ComputeDustEmission(thisISRF, thisGrainModel, thisEmission, DoStochastic);

  for (int i=0;i<nWave;i++) {
   cout << thisWave[i] << " " << thisEmission[0][i] << " ";
    for (int c=0;c<ncomp;++c) {

      cout << thisEmission[2*c+1][i] << " " << thisEmission[2*c+2][i] << " ";
    }
    cout << endl;
  }

//   cout << "EQUILIBRIUM TEMPERATURE AT SIZE " << mysize << "(" 
//        << thisSize[sizeid] << ") IS " << thisTemp[sizeid] << endl; 

//   for (int sz=0;sz<nSize;sz++) 
//     cout << thisSize[sz] << " " << thisTemp[sz] << endl; 
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
