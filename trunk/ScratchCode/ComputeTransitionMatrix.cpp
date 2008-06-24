#include "HeatUtils.h" 

#include <iostream>
#include <vector> 

using namespace std; 

vector <double> ComputeTransitionMatrix(vector <float> & waveE, vector <float> & wave, vector <float> & J, 
					vector <float> & cabs, vector <float> & Temp, vector <float> & Enthalpy, float eqT) 

{

  bool binwidth=false; 
  vector <vector<float> > TM;
  vector <double> _P; 
  int nBins=50; 
  int thisnWave; 
  int nWave=wave.size(); 
  int nWaveLast=nWave-1; // point to last element of wave - less arithmetic in loops.
  float mpe=0.0; 
  float t_abs=HeatUtils::getTauAbs(J,cabs,wave,mpe);
  float Tmean=HeatUtils::getTmean(Temp,Enthalpy,mpe); 
  vector <float>::iterator it,it0,it1,it2,it3,it4,ib,ie; 
  
  // Some initial algorithm to compute search size in constant T
  // Start simple with [0.1,4]*Equilirium Temperature - will probably have to 
  // get more sophisticated.  
  float TMax = 4.0*eqT; 
  TMax = (TMax > Temp[Temp.size()-1])?Temp[Temp.size()-1]:TMax; 
  float TMin = 0.1*eqT;
  //cout << "TMIN, TMEAN, TEQUI, TMAX: " << TMin << " " <<  Tmean << " " << eqT << " " << TMax << endl; 

  // Generate regular T grid with nBins points and interpolate Enthalpy 
  // onto that grid. 
  vector <float> tgrid(nBins+1),enthgrid(nBins+1),enth(nBins),denth(nBins); 
  float DT=(TMax-TMin)/(float)(nBins);
  it=tgrid.begin(); 
  *it = TMin; 
  ib=tgrid.begin()+1; 
  ie=tgrid.end(); 
  for (it=ib;it!=ie;++it) *it = *(it-1)+DT;  
  enthgrid=NumUtils::interpol(Enthalpy,Temp,tgrid,4,1);
  // Compute the enthalpies at grid center
  ib=enthgrid.begin()+1; 
  ie=enthgrid.end(); 
  it1=enth.begin();
  it2=denth.begin(); 
  for (it=ib;it!=ie;++it,++it1,++it2) {
    *it1=(*it+*(it-1))/2.0; 
    *it2=(*it-*(it-1)); 
  }
  
  //for (int i=0;i<nWave;i++) cout << i << " " << wave[i] << " " << waveE[i] << endl; 
  float arg1,arg2,w1,w2,w3,w4,wc,val;
  int idx1,idx2,idx3,idx4,idxc;
  vector <float> thisWave,thisCabs,thisJ; 
  bool lextra=false;
  bool hextra=false;

  // Size the transition matrix.
  TM.resize(nBins); 
  for (int l=0;l<nBins;l++) TM[l].resize(nBins); 

  _P.resize(nBins); 

  if (binwidth) {  // This method accounts for finite bin widths. 

    for (int l=0;l<nBins-1;l++) {
      for (int u=l+1;u<nBins-1;u++) {
	
	// Compute energy integration bounds
	lextra=false;
	hextra=false;
	w1 = enthgrid[u]-enthgrid[l+1];
	arg1 = enthgrid[u]-enthgrid[l];
	arg2 = enthgrid[u+1]-enthgrid[l+1];
	w2 = min(arg1,arg2); 
	w3 = max(arg1,arg2);
	w4 = enthgrid[u+1]-enthgrid[l];
	
	if (w1 > waveE[0] || w4 < waveE[nWave-1]) { 
	  // No photons in our input field can cause transitions between these bins. 
	  // Photon energies between [w1,w4] are necessary
	  TM[u][l] = 0.0;
	} else { 
	  // generate the wavelength grid.  We'll extend beyond nominal bounds
	  // to avoid inserts and push_backs. 
	  idx1=NumUtils::rindex(w1,waveE);
	  if (idx1 < nWave) { ++idx1; hextra=true; }
	  idx4=NumUtils::rindex(w4,waveE);
	  if (idx4 > 0) {  --idx4; lextra=true; } 
	  //cout << "idx4: " << idx4 << endl;
	  
	  //cout << l << " " << u << " " << wave.size() << " " << idx1 << " " << idx4 <<endl; 
	  // Copy necessary sections of input vectors.
	  it = wave.begin()+idx4; 
	  it0 = wave.begin()+idx1; 
	  it1 = cabs.begin()+idx4;
	  it2 = J.begin()+idx4; 
	  it3 = cabs.begin()+idx1; 
	  it4 = J.begin()+idx1; 
	  
	  thisWave.assign(it,it0);
	  thisnWave=thisWave.size(); 
	  thisCabs.assign(it1,it3);
	  thisJ.assign(it2,it4); 
	  
	  if (lextra) { // not at the shortest wavelength end; interpolate
	    //cout << "extrapolating: " << endl; 
	    thisWave[0] = Constant::PLANCKLIGHT/w4;
	    thisCabs[0] = NumUtils::line(*(it),*(it+1),*(it1),*(it1+1),thisWave[0]);
	    thisJ[0] = NumUtils::line(*(it),*(it+1),*(it2),*(it2+1),thisWave[0]);
	  }
	  if (hextra) { // not at the longest wavelength edn; interpolate
	    thisWave[thisnWave-1] = Constant::PLANCKLIGHT/w1;
	    thisCabs[thisnWave-1] = NumUtils::line(*it0,*(it0+1),*it3,*(it3+1),thisWave[thisnWave-1]);
	    thisJ[thisnWave-1] = NumUtils::line(*it0,*(it0+1),*it4,*(it4+1),thisWave[thisnWave-1]);
	  }
	  
	}
      } 
      
      // Do u=M 
      w1 = enthgrid[nBins-1]-enthgrid[l+1];
      wc = enthgrid[nBins-1]-enthgrid[l];
      
      if (w1 > waveE[0] || wc < waveE[nWave-1]) { 
	//cout << "does tm get a value?" << nBins-1 << " " << l << endl; 
	TM[nBins-1][l] = 0.0; 
	//cout << "4" << endl; 
	//cout << "last bin: " <<  l << " " << nBins-1 << " No photons will do this transition" << endl;
      } else {
	//cout << "4" << endl; 
	// generate the wavelength grird
	idx1=NumUtils::rindex(w1,waveE);
	if (idx1 < nWave) { ++idx1; hextra=true; }
	idxc=NumUtils::rindex(wc,waveE);
	if (idxc > 0) { --idxc; lextra=true; } 
	it = wave.begin()+idxc; 
	it0 = wave.begin()+idx1; 
	it1 = cabs.begin()+idxc;
	it2 = J.begin()+idxc; 
	it3 = cabs.begin()+idx1; 
	it4 = J.begin()+idx1; 
	
	thisWave.assign(it,it0);
	thisnWave=thisWave.size(); 
	thisCabs.assign(it1,it3);
	thisJ.assign(it2,it4); 
	
	if (lextra) { // not at the shortest wavelength end; interpolate
	  //cout << "extrapolating: " << endl; 
	  thisWave[0] = Constant::PLANCKLIGHT/wc;
	  thisCabs[0] = NumUtils::line(*(it),*(it+1),*(it1),*(it1+1),thisWave[0]);
	  thisJ[0] = NumUtils::line(*(it),*(it+1),*(it2),*(it2+1),thisWave[0]);
	  
	  // Need to add in the last integral to infinite energy...  
	}
	if (hextra) { // not at the longest wavelength edn; interpolate
	  thisWave[thisnWave-1] = Constant::PLANCKLIGHT/w1;
	  thisCabs[thisnWave-1] = NumUtils::line(*it0,*(it0+1),*it3,*(it3+1),thisWave[thisnWave-1]);
	  thisJ[thisnWave-1] = NumUtils::line(*it0,*(it0+1),*it4,*(it4+1),thisWave[thisnWave-1]);
	}
	//cout << "last bin: " << l << " " << nBins-1 << " " << waveE[idx1] << "(" << wave[idx1] << ") " << waveE[idxc] << "(" << wave[idxc]<<")" << endl; 
	//cout << "index: " << wave[idx1] << " " << wave[idx3] << endl; 
      }
      //wc = Constant::PLANCKLIGHT/(enthgrid[nBins-1]-enthgrid[l]); 
      //cout << l << " " << w1 << " " << wc << endl; 
    }

  } else {  // This one... doesn't

    for (int l=0;l<nBins-1;l++) {
      for (int u=l+1;u<nBins-1;u++) {

      }
    }
  }

  // Compute _P
  
  //
  return _P; 
}
