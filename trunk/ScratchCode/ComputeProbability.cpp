#include "HeatUtils.h" 

#include <iostream>
#include <vector> 

using namespace std; 

vector <double> ComputeProbability(vector <float> & waveE, vector <float> & wave, vector <float> & J, 
				   vector <float> & cabs, vector <float> & Temp, vector <float> & Enthalpy, float eqT) 

{

  bool binwidth=false; 
  vector <vector<double> > TM;
  vector <double> _P; 
  int nBins=500; 
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
  vector <float> tgrid(nBins+1),_temp(nBins),enthgrid(nBins+1),enth(nBins),denth(nBins); 
  float DT=(TMax-TMin)/(float)(nBins);
  it=tgrid.begin(); 
  it1=_temp.begin(); 
  *it = TMin; 
  ib=tgrid.begin()+1; 
  ie=tgrid.end(); 
  for (it=ib;it!=ie;++it,++it1) { *it = *(it-1)+DT; *it1=(*it+*(it-1))/2.0; }
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
  for (int l=0;l<nBins;l++) TM[l].resize(nBins,0.0);
//   for (int l=0;l<nBins;++l) 
//     for (int u=0;u<nBins;++u) 
//       TM[l][u] = 0.0;

  vector <float> integrand; 
  _P.resize(nBins); 

  if (binwidth) {  // This method accounts for finite bin widths. 

//     for (int l=0;l<nBins-1;l++) {
//       for (int u=l+1;u<nBins-1;u++) {
	
// 	// Compute energy integration bounds
// 	lextra=false;
// 	hextra=false;
// 	w1 = enthgrid[u]-enthgrid[l+1];
// 	arg1 = enthgrid[u]-enthgrid[l];
// 	arg2 = enthgrid[u+1]-enthgrid[l+1];
// 	w2 = min(arg1,arg2); 
// 	w3 = max(arg1,arg2);
// 	w4 = enthgrid[u+1]-enthgrid[l];
	
// 	if (w1 > waveE[0] || w4 < waveE[nWave-1]) { 
// 	  // No photons in our input field can cause transitions between these bins. 
// 	  // Photon energies between [w1,w4] are necessary
// 	  TM[u][l] = 0.0;
// 	} else { 
// 	  // generate the wavelength grid.  We'll extend beyond nominal bounds
// 	  // to avoid inserts and push_backs. 
// 	  idx1=NumUtils::rindex(w1,waveE);
// 	  if (idx1 < nWave) { ++idx1; hextra=true; }
// 	  idx4=NumUtils::rindex(w4,waveE);
// 	  if (idx4 > 0) {  --idx4; lextra=true; } 

// 	  // Copy necessary sections of input vectors.
// 	  it = wave.begin()+idx4; 
// 	  it0 = wave.begin()+idx1; 
// 	  it1 = cabs.begin()+idx4;
// 	  it2 = J.begin()+idx4; 
// 	  it3 = cabs.begin()+idx1; 
// 	  it4 = J.begin()+idx1; 
	  
// 	  thisWave.assign(it,it0);
// 	  thisnWave=thisWave.size(); 
// 	  thisCabs.assign(it1,it3);
// 	  thisJ.assign(it2,it4); 
	  
// 	  if (lextra) { // not at the shortest wavelength end; interpolate
// 	    //cout << "extrapolating: " << endl; 
// 	    thisWave[0] = Constant::PLANCKLIGHT/w4;
// 	    thisCabs[0] = NumUtils::line(*(it),*(it+1),*(it1),*(it1+1),thisWave[0]);
// 	    thisJ[0] = NumUtils::line(*(it),*(it+1),*(it2),*(it2+1),thisWave[0]);
// 	  }
// 	  if (hextra) { // not at the longest wavelength edn; interpolate
// 	    thisWave[thisnWave-1] = Constant::PLANCKLIGHT/w1;
// 	    thisCabs[thisnWave-1] = NumUtils::line(*it0,*(it0+1),*it3,*(it3+1),thisWave[thisnWave-1]);
// 	    thisJ[thisnWave-1] = NumUtils::line(*it0,*(it0+1),*it4,*(it4+1),thisWave[thisnWave-1]);
// 	  }
	  
// 	}
//       } 
      
//       // Do u=M 
//       w1 = enthgrid[nBins-1]-enthgrid[l+1];
//       wc = enthgrid[nBins-1]-enthgrid[l];
      
//       if (w1 > waveE[0] || wc < waveE[nWave-1]) { 
// 	//cout << "does tm get a value?" << nBins-1 << " " << l << endl; 
// 	TM[nBins-1][l] = 0.0; 
// 	//cout << "4" << endl; 
// 	//cout << "last bin: " <<  l << " " << nBins-1 << " No photons will do this transition" << endl;
//       } else {
// 	//cout << "4" << endl; 
// 	// generate the wavelength grird
// 	idx1=NumUtils::rindex(w1,waveE);
// 	if (idx1 < nWave) { ++idx1; hextra=true; }
// 	idxc=NumUtils::rindex(wc,waveE);
// 	if (idxc > 0) { --idxc; lextra=true; } 
// 	it = wave.begin()+idxc; 
// 	it0 = wave.begin()+idx1; 
// 	it1 = cabs.begin()+idxc;
// 	it2 = J.begin()+idxc; 
// 	it3 = cabs.begin()+idx1; 
// 	it4 = J.begin()+idx1; 
	
// 	thisWave.assign(it,it0);
// 	thisnWave=thisWave.size(); 
// 	thisCabs.assign(it1,it3);
// 	thisJ.assign(it2,it4); 
	
// 	if (lextra) { // not at the shortest wavelength end; interpolate
// 	  //cout << "extrapolating: " << endl; 
// 	  thisWave[0] = Constant::PLANCKLIGHT/wc;
// 	  thisCabs[0] = NumUtils::line(*(it),*(it+1),*(it1),*(it1+1),thisWave[0]);
// 	  thisJ[0] = NumUtils::line(*(it),*(it+1),*(it2),*(it2+1),thisWave[0]);
	  
// 	  // Need to add in the last integral to infinite energy...  
// 	}
// 	if (hextra) { // not at the longest wavelength edn; interpolate
// 	  thisWave[thisnWave-1] = Constant::PLANCKLIGHT/w1;
// 	  thisCabs[thisnWave-1] = NumUtils::line(*it0,*(it0+1),*it3,*(it3+1),thisWave[thisnWave-1]);
// 	  thisJ[thisnWave-1] = NumUtils::line(*it0,*(it0+1),*it4,*(it4+1),thisWave[thisnWave-1]);
// 	}
// 	//cout << "last bin: " << l << " " << nBins-1 << " " << waveE[idx1] << "(" << wave[idx1] << ") " << waveE[idxc] << "(" << wave[idxc]<<")" << endl; 
// 	//cout << "index: " << wave[idx1] << " " << wave[idx3] << endl; 
//       }
//       //wc = Constant::PLANCKLIGHT/(enthgrid[nBins-1]-enthgrid[l]); 
//       //cout << l << " " << w1 << " " << wc << endl; 
//     }

  } else {  // This one... doesn't

    // Accumulate diagonals concurrently
    float _wT;   // wavelength that will produce the transition. 
    float _cabs; // Absorption cross section at _wT
    float _J;    // Radition field at _wT

//     cout << endl << "*** should be zeros ***" << endl; 
//     for (int l=0;l<nBins;l++) {
//       for (int u=0;u<nBins;u++) { 
// 	cout << TM[l][u] << " ";
// 	TM[l][u] = 0.0;  // reinitialize for next test.
//       }
//       cout << endl; 
//     }

    // Initialization only necessary for diagonal elements as they are -= 
    // whereas remaining elements are =
    for (int i=0;i<nBins;++i) {
      if (i>0) { // Generate cooling elements
	_wT = Constant::PLANCKLIGHT/enth[i];
	idx1 = NumUtils::index(_wT,wave); 
	if (idx1==nWave) { // shortest wavelength is beyond upper defined bound.
	  TM[i][i-1]=0.0; 
	} else { 
	  it0=wave.begin(); it1=wave.end();
	  it2=cabs.begin(); it3=cabs.end(); 
	  if (idx1==0) {    // shortest wavelength falls below lower defined bound; use full grid
	    thisWave.assign(it0,it1); 
	    thisCabs.assign(it2,it3);
	  } else {         // shortest wavelength falls somewhere in the middle
	    it0 += (idx1-1); it2 += (idx1-1);
	    thisWave.assign(it0,it1); 
	    thisCabs.assign(it2,it3);
	    thisWave[0]=_wT; 
	    thisCabs[0]=NumUtils::line(*it0,*(it0+1),*it2,*(it2+1),_wT);
	  }
	  integrand=NumUtils::bbodyCGS(thisWave,_temp[i]);
	  transform(integrand.begin(),integrand.end(),thisCabs.begin(),integrand.begin(),multiplies<float>()); 
	  TM[i][i-1] = Constant::PLANCKLIGHT/(enth[i]-enth[i-1])*NumUtils::integrate(thisWave,integrand); 
	  TM[i-1][i-1] -= TM[i][i-1]; // accumulate diagonal
	}
      }
      for (int f=i+1;f<nBins;++f) {  // Generate heating from i to f=[i+1:nBins-1]
	_wT = Constant::PLANCKLIGHT/(enth[f]-enth[i]);  // Wavelength of photon that can i->f transition.
	idx1 = NumUtils::index(_wT,wave);               // Where that photon falls in our defined grid:
	if ( (idx1==0) || (idx1==nWave) ) {             // no such photon in our defined grid. 
	  TM[i][f] = 0.0; 
	} else {             
	  // Interpolate cabs and J to _wT
	  it0=wave.begin()+idx1; it1=cabs.begin()+idx1; it2=J.begin()+idx1;
	  _cabs=NumUtils::line(*(it0-1),*it0,*(it1-1),*it1,_wT); 
	  _J = NumUtils::line(*(it0-1),*it0,*(it2-1),*it2,_wT);
	  TM[i][f] = Constant::IPLANCKLIGHT*_cabs*_J*denth[f]*pow(_wT,3); 
	  // If photons with energies capable of exciting beyond max bin exist.
 	  // put them into the upper bin.
	  if (idx1 >= 2 && f==(nBins-1)) {
	    thisWave.assign(wave.begin(),it0); 
	    thisnWave = thisWave.size()-1;
	    thisWave[thisnWave] = _wT; 
	    thisCabs.assign(cabs.begin(),it1);
	    thisCabs[thisnWave]=NumUtils::line(*(it0-1),*it0,*(it1-1),*it1,_wT); 
	    thisJ.assign(J.begin(),it2); 
	    thisJ[thisnWave]=NumUtils::line(*(it0-1),*it0,*(it2-1),*it2,_wT);
	    integrand.resize(thisnWave+1); 
	    ib=integrand.begin();
	    ie=integrand.end();
	    it1=thisCabs.begin();
	    it2=thisJ.begin(); 
	    for (it=ib;it!=ie;++it,++it1,++it2) *it = (*it1)*(*it2);
	    TM[i][f] += (_wT*NumUtils::integrate(thisWave,integrand)); 
	  }
	  TM[f][f] -= TM[i][f]; 
	}
      }
    }
    
//     cout << endl << "*** Accumulating diagonal ***" << endl; 
//     for (int f=0;f<nBins;f++) {
//       for (int i=0;i<nBins;i++) { 
// 	cout << TM[i][f] << " ";
// 	TM[i][f] = 0.0;  // reinitialize for next test.
//       }
//       cout << endl; 
//     }

  }
//     for (int u=1;u<nBins;++u) { // Initial levels. 
//       TM[u][u] = 0.0;
//       _wT = Constant::PLANCKLIGHT/enth[u];
//       idx1 = NumUtils::index(_wT,wave); 
//       // Generate the cooling elements.
//       if (idx1==nWave) { // shortest wavelength is beyond upper defined bound.
// 	TM[u][u-1]=0.0; 
//       } else { 
// 	it0=wave.begin(); it1=wave.end();
// 	it2=cabs.begin(); it3=cabs.end(); 
// 	if (idx1==0) {    // shortest wavelength falls below lower defined bound; use full grid
// 	  thisWave.assign(it0,it1); 
// 	  thisCabs.assign(it2,it3);
// 	} else {         // shortest wavelength falls somewhere in the middle
// 	  it0 += (idx1-1); it2 += (idx1-1);
// 	  thisWave.assign(it0,it1); 
// 	  thisCabs.assign(it2,it3);
// 	  thisWave[0]=_wT; 
// 	  thisCabs[0]=NumUtils::line(*it0,*(it0+1),*it2,*(it2+1),_wT);
// 	}
// 	integrand=NumUtils::bbodyCGS(thisWave,_temp[u]);
// 	transform(integrand.begin(),integrand.end(),thisCabs.begin(),integrand.begin(),multiplies<float>()); 
// 	TM[u][u-1] = Constant::PLANCKLIGHT/(enth[u]-enth[u-1])*NumUtils::integrate(thisWave,integrand); 
// 	TM[u-1][u-1] -= TM[u][u-1]; // accumulate diagonal
//       }
//       // Loop over final grain 'states' and comput heating elements. 
//       for (int l=u+1;l<nBins;++l) { 
// 	_wT = Constant::PLANCKLIGHT/(enth[u]-enth[l]);  // Wavelength of photon that can do transition.
// 	idx1 = NumUtils::index(_wT,wave);               // Where that photon falls in our defined grid:
// 	if ( (idx1==0) || (idx1==nWave) ) { // no such photon in our defined grid. 
// 	  TM[u][l] = 0.0; 
// 	} else {
// 	  // Interpolate cabs and J to _wT
// 	  it0=wave.begin()+idx1; it1=cabs.begin()+idx1; it2=J.begin()+idx1;
// 	  _cabs=NumUtils::line(*(it0-1),*it0,*(it1-1),*it1,_wT); 
// 	  _J = NumUtils::line(*(it0-1),*it0,*(it2-1),*it2,_wT);
// 	  TM[u][l] = Constant::IPLANCKLIGHT*_cabs*_J*denth[u]*pow(_wT,3); 
// 	  // If photons with energies capable of exciting beyond max bin exist.
//  	  // put them into the upper bin.
// 	  if (idx1 >= 2 && u==(nBins-1)) { 
// 	    thisWave.assign(wave.begin(),it0); 
// 	    thisnWave = thisWave.size()-1;
// 	    thisWave[thisnWave] = _wT; 
// 	    thisCabs.assign(cabs.begin(),it1);
// 	    thisCabs[thisnWave]=NumUtils::line(*(it0-1),*it0,*(it1-1),*it1,_wT); 
// 	    thisJ.assign(J.begin(),it2); 
// 	    thisJ[thisnWave]=NumUtils::line(*(it0-1),*it0,*(it2-1),*it2,_wT);
// 	    integrand.resize(thisnWave+1); 
// 	    ib=integrand.begin();
// 	    ie=integrand.end();
// 	    it1=thisCabs.begin();
// 	    it2=thisJ.begin(); 
// 	    for (it=ib;it!=ie;++it,++it1,++it2) *it = (*it1)*(*it2);
// 	    TM[u][l] += (_wT*NumUtils::integrate(thisWave,integrand)); 
// 	  }
// 	  TM[u-1][u-1] -= TM[u][l];
// 	}
//       }

//     }
    
 
    
//     // No Diagonals.
//     // Constant 4pi/hc left off of transition matrix to reduce multiplies. 
// //     float _wT;   // wavelength that will produce the transition. 
// //     float _cabs; // Absorption cross section at _wT
// //     float _J;    // Radition field at _wT
//     int uc; 
//     for (int l=0;l<nBins-1;++l) {
//       // Compute cooling elements
//       uc=l+1; 
//       _wT = Constant::PLANCKLIGHT/enth[uc]; 
//       idx1=NumUtils::index(_wT,wave);     
//       if (idx1==nWave) {  // beyond the upper bound of defined wavelength grid
// 	TM[l][uc] = 0.0; 
//       } else { 
// 	it0=wave.begin(); 
// 	it1=wave.end(); 
// 	it2=cabs.begin();
// 	it3=cabs.end(); 
// 	if (idx1!=0) { 
// 	  it0 += (idx1-1);
// 	  it2 += (idx1-1);
// 	  thisWave.assign(it0,it1);
// 	  thisWave[0]=_wT;
// 	  thisCabs.assign(it2,it3);
// 	  thisCabs[0]=NumUtils::line(*it0,*(it0+1),*it2,*(it2+1),_wT); 
// 	} else { 
// 	  thisWave.assign(it0,it1); 
// 	  thisCabs.assign(it2,it3);
// 	}
// 	integrand = NumUtils::bbodyCGS(thisWave,_temp[uc]);
// 	//cout << "size of integrand " << integrand.size() << endl; 
// 	transform(integrand.begin(),integrand.end(),thisCabs.begin(),integrand.begin(),multiplies<float>()); 
// 	TM[l][uc] = Constant::PLANCKLIGHT/(enth[uc]-enth[l])*NumUtils::integrate(thisWave,integrand);
//       }
//       for (int u=l+1;u<nBins;++u) {
// 	_wT = (Constant::PLANCKLIGHT)/(enth[u]-enth[l]); // transition wavelength
// 	idx1=NumUtils::index(_wT,wave); 
// 	if ( (idx1==0) || (idx1==nWave) ) { // Transition wavelength outside of defined field grid
// 	  TM[u][l] = 0.0; 
// 	} else { 
// 	  it0=wave.begin()+idx1; 
// 	  it1=cabs.begin()+idx1; 
// 	  it2=J.begin()+idx1; 
// 	  _cabs = NumUtils::line(*(it0-1),*it0,*(it1-1),*it1,_wT); 
// 	  _J = NumUtils::line(*(it0-1),*it0,*(it2-1),*it2,_wT);
// 	  TM[u][l] = Constant::IPLANCKLIGHT*_cabs*_J*denth[u]*pow(_wT,3); 
// 	}
	

//       }
//       if (idx1>=2) { // photons with wavelengths capable of exciting beyond the max bin exist.
// 	it = wave.begin();
// 	thisWave.assign(it,it0); 
// 	thisnWave = thisWave.size()-1; 
// 	thisWave[thisnWave] = _wT; 
// 	it = cabs.begin();
// 	thisCabs.assign(it,it1); 
// 	thisCabs[thisnWave]=NumUtils::line(*(it0-1),*it0,*(it1-1),*it1,_wT); 
// 	it = J.begin();
// 	thisJ.assign(it,it2); 
// 	thisJ[thisnWave]=NumUtils::line(*(it0-1),*it0,*(it2-1),*it2,_wT);
// 	integrand.resize(thisnWave+1); 
// 	ib = integrand.begin();
// 	ie = integrand.end();
// 	it1 = thisCabs.begin();
// 	it2 = thisJ.begin();
// 	for (it=ib;it!=ie;++it,++it1,++it2) *it = (*it1)*(*it2); 
// 	TM[nBins-1][l] += _wT*NumUtils::integrate(thisWave,integrand); 
//       }
//     }
//   }
//   cout << "*** without accumulating diagonal ***" << endl; 
//   for (int l=0;l<nBins;l++) {
//     for (int u=0;u<nBins;u++) { 
//       cout << TM[l][u] << " "; 
//     }
//     cout << endl; 
//   }
//  cout << nBins << " " << enthgrid.size() << " " << enth.size() << " " << denth.size() << endl; 
  // Compute _P
  
  //
  return _P; 
}
