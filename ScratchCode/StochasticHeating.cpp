#include <iostream>
#include <vector> 

#include "HeatUtils.h" 
#include "NumUtils.h" 

using namespace std; 

int ComputeGrid(vector <float>& enth, vector <float>& denth, vector <float>& temp, 
		vector <float>& tgrid, vector <float> & Temperature, vector <float> & Enthalpy, 
		float & TMax, float & TMin, int & nBins);

vector <double> StochasticHeating(vector <float> & waveE, vector <float> & wave, vector <float> & J, 
				  vector <float> & cabs, vector <float> & Temperature, 
				  vector <float> & Enthalpy, float EAbs, float & TMin, float & TMax) 

{
  
  bool converged=false; 
  bool IncreaseBins=false; 
  int nBins=50; 
  float tol = 0.01; 
  float thistol; 
  double Ptol = 1.0e-14; 
  vector <vector<double> > TM;
  vector <vector<double> > Bij; 
  vector <double> _P; 
  double maxP; 

  int status;
  int thisnWave; 
  int nWave=wave.size(); 
  int nWaveLast=nWave-1; // point to last element of wave - less arithmetic in loops.
  
  int idx,idx1,idxp;
  vector <double> _pofTint; 

  double _pFac; 

  vector <float>::iterator it,it0,it1,it2,it3,it4,ib,ie; 
  vector <double>::iterator idt,idb,ide,idt1,idb1,ide1,idt2; 

  // Some initial algorithm to compute search size in constant T

  //float arg1,arg2,w1,w2,w3,w4,wc,val;
  vector <float> _enth,_denth,_temp,_tgrid;

  vector <float> thisWave,thisCabs,thisJ;

  double _norm; 
  // Size the transition matrix.

  vector <double> integrand; 

  double Eemit; 

  //vector <vector<double> > TM = TransitionMatrix(_enth,_denth,_temp,wave,cabs,J,nBins);
  float _wT;   // wavelength that will produce the transition. 
  float _cabs; // Absorption cross section at _wT
  float _J;    // Radition field at _wT
  
  // Convergence wrapper... 
  while (!converged) {
        
    _enth.resize(nBins,0); _denth.resize(nBins,0); _temp.resize(nBins,0); _tgrid.resize(nBins+1,0); 
    TM.resize(nBins); 
    _P.resize(nBins,0); 
   
    status = ComputeGrid(_enth,_denth,_temp,_tgrid,Temperature,Enthalpy,TMax,TMin,nBins); 
    
    for (int f=0;f<nBins;++f) {  // Loop over available final states 
      
      TM[f].resize(nBins,0);  // Allocate second dim of TM
      
      if (f != nBins-1) { // Cooling transitions - only f+1 to f
	integrand = NumUtils::prod_bbodyCGS<double>(wave,_temp[f+1],cabs); 
	TM[f][f+1] = Constant::PLANCKLIGHT/(_enth[f+1]-_enth[f])*NumUtils::integrate<double>(wave,integrand);
      }
      
      for (int i=0;i<f;++i) { // Heating transitions - all initial < final
	_wT=Constant::PLANCKLIGHT/(_enth[f]-_enth[i]); // wavelength of transition
	if (_wT < wave[0] || _wT >= wave[nWaveLast]) {	    
	  TM[f][i]=0.0; // No photons available to do transition.
	} else { // find out where the photon falls in our grid and interpolate c/J
	  
	  idx=NumUtils::index(_wT,wave);  
	  it0=wave.begin()+idx; it1=cabs.begin()+idx; it2=J.begin()+idx; 
	  _cabs=NumUtils::line(*(it0-1),*it0,*(it1-1),*it1,_wT); 
	  _J=NumUtils::line(*(it0-1),*it0,*(it2-1),*it2,_wT); 
	  TM[f][i] = Constant::IPLANCKLIGHT*_cabs*_J*_denth[f]*pow(_wT,3);  
	  
	  if (f == nBins-1) { // put all transitions to states beyond our defined grid into last bin
	    thisWave.assign(wave.begin(),it0); 
	    thisnWave = thisWave.size()-1; 
	    thisWave[thisnWave] = _wT; 
	    thisCabs.assign(cabs.begin(),it1); 
	    thisCabs[thisnWave] = NumUtils::line(*(it0-1),*it0,*(it1-1),*it1,_wT); 
	    thisJ.assign(J.begin(),it2); 
	    thisJ[thisnWave]=NumUtils::line(*(it0-1),*it0,*(it2-1),*it2,_wT); 
	    integrand.resize(thisnWave+1); 
	    idb=integrand.begin(); 
	    ide=integrand.end(); 
	    it1=thisCabs.begin(); 
	    it2=thisJ.begin(); 
	    for (idt=idb;idt!=ide;++idt,++it1,++it2) *idt=((*it1)*(*it2));
	    TM[f][i] += (_wT*NumUtils::integrate<double>(thisWave,integrand)); 
	  }
	}
      }
    }
    
    // Populate the Diagonal.
    for (int f=1;f<nBins;++f) {
      TM[0][0] -= TM[f][0]; 
      for (int i=f-1;i<nBins;++i) 
	if (f != i) TM[f][f] -= TM[i][f]; 
    }
    
    // Compute B{i,j}
    Bij.resize(nBins); 
    
    for (int i=0;i<nBins;++i) Bij[i].resize(nBins);
    for (int i=0;i<nBins;++i) { 
      Bij[nBins-1][i]=TM[nBins-1][i]; 
      for (int j=nBins-2;j>=0;--j) { 
	Bij[j][i] = Bij[j+1][i]+TM[j][i]; 
      }     
    }
    
    // Compute P
    _P[0]=1.0; 
    _norm=_P[0]; 
    for (int f=1;f<nBins;++f) {
      _P[f] = 0.0; // initialize
      for (int i=0;i<f;++i) { 
	_P[f] += (_P[i]*Bij[f][i]);
      } 
      if (TM[f-1][f] != 0.0) 
	_P[f] /= TM[f-1][f];
      else {
	cout << "Problem with matrix cooling element...." << endl; 
	cout << "Can't solve the transition matrix." << endl; 
	exit(8); 
      }
      _norm += _P[f]; 
    }
    
    // Normalize _P
    idb=_P.begin(); ide=_P.end(); 
    for (idt=idb;idt!=ide;++idt) *idt /= _norm; 

    // Compute emission... 
    // Compute integral over T and return for each lambda.
    integrand.resize(nWave); 
    idb=integrand.begin();
    ide=integrand.end();
    it1=wave.begin(); 
    it=cabs.begin();
    for (idt=idb;idt!=ide;++idt,++it,++it1) { 
      _pofTint = NumUtils::prod_bbodyCGS<double>(*it1,_temp,_P); // P(T)*B(T) at wv,_sz
      _pFac = accumulate(_pofTint.begin(),_pofTint.end(),0.0);   // Sum_T (P(T)*B(T))
      *idt = (*it)*_pFac; 
    }
    Eemit = NumUtils::integrate<double>(wave,integrand); 
    thistol = abs(EAbs-Eemit)/EAbs; 
   
    if (thistol > tol) { // not converged
      if (_P[0] >= 1.) {  // everything is packed into the lowest bin - too broad _tgrid
	TMin = _tgrid[0]; 
	TMax = _tgrid[1]; 
	IncreaseBins = false; 
      } else { 
	idxp = NumUtils::maxID(_P); 
	maxP = _P[idxp]; 
	idx = NumUtils::index(maxP*Ptol,_P); 
	idx1 = NumUtils::rindex(maxP*Ptol,_P); 
	if (idx >= idx1) idx=0; 

	// Adjust upper temperature limit up by 25% or down by probability limi
	if (idx1 == nBins) {
	  TMax = 1.25*TMax - 0.25*TMin; 
	  IncreaseBins = true; 
	} else { 
	  TMax = _tgrid[idx1];
	  IncreaseBins = false;
	}
	
	// Adjust lower temperature limit down by 25% or up by probability limit 
	if (idx == 0) {
	  TMin = 0.75*TMin; 
	  IncreaseBins = true; 
	} else { 
	  TMin = _tgrid[idx]; 
	  IncreaseBins = false; 
	}
      }

      // Increase nbins by 25%
      if (IncreaseBins) 
	nBins = static_cast<int>(static_cast<float>(nBins)*1.25); 
    } else converged=true;  


  }

  return integrand; // this is C(lam)*Sum_T (P(T)*B(T,lam)) ~ stochastic L(lam) 

}


// Generate the temperature and enthalpy grids necessary for calculation.
int ComputeGrid(vector <float>& enth, vector <float>& denth, vector <float>& temp, 
		vector <float>& tgrid, vector <float> & Temperature, vector <float> & Enthalpy, 
		float & TMax, float & TMin, int & nBins)
{

  vector <float>::iterator _itb,_ite,_it,_it1,_it2; 
  float _DT=(TMax-TMin)/static_cast<float>(nBins);  // Linear temperature grid 

  vector <float> _enthgrid(nBins+1); 
 
  // Define the temperature grid, end- and mid-points. 
  _it1 = tgrid.begin(); 
  _it2 = temp.begin();
  *_it1 = TMin; 
  ++_it1; 
  _ite = tgrid.end(); 
  for (_it=_it1;_it!=_ite;++_it,++_it2) { *_it=*(_it-1)+_DT; *_it2=(*_it+*(_it-1))/2.0; }

  // Interpolate Enthalpy onto _tgrid
  _enthgrid = NumUtils::interpol(Enthalpy,Temperature,tgrid,4,1); 

  // Compute enthaply at bin center along with bin width. 
  _itb = _enthgrid.begin()+1;
  _ite = _enthgrid.end(); 
  _it1 = enth.begin(); 
  _it2 = denth.begin(); 
  for (_it=_itb;_it!=_ite;++_it,++_it1,++_it2) { *_it1=(*_it+*(_it-1))/2.0; *_it2 = (*_it-*(_it-1)); }

  return 0; 
}
