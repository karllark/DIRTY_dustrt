#include <iostream>
#include <vector> 

#include "HeatUtils.h" 
#include "NumUtils.h" 

using namespace std; 

int ComputeGrid(vector <float>& enth, vector <float>& denth, vector <float>& temp, 
		vector <float>& tgrid, vector <float> & Temperature, vector <float> & Enthalpy, 
		float & TMax, float & TMin, int & nBins);

vector <double> StochasticHeating(vector <float> & wave, vector <float> & J, 
				  vector <float> & cabs, vector <float> & Temperature, 
				  vector <float> & Enthalpy, float EAbs, float & TMin, float & TMax)

{
  
  int maxBins = 1000; 
  bool converged=false; 
  bool IncreaseBins=false; 
  int nBins=50; 
  int oldnbins=nBins;
  float tol = 0.01; 
  float tol_max_bins = 0.1; 
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
  
  //cout << " in stochastic, initilizing 1.0 " << endl; 
  // Reserve maximum sizes for our vectors.     
  // wavelength
  thisWave.reserve(nWave); 
  thisCabs.reserve(nWave); 
  thisJ.reserve(nWave);
  integrand.reserve(nWave); 
  // bins
  //cout << " in stochastic, initilizing 2.0 " << endl;
  _enth.reserve(maxBins); 
  _denth.reserve(maxBins); 
  _temp.reserve(maxBins); 
  _tgrid.reserve(maxBins); 
  _P.reserve(maxBins); 
  TM.reserve(maxBins);
  Bij.reserve(maxBins); 
  //cout << " in stochastic, initilizing 3.0 " << endl;
  //for (int i=0;i<maxBins;++i) { TM[i].reserve(maxBins); Bij[i].reserve(maxBins); }
  //cout << " in stochastic, initilizing 3.1 " << endl;
  // Convergence wrapper... 
  while (!converged) { // convergence bracket
        
    //cout << "not converged 1" << endl; 
    _enth.resize(nBins,0); _denth.resize(nBins,0); _temp.resize(nBins,0); _tgrid.resize(nBins+1,0); 
    TM.resize(nBins); 
    _P.resize(nBins,0); 
   
    status = ComputeGrid(_enth,_denth,_temp,_tgrid,Temperature,Enthalpy,TMax,TMin,nBins); 
    //cout << "not converged 2 " << endl;
    for (int f=0;f<nBins;++f) {  // Loop over available final states 
      //cout << f <<  " f states, not converged 3.0.0 " << endl; 
      TM[f].resize(nBins,0.0);  // Allocate second dim of TM
      //cout << f <<  " f states, not converged 3.0 " << endl; 
      if (f != nBins-1) { // Cooling transitions - only f+1 to f
	integrand = NumUtils::prod_bbodyCGS<double>(wave,_temp[f+1],cabs); 
	TM[f][f+1] = Constant::PLANCKLIGHT/(_enth[f+1]-_enth[f])*NumUtils::integrate<double>(wave,integrand);
      }
      
      for (int i=0;i<f;++i) { // Heating transitions - all initial < final
	_wT=Constant::PLANCKLIGHT/(_enth[f]-_enth[i]); // wavelength of transition
	if (_wT <= wave[0] || _wT >= wave[nWaveLast]) {	    
	  //cout << i << " f i states, not converged 4 " << endl; 
	  TM[f][i]=0.0; // No photons available to do transition.
	} else { // find out where the photon falls in our grid and interpolate c/J
	  
	  idx=NumUtils::index(_wT,wave);  
	  it0=wave.begin()+idx; it1=cabs.begin()+idx; it2=J.begin()+idx; 
	  _cabs=NumUtils::line(*(it0-1),*it0,*(it1-1),*it1,_wT); 
	  _J=NumUtils::line(*(it0-1),*it0,*(it2-1),*it2,_wT); 
	  TM[f][i] = Constant::IPLANCKLIGHT*_cabs*_J*_denth[f]*pow(_wT,3);  
	  //cout << i << " f i states, not converged 4.1 " << endl; 
	  if (f == nBins-1) { // put all transitions to states beyond our defined grid into last bin
	    thisWave.assign(wave.begin(),it0); 
	    thisnWave = thisWave.size()-1; 
	    thisWave[thisnWave] = _wT; 
	    thisCabs.assign(cabs.begin(),it1); 
	    //cout << i << "  f i states, not converged 4.2 " << idx  << "  " << wave.size() <<  endl; 
	    thisCabs[thisnWave] = NumUtils::line(*(it0-1),*it0,*(it1-1),*it1,_wT);
	    //cout << i << "  f i states, not converged 4.3 " << endl;  
	    thisJ.assign(J.begin(),it2); 
	    thisJ[thisnWave]=NumUtils::line(*(it0-1),*it0,*(it2-1),*it2,_wT); 
	    //cout << i << "  f i states, not converged 4.3 " << endl; 
	    integrand.resize(thisnWave+1); 
	    idb=integrand.begin(); 
	    ide=integrand.end(); 
	    it1=thisCabs.begin(); 
	    it2=thisJ.begin(); 
	    //cout << i << "  f i states, not converged 4.4 " << endl; 
	    for (idt=idb;idt!=ide;++idt,++it1,++it2) *idt=((*it1)*(*it2));
	    //cout << i << "  f i states, not converged 4.5 " << endl; 
	    TM[f][i] += (_wT*NumUtils::integrate<double>(thisWave,integrand)); 
	    //cout << i << "  f i states, not converged 4.6 " << endl; 
	  }
	}
      }
    }
    // this was the last one.
    //cout << "not converged 2 " << endl; 
    // Populate the Diagonal.
    //cout << "  out of states, not converged 5.0 " << endl; 
    for (int f=1;f<nBins;++f) {
      TM[0][0] -= TM[f][0]; 
      for (int i=f-1;i<nBins;++i) 
	if (f != i) TM[f][f] -= TM[i][f]; 
    }
    //cout << "  out of states, not converged 5.1 " << endl; 
    // Compute B{i,j}
    Bij.resize(nBins); 
    
    for (int i=0;i<nBins;++i) Bij[i].resize(nBins);
    for (int i=0;i<nBins;++i) { 
      Bij[nBins-1][i]=TM[nBins-1][i]; 
      for (int j=nBins-2;j>=0;--j) { 
	Bij[j][i] = Bij[j+1][i]+TM[j][i]; 
      }     
    }
    //cout << "  out of states, not converged 5.2 " << endl; 
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
    //cout << "  out of states, not converged 5.3 " << endl; 
    // Normalize _P
    idb=_P.begin(); ide=_P.end(); 
    for (idt=idb;idt!=ide;++idt) *idt /= _norm; 
    //cout << "  out of states, not converged 5.4 " << endl; 
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
    //cout << "  out of states, not converged 5.5 " << endl; 
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
      if (IncreaseBins) {
	oldnbins=nBins; 
	nBins = static_cast<int>(static_cast<float>(nBins)*1.25); 
      }
      
      if (nBins > maxBins) { 
	// accept a lower tolerance for these cases
	cout << "EXCEEDED BIN COUNT IN StochasticHeating()" << endl; 
	if (thistol > tol_max_bins) { // still not converged well enough
	  for (int ii=0;ii<oldnbins;++ii) cout << ii << " " << _temp[ii] << " " << _P[ii] << endl; 
	  for (int ii=0;ii<nWave;++ii) cout << ii << wave[ii] << " " << J[ii] << endl; 
	  cout << EAbs << " " << Eemit << " " << thistol << endl; 
	  exit(8); 
	} else {
	  cout << "but close enough (tol_max_bins = " << tol_max_bins << "; tol = " << thistol << endl;
	  converged=true;
	}
      }
    } else converged=true;  

    //cout << "not converged 2 " << endl; 

  }
  //cout << "  out of states, not converged 6.0 " << endl; 
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
