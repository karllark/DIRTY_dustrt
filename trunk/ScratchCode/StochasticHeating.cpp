#include <iostream>
#include <vector> 

#include "HeatUtils.h" 
#include "NumUtils.h" 

using namespace std; 

int ComputeGrid(vector <float>& enth, vector <float>& denth, vector <float>& temp, 
		vector <float>& tgrid, vector <float> & Temperature, vector <float> & Enthalpy, 
		float & TMax, float & TMin, int & nBins);

void ComputeTransitionMatrix (vector <vector<double> >& TM, vector <float>& wave, vector <float>& temp, 
			      vector <float>& enth, vector <float>& cabs, vector <float>& cJprod, 
			      vector <float>& denth, int & nBins);

void SolveTransitionMatrix ( vector <vector<double> >& TM, int & nBins, vector <double>& P );

vector <double> StochasticHeating(vector <float> & wave, vector <float> & cJprod, 
				  vector <float> & cabs, vector <float> & Temperature, 
				  vector <float> & Enthalpy, float EAbs, float & TMin,
				  float & TMax, float & TEq)

{
  
  int maxBins = 1000; 
  bool converged=false; 
  bool IncreaseBins=false;
  bool lastincrease=true; 
  bool OnePass=false; 
  int nBins=50; 
  //int oldnbins=nBins;
  float oldTMax,oldTMin; 
  float tol = 0.01; 
  float tol_max_bins = 0.1; 
  float thistol,lasttol; 
  double Ptol = 1.0e-15; 
  vector <vector<double> > TM;
  //vector <vector<double> > Bij; 
  vector <double> _P; 
  //double maxP; 
  //
  int status;
  //int thisnWave; 
  int nWave=wave.size(); 
  //int nWaveLast=nWave-1; // point to last element of wave - less arithmetic in loops.
  
  int idx;//,idx1,idxp;
  vector <double> _pofTint; 

  double _pFac; 

  //int crap1=0;
  //int crap2=1;

  vector <float>::iterator it,it0,it1,it2,it3,it4,ib,ie; 
  vector <double>::iterator idt,idb,ide,idt1,idb1,ide1,idt2; 

  // Some initial algorithm to compute search size in constant T

  //float arg1,arg2,w1,w2,w3,w4,wc,val;
  vector <float> _enth,_denth,_temp,_tgrid;

  vector <float> thisWave,thisCabs,thisJ;

  //double _norm; 
  // Size the transition matrix.

  vector <double> integrand; 

  double Eemit; 

  //vector <vector<double> > TM = TransitionMatrix(_enth,_denth,_temp,wave,cabs,J,nBins);
  //float _wT;   // wavelength that will produce the transition. 
  //float _cabs; // Absorption cross section at _wT
  //float _J;    // Radition field at _wT
							
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
  _tgrid.reserve(maxBins+1); 
  _P.resize(maxBins,0); 
  TM.resize(maxBins);
  for (int i=0;i<maxBins;++i) TM[i].resize(maxBins);
  
  //Bij.resize(maxBins); 
  //cout << " in stochastic, initilizing 3.0 " << endl;
  //for (int i=0;i<maxBins;++i) { TM[i].resize(maxBins); Bij[i].resize(maxBins); }
  
  //cout << " in stochastic, initilizing 3.1 " << endl;

  // Isolate temperature region.
  //cout << "Setting up..." << endl; 
  _enth.resize(nBins,0); _denth.resize(nBins,0); _temp.resize(nBins,0); _tgrid.resize(nBins+1,0);
  //cout << "Computing grid 1..." << endl; 
  status = ComputeGrid(_enth,_denth,_temp,_tgrid,Temperature,Enthalpy,TMax,TMin,nBins); 
  idx = 1; 
  float MaxE = Constant::PLANCKLIGHT/wave[0]; 
  while (_enth[idx]-_enth[0] < MaxE) ++idx; 
  TMax=_tgrid[idx]; 
  if (TMax < 1.5*TEq) { 
    if (TEq < 100.) TMax=1.5*TEq; else TMax=TEq+100.; 
  }
  bool _setup=true; 
  while (_setup) { 
    // cout << "computing grid 2 ..." << endl; 
    status = ComputeGrid(_enth,_denth,_temp,_tgrid,Temperature,Enthalpy,TMax,TMin,nBins); 
    //cout << "Computing transition matrix 1 " << endl; 
    //cout << nBins << " " << wave.size() << " " << cabs.size() << " " << cJprod.size() << " "<< _temp.size() << " " << _denth.size() << " " << _enth.size() << endl ;
    ComputeTransitionMatrix(TM,wave,_temp,_enth,cabs,cJprod,_denth,nBins); 
    //cout << "solving transition matrix 2 " << endl; 
    SolveTransitionMatrix(TM,nBins,_P); 
    if (_P[nBins-1] > Ptol) 
      TMax *= 1.5; 
    else 
      _setup=false; 
  }
  //cout << "Starting convergence calculations." << endl; 

  // Convergence wrapper... 
  while (!converged) { // convergence bracket
 
    //cout << "Computing grid " << endl; 
    _enth.resize(nBins,0); _denth.resize(nBins,0); _temp.resize(nBins,0); _tgrid.resize(nBins+1,0); 
    status = ComputeGrid(_enth,_denth,_temp,_tgrid,Temperature,Enthalpy,TMax,TMin,nBins); 
    //cout << "Computing matrix " << endl; 
    ComputeTransitionMatrix(TM,wave,_temp,_enth,cabs,cJprod,_denth,nBins); 
    //cout << "Solving matrix " << endl; 
    SolveTransitionMatrix(TM,nBins,_P); 

    //cout << "Computing emission " << endl; 
    integrand.resize(nWave); 
    idb=integrand.begin();
    ide=integrand.end();
    it1=wave.begin(); 
    it=cabs.begin();
    for (idt=idb;idt!=ide;++idt,++it,++it1) { 
      _pofTint = NumUtils::prod_bbodyCGS<double>(*it1,_temp,_P); // P(T)*B(T) at wv,_sz
      _pFac = 0.0; 
      _pFac = accumulate(_pofTint.begin(),_pofTint.end(),0.0);   // Sum_T (P(T)*B(T))
      *idt = (*it)*_pFac; 
    }
    Eemit = NumUtils::integrate<double>(wave,integrand); 
    thistol = abs(EAbs-Eemit)/EAbs;
    //cout << "adjusting bins "<< endl; 
    //cout << "convergence " << nBins << " " << TMin << " " << TMax << " " <<  EAbs << " " << Eemit << " " << thistol << " "  << lasttol << endl;

    if (thistol > tol) { // not converged
      
      if (thistol > 100.*lasttol && OnePass) {
      	IncreaseBins=true;
	TMin = oldTMin;
	TMax = oldTMax; 
      } else { 
	lasttol = thistol; // save current tolerance.
	oldTMin = TMin; 
	oldTMax = TMax; 
	idx = nBins-1; 
	if (_P[idx] == 0.0) { 
	  while (_P[idx] == 0.0) --idx; 
	  if (idx < 0) { 
	    cout << "Failure 1.0, stochastic heating algorithm" << endl; 
	    exit(8); 
	  }
	  if (_P[idx] < Ptol) { 
	    while (_P[idx] < Ptol) --idx; 
	    if (idx < 0) { 
	      cout << "Failure 1.1, stochastic heating algorithm" << endl; 
	      exit(8); 
	    } 
	    TMax = _tgrid[idx+1]; 
	  } else { 
	    TMax = _tgrid[nBins]; 
	  }
	}
	
	if (_P[idx] < Ptol) { 
	  while (_P[idx] < Ptol) --idx; 
	  if (idx < 0) { 
	    cout << "Failure 2.0, stochastic heating algorithm" << endl; 
	    exit(8); 
	  }
	  TMax=_tgrid[idx+1]; 
	}
	
	if (_P[nBins-1] >= Ptol) { 
	  TMax = 1.5*_tgrid[nBins]; 
	  IncreaseBins=true; 
	}
	
	idx=0; 
	if (_P[idx] < Ptol) { 
	  while (_P[idx] < Ptol) ++idx; 
	  if (idx > nBins-1) { 
	    cout << "Failure 3.0, stochastic heating algorithm" << endl; 
	    exit(8); 
	  }
	  if (idx == 0) 
	    TMin = _tgrid[idx]; 
	  else 
	    TMin = _tgrid[idx-1]; 
	}
      }
      
      if (IncreaseBins) { 
	nBins = static_cast<int>(1.5*static_cast<double>(nBins)); 
      }

      if (nBins > maxBins)  {
	if (lastincrease) { 
	  // accept a lower tolerance for these cases
	  cout << "EXCEEDED BIN COUNT IN StochasticHeating()" << endl; 
	  if (thistol > tol_max_bins) { // still not converged well enough
	    //for (int ii=0;ii<oldnbins;++ii) cout << ii << " " << _temp[ii] << " " << _P[ii] << endl; 
	    //for (int ii=0;ii<nWave;++ii) cout << ii << " " << wave[ii] << " " << J[ii] << " " << integrand[ii] << endl; 
	    cout << EAbs << " " << Eemit << " " << thistol << endl; 
	    exit(8); 
	  } else {
	    cout << "but close enough (tol_max_bins = " << tol_max_bins << "; tol = " << thistol << endl;
	    converged=true;
	  }
	} // else { 
// 	  cout << "Maximum bins will be exceeded - one last pass" << endl; 
// 	  //for (int ii=0;ii<oldnbins;++ii) cout << ii << " " << _temp[ii] << " " << _P[ii] << endl; 
// 	  //for (int ii=0;ii<nWave;++ii) cout << ii << " " << wave[ii] << " " << J[ii] << " " << integrand[ii] << endl; 
// 	  cout << EAbs << " " << Eemit << " " << thistol << endl; 
// 	  lastincrease=true;
// 	  nBins=maxBins; 
// 	  cout << TMin << " " << TMax << endl; 
// 	}	
      }
      OnePass=true; 
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
  for (_it=_it1;_it!=_ite;++_it,++_it2) { *_it=*(_it-1)+_DT;  *_it2=(*_it+*(_it-1))/2.0; }

  // Interpolate Enthalpy onto _tgrid
  //cout << _enthgrid.size() << endl; 
  _enthgrid = NumUtils::interpol(Enthalpy,Temperature,tgrid,4,-99); 
  //NumUtils::interpol(Enthalpy,Temperature,tgrid,_enthgrid,4,-99);
  //cout << _enthgrid.size() << endl;
 //  int output=1; 
//   if (output == 1) for (int i=0;i<_enthgrid.size();++i) cout << i << " " << tgrid[i] << " " << _enthgrid[i] << endl; 
  // Compute enthaply at bin center along with bin width. 
  _itb = _enthgrid.begin()+1;
  _ite = _enthgrid.end(); 
  _it1 = enth.begin(); 
  _it2 = denth.begin(); 
  for (_it=_itb;_it!=_ite;++_it,++_it1,++_it2) { *_it1=(*_it+*(_it-1))/2.0; *_it2 = (*_it-*(_it-1)); }

  return 0; 
}


void ComputeTransitionMatrix (vector <vector<double> >& TM, vector <float>& wave, vector <float>& temp, 
			      vector <float>& enth, vector <float>& cabs, vector <float>& cJprod, 
			      vector <float>& denth, int & nBins) 
{

  float _wT,_slp,_icpt,_thiscjprod; 
  int _nw=wave.size(); 
  int idx; 

  double c1 = 2.0*Constant::PLANCK*pow(Constant::LIGHT,2)*Constant::PLANCKLIGHT;
  double c2 = Constant::PLANCKLIGHT/Constant::BOLTZMAN; 
  
  vector <float> _wave,_cjprod,_integrand; 
  vector <float>::iterator it,itb,ite,it0,it1;

  //cout << "first: " << TM.size() << endl; 
  //for (int f=0;f<nBins;++f) cout << f << " " << TM[f].size() << endl;

  _integrand.resize(_nw);
 
  for (int f=0;f<nBins;++f) { // Loop over final energy states
    
    if (f != nBins-1) { // Cooling transitions - only f+1 to f
      //cout << "cooling " << endl; 
      itb=wave.begin(); ite=wave.end();
      //_integrand.resize(_nw); 
      it0 = _integrand.begin();
      it1 = cabs.begin();
      for (it=itb;it!=ite;++it,++it0,++it1) *it0=(*it1/pow(*it,5))*(1.0/(std::exp(c2/((*it)*temp[f+1]) )-1.0)); 
      TM[f][f+1] = c1/(enth[f+1]-enth[f])*NumUtils::integrate<double>(wave,_integrand); 
    } // Done with cooling transitions. 
    
    for (int i=0;i<f;++i) { // Heating transitions - all initial < final.
      //cout << " heating " << endl; 
      if ((enth[f]-enth[i]) != 0) _wT = Constant::PLANCKLIGHT/(enth[f]-enth[i]); 
      else { cout << "zero denominator " << f << " " << i << endl; exit(8);}// wavelength of transitions
      if (_wT < wave[0] || _wT > wave[_nw-1]) TM[f][i] =0.0; 
      else { // Photons exist that will induce this transition
	idx=0; 
	while (wave[idx] < _wT) idx++;
	if (idx > _nw-1) { cout << "idx error!! " << _wT << " "<< wave[0] << " " << wave[_nw-1] << " " << idx << " " << wave.size() << endl; exit(8); }
	if (idx == 0 || idx == _nw-1) _thiscjprod=cJprod[idx];
	else { 
	  if ((wave[idx]-wave[idx-1]) != 0)
	    _slp=(cJprod[idx]-cJprod[idx-1])/(wave[idx]-wave[idx-1]);
	  else { cout << " zero demon, wave-wave " << idx << endl; exit(8); }
	  _icpt = cJprod[idx]-_slp*wave[idx]; 
	  _thiscjprod = _icpt + _slp*_wT;
	} 
	TM[f][i] = Constant::IPLANCKLIGHT*_thiscjprod*_wT*_wT*_wT*denth[f]; 
	if (f==nBins-1 && idx > 0) { // Put all transitions out of defined states into last state. 
	  it0=wave.begin()+idx; 
	  it1=cJprod.begin()+idx; 
	  _wave.assign(wave.begin(),it0); 
	  _cjprod.assign(cJprod.begin(),it1); 
	  _wave[_wave.size()-1]=_wT; 
	  _cjprod[_wave.size()-1]=_thiscjprod; 
	  TM[f][i] += (_wT*NumUtils::integrate<double>(_wave,_cjprod)); 
	}
	
      } // Done with photons that can produce transition.

    } // Done with heating. 
  } // Done with final statues. 
  //cout << " computing diagonals" << endl; 
  for (int f=0;f<nBins;++f) TM[0][0] -= TM[f][0]; 
  for (int i=1;i<nBins;++i) 
    for (int j=i-1;j<nBins;++j)
      if (i!=j) TM[i][i] -= TM[j][i]; 

}

void SolveTransitionMatrix ( vector <vector<double> >& TM, int & nBins, vector <double>& P )
{
  
  vector <vector<double> > _Bij; 
  _Bij.resize(nBins);  for (int i=0;i<nBins;++i) _Bij[i].resize(nBins); 
  for (int i=0;i<nBins;++i) { 
    _Bij[nBins-1][i]=TM[nBins-1][i]; 
    for (int j=nBins-2;j>=0;--j) { 
      _Bij[j][i] = _Bij[j+1][i]+TM[j][i]; 
    }
  }

  P[0]=1.0; 
  double _norm=1.0; 
  for (int f=1;f<nBins;++f) { 
    P[f] = 0.0; 
    for (int i=0;i<f;++i) { 
      P[f] += (P[i]*_Bij[f][i]); 
    } 
    if (TM[f-1][f] != 0.0)
      P[f] /= TM[f-1][f]; 
    else {
      cout << "Problem with matrix cooling element...." << endl; 
      cout << "Can't solve the transition matrix." << endl; 
      exit(8); 
    }
    _norm += P[f]; 
  }
  
  transform(P.begin(),P.end(),P.begin(),bind2nd(divides<double>(),_norm)); 

}
