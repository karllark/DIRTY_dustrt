// Jun 3-4 2009: Modified to be DirtyFailure aware.  Dust emission failures will pass failure information 
//               back up the call stack. 
#include "HeatUtils.h"
#include "GrainModel.h"
#include "NumUtils.h"

#include "DirtyFailure.h"

//#include "DirtyFlags.h"

extern int StochasticHeating(vector <float> & wave, 
			     vector <float> & J, vector <float> & cabs, 
			     vector <float> & Temp, vector <float> & Enthalpy, 
			     float EAbs, float & TMin, float & TMax, float & TEq,
			     uint _sz, vector <vector<double> > & StochasticLum);

int ComputeDustEmission (vector <float> & J, GrainModel & GrainModel, 
			 vector <vector<double> > & EmmittedEnergy, bool & DoStochastic,
			 float & _FailureSz, int & _FailureComp)
{
 
  // Stochastic state for each size. 
  vector <bool> _dostochastic;
  
  // Hold eq. and st. luminosity as f(size,wave)
  vector <vector<double> > StochasticLum; 
  vector <vector<double> > EquilibriumLum; 
  // Iterators for: 
  // _iteec == emmitted equilibrium for each component at each wavelength. 
  // _itees == emmitted stochastic for each compoenent at each wavelength. 
  // _itet == total emmitted energy for each wavelength.
  vector <double>::iterator _iteec,_itees,_itet; 

  // Some temporary vectors. 
  vector <double> _integrand,_integrand1;
  vector <float> _t;
 
  // Some iterators.
  vector <float>::iterator _itb,_ite;
  vector <float>::iterator _it,_it1,_it2; 

  // Some Grain model derived parameters.
  vector <float> _w = GrainModel.getWave();
  vector <float> _E = GrainModel.getEScale();
  int _ncmp = GrainModel.getNComp(); 
  uint _nw = _w.size();
  vector <float> _Enthalpy; 
  vector <float> _CalTemp;
  vector <float> _cabs,_size,_sizeDist; 
  vector <float> _cJprod(_nw); 
  uint _nsize; 

  int status; 

  // Check J is defined on same wave grid as wave....
  if (J.size() != _nw) { 
    return Flags::FCDE_VECTOR_SIZE_MISMATCH; 
  }
  // Some local variables.  
  // EAbs computed in EqTemp() and passed back for StochasticHeating() to test convergence. 
  float EAbs;
  float _tlo,_thi,_norm,mpe,tau_abs,tau_rad,Tu; 

  // Some default starting values for stochastic caluclation - will return with best guess 
  // from StochasticHeating for next size. 
  int nBins;
  float TMin; 
  float TMax; 

  // Local J that includes the microwave background - Otherwise the really small guys 
  // want to cool to ~0K. 
  //vector <float> _localJ=NumUtils::add_bbodyCGS<float>(_w,2.7,J); 
  //for (int i=0;i<_nw;++i) cout << i << " " << _w[i] <<" " << J[i] << " " << _localJ[i] << endl;
  // Fit power law to J. 
//   vector <float> _LogW(_nw); 
//   vector <float> _LogJ(_nw); 
//   vector <float> _localJ(_nw);
//   transform(_w.begin(),_w.end(),_LogW.begin(),NumUtils::log10<float>());
//   transform(J.begin(),J.end(),_LogJ.begin(),NumUtils::log10<float>());
//   vector <float> _param(2),_sigma(2); 
//   _param=NumUtils::poly_fit(_LogW,_LogJ,1,_localJ,_sigma); 
//   transform(_localJ.begin(),_localJ.end(),_localJ.begin(),NumUtils::pow10<float>()); 
  //for (uint i=0;i<_nw;++i) cout << i << " "<< _w[i] << " " << J[i] << endl; 
  for (int _cmp=0;_cmp<_ncmp;++_cmp) { // Component loop

    _FailureComp = _cmp; 
    // Default starting values.  Works well for size distributions that start small...
    nBins=50; 
    TMin=0.1; 
    TMax=5000.;

    // Populate GrainModel() for each component
    _nsize = GrainModel.nSize(_cmp);
    _cabs = GrainModel.CAbs(_cmp,(int)0); 
    // initilize _t and _dostochastic
    _t.resize(_nsize);
    _dostochastic.resize(_nsize,DoStochastic);  // Populate local with global default. 
    // and first dimension of Luminosities. 
    EquilibriumLum.resize(_nsize);
    StochasticLum.resize(_nsize); 

    // Start values for EqTemp(). 
    _it=_t.begin(); 
    _tlo = 1.0; 
    _thi = 2500.0; 

    _size = GrainModel.Size(_cmp);
    // Loop over all sizes. 
    for (uint _sz=0;_sz<_nsize;++_sz,++_it) { // Size loop

      _FailureSz = _size[_sz]; 
      //cout << "Component " << _cmp << " size " << _sz << " " << _size[_sz] << endl; 
      _cabs = GrainModel.CAbs(_cmp,int(_sz));
      for (uint _wv=0;_wv<_nw;++_wv) _cJprod[_wv] = _cabs[_wv]*J[_wv]; 
      
      // Get equilibrium temperature.  Return emitted energy as well since it's used 
      // in stochastic calculation. 
      status = EqTemp(_w,J,_cabs,EAbs,*_it,((_tlo<0)?1.0:_tlo),((_thi>2500.0)?2500.0:_thi));   
      if (status != Flags::FSUCCESS) return status; // we've had an equilibrium heating failure. 
      
      // Equlibrium luminosity = C_abs*B
      EquilibriumLum[_sz] = NumUtils::prod_bbodyCGS<double>(_w,*_it,_cabs);
      
      _tlo = 0.3*(*_it); 
      _thi = 3.0*(*_it); 

      if (_dostochastic[_sz]) {	
	
	// Since Enthalpy and CalTemp will change as we generate transition matrix,
	// regenerate at each size/cmp pair.
	_Enthalpy = GrainModel.Enthalpy(_sz,_cmp); 
	_CalTemp = GrainModel.CalTemp(_cmp);       	
	// Compute absorptive and radiative timescales to help shut off stochastic. 
	tau_abs = HeatUtils::getTauAbs(J,_cabs,_w,mpe); 
	Tu = HeatUtils::getTmean(_CalTemp,_Enthalpy,mpe); 
	tau_rad = HeatUtils::getTauRad(_cabs,_w,mpe,Tu); 

	if (tau_abs < 10*tau_rad) {   // Turn off stochastic heating for this and all subsequent sizes. 
	  transform(_dostochastic.begin()+_sz,_dostochastic.end(),_dostochastic.begin()+_sz,
		    bind1st(multiplies<bool>(),false));
	} else { // Compute StochasticLum, zero out Eq. 
	  status = StochasticHeating(_w,_cJprod,_cabs,_CalTemp,_Enthalpy,EAbs,TMin,TMax,*_it,_sz,StochasticLum); 
	  if (status != Flags::FSUCCESS) return status; // There's been a stochastic heating failure. 
	  EquilibriumLum[_sz].resize(_nw,0.0);  
	}  
      } else StochasticLum[_sz].resize(_nw,0.0);  // Zero stochastic so we have zero's when we don't do it.  

    }

    // Now we have E_eq(w,a) = C(w,a)*B(w,T(a)) and E_st(w,a) = C(w,a)*SUM(P(T(a))*B(w,T(a))
    // Need to compute integrals over size distribution. 
    _integrand.resize(_nsize); 
    _integrand1.resize(_nsize); 
    _norm=GrainModel.getNormalization(_cmp);
    
    _size = GrainModel.Size(_cmp);
    _sizeDist = GrainModel.getSizeDistribution(_cmp);
    _iteec = EmmittedEnergy[2*_cmp+1].begin(); // component equilibrium
    _itees = EmmittedEnergy[2*_cmp+2].begin(); // component stochastic
    _itet = EmmittedEnergy[0].begin(); 

    for (uint _wv=0;_wv<_nw;++_wv,++_iteec,++_itees,++_itet) { // wavelength loop
    
      // _integrand contains the equilibrium contribution to luminosity
      // _integrand1 contains the stochastic contribution to luminosity
      for (uint _sz=0;_sz<_nsize;++_sz) { // size loop

	if (!_dostochastic[_sz]) { // compute the equilibrium contribution 
	  _integrand[_sz] = _sizeDist[_sz]*EquilibriumLum[_sz][_wv];
	  _integrand1[_sz] = 0.0; 
	} else { // compute the stochastic contribution
	  _integrand[_sz] = 0.0; 
	  _integrand1[_sz] = _sizeDist[_sz]*StochasticLum[_sz][_wv]; 
	}
      }
      
      *_iteec = _norm*Constant::FPI*NumUtils::integrate<double>(_size,_integrand); 
      *_itees = _norm*Constant::FPI*NumUtils::integrate<double>(_size,_integrand1); 
      *_itet += (*_iteec + *_itees); 

    }

  }

  // Successfully completed ComputeDustEmission()
  return Flags::FSUCCESS; 

}  







