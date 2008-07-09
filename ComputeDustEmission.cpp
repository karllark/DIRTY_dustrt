#include "GrainModel.h"
#include "NumUtils.h"

void ComputeDustEmission (vector <float> & J, GrainModel & GrainModel, 
			  vector <vector<double> > & EmmittedEnergy, bool DoStochastic)
{

  // Local copy that can flip when this bin transitions. 
  bool _dostochastic = DoStochastic; 

  vector <float>::iterator _it,_it1,_it2,_itb,_ite;
  // Iterators for: 
  // _iteec == emmitted equilibrium for each component
  // _itet == total emmitted energy. 
  vector <double>::iterator _iteec,_itet; 
  vector <float> _w = GrainModel.getWave();
  vector <float> _E = GrainModel.getEScale();
  vector <float> _Enthalpy,_CalTemp; 
  vector <vector<float> > _TM; 
  vector <float> _cabs,_size,_sizeDist; 
  uint _nsize; 
  uint _nw = _w.size();
  
  if (J.size() != _nw) { cout << "FUCK!!! J/nwave don't match!" << endl; exit(8);}

  vector <float> _integrand; 

  // Do the equilibrium calculation. 
  vector <float> _t;
 
  float _tlo,_thi,_norm; 

  int _ncmp = GrainModel.getNComp(); 

  for (int _cmp=0;_cmp<_ncmp;++_cmp) {

    _nsize = GrainModel.nSize(_cmp);
    _cabs = GrainModel.CAbs(_cmp,(int)0); 
    _t.resize(_nsize); 
    // Compute first temperature to get the ball rolling.
    _it=_t.begin(); 
//     *_it=EqTemp(_w,J,_cabs);
//     ++_it; 
    _tlo = 0.001; 
    _thi = 25000.0; 

    // Compute Equilibrium temperature for each size of this component.
    for (uint _sz=0;_sz<_nsize;++_sz) {
      _cabs = GrainModel.CAbs(_cmp,int(_sz)); 
      *_it = EqTemp(_w,J,_cabs,((_tlo<0)?0.001:_tlo),((_thi>25000.0)?25000.0:_thi));
      ++_it; 
      _tlo = 0.3*(*(_it-1)); 
      _thi = 3.0*(*(_it-1)); 

      if (_dostochastic) { 
	//ComputTransitionMatrix(_E,_w,J,_cabs,_CalTemp,_Enthalpy,_TM); 
      }
    }
    // Have temperature, size distribution for each size 
    
    // Compute luminosity.
    _norm=GrainModel.getNormalization(_cmp);           // Size Dist norm in H^-1

    _size = GrainModel.Size(_cmp);                    // Size in cm
    _sizeDist = GrainModel.getSizeDistribution(_cmp); // Size dist in cm^-1

    _itb = _sizeDist.begin();
    _ite = _sizeDist.end();
   
    _iteec = EmmittedEnergy[2*_cmp+1].begin(); 
    _itet = EmmittedEnergy[0].begin(); 

    for (uint _wv=0;_wv<_nw;++_wv,++_iteec,++_itet) {
      _cabs = GrainModel.wCAbs(_cmp,_wv);             // _cabs at wave[] for all sizes in cm^2.
      _integrand = bbodyCGS(_w[_wv],_t);              // black body at wave[] for all sizes in erg/s/cm^2/cm^1/st
      _it1 = _cabs.begin(); 
      _it2 = _integrand.begin();
      for (_it=_itb;_it!=_ite;++_it,++_it1,++_it2) *_it2 *= ((*_it)*(*_it1));
      
      // now we have n(a)*Cabs(a)*B(T[a]) at wavelength _w[_wv]
      // The luminosity per unit wavelength is given by 
      // 4pi*norm*int[ da n(a)*Cabs(a)*B(T[a]) ]  
      // with units of erg/s/cm/H
      *_iteec = (double)Constant::FPI*_norm*NumUtils::integrate(_size,_integrand); 
      *_itet += *_iteec; 
      if (!finite(*_itet)) { cout << _wv <<"  "<< _w[_wv] << " BLEECH!!!" << endl; exit(8); }
    }
  }

}  
