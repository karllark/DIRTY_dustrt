#include "GrainModel.h"
#include "NumUtils.h"

void ComputeEmission (vector <float> & J, GrainModel & GrainModel)
{

  vector <float> _w = GrainModel.getWave();
  vector <float> _cabs,_size,_sizeDist,_norm,_bb; 
  int _nsize; 
  int _nw = _w.size(); 
  //  _wcabs.resize(_nw); 
  // _bb.resize(_nw);
  
  if (J.size() != _nw) { cout << "FUCK!!! J/nwave don't match!" << endl; exit(8);}

  vector <float> _integrand; 
  // Do the equilibrium calculation. 
  vector <float> _t;
  vector <float>::iterator _it,_it1,_it2,_itb,_ite; 
  float _tlo,_thi; 

  int _ncmp = GrainModel.getNComp(); 
  for (int _cmp=0;_cmp<_ncmp;++_cmp) {

    _size = GrainModel.Size(_cmp); 
    _sizeDist = GrainModel.getSizeDistribution(_cmp); 
    _nsize = GrainModel.nSize(_cmp);

    _integrand.resize(_nsize); 
    _cabs.resize(_nw);
    _cabs = GrainModel.CAbs(_cmp,(int)0); 
    _t.resize(_nsize); 
    _it=_t.begin(); 
    *_it=EqTemp(_w,J,_cabs);
    
    ++_it; 
    
    for (int _sz=1;_sz<_nsize;++_sz,++_it) {
      _cabs = GrainModel.CAbs(_cmp,_sz); 
      _tlo = 0.3*(*(_it-1)); 
      _thi = 4.0*(*(_it-1)); 
      *_it = EqTemp(_w,J,_cabs,((_tlo<0)?1.0:_tlo),((_thi>2500.0)?2500.0:_thi));

      cout << _size[_sz] << " " << *_it << endl; 
    }
    // Have temperature, size distribution for each size 
    
    // Compute luminosity.
    float norm=GrainModel.getNormalization(_cmp); 
    _cabs.resize(_nsize);
    _itb = _sizeDist.begin();
    _ite = _sizeDist.end();
    cout << endl;
    for (int _wv=0;_wv<_nw;++_wv) {
      _cabs = GrainModel.wCAbs(_cmp,_wv); // _cabs at wave[] for all sizes.
      _integrand = bbodyCGS(_w[_wv],_t);       // black body at wave[] for all sizes.  
      _it1 = _cabs.begin(); 
      _it2 = _integrand.begin();
      for (_it=_itb;_it!=_ite;++_it,++_it1,++_it2) *_it2 *= (*_it)*(*_it1);
      cout << _w[_wv] << " " << norm*NumUtils::integrate(_size,_integrand) << endl; 
    }
    
  }

  // Number of components
//   int _nComp = GrainModel.getNComp(); 
//   int _nSize;

//   float _tlo,_thi; 

//   vector <float> _thisCAbs; 
//   vector <float> _thisTemp;
//   vector <float> _thisSize; 
//   vector <float> _thisSizeDist; 
//   vector <float> _thisWave = GrainModel.getWave(); 
  
//   // 
//   vector <float>::iterator _itsz,_ittemp,_itbeg,_itend;  

//   cout << "Starting component analysis" << endl; 
//   // Compute the equilibrium temperature
//   for (int cmp=0;cmp<_nComp;++cmp) { // Component loop
    
//     _nSize = GrainModel.nSize(cmp); 
//     _thisTemp.resize(_nSize); 
//     _thisSize = GrainModel.Size(cmp);
//     _itbeg = _thisSize.begin(); _itend = _thisSize.end();    
//     _ittemp = _thisTemp.begin();
//     _itsz = _thisSize.begin();    
//     _thisCAbs = GrainModel.CAbs(cmp,*_itsz);    
//     *_ittemp = EqTemp(_thisWave,J,_thisCAbs);     
//     ++_ittemp; _itsz=_itbeg+1;

//     // Compute the equilibrium temperature of a grain
//     for (;_itsz!=_itend;++_itsz,++_ittemp) { // Size loop
//       _tlo = 0.3*(*(_ittemp-1)); 
//       _thi = 4.0*(*(_ittemp-1)); 
//       _thisCAbs = GrainModel.CAbs(cmp,*_itsz); 
//       *_ittemp = EqTemp(_thisWave,J,_thisCAbs,((_tlo<0)?1.0:_tlo),((_thi>2500.0)?2500.0:_thi));
//     }

//   }

}
