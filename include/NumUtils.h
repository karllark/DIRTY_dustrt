#ifndef _NUMUTILS_H
#define _NUMUTILS_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector> 

#include "Constants.h"

using namespace std;

namespace NumUtils { // Define a namespace to avoid confusion with other 
                     // classes/libraries. 

  // Note that will you can call these with any type, the meaning is not 
  // intuitively obvious for some.  They were constructed mainly for 
  // float/double types.  

  // ****************************************************************************
  // Template prototypes. 
  // **********************************************************************
  template<class T>
    class sqSum {
  public:
    T operator()(T sum, T val) { return sum+val*val; }
  };
  template<class T> 
    class Square { 
  public: 
    T operator()(T val) { return val*val; }
  }; 
  template<class T> 
    class inverse { 
  public: 
    T operator()(T val) { return 1/val; }
  }; 
  template <class T> 
    class exp { 
  public: 
    T operator()(T val) { return std::exp(val); }
  }; 
  template <class T>  
    class ln { 
  public: 
    T operator()(T val) { return std::log(val); }
  };
  template <class T> 
    class isNan { 
  public: 
    bool operator()(T val) { return val != val; } 
  };
  template <class T> 
    class invmult { 
  public: 
    T operator()(T val1, T val2) { return val1/val2; }
  }; 

  template <class T> void swap( T & a, T & b) { 
    T tmp=a; 
    a=b;
    b=tmp;
  }


  template <typename T> inline T NaN() { return strtod("NaN",NULL); }

  template <typename T> inline T binrepl(T & val1, T & val2)
    { return (val1 <= val2) ? Constant::PLANCKLIGHT:0.0; }
  template <typename T> inline T line(T x1, T x2, T y1, T y2, T x)
    { return ( (y2-y1)/(x2-x1)*(x-x1)+y1 ); }
  template <typename T> inline int index(T val, vector <T>&vect) {
    typename vector <T>::iterator idx;  
    idx=std::lower_bound(vect.begin(),vect.end(),val); 
    return distance(vect.begin(),idx); 
  }
  template <typename T> inline int rindex(T val, vector <T>&vect) {
    typename vector <T>::iterator idx;  
    idx=std::lower_bound(vect.begin(),vect.end(),val,greater<T>()); 
    return distance(vect.begin(),idx); 
  }
  // Some extrapolators. 
  template <typename T> inline T lextra(T & val1, T & val2, T & val3) 
    { return val1*val2/val3; }
  template <typename T> inline T cextra(T & val1, T & val2, T & val3) 
    { T ratio=val2/val3; return val1*ratio*ratio*ratio; } 
  template <typename T> inline T qextra(T & val1, T & val2, T & val3) 
    { T ratio=val2/val3; return val1*ratio*ratio*ratio*ratio; } 
  template <typename T> inline T frootextra(T & val1, T & val2, T & val3) 
    { T ratio=val2/val3; return val1*sqrt(sqrt(ratio)); } 
  
  // Functions with 'source'
  template <typename T> int maxID(vector <T>& vect);
  template <typename T> vector <int> sortID(vector <T>&vect);
  template <typename T> vector <int> sortID_src(vector <T>&vect); 
  template <typename T> void sortIndexOnVector(vector <T> & vect, vector <int> & index, int beg, int end);
  template <typename T, typename T1, typename T2> T integrate (vector <T1>& x, vector <T2>& y);
  template <typename T> vector <T> interpol (vector <T>& v, vector <T>& x, 
					     vector <T>& u, int LoEx=2, int HiEx=2);
  template <typename T> void interpolr (vector <T>& v, vector <T>& x, 
					vector <T>& u, int LoEx=2, int HiEx=2);
  template <typename T, typename T1, typename T2> vector <T> bbodyCGS (T1 wave, vector <T2>& Temperature);
  template <typename T, typename T1, typename T2> vector <T> bbodyCGS (vector <T1>& wave, T2 Temperature);
  template <class T, class T1, class T2, class T3> vector <T> prod_bbodyCGS (T1 wave, vector <T2>& Temperature, vector <T3>& product);
  template <class T, class T1, class T2, class T3> vector <T> prod_bbodyCGS (vector <T1>& wave, T2 Temperature, vector <T3>& product);
  template <class T, class T1, class T2, class T3> vector <T> add_bbodyCGS (T1 wave, vector <T2>& Temperature, vector <T3>& add);
  template <class T, class T1, class T2, class T3> vector <T> add_bbodyCGS (vector <T1>& wave, T2 Temperature, vector <T3>& add);
  template <typename T> string vtos(const T & value);
  template <typename T> class Matrix;
  template <typename T> class Cube; 
  template <typename T> Matrix <T> mInvert(Matrix <T> A); 
  template <typename T> vector <T> poly_fit(vector <T> x, vector <T> y, 
					    int ndegree, vector <T>& yfit, vector <T>& sigma);
  
  
  // ****************************************************************************
  // Function 'source' code. 
  // ****************************************************************************

 // Simple "one-pass" algorithm for mean/variance.
  template <typename T> T stats(vector <T> & vect, T & var)
    {
      T M0=vect[0];
      T Q0=0;
      T M1,Q1;
      for (int i=1;i<vect.size();i++) {
	M1 = M0 + (vect[i]-M0)/(T)(i+1);
	Q1 = Q0 + ((T)i)*pow((vect[i]-M0),2)/(T)(i+1);
	M0 = M1;
	Q0 = Q1;
      }
      var = sqrt(Q1/(T)(vect.size()-1));

      return M1;
    }
  // ****************************************************************************

  template <typename T> T sigma_clip(vector <T> & vect, T & hi_clip, T & lo_clip)
    {
      vect.erase(remove_if(vect.begin(),vect.end(),bind2nd(greater<T>(),hi_clip)),vect.end());
      vect.erase(remove_if(vect.begin(),vect.end(),bind1st(greater<T>(),lo_clip)),vect.end()); 
    }

  // Get max element (index) of vector.
  template <typename T> int maxID(vector <T>& vect)
    {
      T theMax=-1e256;
      int theMaxID; 
      typename vector <T>::iterator _it,_itb,_ite; 
      _itb=vect.begin(); _ite=vect.end();
      for (_it=_itb;_it!=_ite;++_it) {
	if (*_it > theMax) { 
	   theMax = *_it; 
	   theMaxID = distance(_itb,_it); 
        }
      }
      //for (int i=0;i<vect.size();i++) {
//	if (vect[i] > theMax) { 
//	  theMax=vect[i];
//	  theMaxID=i; 
//	}
 //     }
      return theMaxID; 
    }
  // ****************************************************************************
  // ****************************************************************************
  // Returns vector of indices that will sort provided vector
  template <typename T> vector <int> sortID(vector <T>&vect)
    {
      int n=vect.size();
      int i;
      T a; 
      vector <int> retvect;
      retvect.push_back(0); 
      for (int j=1;j<n;j++) {
	retvect.push_back(j); 
	a=vect[j];
	i=j; 
	while (i > 0 && vect[retvect[i-1]] > a) {
	  //retvect[i+1] = retvect[i];
	  swap(retvect[i],retvect[i-1]); 
	  i--; 
	}
      }
      
      return retvect;
    }
  // ****************************************************************************

  template <typename T> void sortIndexOnVector(vector <T> & vect, vector <int> & index, int beg, int end)
    {
      int i; 
      T a; 
      for (int j=beg+1;j<end;j++) {
	//cout << "assigning " << vect[index[j]] << " to a" << endl; 
	a = vect[index[j]];
	i=j; 
	while (i > beg && vect[index[i-1]] > a) { 
	  //cout << "Satisfied that " << vect[index[i-1]] << " is > " << a << endl; 
	  //cout << "Swapping " << vect[index[i]] << " " << vect[index[i-1]] << endl; 
	  swap(index[i],index[i-1]); 
	  i--; 
	}
      }
    }

  // ****************************************************************************
  // returns vector of indices AND Re-orders vector. 
  template <typename T> vector <int> sortID_src(vector <T>&vect)
    {
      
      int n=vect.size();
      int i;
      T a; 
      vector <int> retvect;
      retvect.push_back(0); 
      //for (int j=0;j<n;j++) retvect.push_back(j);
      for (int j=1;j<n;j++) {
	retvect.push_back(j); 
	a=vect[j];
	i=j; 
	while (i > 0 && vect[i-1] > a) {
	  swap(vect[i],vect[i-1]);
	  swap(retvect[i],retvect[i-1]); 
	  i--; 
	}
      }
      
      return retvect;
    }
  // ****************************************************************************

  // ****************************************************************************
  // Simple trapazoidal integration
  template <typename T, typename T1, typename T2> T integrate (vector <T1>& x, vector <T2>& y)
    {
      T answer=0.0;
      typename vector <T2>::iterator iy;
      typename vector <T1>::iterator ix; 
      if (x.size() > 1) {
	iy = y.begin()+1;
	for (ix=x.begin()+1;ix!=x.end();ix++) {
	  answer += 0.5*(*(ix) - *(ix-1))*(*(iy) + *(iy-1));
	  iy++;
	}
      } else { answer=y[0];}
      return answer;
    }
  // ****************************************************************************


 

  // ****************************************************************************
  // Interpolate - should be fast enough for small vectors.
  template <typename T> vector <T> interpol (vector <T>& v, vector <T>& x, 
					     vector <T>& u, int LoEx, int HiEx)
    {

      if (v.size() != x.size()) { 
	cout << "Abscissa and ordinate lengths do not match in interpol" << endl; 
	throw "size mismatch"; 
      } 

      vector <T> r(u.size()); 

      typename vector <T>::iterator iu,iul,iuh,ir,ix,iv; 
      ix = x.begin();
      iv = v.begin(); 
      iul = u.begin(); 
      ir = r.begin(); 

      while ( *iul < *ix ) { // Extrapolate to the left. 
	switch (LoEx) { 
	case -1 :  // set to 0
	  *ir = 0.0; 
	  break; 
	case 0 : // Constant extrapolation
	  *ir = *iv; 
	  break; 
	default : // Power law. 
	  *ir = (*iv)*pow((*iul)/(*ix), LoEx);
	  break; 
	}
	iul++; 
	ir++; 
      }

      iuh = u.end()-1; 
      ir = r.end()-1; 
      ix = x.end()-1; 
      iv = v.end()-1; 
      while ( *iuh > *ix ) { // Extrapolate to the right. 
	switch (HiEx) { 
	case -1 : // set to 0
	  *ir = 0; 
	  break; 
	case 0 : // constant
	  *ir = *iv; 
	  break; 
	default : // Power law. 
	  //cout << *iv << " " << *ix << " " << *iuh << " " << HiEx << endl; 
	  *ir = (*iv)*pow( (*iuh)/(*ix), HiEx); 
	  break; 
	}
	iuh--; 
	ir--;
      }

      ix = x.begin(); 
      iv = v.begin();
      ir = r.begin()+distance(u.begin(),iul); 

      for (iu=iul;iu<=iuh;iu++) {	
	while ( *ix < *iu ) { ix++; iv++; }
	*ir = ((*iu) - (*(ix)))*((*(iv)) - (*(iv-1)))/((*(ix))- (*(ix-1))) + (*iv); 
	ir++;
      }
      
      return r; 
    }
  // ****************************************************************************

  // *************************************************************************
  // Interpolate with replacement of input vector - ie. v comes out with 
  // interpolated vector rather than the input vector. 
  template <typename T> void interpolr (vector <T>& v, vector <T>& x, vector <T>& u, 
					int LoEx, int HiEx) 
    {

      if (v.size() != x.size()) { 
	cout << "Abscissa and ordinate lengths do not match in interpol" << endl; 
	throw "size mismatch"; 
      }  

      // Copy input vector into a temporary vector and use v as the new vector. 
      vector <T> tmpv;
      copy(v.begin(),v.end(),back_inserter(tmpv));
      v.erase(v.begin(),v.end());
      v.resize(u.size());
 
      //T slope,intercept; 

      typename vector <T>::iterator iu,iul,iuh,ir,ix,iv; 
      // Iterators point to begining of vectors.
      ix = x.begin();
      iv = tmpv.begin(); 
      iul = u.begin(); 
      ir = v.begin();
 
      while ( *iul < *ix ) { // Extrapolate to the left. 
	//cout << "EXTRAPOLATING LOW" << endl; 
	switch (LoEx) { 
	case -1 :  // set to 0
	  *ir = 0.0; 
	  break; 
	case 0 : // Constant extrapolation
	  *ir = *iv; 
	  break; 
	default : // Power law. 
	  *ir = (*iv)*pow((*iul)/(*ix), LoEx);
	  break; 
	}
	iul++; 
	ir++; 
      }

      // Iterators point to last element of vectors
      iuh = u.end()-1; 
      ir = v.end()-1; 
      ix = x.end()-1; 
      iv = tmpv.end()-1; 

      while ( *iuh > *ix ) { // Extrapolate to the right. 
	switch (HiEx) { 
	case -1 : // set to 0
	  *ir = 0; 
	  break; 
	case 0 : // constant
	  *ir = *iv; 
	  break; 
	default : // Power law. 
	  *ir = (*iv)*pow( (*ix)/(*iuh), HiEx); 
	  break; 
	}
	iuh--; 
	ir--;
      }

      ix = x.begin(); 
      iv = tmpv.begin();
      ir = v.begin()+distance(u.begin(),iul); 
  
      for (iu=iul;iu<=iuh;iu++) {      
	while ( *ix < *iu ) { ix++; iv++; }
	*ir = ((*iu) - (*(ix)))*((*(iv)) - (*(iv-1)))/((*(ix))- (*(ix-1))) + (*iv); 
	ir++;
      }
 
    }
  // ****************************************************************************

    // ****************************************************************************
  // BB inputs and outputs all in cgs - cm and erg/s/cm^2/st/cm
  // Returns BB as f(T) for a single wavelength
  template <typename T, typename T1, typename T2> vector <T> bbodyCGS (T1 wave, vector <T2>& Temperature)
    {   
      // Wave is in CGS (cm), return in CGS: erg/cm^2/s/cm/st
      T c1 = 2.0*Constant::PLANCK*pow(Constant::LIGHT,2);
      T c2 = Constant::PLANCKLIGHT/(Constant::BOLTZMAN*static_cast<T>(wave));
/*       cout << Constant::PLANCKLIGHT << " " << Constant::BOLTZMAN << " " << wave << endl;  */
/*       cout << Constant::BOLTZMAN*wave << endl;  */
/*       cout << Constant::PLANCKLIGHT/(Constant::BOLTZMAN*wave) << endl; */
/*       cout << "C2 in bb: " << c2 << endl;  */
      typename vector <T2>::iterator _itbeg,_itend,_it;
      //cout << "comp bb: " <<  c1 << " " << c2 << endl; 
      vector <T> theBB;
      _itbeg = Temperature.begin();
      _itend = Temperature.end();
      int i=0; 
      for (_it=_itbeg;_it!=_itend;++_it,++i) {
	//cout << "INBB 1: " << static_cast<T1>(wave) << " " << c1/pow((static_cast<T1>(wave)),5) << endl; 
	//cout << (std::exp((c2/static_cast<T1>((*_it))))) << endl; 
	theBB.push_back(c1/pow((wave),5)*(1.0/(std::exp(c2/(*_it)) - 1.0)));
	//cout << "IN BB: " << wave << " " << *_it << " "<<  theBB[i] << endl; 
      }
      return theBB;
    }
  // ****************************************************************************

  // ****************************************************************************
  // BB inputs and outputs all in cgs - cm and erg/s/cm^2/st/cm
  // Returns BB multiplied by another vector as f(T) for a single wavelength
  template <class T, class T1, class T2, class T3> vector <T> prod_bbodyCGS (T1 wave, vector <T2>& Temperature, vector <T3>& product)
    {   
      // Wave is in CGS (cm), return in CGS: erg/cm^2/s/cm/st
      T c1 = 2.0*Constant::PLANCK*pow(Constant::LIGHT,2);
      T c2 = Constant::PLANCKLIGHT/(Constant::BOLTZMAN*static_cast<T>(wave));
      typename vector <T2>::iterator _itbeg,_itend,_it;
      typename vector <T3>::iterator _itp;

      vector <T> theBB;
      _itbeg = Temperature.begin();
      _itend = Temperature.end();
      _itp = product.begin(); 
      for (_it=_itbeg;_it!=_itend;++_it,++_itp) {
	theBB.push_back((*_itp)*c1/pow((wave),5)*(1.0/(std::exp(c2/(*_it)) - 1.0)));	
      }
      return theBB;
    }
  // ****************************************************************************
  // ****************************************************************************
  // BB inputs and outputs all in cgs - cm and erg/s/cm^2/st/cm
  // Returns BB summed with another vector as f(T) for a single wavelength
  template <class T, class T1, class T2, class T3> vector <T> add_bbodyCGS (T1 wave, vector <T2>& Temperature, vector <T3>& add)
    {   
      // Wave is in CGS (cm), return in CGS: erg/cm^2/s/cm/st
      T c1 = 2.0*Constant::PLANCK*pow(Constant::LIGHT,2);
      T c2 = Constant::PLANCKLIGHT/(Constant::BOLTZMAN*static_cast<T>(wave));
      typename vector <T2>::iterator _itbeg,_itend,_it;
      typename vector <T3>::iterator _ita;

      vector <T> theBB;
      _itbeg = Temperature.begin();
      _itend = Temperature.end();
      _ita = add.begin(); 
      for (_it=_itbeg;_it!=_itend;++_it,++_ita) {
	theBB.push_back( (c1/pow((wave),5)*(1.0/(std::exp(c2/(*_it)) - 1.0))) + *_ita);	
      }
      return theBB;
    }
  // ****************************************************************************

  // ****************************************************************************
  // BB inputs and outputs all in cgs - cm and erg/s/cm^2/st/cm
  // Returns BB as f(wave) for a single T
  template <typename T, typename T1, typename T2> vector <T> bbodyCGS (vector <T1>& wave, T2 Temperature)
    {   
      // Wave is in CGS (cm), return in CGS: erg/cm^2/s/cm/st
      T c1 = 2.0*Constant::PLANCK*pow(Constant::LIGHT,2);
      T c2 = Constant::PLANCKLIGHT/(Constant::BOLTZMAN*Temperature);
      class vector <T1>::iterator _itbeg,_itend,_it;

      vector <T> theBB;
      _itbeg = wave.begin();
      _itend = wave.end();
      for (_it=_itbeg;_it!=_itend;++_it) 
	theBB.push_back(c1/pow((*_it),5)*(1.0/(std::exp(c2/(*_it)) - 1.0)));
      
      return theBB;
    }
  // ****************************************************************************

  // ****************************************************************************
  // BB inputs and outputs all in cgs - cm and erg/s/cm^2/st/cm
  // Returns BB multiplied by another vector as f(wave) for a single T
  template <class T, class T1, class T2, class T3> vector <T> prod_bbodyCGS (vector <T1>& wave, T2 Temperature, vector <T3>& product)
    {   
      // Wave is in CGS (cm), return in CGS: erg/cm^2/s/cm/st
      T c1 = 2.0*Constant::PLANCK*pow(Constant::LIGHT,2);
      T c2 = Constant::PLANCKLIGHT/(Constant::BOLTZMAN*Temperature);
      class vector <T1>::iterator _itbeg,_itend,_it;
      class vector <T3>::iterator _itp; 

      vector <T> theBB;
      _itbeg = wave.begin();
      _itend = wave.end();
      _itp = product.begin();
      for (_it=_itbeg;_it!=_itend;++_it,++_itp) 
	theBB.push_back((*_itp)*c1/pow((*_it),5)*(1.0/(std::exp(c2/(*_it)) - 1.0)));
      
      return theBB;
    }
  // ****************************************************************************
  // ****************************************************************************
  // BB inputs and outputs all in cgs - cm and erg/s/cm^2/st/cm
  // Returns BB multiplied by another vector as f(wave) for a single T
  template <class T, class T1, class T2, class T3> vector <T> add_bbodyCGS (vector <T1>& wave, T2 Temperature, vector <T3>& add)
    {   
      // Wave is in CGS (cm), return in CGS: erg/cm^2/s/cm/st
      T c1 = 2.0*Constant::PLANCK*pow(Constant::LIGHT,2);
      T c2 = Constant::PLANCKLIGHT/(Constant::BOLTZMAN*Temperature);
      class vector <T1>::iterator _itbeg,_itend,_it;
      class vector <T3>::iterator _ita; 

      vector <T> theBB;
      _itbeg = wave.begin();
      _itend = wave.end();
      _ita = add.begin();
      for (_it=_itbeg;_it!=_itend;++_it,++_ita) 
	theBB.push_back( (c1/pow((*_it),5)*(1.0/(std::exp(c2/(*_it)) - 1.0))) + *_ita);
      
      return theBB;
    }
  // ****************************************************************************

  // ****************************************************************************
  // Convert value of type T to a string - useful for constructing filenames
  // out of numbers. 
  template <typename T> string vtos(const T & value)
    {
      stringstream ss;
      ss << value;
      return ss.str();
    }
  // ****************************************************************************

  // ****************************************************************************
  // Do a polynomial fit to a set of data. Simple re-write of IDL poly_fit.pro
  // TODO: add HAVEMEASUREERRORS.
  template <typename T> vector <T> poly_fit(vector <T> x, 
					    vector <T> y, 
					    int ndegree, 
					    vector <T>& yfit,
					    vector <T>& sigma)
    {

      int m=ndegree+1;
      int n=x.size();
      //cout << m << " " << ndegree << endl; 
      vector <T> param(m);
      NumUtils::Matrix <T> covar(m,m);
      
      vector <T> b(m); 
      vector <T> z(n,1.0);
      T sdev=1.0,sdev2=1.0;
      vector <T> wy=y; 
      //cout << "in polyfit 1: " << param[0] << " " << param[1] << endl; 
      if (yfit.size() != n) { 
	cout << "Size of fit return not equal input vector sizes; poly_fit()" << endl; 
	exit(8);
      }
      if (sigma.size() != m) { 
	cout << "Size of parameter error array not equal to order+1; poly_fit()" << endl; 
	exit(8); 
      }

      T sum; 
      covar(0,0)=(T)n; 
      b[0] = accumulate(wy.begin(),wy.end(),0.0); 
      for (long p=1;p<=2*ndegree;p++) { 
	// z=z*x
	transform(z.begin(),z.end(),x.begin(),z.begin(),multiplies<T>());
	// b[p] = wy*z
	if (p < m) b[p]=inner_product(wy.begin(),wy.end(),z.begin(),0.0);
	// sum = total(z)
	T sum=accumulate(z.begin(),z.end(),0.0);
	for (int j=( (0>(p-ndegree)) ? 0:(p-ndegree) );j<=( (ndegree<p) ? ndegree : p);j++) covar(j,p-j)=sum;
      }
      covar = mInvert(covar);
      for (int j=0;j<m;j++) 
	for (int i=0;i<m;i++) param[j] += b[i]*covar(i,j);
      //vector <double> yfit(n,param[ndegree]);
      for (int i=0;i<n;++i) yfit[i]=param[ndegree]; 
      //for (int i=0;i<n;i++) yfit.push_back(param[ndegree]); 
      for (int k=ndegree-1;k>=0;k--) {
	transform(yfit.begin(),yfit.end(),x.begin(),yfit.begin(),multiplies<T>());
	transform(yfit.begin(),yfit.end(),yfit.begin(), 
		  bind2nd(plus<T>(),param[k]));
      }
      // Sigma == diagonal 
      for (int i=0;i<m;i++) sigma[i] = (sqrt(abs(covar(i,i))));
      // diff = yfit - y; chisq=SUM[(diff)^2]
      vector <T> diff;
      transform(yfit.begin(),yfit.end(),y.begin(),back_inserter(diff),minus<T>());
      T chisq=accumulate(diff.begin(),diff.end(),0.0,NumUtils::sqSum<T>());
      T var = (n > m) ? chisq/T(n-m) : 0.0; 
      transform(sigma.begin(),sigma.end(),sigma.begin(),
		bind2nd(plus<T>(),sqrt(chisq/(T)(n-m)))); 
      T yerr = sqrt(var);

      //cout << "in polyfit 2: " << param[0] << " " << param[1] << endl; 
      return param; 
    }
  // ****************************************************************************
  
  // ****************************************************************************
  // Compute inverse of a matrix; stolen from the web... Forgot where...
  template <typename T> NumUtils::Matrix <T> mInvert(NumUtils::Matrix <T> A) 
    {
      int N = A.nRow(); 
      NumUtils::Matrix <T> Ainv=A; 
      NumUtils::Matrix <T> b(N,N,0.0);
      vector <T> scale(N);
      vector <int> index(N+1);     
      for (int i=0;i<N;i++) { 
	b(i,i)=1.0;
	index[i]=i;
	T scalemax=0.0;
	for (int j=0;j<N;j++) scalemax = (scalemax > fabs(A(i,j))) ? scalemax : fabs(A(i,j));
	scale[i] = scalemax; 
      }
      int signDet=1;
      for (int k=0;k<N-1;k++) { 
	T ratiomax=0.0;
	int jPivot=k;
	for (int i=k;i<N;i++) { 
	  T ratio=fabs(A(index[i],k))/scale[index[i]]; 
	  if (ratio > ratiomax) { 
	    jPivot=i; 
	    ratiomax=ratio; 
	  }
	}
	int indexJ=index[k]; 
	if (jPivot != k) { 
	  indexJ = index[jPivot];
	  index[jPivot] = index[k]; 
	  index[k] = indexJ;
	  signDet *= -1; 
	}
	for (int i=k+1;i<N;i++) {
	  T coeff=A(index[i],k)/A(indexJ,k);
	  for (int j=k+1;j<N;j++) A(index[i],j) -= coeff*A(indexJ,j);
	  A(index[i],k)=coeff; 
	  for (int j=0;j<N;j++) b(index[i],j) -= A(index[i],k)*b(indexJ,j);
	}
      }
      for (int k=0;k<N;k++) { 
	Ainv(N-1,k) = b(index[N-1],k)/A(index[N-1],N-1);
	for (int i=N-2;i>=0;i--) { 
	  T sum=b(index[i],k);
	  for (int j=i+1;j<N;j++) sum -= A(index[i],j)*Ainv(j,k);
	  Ainv(i,k) = sum/A(index[i],i); 
	}
      }
      return Ainv; 
    }
  // ****************************************************************************



  // ****************************************************************************
  //  Separating out header and source doesn't work so well...
  template <typename T> class  Matrix : public std::vector<T> {
    
  public:
    // Constructors/destructors.
    Matrix() : std::vector<T>() {}
    Matrix(int n1, int n2, const T& ival) : std::vector<T>(n1*n2,ival), _n1 (n1), _n2 (n2) {}
    explicit Matrix(int n1, int n2) : std::vector<T>(n1*n2), _n1 (n1), _n2 (n2) {}
    ~Matrix() {}
    // Return reference to correct element in vector; All vector type assignments
    //   should work.
    T& operator() (int n1_id, int n2_id);
    int nRow ();
    int nCol ();
    void MSize(int n1_id, int n2_id); 
    
  private:
    // Keep track of how many elements we're allowed to have in each dimension..
    int _n1,_n2;
  };
  
  //  Overload of operator()
  template <typename T> inline T& Matrix<T>::operator() (int n1_id, int n2_id)
    {
      if (n1_id < 0 || n1_id >= _n1) {
	cout << "out_of_range Matrix::operator(), element 1" << endl;
	std::string ExceptionObject = "out_of_range Matrix::operator(), element 1";
	throw std::out_of_range(ExceptionObject);
      }
      if (n2_id < 0 || n2_id >= _n2) {
	cout << "out_of_range Matrix::operator(), element 2" << endl;
	std::string ExceptionObject = "out_of_range Matrix::operator(), element 2";
	throw std::out_of_range(ExceptionObject);
      }
      return *( this->begin()+ n2_id*_n1 + n1_id );
    }
  
  // Size the Matrix after instatiation.
  template <typename T> void Matrix<T>::MSize(int n1_id, int n2_id)
    {
      _n1 = n1_id;
      _n2 = n2_id;
      Matrix::clear();
      Matrix::resize(n1_id*n2_id);
    }
  
  template <typename T> int Matrix<T>::nRow() { return _n1; }
  template <typename T> int Matrix<T>::nCol() { return _n2; }
  // ****************************************************************************

  // ****************************************************************************
  // Cube class.
  template <typename T> class Cube : public std::vector<T> {
    
  public:
    // Constructors/destructors.
    Cube() : std::vector<T>() {} 
    Cube(int n1, int n2, int n3, const T& ival) : std::vector<T>(n1*n2*n3,ival), _n1 (n1),_n2 (n2), _n3 (n3) {}
    explicit Cube(int n1, int n2, int n3) : std::vector<T>(n1*n2*n3), _n1 (n1), _n2 (n2), _n3 (n3) {}
    ~Cube() {}
    // Return reference to correct element in vector; All vector type assignments
    //   should work.
    T& operator() (int n1_id, int n2_id, int n3_id);
    void CSize (int n1_id, int n2_id, int n3_id);
    
    int nRow ();
    int nCol ();
    int n3rd ();
    
  private: 
    // Keep track of how many elements we're allowed to have in each dimension. 
    int _n1,_n2,_n3;
    
  };
  
  // Overload of operator() 
  template <typename T> inline T& Cube<T>::operator() (int n1_id, int n2_id, int n3_id)
    {
      if (n1_id < 0 || n1_id >= _n1) { 
	cout << "out_of_range Cube::operator(), element 1 " << n1_id << "," << n2_id << "," << n3_id << " " << _n1 << 
endl;
	std::string ExceptionObject = "out_of_range Cube::operator(), element 1"; 
	throw std::out_of_range(ExceptionObject);
      }
      if (n2_id < 0 || n2_id >= _n2) { 
	cout << "out_of_range Cube::operator(), element 2 " << n1_id << "," << n2_id << "," << n3_id << " " << _n2 << 
endl;
	std::string ExceptionObject = "out_of_range Cube::operator(), element 2"; 
	throw std::out_of_range(ExceptionObject);
      }
      if (n3_id < 0 || n3_id >= _n3) {
	cout << "out_of_range Cube::operator(), element 3 " << n1_id << "," << n2_id << "," << n3_id << " " << _n3 << 
endl; 
	std::string ExceptionObject = "out_of_range Cube::operator(), element 3"; 
	throw std::out_of_range(ExceptionObject);
      }
      return *( this->begin() + n3_id*_n1*_n2 + n2_id*_n1 +  n1_id);
      
    }

  // Size the Cube after instatiation. 
  template <typename T> void Cube<T>::CSize(int n1_id, int n2_id, int n3_id)
    {
      _n1 = n1_id;
      _n2 = n2_id; 
      _n3 = n3_id;
      Cube::clear();
      Cube::resize(n1_id*n2_id*n3_id);
    }
  
  template <typename T> int Cube<T>::nRow() { return _n1; }
  template <typename T> int Cube<T>::nCol() { return _n2; }
  template <typename T> int Cube<T>::n3rd() { return _n3; }
  
  // ****************************************************************************

}

#endif
