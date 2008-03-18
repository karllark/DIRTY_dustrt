#include "EqTemp.h"

// Given a wavelength in cm (wave) and the mean intensity J 
// (=1/(4pi)*int(I domega)) in erg/cm^2/s/cm and a grain with Q(a,wave) of Q, 
// return the equilibrium temperature of  the grain. NOTE: Q can be Q, C, sigma,
// etc....
//
// Restrictions: assumption is that J and Q are on the SAME WAVELENGTH scale. 
//               No checks other than a <vector>.size() check are made. 
//               Why? because it's not really possible to check...


double EqTemp(vector <double> wave, vector <double> J, vector <double> Q, double Tl, double Th)

{

  double Tk;                   // Equilibrium temperature of grain, returned value
  double LHS,RHS;              // Left and Right hand sides of equilibrium equation
  vector <double> integrand;   // Um, the integrands? 
  double Convergence=0.001;    // Fractional convergence in (RHS/LHS) - 0.1%
                               // Generally corresponsds to < 0.1K.  This is overkill 
                               // in the situations tested (one can get away with 1%),
                               // but the hit in speed is minimal. 
  double TestConvergence=1.0;  // Has our algorithm converged? 
  double Delta;                // Computed difference between LHS and RHS

  
  if (J.size() != Q.size()) { 
    cout << "Input energy density arrays and grain arrays do not match in size" << endl;
    cout << "IN EqTemp(), J.size()==" << J.size() << " Q.size()=="<< Q.size() << endl; 
    return -1; 
  } 

  // make the integrand equal to Q*J
  copy(J.begin(),J.end(),back_inserter(integrand)); 
  transform(integrand.begin(),integrand.end(),Q.begin(),integrand.begin(),multiplies<double>());

  // compute the left hand side of the equilibrium equation - this doesn't change anymore!
  LHS=NumUtils::integrate(wave,integrand); 
  
  // Clear integrand as we will re-use it to compute RHS. 
  //   bbody uses a push_back to populate the vector, so we want it empty.
  integrand.clear(); 

  // Terminate calculation if we've converged 
  while (TestConvergence > Convergence) { 
    
    Tk = (Tl+Th)/2.0; 
    
    // Get black body at Tk, compute RHS integrand. 
    integrand = NumUtils::bbodyCGS(wave,Tk); 
    transform(integrand.begin(),integrand.end(),Q.begin(),integrand.begin(),multiplies<double>()); 
    
    RHS = NumUtils::integrate(wave,integrand);
    Delta = (LHS-RHS)/LHS; 

    // If we've converged, Tk is the correct temperature.
    if (fabs(Delta) < Convergence) 
      TestConvergence=Convergence; 
    else { // Reset bounds, iterate again. 
      if (Delta < 0.0) 
	Th = Tk; 
      else 
	Tl = Tk;
    }
      
    // .clear() integrand as bbody uses a push_back()
    integrand.clear();

  }

  return Tk;

}

float EqTemp(vector <float> wave, vector <float> J, vector <float> Q, float Tl, float Th)

{

  float Tk;                   // Equilibrium temperature of grain, returned value
  float LHS,RHS;              // Left and Right hand sides of equilibrium equation
  vector <float> integrand;   // Um, the integrands? 
  float Convergence=0.001;    // Fractional convergence in (RHS/LHS) - 0.1%
                              // Generally corresponsds to < 0.1K.  This is overkill 
                              // in the situations tested (one can get away with 1%),
                              // but the hit in speed is minimal. 
  float TestConvergence=1.0;  // Has our algorithm converged? 
  float Delta;                // Computed difference between LHS and RHS

  if (J.size() != Q.size()) { 
    cout << "Input energy density arrays and grain arrays do not match in size" << endl;
    cout << "IN EqTemp(), J.size()==" << J.size() << " Q.size()=="<< Q.size() << endl; 
    return -1; 
  } 

  // make the integrand equal to  Q*J
  copy(J.begin(),J.end(),back_inserter(integrand)); 
  transform(integrand.begin(),integrand.end(),Q.begin(),integrand.begin(),multiplies<float>());

  // compute the left hand side of the equilibrium equation - this doesn't change anymore!
  LHS=NumUtils::integrate(wave,integrand); 

  // Clear integrand as we will re-use it to compute RHS. 
  //   bbody uses a push_back to populate the vector, so we want it empty.
  integrand.clear(); 

  // Terminate calculation if we've converged
  while (TestConvergence > Convergence) { 
    
    Tk = (Tl+Th)/2.0; 

    // Get black body at Tk, compute RHS integrand. 
    integrand = NumUtils::bbodyCGS(wave,Tk); 
    transform(integrand.begin(),integrand.end(),Q.begin(),integrand.begin(),multiplies<float>()); 
    
    RHS = NumUtils::integrate(wave,integrand);
    Delta = (LHS-RHS)/LHS; 

    // If we've converged, Tk is the correct temperature.
    if (fabs(Delta) < Convergence) 
      TestConvergence=Convergence; 
    else {  // Reset bounds, iterate again. 
      if (Delta < 0.0) 
	Th = Tk; 
      else 
	Tl = Tk;
    }

    // .clear() integrand as bbody uses a push_back()
    integrand.clear();
 
  }

  return Tk;

}
