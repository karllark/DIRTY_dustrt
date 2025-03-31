// Jun 3-4 2009: Modified to be DirtyFailure aware, eg. Teq is return in parameter list and return
// value
//               is an int holding failure status from DirtyFlags.
#ifndef _EQTMP_H
#define _EQTMP_H

#include <iostream>
#include <vector>

#include "DirtyFlags.h"
#include "NumUtils.h"

using namespace std;

template <class T1, class T2, class T3>
int
EqTemp (vector<T2> &wave, vector<T2> &J, vector<T2> &Q, T1 &EAbs, T1 &_Tk, T3 Tlo = 1,
        T3 Thi = 2500)
{
  // Local variables
  // T1 _Tk;                                  // Equilibrium temperature of grain, returned value
  T1 _LHS, _RHS;                          // Left and Right hand sides of equilibrium equation
  vector<T1> _integrand;                  // Um, the integrands?
  T1 _Convergence = 0.001;                // Fractional convergence in (RHS/LHS) - 0.1%
                                          // Generally corresponsds to < 0.1K.  This is overkill
                                          // in the situations tested (one can get away with 1%),
  T1 _TestConvergence = 1.0;              // Has our algorithm converged?
  T1 _Delta;                              // Computed difference between LHS and RHS
  typename vector<T2>::iterator _iJ, _iQ; // Iterators for input vectors
  typename vector<T2>::iterator _ib, _ie, _it; // Generic iterators

  T3 lolim = Tlo;
  T3 hilim = Thi;

  if (J.size () != Q.size ())
    {
      cout << "Input energy density arrays and grain arrays do not match in size" << endl;
      cout << "IN EqTemp(), J.size()==" << J.size () << " Q.size()==" << Q.size () << endl;
      return static_cast<T1> (-1);
    }

  _integrand.resize (J.size ());
  _ib = J.begin ();
  _ie = J.end ();
  _iQ = Q.begin ();
  _it = _integrand.begin ();
  for (_iJ = _ib; _iJ != _ie; ++_iJ, ++_iQ, ++_it)
    *_it = (*_iJ) * (*_iQ);

  // LHS = energy absorbed by the grain.
  _LHS = NumUtils::integrate<T1> (wave, _integrand);
  EAbs = _LHS;

  while (_TestConvergence > _Convergence)
    {

      _Tk = (Tlo + Thi) / 2.0;
      _integrand = NumUtils::prod_bbodyCGS<T1> (wave, _Tk, Q);
      _RHS = NumUtils::integrate<T1> (wave, _integrand);
      _Delta = (_LHS - _RHS) / _LHS;

      if (fabs (_Delta) < _Convergence)
        {
          _TestConvergence = _Convergence;
        }
      else
        {
          if (_Delta < 0.0)
            Thi = _Tk;
          else
            Tlo = _Tk;
        }
      if (Tlo == Thi)
        {
          // cout << "Equilibrium heating algorithm failure. Bye-bye." << endl;
          // exit(8);
          if (Tlo == lolim)
            return Flags::FEQ_ZEROBOUND_LOLIM;
          if (Thi == hilim)
            return Flags::FEQ_ZEROBOUND_HILIM;
          return Flags::FEQ_ZEROBOUND;
        }
    }

  return Flags::FSUCCESS;
}

#endif
