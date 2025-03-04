#ifndef _HEATUTILS_H
#define _HEATUTILS_H

#include <iostream>
#include <vector>

#include "Constants.h"
#include "NumUtils.h"

using namespace std;

namespace HeatUtils
{ // Define a namespace to avoid confusion with other
  // classes/libraries.

// Note that you can call these with any type, the meaning is not
// intuitively obvious for some.  They were constructed mainly for
// float/double types.

// ****************************************************************************
// Inline Templates
// ****************************************************************************

// Get T at mean_photon_energy of internal energy.
// If mean energy is below threshold of tabulated enthalpy (VERY unlikely!),
//  extrapolate with E^(1/4) to get the corresponding temperature.
// If mean energe is above threshold (destroy the grain?), extrapolate with E^1
// Otherwise, simply find bracketing values and interpolate.
template <typename T> T getTmean(vector<T> &Temp, vector<T> &Enth, T &mpe)
{
    int nEnth = Enth.size();
    if (mpe < Enth[0])
        return Temp[0] * sqrt(sqrt(mpe / Enth[0]));
    else if (mpe > Enth[nEnth - 1])
        return Temp[nEnth - 1] * (mpe / Enth[nEnth - 1]);
    else
    {
        int Tid = NumUtils::index(mpe, Enth);
        return NumUtils::line(Enth[Tid], Enth[Tid + 1], Temp[Tid], Temp[Tid + 1], mpe);
    }
}
// ****************************************************************************

// ****************************************************************************
// Template prototypes.
// ****************************************************************************
template <typename T> T getTauAbs(vector<T> &J, vector<T> &C, vector<T> &wave, T &mpe);
template <typename T> T getTauRad(vector<T> &C, vector<T> &w, T mpe, T Tu);
// ****************************************************************************
// Function 'source' code.
// ****************************************************************************

// ****************************************************************************
// Compute absorption timescale and mean absorbed photon energy
template <typename T> T getTauAbs(vector<T> &J, vector<T> &C, vector<T> &w, T &mpe)
{
    vector<T> integrand1(J.size());
    vector<T> integrand2(J.size());
    // Cannot template a typedef...
    // typedef vector<T>::iterator iter;
    typename vector<T>::iterator iJ;
    typename vector<T>::iterator iC = C.begin();
    typename vector<T>::iterator iw = w.begin();
    typename vector<T>::iterator ii1 = integrand1.begin();
    typename vector<T>::iterator ii2 = integrand2.begin();
    for (iJ = J.begin(); iJ != J.end(); iJ++)
    {
        *ii1 = (*iJ) * (*iC);
        *ii2 = (*ii1) * (*iw);
        ii1++;
        ii2++;
        iC++;
        iw++;
    }
    mpe = NumUtils::integrate<T>(w, integrand1);
    T TauAbs = NumUtils::integrate<T>(w, integrand2);
    TauAbs = (Constant::PLANCKLIGHT) / TauAbs;
    mpe *= TauAbs;
    TauAbs /= (4.0 * (Constant::PI));

    return TauAbs;
}
// ****************************************************************************

// ****************************************************************************
// Compute radiative cooling time from Temperature Tu.
template <typename T> T getTauRad(vector<T> &C, vector<T> &w, T mpe, T Tu)
{
    T ret;

    T mw = (Constant::PLANCKLIGHT / mpe);
    int uil = NumUtils::index(mw, w);
    int tsize = w.size() - uil;

    vector<T> tw(tsize);
    vector<T> integrand(tsize);
    typename vector<T>::iterator iC = C.begin() + uil;
    typename vector<T>::iterator iw = w.begin() + uil;
    typename vector<T>::iterator itw = tw.begin();

    *(itw) = mw;
    T midC = NumUtils::line(*(iw), *(iw + 1), *(iC), *(iC + 1), mw);
    iw++;
    itw++;

    for (; itw != tw.end(); itw++)
    {
        *(itw) = *(iw);
        iw++;
    }

    integrand = NumUtils::bbodyCGS<T>(tw, Tu);

    iC++; // Now at first CAbs after interpolated postion.
    typename vector<T>::iterator ii = integrand.begin();
    *(ii) *= midC;
    for (ii = integrand.begin() + 1; ii != integrand.end(); ii++)
    {
        *(ii) *= *(iC);
        iC++;
    }

    ret = NumUtils::integrate<T>(tw, integrand);
    return (mpe / (ret * 4.0 * (Constant::PI)));
}
// ****************************************************************************

} // namespace HeatUtils

#endif
