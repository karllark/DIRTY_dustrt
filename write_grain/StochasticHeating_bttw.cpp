// Jun 3-4 2009: Changed return value to int indicating status.  StochasticLum
// is now returned in the parameter list.
//               Added DirtyFlags.h dependency.
//               Modified return values of all sub calls to be int - will return
//               status to be passed back up the line -KAM
#include <iostream>
#include <vector>

#include "DirtyFlags.h"
#include "HeatUtils.h"
#include "NumUtils.h"

using namespace std;

int ComputeGrid (vector<float> &enth, vector<float> &denth, vector<float> &temp,
                 vector<float> &tgrid, vector<float> &Temperature, vector<float> &Enthalpy,
                 float &TMax, float &TMin, int &nBins);

int ComputeTransitionMatrix (vector<vector<double> > &TM, vector<float> &wave, vector<float> &temp,
                             vector<float> &enth, vector<float> &cabs, vector<float> &cJprod,
                             vector<float> &denth, int &nBins);

int SolveTransitionMatrix (vector<vector<double> > &TM, int &nBins, vector<double> &P);

int
StochasticHeating_bttw (vector<float> &wave, vector<float> &cJprod, vector<float> &cabs,
                        vector<float> &Temperature, vector<float> &Enthalpy, float &EAbs,
                        float &TMin, float &TMax, vector<double> &Emission,
                        vector<double> &ProbabilityDistribution, vector<float> &temp)
{
  int status;
  int idx;
  int nbins;
  int nwave;
  float PTol;
  float ETol;
  vector<float> enth;
  vector<float> denth;
  // vector <float>  temp;
  vector<float> tgrid;
  vector<double> integrand;
  vector<double> Probability_BBody; // Probability(T) * Black Body(T)
  double ProbabilitySum;            // Sum_T (Probability_BBody)
  double EnergyEmitted;

  vector<vector<double> > TransitionMatrix;
  // vector <double> ProbabilityDistribution;

  nbins = 2500; // Use 2500 bins for each size.  We no care about time here.
  PTol = 1.0e-15;
  ETol = 0.0;
  nwave = wave.size ();
  enth.resize (nbins, 0);
  denth.resize (nbins, 0);
  temp.resize (nbins, 0);
  tgrid.resize (nbins + 1, 0);
  integrand.resize (nwave, 0);
  ProbabilityDistribution.resize (nbins, 0);
  TransitionMatrix.resize (nbins);
  for (int i = 0; i < nbins; ++i)
    TransitionMatrix[i].resize (nbins, 0);

  // Initial setup
  status = ComputeGrid (enth, denth, temp, tgrid, Temperature, Enthalpy, TMax, TMin, nbins);
  // cout << "First grid" << TMin << " " << TMax << endl;
  // for (uint i=0;i<temp.size();++i) cout << i << " " << temp[i] << " " <<
  // enth[i] << " " << denth[i] << endl;
  if (status != Flags::FSUCCESS)
    return status;
  // refine the grid limits using maximum photon energy to truncate the enthalpy
  // grid
  //   idx = 1;
  //   MaxE = Constant::PLANCKLIGHT/wave[0];
  //   while (enth[idx]-enth[0] < MaxE && idx < (nbins-1)) ++idx;
  //   if (tgrid[idx] < TMax) {
  //     TMax = tgrid[idx];
  //     status =
  //     ComputeGrid(enth,denth,temp,tgrid,Temperature,Enthalpy,TMax,TMin,nbins);
  //     if (status != Flags::FSUCCESS) return status;
  //   }
  // cout << "Second grid " << TMin << " " << TMax << endl;
  // for (uint i=0;i<temp.size();++i) cout << i << " " << tgrid[i] << " " <<
  // temp[i] << " " << tgrid[i+1] << endl;
  status = ComputeTransitionMatrix (TransitionMatrix, wave, temp, enth, cabs, cJprod, denth, nbins);
  if (status != Flags::FSUCCESS)
    return status;
  status = SolveTransitionMatrix (TransitionMatrix, nbins, ProbabilityDistribution);
  if (status != Flags::FSUCCESS)
    return status;

  // cout << "first pass" << endl;
  // for (int i=0;i<nbins;++i) cout << i << " " << temp[i] << " " <<
  // ProbabilityDistribution[i] << endl;

  // Trim the grid, from upper limit
  idx = nbins - 1;
  while (ProbabilityDistribution[idx] <= PTol && idx > 0)
    --idx;
  if (idx < 0)
    return Flags::FST_SMALL_PROBABILITY_HI;
  TMax = tgrid[idx + 1];

  // cout << "Upper trim " << idx <<  " " << TMax << endl;
  // for (int i=0;i<nbins;++i) cout << i << " " << temp[i] << " " <<
  // ProbabilityDistribution[i] << endl;
  //  Trim the grid, lower limit
  idx = 0;
  if (ProbabilityDistribution[idx] < PTol)
    {
      while (ProbabilityDistribution[idx] < PTol)
        ++idx;
      if (idx > nbins - 1)
        return Flags::FST_SMALL_PROBABILITY_LO;
      if (idx == 0)
        TMin = tgrid[idx];
      else
        TMin = tgrid[idx - 1];
    }
  // cout << "Lower trim " << idx <<  " " << TMin << endl;

  // Refine nbins so we don't have too many elements at essentially the same
  // enthalpy.
  if ((TMax - TMin) / ((float)nbins) < 0.05)
    {
      nbins = static_cast<int> ((TMax - TMin) / 0.05);
      enth.resize (nbins, 0);
      denth.resize (nbins, 0);
      temp.resize (nbins, 0);
      tgrid.resize (nbins + 1, 0);
      ProbabilityDistribution.resize (nbins, 0);
      TransitionMatrix.resize (nbins);
      for (int i = 0; i < nbins; ++i)
        TransitionMatrix[i].resize (nbins, 0);
    }

  status = ComputeGrid (enth, denth, temp, tgrid, Temperature, Enthalpy, TMax, TMin, nbins);
  // for (int i=0;i<nbins;++i) cout << i << " " << temp[i] << " " << enth[i] <<
  // " " << denth[i] << endl;
  if (status != Flags::FSUCCESS)
    return status;
  status = ComputeTransitionMatrix (TransitionMatrix, wave, temp, enth, cabs, cJprod, denth, nbins);
  if (status != Flags::FSUCCESS)
    return status;
  status = SolveTransitionMatrix (TransitionMatrix, nbins, ProbabilityDistribution);
  if (status != Flags::FSUCCESS)
    return status;

  // Final pass to refine temperature limits
  // Trim the grid, from upper limit
  idx = nbins - 1;
  if (ProbabilityDistribution[idx] > PTol)
    TMax *= 1.5;
  idx = 0;
  if (ProbabilityDistribution[idx] < PTol)
    {
      while (ProbabilityDistribution[idx] < PTol)
        ++idx;
      if (idx > nbins - 1)
        return Flags::FST_SMALL_PROBABILITY_LO;
      if (idx == 0)
        TMin = tgrid[idx];
      else
        TMin = tgrid[idx - 1];
    }
  status = ComputeGrid (enth, denth, temp, tgrid, Temperature, Enthalpy, TMax, TMin, nbins);
  // for (int i=0;i<nbins;++i) cout << i << " " << temp[i] << " " << enth[i] <<
  // " " << denth[i] << endl;
  if (status != Flags::FSUCCESS)
    return status;
  status = ComputeTransitionMatrix (TransitionMatrix, wave, temp, enth, cabs, cJprod, denth, nbins);
  if (status != Flags::FSUCCESS)
    return status;
  status = SolveTransitionMatrix (TransitionMatrix, nbins, ProbabilityDistribution);
  if (status != Flags::FSUCCESS)
    return status;
  // cout << "Final: " << endl;
  // for (int i=0;i<nbins;++i) cout << i << " " << temp[i] << " " <<
  // ProbabilityDistribution[i] << endl;
  //  Compute energy coservation.
  for (int i = 0; i < nwave; ++i)
    {
      Probability_BBody
          = NumUtils::prod_bbodyCGS<double> (wave[i], temp, ProbabilityDistribution); // P(T)*B(T)
      ProbabilitySum = 0.0; // Redundant...
      ProbabilitySum = accumulate (Probability_BBody.begin (), Probability_BBody.end (),
                                   0.0); // Sum_T (P(T)*B(T))
      integrand[i] = cabs[i] * ProbabilitySum;
      Emission[i] = integrand[i]; // this is C(lam)*Sum_T (P(T)*B(T,lam)) ~
                                  // stochastic L(lam)
    }
  EnergyEmitted = NumUtils::integrate<double> (wave, integrand);
  ETol = abs (EAbs - EnergyEmitted) / EAbs;
  // cout << "Energy Balance: " << EnergyEmitted << " " << EAbs << " " << ETol
  // << " " << TMax << " " << TMin << endl; for (uint i=0;i<temp.size();i++)
  // cout << i << " " << temp[i] << " " << ProbabilityDistribution[i] << endl;

  return Flags::FSUCCESS;
}

// Generate the temperature and enthalpy grids necessary for calculation.
int
ComputeGrid (vector<float> &enth, vector<float> &denth, vector<float> &temp, vector<float> &tgrid,
             vector<float> &Temperature, vector<float> &Enthalpy, float &TMax, float &TMin,
             int &nBins)
{
  vector<float>::iterator _itb, _ite, _it, _it1, _it2;
  float _DT = (TMax - TMin) / static_cast<float> (nBins); // Linear temperature grid

  vector<float> _enthgrid (nBins + 1);

  // Define the temperature grid, end- and mid-points.
  _it1 = tgrid.begin ();
  _it2 = temp.begin ();
  *_it1 = TMin;
  ++_it1;
  _ite = tgrid.end ();
  for (_it = _it1; _it != _ite; ++_it, ++_it2)
    {
      *_it = *(_it - 1) + _DT;
      *_it2 = (*_it + *(_it - 1)) / 2.0;
    }

  // Interpolate Enthalpy onto _tgrid
  _enthgrid = NumUtils::interpol (Enthalpy, Temperature, tgrid, 4, -99);

  // Compute enthaply at bin center along with bin width.
  _itb = _enthgrid.begin () + 1;
  _ite = _enthgrid.end ();
  _it1 = enth.begin ();
  _it2 = denth.begin ();
  for (_it = _itb; _it != _ite; ++_it, ++_it1, ++_it2)
    {
      *_it1 = (*_it + *(_it - 1)) / 2.0;
      *_it2 = (*_it - *(_it - 1));
    }

  // No Failure modes...
  return Flags::FSUCCESS;
}

int
ComputeTransitionMatrix (vector<vector<double> > &TM, vector<float> &wave, vector<float> &temp,
                         vector<float> &enth, vector<float> &cabs, vector<float> &cJprod,
                         vector<float> &denth, int &nBins)
{
  float _wT, _slp, _icpt, _thiscjprod;
  int _nw = wave.size ();
  int idx;

  double c1 = 2.0 * Constant::PLANCK * pow (Constant::LIGHT, 2) * Constant::PLANCKLIGHT;
  double c2 = Constant::PLANCKLIGHT / Constant::BOLTZMAN;

  vector<float> _wave, _cjprod, _integrand;
  vector<float>::iterator it, itb, ite, it0, it1;

  _integrand.resize (_nw);

  for (int f = 0; f < nBins; ++f)
    { // Loop over final energy states

      if (f != nBins - 1)
        { // Cooling transitions - only f+1 to f

          itb = wave.begin ();
          ite = wave.end ();
          it0 = _integrand.begin ();
          it1 = cabs.begin ();
          for (it = itb; it != ite; ++it, ++it0, ++it1)
            *it0 = (*it1 / pow (*it, 5)) * (1.0 / (std::exp (c2 / ((*it) * temp[f + 1])) - 1.0));
          TM[f][f + 1]
              = c1 / (enth[f + 1] - enth[f]) * NumUtils::integrate<double> (wave, _integrand);
        } // Done with cooling transitions.

      for (int i = 0; i < f; ++i)
        { // Heating transitions - all initial < final.
          if ((enth[f] - enth[i]) != 0)
            _wT = Constant::PLANCKLIGHT / (enth[f] - enth[i]);
          else
            return Flags::FST_CTM_ZERO_ENTH_DEN;
          if (_wT < wave[0] || _wT > wave[_nw - 1])
            TM[f][i] = 0.0;
          else
            { // Photons exist that will induce this transition
              idx = 0;
              while (wave[idx] < _wT)
                idx++;
              if (idx > _nw - 1)
                return Flags::FST_CTM_TRANSITIONID_EXCEEDS_WAVEID;
              if (idx == 0 || idx == _nw - 1)
                _thiscjprod = cJprod[idx];
              else
                {
                  if ((wave[idx] - wave[idx - 1]) != 0)
                    _slp = (cJprod[idx] - cJprod[idx - 1]) / (wave[idx] - wave[idx - 1]);
                  else
                    return Flags::FST_CTM_WAVE_ZERO;
                  _icpt = cJprod[idx] - _slp * wave[idx];
                  _thiscjprod = _icpt + _slp * _wT;
                }
              TM[f][i] = Constant::IPLANCKLIGHT * _thiscjprod * _wT * _wT * _wT * denth[f];
              if (f == nBins - 1 && idx > 0)
                { // Put all transitions out of defined
                  // states into last state.
                  it0 = wave.begin () + idx;
                  it1 = cJprod.begin () + idx;
                  _wave.assign (wave.begin (), it0);
                  _cjprod.assign (cJprod.begin (), it1);
                  _wave[_wave.size () - 1] = _wT;
                  _cjprod[_wave.size () - 1] = _thiscjprod;
                  TM[f][i] += (_wT * NumUtils::integrate<double> (_wave, _cjprod));
                }

            } // Done with photons that can produce transition.

        } // Done with heating.
    } // Done with final statues.

  for (int f = 0; f < nBins; ++f)
    TM[0][0] -= TM[f][0];
  for (int i = 1; i < nBins; ++i)
    for (int j = i - 1; j < nBins; ++j)
      if (i != j)
        TM[i][i] -= TM[j][i];

  return Flags::FSUCCESS;
}

int
SolveTransitionMatrix (vector<vector<double> > &TM, int &nBins, vector<double> &P)
{
  vector<vector<double> > _Bij;
  _Bij.resize (nBins);
  for (int i = 0; i < nBins; ++i)
    _Bij[i].resize (nBins);
  for (int i = 0; i < nBins; ++i)
    {
      _Bij[nBins - 1][i] = TM[nBins - 1][i];
      for (int j = nBins - 2; j >= 0; --j)
        {
          _Bij[j][i] = _Bij[j + 1][i] + TM[j][i];
        }
    }

  P[0] = 1.0;
  double _norm = 1.0;
  for (int f = 1; f < nBins; ++f)
    {
      P[f] = 0.0;
      for (int i = 0; i < f; ++i)
        {
          P[f] += (P[i] * _Bij[f][i]);
        }
      if (TM[f - 1][f] != 0.0)
        P[f] /= TM[f - 1][f];
      else
        return Flags::FST_STM_ZERO_COOLING;
      _norm += P[f];
    }

  transform (P.begin (), P.end (), P.begin (), bind2nd (divides<double> (), _norm));

  return Flags::FSUCCESS;
}
