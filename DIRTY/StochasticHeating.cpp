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

int ComputeGrid(vector<float> &enth, vector<float> &denth, vector<float> &temp, vector<float> &tgrid,
                vector<float> &Temperature, vector<float> &Enthalpy, float &TMax, float &TMin, int &nBins);

int ComputeTransitionMatrix(vector<vector<double>> &TM, vector<float> &wave, vector<float> &temp, vector<float> &enth,
                            vector<float> &cabs, vector<float> &cJprod, vector<float> &denth, int &nBins);

int SolveTransitionMatrix(vector<vector<double>> &TM, int &nBins, vector<double> &P);

int StochasticHeating(vector<float> &wave, vector<float> &cJprod, vector<float> &cabs, vector<float> &Temperature,
                      vector<float> &Enthalpy, float EAbs, float &TMin, float &TMax, float &TEq, uint _sz,
                      vector<vector<double>> &StochasticLum)

{

    int maxBins = 1000;
    bool converged = false;
    bool IncreaseBins = false;
    bool lastincrease = true;
    bool OnePass = false;
    int nBins = 50;

    float oldTMax = 0.0;
    float oldTMin = 0.0;
    float tol = 0.01;
    float tol_max_bins = 0.1;
    float thistol, lasttol = 0.0;
    double Ptol = 1.0e-15;
    vector<vector<double>> TM;

    vector<double> _P;

    int status;
    int nWave = wave.size();

    int idx; //,idx1,idxp;
    vector<double> _pofTint;

    double _pFac;

    vector<float>::iterator it, it0, it1, it2, it3, it4, ib, ie;
    vector<double>::iterator idt, idb, ide, idt1, idb1, ide1, idt2;

    // Some initial algorithm to compute search size in constant T

    vector<float> _enth, _denth, _temp, _tgrid;

    vector<double> integrand(nWave);

    StochasticLum[_sz].resize(nWave);

    double Eemit;

    // bins
    _enth.reserve(maxBins);
    _denth.reserve(maxBins);
    _temp.reserve(maxBins);
    _tgrid.reserve(maxBins + 1);
    // Size to max - memory inefficient, cpu efficient(?)
    _P.resize(maxBins, 0);
    TM.resize(maxBins);
    for (int i = 0; i < maxBins; ++i)
        TM[i].resize(maxBins);

    // Isolate temperature region.
    _enth.resize(nBins, 0);
    _denth.resize(nBins, 0);
    _temp.resize(nBins, 0);
    _tgrid.resize(nBins + 1, 0);
    status = ComputeGrid(_enth, _denth, _temp, _tgrid, Temperature, Enthalpy, TMax, TMin, nBins);
    idx = 1;
    float MaxE = Constant::PLANCKLIGHT / wave[0];
    int idx_max = _tgrid.size() - 1;
    while (_enth[idx] - _enth[0] < MaxE && idx < idx_max)
        ++idx;
    TMax = _tgrid[idx];
    if (TMax < 1.5 * TEq)
    {
        if (TEq < 100.)
            TMax = 1.5 * TEq;
        else
            TMax = TEq + 100.;
    }

    bool _setup = true;
    while (_setup)
    {
        // Setup Grid
        status = ComputeGrid(_enth, _denth, _temp, _tgrid, Temperature, Enthalpy, TMax, TMin, nBins);
        // Setup transition matrix
        status = ComputeTransitionMatrix(TM, wave, _temp, _enth, cabs, cJprod, _denth, nBins);
        // Setup probability distribution
        status = SolveTransitionMatrix(TM, nBins, _P);
        if (status != Flags::FSUCCESS)
            return status;
        // Adjust temperature
        if (_P[nBins - 1] > Ptol)
            TMax *= 1.5;
        else
            _setup = false;
    }

    int loop_count = 0;
    // Convergence wrapper...
    while (!converged)
    { // convergence bracket

        // Compute grid
        _enth.resize(nBins, 0);
        _denth.resize(nBins, 0);
        _temp.resize(nBins, 0);
        _tgrid.resize(nBins + 1, 0);
        status = ComputeGrid(_enth, _denth, _temp, _tgrid, Temperature, Enthalpy, TMax, TMin, nBins);
        if (status != Flags::FSUCCESS)
            return status;
        // Compute transition matrix
        status = ComputeTransitionMatrix(TM, wave, _temp, _enth, cabs, cJprod, _denth, nBins);
        if (status != Flags::FSUCCESS)
            return status;
        // Solve transition matrix
        status = SolveTransitionMatrix(TM, nBins, _P);
        if (status != Flags::FSUCCESS)
            return status;

        // integrand.resize(nWave);
        idb = integrand.begin();
        it1 = wave.begin();
        it = cabs.begin();
        for (int i = 0; i < nWave; i++, idb++, it1++, it++)
        {
            _pofTint = NumUtils::prod_bbodyCGS<double>(*it1, _temp,
                                                       _P); // P(T)*B(T) at wv,_sz
            _pFac = 0.0;
            _pFac = accumulate(_pofTint.begin(), _pofTint.end(),
                               0.0); // Sum_T (P(T)*B(T))
            *idb = (*it) * _pFac;
            StochasticLum[_sz][i] = *idb; // this is C(lam)*Sum_T (P(T)*B(T,lam)) ~ stochastic L(lam)
        }
        Eemit = NumUtils::integrate<double>(wave, integrand);
        thistol = abs(EAbs - Eemit) / EAbs;

        if (thistol > tol)
        { // not converged

            if (thistol > 100. * lasttol && OnePass)
            {
                IncreaseBins = true;
                TMin = oldTMin;
                TMax = oldTMax;
            }
            else
            {
                lasttol = thistol; // save current tolerance.
                oldTMin = TMin;
                oldTMax = TMax;
                idx = nBins - 1;

                // Work on upper bounds
                if (_P[idx] == 0.0)
                {
                    while (_P[idx] == 0.0)
                        --idx;
                    if (idx < 0)
                        return Flags::FST_ZERO_PROBABILITY;
                    if (_P[idx] < Ptol)
                    {
                        while (_P[idx] < Ptol)
                            --idx;
                        if (idx < 0)
                            return Flags::FST_SMALL_PROBABILITY_HI;
                        TMax = _tgrid[idx + 1];
                    }
                    else
                    {
                        TMax = _tgrid[nBins];
                    }
                }

                if (_P[idx] < Ptol)
                {
                    while (_P[idx] < Ptol)
                        --idx;
                    if (idx < 0)
                        return Flags::FST_SMALL_PROBABILITY_HI;
                    TMax = _tgrid[idx + 1];
                }

                if (_P[nBins - 1] >= Ptol)
                {
                    TMax = 1.5 * _tgrid[nBins];
                    IncreaseBins = true;
                }

                // set lower bounds
                idx = 0;
                if (_P[idx] < Ptol)
                {
                    while (_P[idx] < Ptol)
                        ++idx;
                    if (idx > nBins - 1)
                        return Flags::FST_SMALL_PROBABILITY_LO;
                    if (idx == 0)
                        TMin = _tgrid[idx];
                    else
                        TMin = _tgrid[idx - 1];
                }
            }

            if (IncreaseBins)
            {
                nBins = static_cast<int>(1.5 * static_cast<double>(nBins));
            }

            if (nBins > maxBins)
            {
                if (lastincrease)
                {
                    // accept a lower tolerance for these cases
                    if (thistol > tol_max_bins)
                    { // still not converged well enough
                        return Flags::FST_EXCEED_MAXBINS;
                    }
                    else
                    { // close enough for DIRTY work.
                        converged = true;
                    }
                }
            }
            OnePass = true;
        }
        else
            converged = true;

        if (loop_count >= 100)
        {
            cout << "** Breaking the convergence loop in StochasticHeating(): loop " << loop_count << endl;
            cout << "** Input Parameters:" << endl;
            cout << "wave:\t";
            for (uint i = 0; i < wave.size(); i++)
                cout << wave[i] << ' ';
            cout << endl;
            cout << "cJprod:\t";
            for (uint i = 0; i < cJprod.size(); i++)
                cout << cJprod[i] << ' ';
            cout << endl;
            cout << "cabs:\t";
            for (uint i = 0; i < cabs.size(); i++)
                cout << cabs[i] << ' ';
            cout << endl;
            cout << "Temperature:\t";
            for (uint i = 0; i < Temperature.size(); i++)
                cout << Temperature[i] << ' ';
            cout << endl;
            cout << "Enthalpy:\t";
            for (uint i = 0; i < Enthalpy.size(); i++)
                cout << Enthalpy[i] << ' ';
            cout << endl;
            cout << "EAbs:\t" << EAbs << endl;
            cout << "TMin:\t" << TMin << endl;
            cout << "TMax:\t" << TMax << endl;
            cout << "TEq:\t" << TEq << endl;
            cout << "** Calculation variables:" << endl;
            cout << "integrand:\t";
            for (uint i = 0; i < integrand.size(); i++)
                cout << integrand[i] << ' ';
            cout << endl;
            cout << "Eemit:\t" << Eemit << endl;
            cout << "tol:\t" << tol << endl;
            cout << "thistol:\t" << thistol << endl;
            cout << "lasttol:\t" << lasttol << endl;
            cout << "nBins:\t" << nBins << endl;
            cout << "maxBins:\t" << maxBins << endl;
            return Flags::FST_EXCEED_CONVERGENCE_LOOP_COUNT;
        }
        ++loop_count;
    }

    return Flags::FSUCCESS;
}

// Generate the temperature and enthalpy grids necessary for calculation.
int ComputeGrid(vector<float> &enth, vector<float> &denth, vector<float> &temp, vector<float> &tgrid,
                vector<float> &Temperature, vector<float> &Enthalpy, float &TMax, float &TMin, int &nBins)
{

    vector<float>::iterator _itb, _ite, _it, _it1, _it2;
    float _DT = (TMax - TMin) / static_cast<float>(nBins); // Linear temperature grid

    vector<float> _enthgrid(nBins + 1);

    // Define the temperature grid, end- and mid-points.
    _it1 = tgrid.begin();
    _it2 = temp.begin();
    *_it1 = TMin;
    ++_it1;
    _ite = tgrid.end();
    for (_it = _it1; _it != _ite; ++_it, ++_it2)
    {
        *_it = *(_it - 1) + _DT;
        *_it2 = (*_it + *(_it - 1)) / 2.0;
    }

    // Interpolate Enthalpy onto _tgrid
    _enthgrid = NumUtils::interpol(Enthalpy, Temperature, tgrid, 4, -99);

    // Compute enthaply at bin center along with bin width.
    _itb = _enthgrid.begin() + 1;
    _ite = _enthgrid.end();
    _it1 = enth.begin();
    _it2 = denth.begin();
    for (_it = _itb; _it != _ite; ++_it, ++_it1, ++_it2)
    {
        *_it1 = (*_it + *(_it - 1)) / 2.0;
        *_it2 = (*_it - *(_it - 1));
    }

    // No Failure modes...
    return Flags::FSUCCESS;
}

int ComputeTransitionMatrix(vector<vector<double>> &TM, vector<float> &wave, vector<float> &temp, vector<float> &enth,
                            vector<float> &cabs, vector<float> &cJprod, vector<float> &denth, int &nBins)
{

    float _wT, _slp, _icpt, _thiscjprod;
    int _nw = wave.size();
    int idx;

    double c1 = 2.0 * Constant::PLANCK * pow(Constant::LIGHT, 2) * Constant::PLANCKLIGHT;
    double c2 = Constant::PLANCKLIGHT / Constant::BOLTZMAN;

    vector<float> _wave, _cjprod, _integrand;
    vector<float>::iterator it, itb, ite, it0, it1;

    _integrand.resize(_nw);

    for (int f = 0; f < nBins; ++f)
    { // Loop over final energy states

        if (f != nBins - 1)
        { // Cooling transitions - only f+1 to f

            itb = wave.begin();
            ite = wave.end();
            it0 = _integrand.begin();
            it1 = cabs.begin();
            for (it = itb; it != ite; ++it, ++it0, ++it1)
                *it0 = (*it1 / pow(*it, 5)) * (1.0 / (std::exp(c2 / ((*it) * temp[f + 1])) - 1.0));
            TM[f][f + 1] = c1 / (enth[f + 1] - enth[f]) * NumUtils::integrate<double>(wave, _integrand);
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
                    it0 = wave.begin() + idx;
                    it1 = cJprod.begin() + idx;
                    _wave.assign(wave.begin(), it0);
                    _cjprod.assign(cJprod.begin(), it1);
                    _wave[_wave.size() - 1] = _wT;
                    _cjprod[_wave.size() - 1] = _thiscjprod;
                    TM[f][i] += (_wT * NumUtils::integrate<double>(_wave, _cjprod));
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

int SolveTransitionMatrix(vector<vector<double>> &TM, int &nBins, vector<double> &P)
{

    vector<vector<double>> _Bij;
    _Bij.resize(nBins);
    for (int i = 0; i < nBins; ++i)
        _Bij[i].resize(nBins);
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

    transform(P.begin(), P.end(), P.begin(), bind2nd(divides<double>(), _norm));

    return Flags::FSUCCESS;
}
