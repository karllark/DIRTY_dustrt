#include "ISRF.h"

ISRF::ISRF (vector<float> in_wave, float in_XMMP)
{
    // cout << "Instantiating ISRF Object with scaling " <<  in_XMMP << endl;

    float _a[] = { -6.46657e-09, 8.80371e-06, 1.29912e-06, 1.84887e-08, -0.00106653, 13.1090 };
    float _b[] = { 1.0e-14, 1.0e-13, 4.0e-13 };
    float _c[] = { 7500., 4000., 3000. };
    vector<float> _UVP (_a, _a + 6);
    vector<float> _W (_b, _b + 3);
    vector<float> _T (_c, _c + 3);
    vector<float> _TempISRF;

    float JtoU = 4.0 * Constant::PI / Constant::LIGHT;

    // Populate class members.
    XMMP = in_XMMP;
    wave.assign (in_wave.begin (), in_wave.end ());
    nWave = wave.size ();
    theISRF.resize (nWave, 0);

    vector<float>::iterator iW;
    vector<float>::iterator iT = _T.begin ();
    vector<float>::iterator iTempISRF;
    vector<float>::iterator iISRF;

    for (iW = _W.begin (); iW != _W.end (); iW++, ++iT)
        {
            // _TempISRF = black body at T_i
            _TempISRF = NumUtils::bbodyCGS<float> (wave, *iT);
            iISRF = theISRF.begin ();

            // _TempISRF = bbody at T_i scaled by W_i
            // ISRF = sum of scaled bb's at T converted to U (UV fitting done in U...)
            for (iTempISRF = _TempISRF.begin (); iTempISRF != _TempISRF.end (); iTempISRF++, ++iISRF)
                {
                    *iTempISRF *= *iW;
                    *iISRF += (*iTempISRF) * JtoU;
                }

            // transform(TempISRF.begin(),TempISRF.end(),TempISRF.begin(),bind2nd(multiplies<double>(),(*iW)));
            // transform(ISRF.begin(),ISRF.end(),TempISRF.begin(),ISRF.begin(),plus<double>());
            //  bbodyCGS uses push_back
            _TempISRF.clear ();
        }

    // For UV wavelengths, add the UV portion of the ISRF
    // Wavelength range over which the UV fit is valid.
    float UVMin = 0.09e-4;
    float UVMax = 0.25e-4;

    vector<float>::iterator loIter, hiIter;
    // float UVMod;
    loIter = find_if (wave.begin (), wave.end (), bind2nd (greater_equal<float> (), UVMin));
    if (loIter != wave.end ())
        {
            hiIter = find_if (wave.begin (), wave.end (), bind2nd (greater<float> (), UVMax));
            int loIdx = distance (wave.begin (), loIter);
            int hiIdx = distance (wave.begin (), hiIter);
            for (int i = loIdx; i < hiIdx; i++)
                theISRF[i] += _UVP[0] * exp (-pow ((wave[i] - _UVP[1]) / _UVP[2], 2) / 2.0) + _UVP[3]
                              + _UVP[4] * wave[i] + _UVP[5] * wave[i] * wave[i];

            // Make everything below low cutoff equal to 0.0
            if (loIdx >= 0)
                for (int i = 0; i < loIdx; i++)
                    theISRF[i] = 0.0;
        }

    // If we happen to have snuck some negative numbers in (shouldn't happen)...
    replace_if (theISRF.begin (), theISRF.end (), bind2nd (less<float> (), 0), 0);

    // Since we are returning J, transform U to J.  Convulted, I know...
    // One division as opposed to ISRF.size() divisions.  Scale the radiation
    // field by X_MMP
    float UtoJ = 1.0 / JtoU * XMMP;
    transform (theISRF.begin (), theISRF.end (), theISRF.begin (), bind2nd (multiplies<float> (), UtoJ));
    // cout << "DONE!" << endl;
}
