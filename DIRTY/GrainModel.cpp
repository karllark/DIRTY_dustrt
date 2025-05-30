//*******************************************************************************
// GrainModel -- GrainModel Class
//
// Requires:
//     ConfigFile
//     Grain
//     StringManip
//
// History:
//     Add THEMIS2 model: Jan 2024
//     Add HD2023 model: Jan 2024
//     Written by:   Putrid, Winter 2006/2007
//     Contact: misselt@as.arizona.edu
//
//*******************************************************************************

#include "GrainModel.h"

GrainModel::GrainModel()
{
}

void GrainModel::MakeGrainModel(ConfigFile &_mbkcf, vector<float> &MasterWave)
{
    string _thisSection;
    string _line;

    string _thisCrossSectionFile;
    string _thisCalorimetryFile;
    string _thisSizeDistribution;
    string _thisSizeType;
    string _thisSizeScale;
    string _thisSizeFile;

    string _TopLevelPath;
    string _CrossSectionsPath;
    string _CalorimetryPath;
    string _ModelDefinitionPath;
    string _FullModelDefinitionFile;

    int _thisSizeIN;
    float _thisA_min;
    float _thisA_max;
    float _thisClip;
    float _delta;
    float _totalMassFraction = 0;

    vector<float> _thisSize;
    vector<float> _thisCAbs;
    vector<float> _thisCSca;
    vector<float> _thisphFunc;
    vector<float> _integrand;

    float delta, _CAbsIntegral, _CScaIntegral, _phFuncIntegral; //,_MassConstants;

    typedef vector<float>::iterator iter;
    iter iCAbsEff, iCScaEff, iphFuncEff, iTau, iAlbedo;
    iter ithisCAbs, ithisCSca, ithisphFunc, ithisSize, ithisSizeDist;

    // int _thisNSize;

    // Constants that will multiply the mass integrals
    //_MassConstants = (4.0/3.0)*(Constant::PI);

    // Make the wavelength vector for the GrainModel the input MasterWave.
    // Generate an energy scale at the same time.
    copy(MasterWave.begin(), MasterWave.end(), back_inserter(Wave));
    nWave = Wave.size();
    vector<float>::iterator _itb, _ite, _it;
    _itb = Wave.begin();
    _ite = Wave.end();
    for (_it = _itb; _it != _ite; ++_it)
        EScale.push_back(Constant::PLANCKLIGHT / (*_it));

    // Setup the model string->id map for case statements.
    ModelMapping();
    // Setup the size distribution string->id map for case statements.
    SizeDistMapping();
    // Setup the size type string->id map for case statements.
    SizeTypeMapping();

    // Take input size tabulation as gospel for each grain material.
    vector<float> _MasterSize;
    //_MasterSize.push_back(-1);

    // Parse the Grain config file.
    // ConfigFile _mbkcf(ModelBookKeeping);

    // Model Book Keeping.  Paths, Model definition configuration file.
    // Top Level:
    _TopLevelPath = _mbkcf.SValue("Model Book Keeping", "Path to Dust Properties");
    // CrossSections:
    _CrossSectionsPath = _mbkcf.SValue("Model Book Keeping", "Cross Section SubDir");
    // Calorimetry:
    _CalorimetryPath = _mbkcf.SValue("Model Book Keeping", "Calorimetry SubDir");
    // Model Definiiton:
    _ModelDefinitionPath = _mbkcf.SValue("Model Book Keeping", "Model SubDir");
    // Model Name:
    ModelName = _mbkcf.SValue("Model Book Keeping", "Model Name");

    // Check that model name exists in pre-defined map.
    if (ModelID.find(ModelName) == ModelID.end())
    { // NO!
        cout << "****************************************************" << endl;
        cout << "You are not using a predefined grain model.         " << endl;
        cout << endl << "Pre-deined models are: " << endl;
        for (iMap = ModelID.begin(); iMap != ModelID.end(); iMap++)
            cout << iMap->first << endl;
        cout << "****************************************************" << endl;
        // exit(8);
    }
    // YES! Check and make sure model definition file exists.
    _FullModelDefinitionFile = _TopLevelPath + _ModelDefinitionPath + ModelName;
    if (!StringManip::FileExists(_FullModelDefinitionFile))
    { // NO!
        cout << "****************************************************" << endl;
        cout << "Model Definition File " << _FullModelDefinitionFile << " not found." << endl;
        cout << "****************************************************" << endl;
        exit(8);
    }
    // YES! Proceed with model definition.

    // Create a ConfigFile Object for the actual model definition.
    ConfigFile _mcf(_FullModelDefinitionFile);

    // Get number of components and resize member vectors as necessary.
    nComp = _mcf.IValue("Model", "Number of Components");
    Component.resize(nComp);
    Normalization.resize(nComp);
    SizeDistribution.resize(nComp);
    SizeDistributionNorm.resize(nComp);
    Temperature.resize(nComp);
    DustMass.resize(nComp);
    DustVolume.resize(nComp);
    TotalDustMass = 0.0;

    // Size the size distribution/composition averaged quantities to masterwave
    // size
    CAbsEff.resize(nWave, 0);
    CScaEff.resize(nWave, 0);
    phFuncEff.resize(nWave, 0);
    Tau.resize(nWave, 0);
    Albedo.resize(nWave, 0);

    // Get each component.
    //  - Wavelength grid will be set by MasterWave.
    //  - Size grid will be set by the Cross Section file for each component.

    totalnumber = 0;
    number.resize(nComp);
    // Construct the Size distribution for each component.
    for (int cmp = 0; cmp < nComp; cmp++)
    { // Loop over all components

        // Construct section key for this component.
        // Assume key is of the form "Comp i" where i=1,2,...,nComp
        _thisSection = "Component " + StringManip::vtos(cmp + 1);

        // Get the Size Distribution information.  Map in GrainModel.h is in
        // upper case, so convert to upper case.
        _thisSizeDistribution = _mcf.SValue(_thisSection, "Size Distribution");
        transform(_thisSizeDistribution.begin(), _thisSizeDistribution.end(), _thisSizeDistribution.begin(),
                  (int (*)(int))toupper);

        // Clean out _MasterSize
        _MasterSize.erase(_MasterSize.begin(), _MasterSize.end());

        // Get grain property definitions.
        _thisCrossSectionFile = _mcf.SValue(_thisSection, "Cross Sections");
        _thisCalorimetryFile = _mcf.SValue(_thisSection, "Calorimetry");

        _thisA_min = _mcf.FValue(_thisSection, "a_min") * Constant::UM_CM;
        _thisA_max = _mcf.FValue(_thisSection, "a_max") * Constant::UM_CM;

        // For a Gaussian (or arbitrary in the future?) size distribution, check
        // that a_min and _amax are correctly defined.
        if (SizeDistID.find(_thisSizeDistribution) != SizeDistID.end())
        {

            if (_thisSizeDistribution == "GAUSS")
            { // Gaussian - add other types as necessary

                float _a0 = _mcf.FValue(_thisSection, "a0");
                if (_mcf.isBadFloat(_a0))
                {
                    cout << endl;
                    cout << "************************************************" << endl;
                    cout << "Did not find key a0 in section " << _thisSection << endl;
                    cout << "a0 is required for GAUSSIAN size distribution." << endl;
                    cout << "************************************************" << endl;
                    exit(8);
                }

                float _as = _mcf.FValue(_thisSection, "as");
                if (_mcf.isBadFloat(_as))
                {
                    cout << endl;
                    cout << "************************************************" << endl;
                    cout << "Did not find key as in section " << _thisSection << endl;
                    cout << "as is required for GAUSSIAN size distribution." << endl;
                    cout << "************************************************" << endl;
                    exit(8);
                }

                _thisClip = _mcf.FValue(_thisSection, "clip");
                if (_mcf.isBadFloat(_thisClip))
                    _thisClip = 3;

                if (_mcf.isBadFloat(_thisA_min))
                    _thisA_min = _a0 - _thisClip * _as;
                if (_thisA_min < 0)
                    _thisA_min = 0;

                if (_mcf.isBadFloat(_thisA_max))
                    _thisA_max = _a0 + _thisClip * _as;

                _thisA_min *= Constant::UM_CM;
                _thisA_max *= Constant::UM_CM;

            } // End Gaussian size distribution
        }
        else
        { // Encountered a SizeDist key that is not defined in the map.

            cout << endl;
            cout << "****************************************************" << endl;
            cout << _thisSizeDistribution << " is not a defined size distribution type." << endl;
            cout << "Allowed size distribution types are: " << endl;
            for (iMap = SizeDistID.begin(); iMap != SizeDistID.end(); iMap++)
                //	cout << iMap->first << " " << iMap->second << endl;
                cout << "****************************************************" << endl;
            exit(8);
        }

        _thisSizeType = _mcf.SValue(_thisSection, "Size Definition Type");
        transform(_thisSizeType.begin(), _thisSizeType.end(), _thisSizeType.begin(), (int (*)(int))toupper);
        // SizeTypeID is a STRING (_thisSizeType) to INTEGER map defined in
        // GrainModel.h NULL  -> 0  ; size grid will be taken from the optical
        // constants definition file DEF   -> 0  ; size grid will be taken from the
        // optical constants definition file FILE  -> 1  ; size grid will be read
        // from an on disk tabulation NSIZE -> 2  ; size grid will be constructed
        // from amin,amax,nsize for each component
        switch (SizeTypeID[_thisSizeType])
        {
        case 0: // Explicit DEF definition as well as key not found
        {
            _MasterSize.push_back(-1);
            break;
        }
        case 1: {
            // Filename
            // Assumptions:
            //   - Model configuration file provides the full path to file.
            //     eg. if the file is in the current directory, the config
            //     file should have "Size Definition File=MySizeGrid.dat"
            //     If the file lives in /home/dirty, the config file should have
            //     "Size Definition File=/home/dirty/MySizeGrid.dat"
            //   - The units of the size grid are assumed to be in um and are
            //   converted
            //     to cm via Constant::UM_CM
            _thisSizeFile = _mcf.SValue(_thisSection, "Size Definition File");
            if (_thisSizeFile == _mcf.BadString())
            {
                cout << "Unable to retrieve size grid file name, resorting to default." << endl;
                _MasterSize.push_back(-1);
            }
            else
            {
                // ifstream
                // _file((_TopLevelPath+_ModelDefinitionPath+_thisSizeFile).c_str());
                ifstream _file((_thisSizeFile).c_str());
                if (!_file.is_open())
                {
                    cout << "Unable to open size grid definition file " << _thisSizeFile << endl;
                    exit(8);
                }
                while (getline(_file, _line))
                {
                    if (_line[0] == '#')
                        continue;
                    else
                    {
                        _MasterSize.push_back(atof(_line.c_str()) * Constant::UM_CM);
                    }
                }
                // Trim the master size to our min/max
                _MasterSize.erase(remove_if(_MasterSize.begin(), _MasterSize.end(), bind2nd(less<float>(), _thisA_min)),
                                  _MasterSize.end());
                _MasterSize.erase(
                    remove_if(_MasterSize.begin(), _MasterSize.end(), bind2nd(greater<float>(), _thisA_max)),
                    _MasterSize.end());
                // And make sure end points match the min/max
                if (_MasterSize[0] > _thisA_min)
                    _MasterSize.insert(_MasterSize.begin(), _thisA_min);
                if (_MasterSize[_MasterSize.size() - 1] < _thisA_max)
                    _MasterSize.push_back(_thisA_max);
            }
            break;
        }
        case 2: {
            // Get the number of sizes
            _thisSizeIN = _mcf.IValue(_thisSection, "Size Definition Number");
            cout << "Generating size grid internally" << endl;
            if (_thisSizeIN <= 0)
            {
                cout << "Unable to define size grid manually.  Resorting to default." << endl;
                _thisSizeType = "DEF";
                _MasterSize.push_back(-1);
            }
            else
            {
                _MasterSize.resize(_thisSizeIN);
                // log/linear
                _thisSizeScale = _mcf.SValue(_thisSection, "Size Definition Scale");
                if (_thisSizeScale == _mcf.BadString())
                {
                    cout << "Did not find size gridding definition - assuming linear" << endl;
                    _thisSizeScale = "LINEAR";
                }
                else
                {
                    transform(_thisSizeScale.begin(), _thisSizeScale.end(), _thisSizeScale.begin(),
                              (int (*)(int))toupper);
                }
                if (_thisSizeScale == "LINEAR")
                { // Generate Linear scale
                    _MasterSize[0] = _thisA_min;
                    _delta = (_thisA_max - _thisA_min) / static_cast<float>(_thisSizeIN - 1);
                    for (int _sz = 1; _sz < _thisSizeIN; ++_sz)
                        _MasterSize[_sz] = _MasterSize[_sz - 1] + _delta;
                }
                if (_thisSizeScale == "LOG")
                { // Generate Log scale.
                    _MasterSize[0] = _thisA_min;
                    _delta = (1.0 / static_cast<float>(_thisSizeIN - 1)) * log10(_thisA_max / _thisA_min);
                    for (int _sz = 1; _sz < _thisSizeIN; ++_sz)
                        _MasterSize[_sz] = _MasterSize[_sz - 1] * pow(10, _delta);
                }
            }
            break;
        }
        default: // Unable to identify size grid, use the default.
        {
            _thisSizeType = "DEF";
            _MasterSize.push_back(-1);
            break;
        }
        }

        // Populate Component[i] with the appropriate Grain object.
        Component[cmp].MakeGrain(_CrossSectionsPath + _thisCrossSectionFile, _CalorimetryPath + _thisCalorimetryFile,
                                 MasterWave, _MasterSize, _TopLevelPath, _thisA_min, _thisA_max);

        SizeDistribution[cmp].resize(Component[cmp].nsize);

        // Wrap in a .find since SizeDist[key] will add members.
        if (SizeDistID.find(_thisSizeDistribution) != SizeDistID.end())
        { // Determine which sd to use

            // SizeDistID is a STRING -> INTEGER map defined in GrainModel.h
            // ZDA    = 0;    // Zubko et al. Size distribution function
            // WD01   = 1;    // Wiengartner+Draine Size distribution function
            // PWRLAW = 2;    // Power law size distribution - e.g. THEMIS
            // GAUSS  = 3;    // Appox. single size grain distribution - Guassian
            // HD23   = 4;    // Hensley & Draine 2023
            // LOGNRM = 5;    // Log-Normal size distribution - e.g. THEMIS
            // NULL   = 99;
            switch (SizeDistID[_thisSizeDistribution])
            {

            case 0: // Zubko, Dwek, & Arendt 2004, ApJS, 152, 211
            {
                // So we need to get the co-efficients for a ZDA fit.
                vector<float> _zda_coeff(15);
                _zda_coeff[0] = _mcf.FValue(_thisSection, "A");
                Normalization[cmp] = _zda_coeff[0];
                _zda_coeff[1] = _mcf.FValue(_thisSection, "c0");
                _zda_coeff[2] = _mcf.FValue(_thisSection, "b0");
                _zda_coeff[3] = _mcf.FValue(_thisSection, "b1");
                _zda_coeff[4] = _mcf.FValue(_thisSection, "a1");
                _zda_coeff[5] = _mcf.FValue(_thisSection, "m1");
                _zda_coeff[6] = _mcf.FValue(_thisSection, "b2");
                _zda_coeff[7] = _mcf.FValue(_thisSection, "a2");
                _zda_coeff[8] = _mcf.FValue(_thisSection, "m2");
                _zda_coeff[9] = _mcf.FValue(_thisSection, "b3");
                _zda_coeff[10] = _mcf.FValue(_thisSection, "a3");
                _zda_coeff[11] = _mcf.FValue(_thisSection, "m3");
                _zda_coeff[12] = _mcf.FValue(_thisSection, "b4");
                _zda_coeff[13] = _mcf.FValue(_thisSection, "a4");
                _zda_coeff[14] = _mcf.FValue(_thisSection, "m4");
                SizeDistribution[cmp] = getZDA_sdist(_zda_coeff, cmp);
            }
            break;

            case 1: // Wiengartner & Draine 2001, ApJ, 548, 296
            {
                vector<float> _wd_coeff(8);
                Normalization[cmp] = 1.0; // Defined functions already return per H
                _wd_coeff[0] = _mcf.FValue(_thisSection, "BC");
                _wd_coeff[1] = _mcf.FValue(_thisSection, "B1");
                _wd_coeff[2] = _mcf.FValue(_thisSection, "B2");
                _wd_coeff[3] = _mcf.FValue(_thisSection, "C");
                _wd_coeff[4] = _mcf.FValue(_thisSection, "at");
                _wd_coeff[5] = _mcf.FValue(_thisSection, "ac");
                _wd_coeff[6] = _mcf.FValue(_thisSection, "alpha");
                _wd_coeff[7] = _mcf.FValue(_thisSection, "beta");
                SizeDistribution[cmp] = getWD_sdist(_wd_coeff, cmp);
            }
            break;

            case 2: // Power law e.g. THEMIS-2 implement only necessary coefficients
                    // for THEMIS (pl + exp cutoff)
                // - add terms as necessary for other forms.  We will need to update
                // model definition files for e.g. 'incomplete' (less complex) power
                // laws when we add terms POWER LAW definition:  A (normalization),
                // alpha (slope of power law), at, ac, gamma (for a cutoff term, and au,
                // zeta, zxp, gama (curvature) cutoff: exp(-((a-at)/aC)^gamma)
                // curvature: (1 + zeta*(a/au)^gamma)^zxp  - NOT IMPLEMENTED!
                {
                    vector<float> _pwrlaw_coeff(5);
                    Normalization[cmp] = 1.0; // Defined functions already return per H
                    _pwrlaw_coeff[0] = _mcf.FValue(_thisSection, "A");
                    _pwrlaw_coeff[1] = _mcf.FValue(_thisSection, "alpha");
                    _pwrlaw_coeff[2] =
                        _mcf.FValue(_thisSection, "at") * Constant::UM_CM; // at and ac are defined in um in input file.
                    _pwrlaw_coeff[3] = _mcf.FValue(_thisSection, "ac") * Constant::UM_CM;
                    _pwrlaw_coeff[4] = _mcf.FValue(_thisSection, "gamma");
                    SizeDistribution[cmp] = getPwrLaw_sdist(_pwrlaw_coeff, cmp);
                }
                break;

            case 3: // Gaussian
            {
                vector<float> _gauss_coeff(3);
                // Get input dust to gas ratio and mean molecular weight
                // Shouldn't have to do this for every component, but I'm lazy...
                DustToGasMassRatio = _mcf.FValue("Model", "Dust to Gas Mass Ratio");
                if (_mcf.isBadFloat(DustToGasMassRatio))
                {
                    cout << "Dust to gas mass ratio not defined - required for this size "
                            "distribution"
                         << endl;
                    cout << "Looked under section 'Model' for key 'Dust to Gas Mass Ratio'" << endl;
                    exit(8);
                }
                MeanMolecularWeight = _mcf.FValue("Model", "Mean Molecular Weight");
                if (_mcf.isBadFloat(MeanMolecularWeight))
                {
                    cout << "Mean molecular weight not defined - required for this size "
                            "distribution"
                         << endl;
                    cout << "Looked under section 'Model' for key 'Mean Molecular Weight'" << endl;
                    exit(8);
                }

                _gauss_coeff[0] = 1; // Force amplitude to be 1 for now

                // Convert all length units in gaussian parameters to CM - that way, the
                // size distribution will be returned in cm^(-1) H^(-1)
                _gauss_coeff[1] = _mcf.FValue(_thisSection, "a0") * Constant::UM_CM;
                _gauss_coeff[2] = _mcf.FValue(_thisSection, "as") * Constant::UM_CM;
                SizeDistribution[cmp] = getGauss_sdist(_gauss_coeff, cmp);

                // Now renormalize to assure proper dust mass/tau/etc computations
                _integrand.resize(Component[cmp].nsize);

                // Need to keep track so that we don't exceed 1 for mass frac.
                float _massfrac = _mcf.FValue(_thisSection, "X");
                _totalMassFraction += _massfrac;
                for (int sz = 0; sz < Component[cmp].nsize; sz++)
                    _integrand[sz] = Component[cmp].mass[sz] * SizeDistribution[cmp][sz];

                // size is in cm, integrand is in grams, _mmw is dimensionless, gdmris
                // dimensionless, AMU_CGS is gm/H, the dimensions os _newAmp are cm^-1
                // H^-1 ... float _newAmp =
                // (DustToGasMassRatio*MeanMolecularWeight*Constant::AMU_CGS)/NumUtils::integrate<float>(Component[cmp].size,_integrand);
                // ... as are the dim. of SizeDistribution  since SizeDistribution was
                // previously dimensionless = 1.0 * exp(-z^2/2).
                Normalization[cmp] = (DustToGasMassRatio * MeanMolecularWeight * Constant::AMU_CGS * _massfrac) /
                                     NumUtils::integrate<float>(Component[cmp].size, _integrand);
            }
            break;

            case 4: // Hensley & Draine 2023, ApJ, 948, 55
            {
                vector<float> _hd_coeff(11);
                Normalization[cmp] = 1.0; // Defined functions already return per H
                _hd_coeff[0] = _mcf.FValue(_thisSection, "B1");
                _hd_coeff[1] = _mcf.FValue(_thisSection, "B2");
                _hd_coeff[2] = _mcf.FValue(_thisSection, "a1");
                _hd_coeff[3] = _mcf.FValue(_thisSection, "a2");
                _hd_coeff[4] = _mcf.FValue(_thisSection, "sig");
                _hd_coeff[5] = _mcf.FValue(_thisSection, "A0");
                _hd_coeff[6] = _mcf.FValue(_thisSection, "A1");
                _hd_coeff[7] = _mcf.FValue(_thisSection, "A2");
                _hd_coeff[8] = _mcf.FValue(_thisSection, "A3");
                _hd_coeff[9] = _mcf.FValue(_thisSection, "A4");
                _hd_coeff[10] = _mcf.FValue(_thisSection, "A5");
                SizeDistribution[cmp] = getHD_sdist(_hd_coeff, cmp);
            }
            break;

            case 5: // Log-Normal distribution e.g. THEMIS-2 implement only necessary
                    // coefficients for THEMIS
                    // - add terms as necessary for other forms.  We will need to
                    // update model definition files for e.g. 'incomplete' (less
                    // complex) log-normal when we add terms A*(1/size)*exp(-1/2
                    // *[ln(size/a0)/sigma]^2)
            {
                vector<float> _logn_coeff(3);
                Normalization[cmp] = 1.0; // Defined functions already return per H
                _logn_coeff[0] = _mcf.FValue(_thisSection, "A");
                _logn_coeff[1] = _mcf.FValue(_thisSection, "a0") * Constant::UM_CM; // a0 is in um in the config file.
                _logn_coeff[2] = _mcf.FValue(_thisSection, "sigma");
                SizeDistribution[cmp] = getLogN_sdist(_logn_coeff, cmp);
            }
            break;

            default: {
                cout << endl;
                cout << "****************************************************" << endl;
                cout << "Ooops! Something BROKE defining size distribution   " << endl;
                cout << "trying to set up " << _thisCrossSectionFile << "!   " << endl;
                cout << "Key " << _thisSizeDistribution << " exists, but has no corresponding case." << endl;
                cout << "(Looking for " << SizeDistID[_thisSizeDistribution] << ")" << endl;
                cout << "This is PROGRAMMER ERROR." << endl;
                cout << "****************************************************" << endl;
                exit(8);
            }
            break;
            }
        }
        else
        { // Encountered a SizeDist key that is not defined in the map.

            cout << endl;
            cout << "****************************************************" << endl;
            cout << _thisSizeDistribution << " is not a defined size distribution type." << endl;
            cout << "Allowed size distribution types are: " << endl;
            for (iMap = SizeDistID.begin(); iMap != SizeDistID.end(); iMap++)
                cout << iMap->first << " " << iMap->second << endl;
            cout << "****************************************************" << endl;
            exit(8);
        }

        // Compute this components contribution to the total mass
        _integrand.resize(Component[cmp].nsize);
        for (int sz = 0; sz < Component[cmp].nsize; sz++)
            _integrand[sz] = Component[cmp].mass[sz] * SizeDistribution[cmp][sz];
        DustMass[cmp] = Normalization[cmp] * NumUtils::integrate<float>(Component[cmp].size, _integrand);
        TotalDustMass += DustMass[cmp];
        _integrand.erase(_integrand.begin(), _integrand.end());

        // Compute the 'volume' in this component; we could roll into the loop above
        // to save computation at the cost of additional RAM, but computational cost
        // at this stage of the process is not a large factor, RAM may be.
        for (int sz = 0; sz < Component[cmp].nsize; sz++)
            _integrand[sz] = Constant::VFAC * pow(Component[cmp].size[sz], 3) * SizeDistribution[cmp][sz];
        DustVolume[cmp] = Normalization[cmp] * NumUtils::integrate<float>(Component[cmp].size, _integrand);

        number[cmp] = Normalization[cmp] * NumUtils::integrate<float>(Component[cmp].size, SizeDistribution[cmp]);
        totalnumber += number[cmp];

        // Generate the normalized size distribution to avoid computational loops
        // everytime it's called
        transform(SizeDistribution[cmp].begin(), SizeDistribution[cmp].end(), back_inserter(SizeDistributionNorm[cmp]),
                  bind2nd(multiplies<float>(), Normalization[cmp]));

    } // End Component loop.

    // Total normalization; this is the denominator for C_abs,sca size integrated,
    // comp summed quantities.
    TotalNormalization = accumulate(Normalization.begin(), Normalization.end(), 0.0);

    // Now construct the Effective grain properties.  Separate out from the size
    // distribution definition so that the total normalization can be applied
    // during the effective property construction rather than as a new  nwave
    // loop.

    iTau = Tau.begin();
    iAlbedo = Albedo.begin();
    iCAbsEff = CAbsEff.begin();
    iCScaEff = CScaEff.begin();
    iphFuncEff = phFuncEff.begin();

    for (int wv = 0; wv < nWave; wv++)
    { // Loop over wavelength

        for (int cmp = 0; cmp < nComp; cmp++)
        { // Loop over Chemical composition

            // Get size as well as Cs and g as f(a)
            //_thisSize = Component[cmp].getSize();         // size grid.
            _thisCAbs = Component[cmp].w_getCAbs(wv);     // CAbs as f(a)
            _thisCSca = Component[cmp].w_getCSca(wv);     // CSca as f(a)
            _thisphFunc = Component[cmp].w_getphFunc(wv); // g as f(a)

            // Initialize integrals
            _CAbsIntegral = 0.0;
            _CScaIntegral = 0.0;
            _phFuncIntegral = 0.0;

            // Point iterators to begining of vectors + 1 for integration over a.
            // Rather than call integrate 3 times and do the size loop 3 times,
            // roll the integration into the loop.  Doesn't allow flexibility in
            // specifying the integrator, but...
            ithisCAbs = _thisCAbs.begin() + 1;                 // begining of CAbs+1
            ithisCSca = _thisCSca.begin() + 1;                 // begining of Csca+1
            ithisphFunc = _thisphFunc.begin() + 1;             // begining of g+1
            ithisSizeDist = SizeDistribution[cmp].begin() + 1; // beginning of SizeDistribution+1

            // for (ithisSize =
            // _thisSize.begin()+1;ithisSize!=_thisSize.end();ithisSize++) { //
            // Integration loop.
            for (ithisSize = Component[cmp].size.begin() + 1; ithisSize != Component[cmp].size.end(); ithisSize++)
            { // Integration loop.
                delta = 0.5 * ((*ithisSize) - *(ithisSize - 1));
                // C(abs,sca) = int(C(a)*f(a)da)
                _CAbsIntegral += delta * ((*ithisCAbs) * (*ithisSizeDist) + *(ithisCAbs - 1) * (*(ithisSizeDist - 1)));
                _CScaIntegral += delta * ((*ithisCSca) * (*ithisSizeDist) + *(ithisCSca - 1) * (*(ithisSizeDist - 1)));
                // g = int(g(a)*Csca(a)*f(a)da
                _phFuncIntegral += delta * ((*ithisphFunc) * (*ithisCSca) * (*ithisSizeDist) +
                                            (*(ithisphFunc - 1)) * (*(ithisCSca - 1)) * (*(ithisSizeDist - 1)));

                // Increment non-loop iterators to next position.
                ithisCAbs++;
                ithisCSca++;
                ithisphFunc++;
                ithisSizeDist++;
            }

            // Multiply by the appropriate normalization integral (eg. now we have
            // []/Hatom
            *(iCAbsEff) += (_CAbsIntegral * Normalization[cmp]);
            *(iCScaEff) += (_CScaIntegral * Normalization[cmp]);
            *(iphFuncEff) += (_phFuncIntegral * Normalization[cmp]);
        }

        // Overall normalization - this is the normalization that comes from the
        // size integration.
        *(iTau) = (*(iCAbsEff)) + (*(iCScaEff));
        *(iphFuncEff) /= *(iCScaEff);
        // Do the normalization first to make more... reasonable... numbers for
        // albedo computation
        *(iCAbsEff) /= TotalNormalization;
        *(iCScaEff) /= TotalNormalization;
        *(iAlbedo) = *(iCScaEff) / (*(iCAbsEff) + *(iCScaEff));

        // Increment iterator to next member of effective vectors
        iTau++;
        iphFuncEff++;
        iCAbsEff++;
        iCScaEff++;
        iAlbedo++;
    }

    cout << "Done!" << endl;
}

// Return the size distribution in cm^-1 for a ZDA definition.
//  Zubko, Dwek, & Arendt 2004, ApJS, 152, 211
//  Bare grains, PAHs, Graphite and Silicates.
//  Size distributions defined by:
//  log g(a) = c0 + b0*log(a)
//                - b1*|log(a/a1)|^m1
//                - b2*|log(a/a2)|^m2
//                - b3*|a-a3|^m3
//                - b4*|a-a4|^m4
//  with f(a) = A*g(a) and Int_(a_min)^(a_max) g(a) da = 1

// vector <float> GrainModel::getZDA_sdist(vector <float> coeff, vector <float>
// sz)
vector<float> GrainModel::getZDA_sdist(vector<float> coeff, int cmp)
{

    float thissz;                                     // Hold um copy of size at each point.
    vector<float> retvec(Component[cmp].size.size()); // The size distribution [cm^-1]
    vector<float>::iterator isz;                      // Iterator through size vector
    vector<float>::iterator iretvec;                  // Iterator through size dist vector

    // Set iterator to beginning of size dist vector.
    iretvec = retvec.begin();

    // for (isz=sz.begin();isz!=sz.end();isz++) {
    for (isz = Component[cmp].size.begin(); isz != Component[cmp].size.end(); isz++)
    {
        thissz = (*isz) * Constant::CM_UM; // size[i] in [um]
        *iretvec = coeff[1];               // log(g) = c0
        if (!isnan(coeff[2]))
        { // b0 is defined
            *iretvec += coeff[2] * log10(thissz);
            // cout << "b0: " << cmp << " " << thissz << " " << coeff[2]*log10(thissz)
            // << endl;
        }
        if (!isnan(coeff[3]))
        { // b1 is defined
            *iretvec -= coeff[3] * pow(fabs(log10(thissz / coeff[4])), coeff[5]);
            // cout << "b1: " << cmp << " " << thissz << " " << coeff[3] << " " <<
            // coeff[4] << " " << coeff[5] << endl; cout << "    " << cmp << " " <<
            // thissz << " " << pow(fabs(log10(thissz/coeff[4])),coeff[5]) << endl;
        }
        if (!isnan(coeff[6]))
        { // b2 is defined
            *iretvec -= coeff[6] * pow(fabs(log10(thissz / coeff[7])), coeff[8]);
            // cout << "b2: " << cmp << " " << thissz << " " << coeff[6] << " " <<
            // coeff[7] << " " << coeff[8] << endl; cout << "    " << cmp << " " <<
            // thissz << " " << pow(fabs(log10(thissz/coeff[7])),coeff[8]) << endl;
        }
        if (!isnan(coeff[9]))
        { // b3 is defined
            *iretvec -= coeff[9] * pow(fabs((thissz)-coeff[10]), coeff[11]);
            // cout << "b3: " << cmp << " " << thissz << " " <<  coeff[9]*pow(fabs(
            // (thissz) - coeff[10]),coeff[11]) << endl;
        }
        if (!isnan(coeff[12]))
        { // b4 is defined
            *iretvec -= coeff[12] * pow(fabs((thissz)-coeff[13]), coeff[14]);
            // cout << "b4: " << cmp << " " << thissz << " " <<  coeff[12]*pow(fabs(
            // (thissz) - coeff[13]),coeff[14]) << endl;
        }
        // cout << "LogG: " << cmp << " " << thissz << " " << *iretvec << endl;
        *iretvec = pow(10, (*iretvec)); // g(a) in um^-1
        // cout << "G:    " << cmp << " " << thissz << " " << *iretvec << endl;
        //  Zubko et al. functional definitions give g(a) in um^-1, convert to cm.
        *iretvec *= Constant::CM_UM; // g(a) in cm^-1
        iretvec++;
    }

    // Renormalize retvec to 1; ~g(a) is not correctly normalized like g(a) as it
    // is a numerical approximation to the normalized size distribution. This
    // normalization will depend slightly on the choice of size grid, so we can't
    // simply incorporate it into Normalization[]
    // cout << "re-normalizing" <<endl;
    // for (int i=0;i<retvec.size();++i) cout << Component[cmp].size[i] << " " <<
    // retvec[i] << endl;
    float _Nfac = 1.0 / NumUtils::integrate<float>(Component[cmp].size, retvec);
    // for (iretvec=retvec.begin();iretvec!=retvec.end();++iretvec) cout <<
    // *iretvec << endl; cout << "The normalization was: " << _Nfac << endl <<
    // endl;
    transform(retvec.begin(), retvec.end(), retvec.begin(), bind2nd(multiplies<float>(), _Nfac));

    return retvec;
}

// Return the size distribution in cm^-1 H^-1
//  Wiengartner & Draine 2001, ApJ 548, 296
vector<float> GrainModel::getWD_sdist(vector<float> coeff, int cmp)
{

    // float thissz;                     // Hold um copy of size at each point.
    vector<float> retvec(Component[cmp].size.size()); // The size distribution [cm^-1]
    vector<float>::iterator isz;                      // Iterator through size vector
    vector<float>::iterator iretvec;                  // Iterator through size dist vector

    // Set iterator to beginning of size dist vector.
    iretvec = retvec.begin();

    // Wiengartner & Draine size distribution functions have mixed units - convert
    // where necessary to ensure out put in cm^-1 H^-1
    for (isz = Component[cmp].size.begin(); isz != Component[cmp].size.end(); isz++)
    {
        // Compute the common part of the size distributions.
        *iretvec = coeff[3] / (*isz) * pow((*isz / (coeff[4] * static_cast<float>(Constant::UM_CM))), coeff[6]);
        if (coeff[7] >= 0.0)
            *iretvec *= (1.0 + coeff[7] * (*isz) / (coeff[4] * Constant::UM_CM));
        else
            *iretvec /= (1.0 - coeff[7] * (*isz) / (coeff[4] * Constant::UM_CM));
        if ((*isz) > (coeff[4] * Constant::UM_CM))
        {
            *iretvec *= exp(-1 * pow(((*isz) - coeff[4] * Constant::UM_CM) / (coeff[5] * Constant::UM_CM), 3));
        }
        if (coeff[0] > 0 && !isnan(coeff[0]))
        {
            *iretvec +=
                coeff[1] * coeff[0] / (*isz) * exp(-0.5 * (pow((log((*isz) / (3.5 * Constant::ANG_CM)) / 0.4), 2)));
            *iretvec +=
                coeff[2] * coeff[0] / (*isz) * exp(-0.5 * (pow((log((*isz) / (30.0 * Constant::ANG_CM)) / 0.4), 2)));
        }
        iretvec++;
    }

    return retvec;
}

// Return the size distribution in cm^-1 H^-1
//  Hensley & Draine 2023, ApJ, 948, 55
vector<float> GrainModel::getHD_sdist(vector<float> coeff, int cmp)
{

    vector<float> retvec(Component[cmp].size.size()); // The size distribution [cm^-1]
    vector<float>::iterator isz;                      // Iterator through size vector
    vector<float>::iterator iretvec;                  // Iterator through size dist vector

    // the size grid Component[cmp].size() should be in CM, so just make sure
    // anywere we ratio to other sizes, we are consistent.  Since these
    // distributions are of the form FUNCTION/size and size is in CM, we are
    // returning cm^-1 H^-1
    iretvec = retvec.begin();

    // Hensley & Draine size distribution functions have mixed units (because of
    // course they do)
    // - convert where necessary to ensure out put in cm^-1 H^-1.  Since some of
    // the random selection of units occur in powers of natural logs in sum of
    // expoentials, not trivial to absord into a redefinition of constants.
    float _sig = 2. * pow(coeff[4], 2);
    // cout << cmp << ": " << Component[cmp].size.size() << endl;
    for (isz = Component[cmp].size.begin(); isz != Component[cmp].size.end(); isz++, iretvec++)
    {
        *iretvec = (coeff[0] / (*isz)) * exp((-1) * (pow(log((*isz) / (coeff[2] * Constant::UM_CM)), 2) / _sig));
        // cout << *isz << " " << *iretvec << endl;
        if (coeff[1] != 0)
            *iretvec += (coeff[1] / (*isz)) * exp((-1) * (pow(log((*isz) / (coeff[3] * Constant::UM_CM)), 2) / _sig));
        // cout << *isz << " " << *iretvec << endl;
        if (coeff[5] != 0)
        {
            // cout << "ADDING SUM OF EXP" << endl;
            float _sumofexp = 0;
            // cout << *isz << " " << *iretvec << endl;
            for (int _i = 6; _i < 11; ++_i)
            {
                // -5 on the power so we have powers of 1-5 (not 6-11)
                // cout << _i << " " << coeff[_i] << " " << _i-5 << endl;
                _sumofexp += coeff[_i] * pow(log((*isz) / Constant::ANG_CM), _i - 5);
            }
            // cout << coeff[5] << " " << exp(_sumofexp) << endl << endl;
            *iretvec += (coeff[5] / (*isz)) * exp(_sumofexp);
            // cout << *isz << " " << *iretvec << endl;
        }
    }

    return retvec;
}

// Return size distribution in cm^-1 H^-1
//   Power-law type distribution (e.g. specifically added for THEMIS2, Ysard et
//   al 2024 preprint
//                // POWER LAW definition:  A (normalization), alpha (slope of
//                power law), at, ac, gamma (for a 1/n_H dn/da = A * sz^alpha *
//                cutoff * curvature
//                // cutoff: exp(-((a-at)/aC)^gamma)
//                // curvature: (1 + zeta*(a/au)^gamma)^zxp  - NOT IMPLEMENTED!
//          _pwrlaw_coeff[0] = _mcf.FValue(_thisSection,"A");
//          _pwrlaw_coeff[1] = _mcf.FValue(_thisSection,"alpha");
//          _pwrlaw_coeff[2] = _mcf.FValue(_thisSection,"at");
//          _pwrlaw_coeff[3] = _mcf.FValue(_thisSection,"ac");
//          _pwrlaw_coeff[4] = _mcf.FValue(_thisSection,"gamma");

vector<float> GrainModel::getPwrLaw_sdist(vector<float> coeff, int cmp)
{
    vector<float> retvec(Component[cmp].size.size()); // The size distribution [cm^-1]
    vector<float>::iterator isz;                      // Iterator through size vector
    vector<float>::iterator iretvec;                  // Iterator through size dist vector

    // Set iterator to beginning of size dist vector.
    iretvec = retvec.begin();

    // We assume coefficients are in cm where appropriate - either defined as such
    // in the model definition file or converted as necessary in the read code.
    for (isz = Component[cmp].size.begin(); isz != Component[cmp].size.end(); isz++, iretvec++)
    {
        *iretvec = coeff[0] * pow(*isz, coeff[1]); // this is the power law
        if (*isz >= coeff[3])
        { // this is the exponential cutoff
            *iretvec *= (exp((-1.0) * pow((*isz - coeff[2]) / coeff[3], coeff[4])));
        }
    }
    // Test this! But we don't need an extra 'loop' even if 'optimized
    // transform(retvec.begin(),retvec.end(),retvec.begin(), [&](float elm)
    // {return elm *= (coeff[0]); });
    return retvec;
}

// Retrun size distribution in cm^-1 H^-1 (or whatever you normalization or lack
// thereof isj from "A" constant Log-Normal distribution e.g. THEMIS-2 implement
// only necessary coefficients for THEMIS, Ysard et al 2024 preprint
//         Log-Normal defintion: 1/n_H dn/da = A*(1/size)*exp(-1/2
//         *[ln(size/a0)/sigma]^2)
//             _logn_coeff[0] = _mcf.FValue(_thisSection,"A");
//             _logn_coeff[1] = _mcf.FValue(_thisSection,"a0");
//             _logn_coeff[2] = _mcf.FValue(_thisSection,"sigma");
//
// - add terms as necessary for other forms.  Will need to update Model
// definition files
vector<float> GrainModel::getLogN_sdist(vector<float> coeff, int cmp)
{
    vector<float> retvec(Component[cmp].size.size()); // The size distribution [cm^-1]
    vector<float>::iterator isz;                      // Iterator through size vector
    vector<float>::iterator iretvec;                  // Iterator through size dist vector

    float _temp;

    // Set iterator to beginning of size dist vector.
    iretvec = retvec.begin();

    // We assume coefficients are in cm where appropriate - either defined as such
    // in the model definition file or converted as necessary in the read code.
    for (isz = Component[cmp].size.begin(); isz != Component[cmp].size.end(); isz++, iretvec++)
    {
        *iretvec = coeff[0] / (*isz) * exp((-0.5) * pow(log(*isz / coeff[1]) / coeff[2], 2));
    }
    // Test this! But we don't need an extra 'loop' even if 'optimized
    // transform(retvec.begin(),retvec.end(),retvec.begin(), [&](float elm)
    // {return elm *= (coeff[0]); });
    return retvec;
}

// Return Gaussian size distribution in cm^-1
vector<float> GrainModel::getGauss_sdist(vector<float> coeff, int cmp)
{

    vector<float> retvec(Component[cmp].size.size()); // The size distribution [cm^-1]
    vector<float>::iterator isz;                      // Iterator through size vector
    vector<float>::iterator isze;                     // end of size vector
    vector<float>::iterator iretvec;                  // Iterator through size dist vector

    float z;

    isz = Component[cmp].size.begin();
    isze = Component[cmp].size.end();
    iretvec = retvec.begin();

    for (; isz != isze; ++isz, ++iretvec)
    {
        z = (*isz - coeff[1]) / coeff[2];
        *iretvec = coeff[0] * exp(-pow(z, 2) / 2);
    }

    return retvec;
}

float GrainModel::getTau(float a_wave)
{

    if (a_wave < Wave[0] || a_wave > Wave[nWave - 1])
    {
        cout << "You are attempting to retrieve \tau_(lam) outside of the defined" << endl;
        cout << "wavelength range. " << endl;
        cout << "Error in Grain::getCAbs(float)." << endl;
        cout << "Defined wavelength range is " << Wave[0] << "--" << Wave[nWave - 1] << endl;
        cout << "Requested wavelength: " << a_wave << endl;
        exit(8);
    }

    typedef vector<float>::iterator wvIter;
    wvIter it_beg = Wave.begin();
    wvIter it_loc = find(Wave.begin(), Wave.end(), a_wave);
    int idx = it_loc - it_beg;

    if (it_loc != Wave.end())
    { // Found an exact match.
        return Tau[idx];
    }
    else
    { // weighted average

        it_loc = lower_bound(Wave.begin(), Wave.end(), a_wave);
        idx = it_loc - it_beg;

        // cout << Wave[idx-1] << " " << a_wave << " " << Wave[idx] << endl;
        // float slope=(Tau[idx]-Tau[idx-1])/(Wave[idx]-Wave[idx-1]);
        // float D1 = 1.0/abs(Wave[idx] - a_wave);
        float D1 = 1.0 / (pow((Wave[idx] - a_wave), 2));
        float D2 = 1.0 / (pow((a_wave - Wave[idx - 1]), 2));
        // float D2 = 1.0/abs(a_wave - Wave[idx-1]);
        // cout << D1 << " " << " " << D2 << endl;
        // cout << Tau[idx-1] << " " << Tau[idx] << endl;
        // cout << D1+D2 << endl;
        // cout << Tau[idx]*D1 << " " << Tau[idx-1]*D2 << endl;
        return ((Tau[idx] * D1 + Tau[idx - 1] * D2) / (D1 + D2));
        // cout << Tau[idx]+slope*(a_wave-Wave[idx]) << " " <<
        // Tau[idx-1]+slope*(a_wave-Wave[idx-1]) << endl; return
        // Tau[idx]+slope*(a_wave-Wave[idx]);
    }
}

// void GrainModel::ComputeEmission ( vector <float> rField )
// {

//   vector <float> _t;
//   vector <float>::iterator _it;
//   float _tlo,_thi;

//   for (int _cmp=0;_cmp<nComp;++_cmp) { // Component loop

//     Temperature[_cmp].resize(Component[_cmp].nsize);
//     _it=Temperature[_cmp].begin();
//     *_it = EqTemp(Wave, rField, Component[_cmp].CAbs[0]);
//     ++_it;
//     for (int _sz=1;_sz<Component[_cmp].nsize;++_sz,++_it) {
//       _tlo = 0.3*(*(_it-1));
//       _thi = 4.0*(*(_it-1));
//       *_it = EqTemp(Wave, rField,
//       Component[_cmp].CAbs[_sz],((_tlo<0)?1.0:_tlo),((_thi>2500.0)?2500.0:_thi));
//     }

//     // Logic for Equilibrium vs. non-Equilibrium computation.
// //     for (int sz=0;sz<Component[_cmp].nsize;++sz) {
// //       cout << Component[_cmp].size[sz] << " " << Temperature[_cmp][sz] <<
// endl;
// //     }
//   }

// }
