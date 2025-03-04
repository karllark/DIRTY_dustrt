//*******************************************************************************
// GrainModel -- GrainModel Class
//
// Requires:
//     ConfigFile
//     Grain
//     StringManip
//
// History:
//     Written by:   Putrid, Winter 2006/2007
//     Contact: misselt@as.arizona.edu
//
//    Added a few informational functions/capabilities (KDG - Oct 2014)
//       like _thisComponentName
//
//*******************************************************************************

#include "GrainModel.h"

GrainModel::GrainModel() {}

void GrainModel::MakeGrainModel(ConfigFile &_mbkcf, vector<float> &MasterWave) {

  string _thisSection;
  string _line;

  string _thisComponentName;
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
  _TopLevelPath =
      _mbkcf.SValue("Model Book Keeping", "Path to Dust Properties");
  // CrossSections:
  _CrossSectionsPath =
      _mbkcf.SValue("Model Book Keeping", "Cross Section SubDir");
  // Calorimetry:
  _CalorimetryPath = _mbkcf.SValue("Model Book Keeping", "Calorimetry SubDir");
  // Model Definiiton:
  _ModelDefinitionPath = _mbkcf.SValue("Model Book Keeping", "Model SubDir");
  // Model Name:
  ModelName = _mbkcf.SValue("Model Book Keeping", "Model Name");

  // Check that model name exists in pre-defined map.
  if (ModelID.find(ModelName) == ModelID.end()) { // NO!
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
  if (!StringManip::FileExists(_FullModelDefinitionFile)) { // NO!
    cout << "****************************************************" << endl;
    cout << "Model Definition File " << _FullModelDefinitionFile
         << " not found." << endl;
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
  for (int cmp = 0; cmp < nComp; cmp++) { // Loop over all components

    // Construct section key for this component.
    // Assume key is of the form "Comp i" where i=1,2,...,nComp
    _thisSection = "Component " + StringManip::vtos(cmp + 1);

    // Get the name of this component
    _thisComponentName = _mcf.SValue(_thisSection, "Component Name");

    // Get the Size Distribution information.  Map in GrainModel.h is in
    // upper case, so convert to upper case.
    _thisSizeDistribution = _mcf.SValue(_thisSection, "Size Distribution");
    transform(_thisSizeDistribution.begin(), _thisSizeDistribution.end(),
              _thisSizeDistribution.begin(), (int (*)(int))toupper);

    // Clean out _MasterSize
    _MasterSize.erase(_MasterSize.begin(), _MasterSize.end());

    // Get grain property definitions.
    _thisCrossSectionFile = _mcf.SValue(_thisSection, "Cross Sections");
    _thisCalorimetryFile = _mcf.SValue(_thisSection, "Calorimetry");

    _thisA_min = _mcf.FValue(_thisSection, "a_min") * Constant::UM_CM;
    _thisA_max = _mcf.FValue(_thisSection, "a_max") * Constant::UM_CM;

    // For a Gaussian (or arbitrary in the future?) size distribution, check
    // that a_min and _amax are correctly defined.
    if (SizeDistID.find(_thisSizeDistribution) != SizeDistID.end()) {

      if (_thisSizeDistribution ==
          "GAUSS") { // Gaussian - add other types as necessary

        float _a0 = _mcf.FValue(_thisSection, "a0");
        if (_mcf.isBadFloat(_a0)) {
          cout << endl;
          cout << "************************************************" << endl;
          cout << "Did not find key a0 in section " << _thisSection << endl;
          cout << "a0 is required for GAUSSIAN size distribution." << endl;
          cout << "************************************************" << endl;
          exit(8);
        }
        float _as = _mcf.FValue(_thisSection, "as");
        if (_mcf.isBadFloat(_as)) {
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
      }
    } else { // Encountered a SizeDist key that is not defined in the map.

      cout << endl;
      cout << "****************************************************" << endl;
      cout << _thisSizeDistribution
           << " is not a defined size distribution type." << endl;
      cout << "Allowed size distribution types are: " << endl;
      for (iMap = SizeDistID.begin(); iMap != SizeDistID.end(); iMap++)
        cout << iMap->first << " " << iMap->second << endl;
      cout << "****************************************************" << endl;
      exit(8);
    }
    _thisSizeType = _mcf.SValue(_thisSection, "Size Definition Type");
    transform(_thisSizeType.begin(), _thisSizeType.end(), _thisSizeType.begin(),
              (int (*)(int))toupper);
    switch (SizeTypeID[_thisSizeType]) {
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
      if (_thisSizeFile == _mcf.BadString()) {
        cout << "Unable to retrieve size grid file name, resorting to default."
             << endl;
        _MasterSize.push_back(-1);
      } else {
        // ifstream
        // _file((_TopLevelPath+_ModelDefinitionPath+_thisSizeFile).c_str());
        ifstream _file((_thisSizeFile).c_str());
        if (!_file.is_open()) {
          cout << "Unable to open size grid definition file " << _thisSizeFile
               << endl;
          exit(8);
        }
        while (getline(_file, _line)) {
          if (_line[0] == '#')
            continue;
          else {
            _MasterSize.push_back(atof(_line.c_str()) * Constant::UM_CM);
          }
        }
        _MasterSize.erase(remove_if(_MasterSize.begin(), _MasterSize.end(),
                                    bind2nd(less<float>(), _thisA_min)),
                          _MasterSize.end());
        _MasterSize.erase(remove_if(_MasterSize.begin(), _MasterSize.end(),
                                    bind2nd(greater<float>(), _thisA_max)),
                          _MasterSize.end());
        if (_MasterSize[0] > _thisA_min)
          _MasterSize.insert(_MasterSize.begin(), _thisA_min);
        if (_MasterSize[_MasterSize.size() - 1] < _thisA_max)
          _MasterSize.push_back(_thisA_max);
        // for (int _sz=0;_sz<_MasterSize.size();++_sz) cout << _MasterSize[_sz]
        // << endl;
      }
      break;
    }
    case 2: {
      // Get the number of sizes
      _thisSizeIN = _mcf.IValue(_thisSection, "Size Definition Number");
      // cout << _thisSizeIN << endl;
      if (_thisSizeIN <= 0) {
        cout << "Unable to define size grid manually.  Resorting to default."
             << endl;
        _thisSizeType = "DEF";
        _MasterSize.push_back(-1);
      } else {
        _MasterSize.resize(_thisSizeIN);
        // log/linear
        _thisSizeScale = _mcf.SValue(_thisSection, "Size Definition Scale");
        if (_thisSizeScale == _mcf.BadString()) {
          cout << "Did not find size gridding definition - assuming linear"
               << endl;
          _thisSizeScale = "LINEAR";
        } else {
          transform(_thisSizeScale.begin(), _thisSizeScale.end(),
                    _thisSizeScale.begin(), (int (*)(int))toupper);
        }
        if (_thisSizeScale == "LINEAR") { // Generate Linear scale
          _MasterSize[0] = _thisA_min;
          _delta =
              (_thisA_max - _thisA_min) / static_cast<float>(_thisSizeIN - 1);

          for (int _sz = 1; _sz < _thisSizeIN; ++_sz)
            _MasterSize[_sz] = _MasterSize[_sz - 1] + _delta;
        }
        if (_thisSizeScale == "LOG") { // Generate Log scale.
          _MasterSize[0] = _thisA_min;
          _delta = (1.0 / static_cast<float>(_thisSizeIN - 1)) *
                   log10(_thisA_max / _thisA_min);
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
    Component[cmp].MakeGrain(
        _thisComponentName, _CrossSectionsPath + _thisCrossSectionFile,
        _CalorimetryPath + _thisCalorimetryFile, MasterWave, _MasterSize,
        _TopLevelPath, _thisA_min, _thisA_max);

    SizeDistribution[cmp].resize(Component[cmp].nsize);

    // Wrap in a .find since SizeDist[key] will add members.
    if (SizeDistID.find(_thisSizeDistribution) !=
        SizeDistID.end()) { // Determine which sd to use

      switch (SizeDistID[_thisSizeDistribution]) {

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
        // SizeDistribution[cmp].resize(_thisSize.size());
        SizeDistribution[cmp] = getZDA_sdist(_zda_coeff, cmp);
      } break;

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
      } break;

      case 2: // Power law
      {
        vector<float> _pl_coeff(2);
      }

      case 3: // Gaussian
      {
        vector<float> _gauss_coeff(3);
        // Get input dust to gas ratio and mean molecular weight
        // Shouldn't have to do this for every component, but I'm lazy...
        DustToGasMassRatio = _mcf.FValue("Model", "Dust to Gas Mass Ratio");
        if (_mcf.isBadFloat(DustToGasMassRatio)) {
          cout << "Dust to gas mass ratio not defined - required for this size "
                  "distribution"
               << endl;
          cout
              << "Looked under section 'Model' for key 'Dust to Gas Mass Ratio'"
              << endl;
          exit(8);
        }
        MeanMolecularWeight = _mcf.FValue("Model", "Mean Molecular Weight");
        if (_mcf.isBadFloat(MeanMolecularWeight)) {
          cout << "Mean molecular weight not defined - required for this size "
                  "distribution"
               << endl;
          cout << "Looked under section 'Model' for key 'Mean Molecular Weight'"
               << endl;
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
        Normalization[cmp] =
            (DustToGasMassRatio * MeanMolecularWeight * Constant::AMU_CGS *
             _massfrac) /
            NumUtils::integrate<float>(Component[cmp].size, _integrand);
        // cout << Normalization[cmp] << endl;
        //  transform(SizeDistribution[cmp].begin(),SizeDistribution[cmp].end(),SizeDistribution[cmp].begin(),
        // 		      bind2nd(multiplies<float>(),_newAmp));
      } break;

      default: {
        cout << endl;
        cout << "****************************************************" << endl;
        cout << "Ooops! Something BROKE defining size distribution   " << endl;
        cout << "trying to set up " << _thisCrossSectionFile << "!   " << endl;
        cout << "Key " << _thisSizeDistribution
             << " exists, but has no corresponding case." << endl;
        cout << "(Looking for " << SizeDistID[_thisSizeDistribution] << ")"
             << endl;
        cout << "This is PROGRAMMER ERROR." << endl;
        cout << "****************************************************" << endl;
        exit(8);
      } break;
      }

    } else { // Encountered a SizeDist key that is not defined in the map.

      cout << endl;
      cout << "****************************************************" << endl;
      cout << _thisSizeDistribution
           << " is not a defined size distribution type." << endl;
      cout << "Allowed size distribution types are: " << endl;
      for (iMap = SizeDistID.begin(); iMap != SizeDistID.end(); iMap++)
        cout << iMap->first << " " << iMap->second << endl;
      cout << "****************************************************" << endl;
      exit(8);
    }

    // Compute this components contribution to the total mass.
    _integrand.resize(Component[cmp].nsize);
    for (int sz = 0; sz < Component[cmp].nsize; sz++)
      _integrand[sz] = Component[cmp].mass[sz] * SizeDistribution[cmp][sz];
    // cout << Normalization[cmp] << endl;
    DustMass[cmp] = Normalization[cmp] *
                    NumUtils::integrate<float>(Component[cmp].size, _integrand);
    TotalDustMass += DustMass[cmp];
    _integrand.erase(_integrand.begin(), _integrand.end());

    number[cmp] =
        Normalization[cmp] *
        NumUtils::integrate<float>(Component[cmp].size, SizeDistribution[cmp]);
    totalnumber += number[cmp];

    // for (int sz=0;sz<Component[cmp].nsize;++sz)
    //   cout << sz << " " << Component[cmp].size[sz] << " " <<
    //   SizeDistribution[cmp][sz] << endl;

    // Generate the normalized size distribution to avoid computational loops
    // everytime it's called
    transform(SizeDistribution[cmp].begin(), SizeDistribution[cmp].end(),
              back_inserter(SizeDistributionNorm[cmp]),
              bind2nd(multiplies<float>(), Normalization[cmp]));

  } // End Component loop.

  // Total normalization; this is the denominator for C_abs,sca size integrated,
  // comp summed quantities.
  TotalNormalization =
      accumulate(Normalization.begin(), Normalization.end(), 0.0);

  // Now construct the Effective grain properties.  Separate out from the size
  // distribution definition so that the total normalization can be applied
  // during the effective property construction rather than as a new  nwave
  // loop.
  // cout << "Computing size averaged quantities....";

  iTau = Tau.begin();
  iAlbedo = Albedo.begin();
  iCAbsEff = CAbsEff.begin();
  iCScaEff = CScaEff.begin();
  iphFuncEff = phFuncEff.begin();

  for (int wv = 0; wv < nWave; wv++) { // Loop over wavelength

    for (int cmp = 0; cmp < nComp; cmp++) { // Loop over Chemical composition

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
      ithisCAbs = _thisCAbs.begin() + 1;     // begining of CAbs+1
      ithisCSca = _thisCSca.begin() + 1;     // begining of Csca+1
      ithisphFunc = _thisphFunc.begin() + 1; // begining of g+1
      ithisSizeDist =
          SizeDistribution[cmp].begin() + 1; // beginning of SizeDistribution+1

      // for (ithisSize =
      // _thisSize.begin()+1;ithisSize!=_thisSize.end();ithisSize++) { //
      // Integration loop.
      for (ithisSize = Component[cmp].size.begin() + 1;
           ithisSize != Component[cmp].size.end();
           ithisSize++) { // Integration loop.
        // 1/2 dx
        delta = 0.5 * ((*ithisSize) - *(ithisSize - 1));
        // C(abs,sca) = int(C(a)*f(a)da)
        _CAbsIntegral += delta * ((*ithisCAbs) * (*ithisSizeDist) +
                                  *(ithisCAbs - 1) * (*(ithisSizeDist - 1)));
        _CScaIntegral += delta * ((*ithisCSca) * (*ithisSizeDist) +
                                  *(ithisCSca - 1) * (*(ithisSizeDist - 1)));
        // g = int(g(a)*Csca(a)*f(a)da
        _phFuncIntegral +=
            delta * ((*ithisphFunc) * (*ithisCSca) * (*ithisSizeDist) +
                     (*(ithisphFunc - 1)) * (*(ithisCSca - 1)) *
                         (*(ithisSizeDist - 1)));

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
vector<float> GrainModel::getZDA_sdist(vector<float> coeff, int cmp) {

  float thissz; // Hold um copy of size at each point.
  vector<float> retvec(
      Component[cmp].size.size()); // The size distribution [cm^-1]
  vector<float>::iterator isz;     // Iterator through size vector
  vector<float>::iterator iretvec; // Iterator through size dist vector

  // Set iterator to beginning of size dist vector.
  iretvec = retvec.begin();

  // for (isz=sz.begin();isz!=sz.end();isz++) {
  for (isz = Component[cmp].size.begin(); isz != Component[cmp].size.end();
       isz++) {
    thissz = (*isz) * Constant::CM_UM; // size[i] in [um]
    *iretvec = coeff[1];               // log(g) = c0
    if (!isnan(coeff[2]))              // b0 is defined
      *iretvec += coeff[2] * log10(thissz);
    if (!isnan(coeff[3])) // b1 is defined
      *iretvec -= coeff[3] * pow(fabs(log10(thissz / coeff[4])), coeff[5]);
    if (!isnan(coeff[6])) // b2 is defined
      *iretvec -= coeff[6] * pow(fabs(log10(thissz / coeff[7])), coeff[8]);
    if (!isnan(coeff[9])) // b3 is defined
      *iretvec -= coeff[9] * pow(fabs((thissz)-coeff[10]), coeff[11]);
    if (!isnan(coeff[12])) // b4 is defined
      *iretvec -= coeff[12] * pow(fabs((thissz)-coeff[13]), coeff[14]);
    *iretvec = pow(10, (*iretvec)); // g(a) in um^-1
    // Zubko et al. functional definitions give g(a) in um^-1, convert to cm.
    *iretvec *= Constant::CM_UM; // g(a) in cm^-1
    iretvec++;
  }

  // Renormalize retvec to 1; ~g(a) is not correctly normalized like g(a) as it
  // is a numerical approximation to the normalized size distribution. This
  // normalization will depend slightly on the choice of size grid, so we can't
  // simply incorporate it into Normalization[]
  float _Nfac = 1.0 / NumUtils::integrate<float>(Component[cmp].size, retvec);
  transform(retvec.begin(), retvec.end(), retvec.begin(),
            bind2nd(multiplies<float>(), _Nfac));

  return retvec;
}

// Return the size distribution in cm^-1 H^-1
//  Wiengartner & Draine 2001, ApJ 548, 296
vector<float> GrainModel::getWD_sdist(vector<float> coeff, int cmp) {

  // float thissz;                     // Hold um copy of size at each point.
  vector<float> retvec(
      Component[cmp].size.size()); // The size distribution [cm^-1]
  vector<float>::iterator isz;     // Iterator through size vector
  vector<float>::iterator iretvec; // Iterator through size dist vector

  // Set iterator to beginning of size dist vector.
  iretvec = retvec.begin();

  // Wiengartner & Draine size distribution functions have mixed units - convert
  // where necessary to ensure out put in cm^-1 H^-1
  for (isz = Component[cmp].size.begin(); isz != Component[cmp].size.end();
       isz++) {
    // Compute the common part of the size distributions.
    *iretvec = coeff[3] / (*isz) *
               pow((*isz / (coeff[4] * static_cast<float>(Constant::UM_CM))),
                   coeff[6]);
    if (coeff[7] >= 0.0)
      *iretvec *= (1.0 + coeff[7] * (*isz) / (coeff[4] * Constant::UM_CM));
    else
      *iretvec /= (1.0 - coeff[7] * (*isz) / (coeff[4] * Constant::UM_CM));
    if ((*isz) > (coeff[4] * Constant::UM_CM)) {
      *iretvec *= exp(-1 * pow(((*isz) - coeff[4] * Constant::UM_CM) /
                                   (coeff[5] * Constant::UM_CM),
                               3));
    }
    if (coeff[0] > 0 && !isnan(coeff[0])) {
      *iretvec +=
          coeff[1] * coeff[0] / (*isz) *
          exp(-0.5 * (pow((log((*isz) / (3.5 * Constant::ANG_CM)) / 0.4), 2)));
      *iretvec +=
          coeff[2] * coeff[0] / (*isz) *
          exp(-0.5 * (pow((log((*isz) / (30.0 * Constant::ANG_CM)) / 0.4), 2)));
    }
    iretvec++;
  }

  return retvec;
}

// Return a power law size distribution in cm^-1
// vector <float> GrainModel::getPL_sdist(vector <float> coeff, int cmp)
// {
//   vector <float> retvec(Component[cmp].size.size()); // The size distribution
//   [cm^-1] vector <float>::iterator isz;     // Iterator through size vector
//   vector <float>::iterator iretvec; // Iterator through size dist vector

//   // Set iterator to beginning of size dist vector.
//   iretvec = retvec.begin();

//   return retvec;

// }

// Return Gaussian size distribution in cm^-1
vector<float> GrainModel::getGauss_sdist(vector<float> coeff, int cmp) {

  vector<float> retvec(
      Component[cmp].size.size()); // The size distribution [cm^-1]
  vector<float>::iterator isz;     // Iterator through size vector
  vector<float>::iterator isze;    // end of size vector
  vector<float>::iterator iretvec; // Iterator through size dist vector

  float z;

  isz = Component[cmp].size.begin();
  isze = Component[cmp].size.end();
  iretvec = retvec.begin();

  for (; isz != isze; ++isz, ++iretvec) {
    z = (*isz - coeff[1]) / coeff[2];
    *iretvec = coeff[0] * exp(-pow(z, 2) / 2);
  }

  return retvec;
}

float GrainModel::getTau(float a_wave) {

  if (a_wave < Wave[0] || a_wave > Wave[nWave - 1]) {
    cout << "You are attempting to retrieve \tau_(lam) outside of the defined"
         << endl;
    cout << "wavelength range. " << endl;
    cout << "Error in Grain::getCAbs(float)." << endl;
    cout << "Defined wavelength range is " << Wave[0] << "--" << Wave[nWave - 1]
         << endl;
    cout << "Requested wavelength: " << a_wave << endl;
    exit(8);
  }

  typedef vector<float>::iterator wvIter;
  wvIter it_beg = Wave.begin();
  wvIter it_loc = find(Wave.begin(), Wave.end(), a_wave);
  int idx = it_loc - it_beg;

  if (it_loc != Wave.end()) { // Found an exact match.
    return Tau[idx];
  } else { // weighted average

    it_loc = lower_bound(Wave.begin(), Wave.end(), a_wave);
    idx = it_loc - it_beg;

    float D1 = 1.0 / (pow((Wave[idx] - a_wave), 2));
    float D2 = 1.0 / (pow((a_wave - Wave[idx - 1]), 2));

    return ((Tau[idx] * D1 + Tau[idx - 1] * D2) / (D1 + D2));
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
