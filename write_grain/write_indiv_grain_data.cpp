// ======================================================================
//   Program to write out the grain properties for individual grains.
// This is to create the input files for python code to determine
// the grain size and composition distributions that fit observed
// extinction, abundances, scattering properties, infrared emission,
// and whatever else can be found to constrain dust grain models.
// Using the DIRTY dust grain code to avoid redevelopment and to
// take advantage of the extensive work that has already been on
// in the DIRTY coding effort.
//
// 2014 Oct/KDG - started
// ======================================================================
#include "write_indiv_grain_data.h"

int main(int argc, char* argv[])

{
  vector<float> RadiationField;
  int status;
  ofstream myfile;
  string myfilename;
  bool StochasticallyHeated;

  vector<float> RadiationFieldScale {0.5, 1.0, 2.0, 5.0, 10.0};
  int nISRFs;
  nISRFs = RadiationFieldScale.size();
  string RadiationFieldType;
  float RadiationFieldTemperature;

  float TauHeating;
  float TauCooling;
  float TauScaling;

  int nSize;
  vector<float> sizedist;
  vector<float> sizevals;
  vector<float> cabs;
  vector<float> csca;
  vector<float> phfunc;
  vector<float> Enthalpy;
  vector<float> CalorimetryTGrid;
  vector<double> EquilibriumEmission;
  vector<double> StochasticEmission;

  // parse the command line
  string param_filename;
  if (argc == 2) {
    param_filename = argv[1];
    // check that the file exists
    ifstream test_file(param_filename.c_str());
    if (test_file.fail()) {
      cout << "Parameter file (" << param_filename << ") does not exist."
           << endl;
      exit(8);
    }
    test_file.close();
  } else {
    cout << "Usage: write_indiv_grain_data parameter_file" << endl;
    exit(8);
  }

  // read parameter file
  ConfigFile param_data(param_filename);
  // get info on wavelength grid type
  int n_waves;
  std::vector<float> wavelength;

  int i = 0;

  string grid_type = param_data.SValue("Dust Grains", "wave_type");
  check_input_param("wave_type", grid_type, "res");
  if (grid_type == "res") {
    float wave_min = param_data.FValue("Dust Grains", "wave_min");
    check_input_param("dust_model: wave_min", wave_min, 0.05, 1e4);
    float wave_max = param_data.FValue("Dust Grains", "wave_max");
    check_input_param("dust_model: wave_min", wave_min, 0.05, 1e4);
    float wave_res = param_data.FValue("Dust Grains", "wave_resolution");
    check_input_param("dust_model: wave_resolution", wave_res, 0.05, 1e4);

    // determine the number of wavelength points needed
    n_waves = int(log10(wave_max / wave_min) /
                      log10((1.0 + 2.0 * wave_res) / (2.0 * wave_res - 1.0)) +
                  1.0);

    // cout << "n_waves = " << n_waves << endl;

    // construct the wavelength grid
    float log_wave_min = log10(wave_min);
    float delta_log_wave = (log10(wave_max) - log_wave_min) / (n_waves - 1);
    float cur_wave = 0.0;
    for (i = 0; i < n_waves; i++) {
      cur_wave = log_wave_min + delta_log_wave * float(i);
      cur_wave = pow(double(10.0), double(cur_wave));
      // convert wavelengths from micron to cm
      cur_wave *= Constant::UM_CM;
      wavelength.push_back(cur_wave);
    }
  } else if (grid_type == "file") {
    // get the filename
    string wave_filename = param_data.SValue("Dust Grains", "wave_file");
    // check that the file exists
    ifstream wave_file(wave_filename.c_str());
    if (wave_file.fail()) {
      cout << "Empirical wave file (" << wave_filename << ") does not exist."
           << endl;
      exit(8);
    }
    wave_file.close();

    // get the wavelength grid
    vector<double> in_wavelength;
    DataFile(wave_filename, in_wavelength);

    // go through the values and make sure they are within bounds
    double check_wave = 0.005;
    n_waves = wavelength.size();
    for (i = 0; i < n_waves; i++) {
      check_input_param("wavelength from file", in_wavelength[i], check_wave,
                        1e4);
      check_wave = in_wavelength[i];
      wavelength.push_back(wavelength[i] * Constant::UM_CM);
    }

  } else {
    cout << "wavelength grid wave_type = " << grid_type << " not known" << endl;
    cout << "allowed types = (res)" << endl;
    exit(8);
  }

  // initialize the 2D vector needed to store the emission from different ISRFs
  vector<float> row(nISRFs, 0.0);
  vector<vector<float>> EquilEmiss_All(n_waves, row);
  vector<vector<float>> StochEmiss_All(n_waves, row);

  // Grain definitions.
  FitGrains FitGrainDefines;
  FitGrainDefines.MakeFitGrains(param_data, wavelength);

  // Generate
  // Define a radiation field.
  // The local ISRF.
  RadiationFieldType = "IRSF";
  RadiationFieldTemperature = -1.0;  // for e.g. Black body field
  TauScaling = 100.0;

  // now get the individual grain information
  EquilibriumEmission.resize(n_waves);
  StochasticEmission.resize(n_waves);
  for (int i = 0; i < FitGrainDefines.getNComp(); ++i) {
    // for (int i=FitGrainDefines.getNComp()-1;i<FitGrainDefines.getNComp();++i)
    // {
    cout << "Working on " << FitGrainDefines.ComponentName(i) << " (" << i + 1
         << " of " << FitGrainDefines.getNComp() << ")" << endl;
    nSize = FitGrainDefines.nSize(i);
    sizevals = FitGrainDefines.Size(i);
    myfile << endl << i << " " << nSize << endl;

    for (int j = 0; j < nSize; ++j) {
      // for (int j=40;j<41;++j) {
      cout << "Size: " << j+1 << " " << sizevals[j] << endl;
      stringstream ss;
      ss << FitGrainDefines.getModelName() << "_c_" << (i + 1) << "_"
         << FitGrainDefines.ComponentName(i);
      ss << "_s_" << (j + 1) << ".dat";
      myfilename = ss.str();

      cabs = FitGrainDefines.CAbs(i, j);
      csca = FitGrainDefines.CSca(i, j);
      phfunc = FitGrainDefines.phFunc(i, j);
      Enthalpy = FitGrainDefines.Enthalpy(j, i);
      CalorimetryTGrid = FitGrainDefines.CalTemp(i);

      // loop over the desired radiation fields
      for (int jj = 0; jj < nISRFs; jj++) {
        cout << "RadField = " << RadiationFieldScale[jj] << endl;
        ISRF* myISRF = new ISRF(wavelength, RadiationFieldScale[jj]);
        RadiationField = myISRF->getISRF();

        // pass in size index, component index, component name, so I can write
        // out the pdist
        status = ComputeDustEmission_igrain(
            RadiationField, cabs, wavelength, Enthalpy, CalorimetryTGrid,
            EquilibriumEmission, StochasticEmission, TauHeating, TauCooling,
            TauScaling, StochasticallyHeated, FitGrainDefines.getModelName(),
            FitGrainDefines.ComponentName(i), i, sizevals[j], j,
            RadiationFieldType, RadiationFieldScale[jj], RadiationFieldTemperature);
        if (status != Flags::FSUCCESS) cout << "FUCK" << endl;

        // save the radiation field for output
        for (int k = 0; k < n_waves; k++) {
          EquilEmiss_All[k][jj] = EquilibriumEmission[k];
          StochEmiss_All[k][jj] = StochasticEmission[k];
        }

      }

      // exit(8);
      myfile.open(myfilename.c_str());
      myfile << "# " << sizevals[j] << " " << i << " " << j
             << " : grain size, compostion index, size index" << endl;
      if (StochasticallyHeated)
        myfile << "# This grain was stochastically heated: Yes" << endl;
      else
        myfile << "# This grain was stochastically heated: No" << endl;
      myfile << "# Heating and cooling timescale: " << TauHeating << " "
             << TauCooling << endl;
      myfile << "# Ratio of heating to cooling to determine stochastic cutoff: "
             << TauScaling << endl;
      myfile << "# Equilibrium emission is recorded even for those grains in "
                "the stochastic regime."
             << endl;
      myfile << "# Emission spectrum was computed for: " << endl;
      myfile << "# Type: " << RadiationFieldType << endl;
      myfile << "# Scales: ";
      for (int k = 0; k < nISRFs; k++)
        myfile << " " << RadiationFieldScale[k];
      myfile << endl;
      myfile << "# Temperature: " << RadiationFieldTemperature << endl;
      myfile << "# Wavelength  CExt  CAbs  CSca  Albedo G";
      for (int k = 0; k < nISRFs; k++)
        myfile << " EqEm" << (k+1) << " StEm" << (k+1);
      myfile << endl;

      for (int k = 0; k < n_waves; k++) {
        myfile << wavelength[k] * Constant::CM_UM;
        myfile << " " << cabs[k] + csca[k];
        myfile << " " << cabs[k];
        myfile << " " << csca[k];
        myfile << " " << csca[k] / (cabs[k] + csca[k]);
        myfile << " " << phfunc[k];
        for (int kk = 0; kk < nISRFs; kk++) {
          myfile << " " << EquilEmiss_All[k][kk];
          myfile << " " << StochEmiss_All[k][kk];
        }
        myfile << endl;
      }
      myfile.close();
    }
    // cout << endl;
  }

  // using namespace std;

  //   int j,k;
  //   for (i = 0; i < CurGrainModel.getNComp(); i++) {
  //     sizedist = CurGrainModel.getSizeDistribution(i);
  //     sizevals = CurGrainModel.Size(i);
  //     myfile << endl << i << " " << sizedist.size() << endl;

  //     cout << sizedist.size() << endl;

  //     for (j = 0; j < int(sizedist.size()); j++) {

  //       status=0;
  //       stringstream ss;
  //       ss << CurGrainModel.getModelName() << "_c_" << (i+1) << "_" <<
  //       CurGrainModel.getComponentName(i); ss << "_s_" << (j+1) << ".dat";
  //       myfilename = ss.str();

  //       cout << myfilename << endl;

  //       myfile.open(myfilename.c_str());
  //       myfile << "# " << sizevals[j] << " " << i << " " << j << " : grain
  //       size, compostion index, size index" << endl; myfile << "# Wavelength
  //       CExt  CAbs  CSca  Albedo" << endl;

  //       cabs = CurGrainModel.CAbs(i,j);
  //       csca = CurGrainModel.CSca(i,j);

  //       cout << "Size: " << sizevals[j]*Constant::CM_UM << endl;
  //       status = ComputeDustEmission_igrain(RadiationField, cabs,
  //       wavelength);

  //       cout << cabs.size() << endl;
  //       for (k = 0; k < int(cabs.size()); k++) {
  // 	myfile << wavelength[k]*Constant::CM_UM;
  // 	myfile << " " << cabs[k] + csca[k];
  // 	myfile << " " << cabs[k];
  // 	myfile << " " << csca[k];
  // 	myfile << " " << csca[k]/(cabs[k]+csca[k]);
  // 	myfile << endl;
  //       }
  //       myfile.close();
  //     }
  //   }
}
