// ******************************************************************************
// Grain -- Grain Class
// History:
//     Written by:   Putrid, August 2004
//                   Contact: misselt@as.arizona.edu
//                   Re-written completely, Sept 2006.
//
//     RESTRICTIONS: MasterWave and MasterSize MUST BE PASSED IN CGS!
//                   ie. Both should be in 'cm'.
//                   For calorimetry, we expect specific heat capacity and
//                   specific enthalpy files to have the form:
//                       T      E      C_v
//                   # are interpreted as comment lines.
//                   The first line after the comments and before the data sould
//                   have the form:
//                   i   TMIN,TMAX NT:  1.0,2500.0 50
// ******************************************************************************
#include "Grain.h"

Grain::Grain ()
{
  // Empty Constructor
}

// Pass copies; faster and we can return through the call.

void
Grain::MakeGrain (string const &inComponentName, string const &fOpticalConstants,
                  string const &fCalorimetry, vector<float> &MasterWave, vector<float> &MasterSize,
                  string const &Path, float a_min, float a_max)
{
  cout << "Instantiating Grain Object for " << fOpticalConstants << endl;

  // set the name of this component
  ComponentName = inComponentName;

  // ***************************************************************************
  // Check if we are going to interpolate onto new wavelength and/or size grids
  bool interpWave = true;
  bool interpSize = true;
  if (MasterWave[0] == -1)
    interpWave = false;
  if (MasterSize[0] == -1)
    interpSize = false;
  // ***************************************************************************

  // ***************************************************************************
  // Open grain "constants" file.
  OpticalFile = fOpticalConstants;
  ifstream _file ((Path + fOpticalConstants).c_str ());
  if (!_file.is_open ())
    {
      cout << "Unable to open optical constants file " << fOpticalConstants << endl;
      exit (8);
    }
  // ***************************************************************************

  // ***************************************************************************
  // Local Variable defines.
  // Local variables for processing the input file.
  string _line, _filename;
  vector<string> _parts;
  bool _StartData = false;
  int _offset;

  // Local variables to hold grain info for reading/interpolation/etc
  int _nsize, _nwave;
  vector<float> _Qabs;
  vector<float> _Qsca;
  vector<float> _g;
  float _Const;
  vector<float> _wave;
  vector<float> _size;
  vector<float> Crap;
  vector<vector<float> > _Qabs2d;
  vector<vector<float> > _Qsca2d;
  vector<vector<float> > _g2d;
  // ***************************************************************************

  // ***************************************************************************
  // Read over any header information

  // Comments:
  while (!_StartData)
    {
      getline (_file, _line);
      if (_line[0] == '#')
        continue;
      _StartData = true;
    }

  // General description of material in file.
  // First line after comments is a line describing the grain sizes -
  // We've just read that one in.
  _parts = StringManip::strsplit (_line);
  _nsize = atoi (_parts[1].c_str ());
  // Next line is a line describing the wavelengths tabulated
  getline (_file, _line);
  _parts = StringManip::strsplit (_line);
  _nwave = atoi (_parts[1].c_str ());

  // Stoichiometric Info.
  getline (_file, _line);
  _parts = StringManip::strsplit (_line);
  nelements = atoi (_parts[0].c_str ());
  density = atof (_parts[1].c_str ());
  for (int _elem = 0; _elem < nelements; _elem++)
    {
      getline (_file, _line);
      _parts = StringManip::strsplit (_line);
      amu.push_back (atof (_parts[0].c_str ()));
      stoich.push_back (atof (_parts[1].c_str ()));
      atom.push_back (_parts[2]);
    }

  // Set the grain sublimation temperature.
  SubT = 1500.0; // Will want to have this read from input file.
  // ***************************************************************************

  // ***************************************************************************
  // Get tabulated optical constants.

  // Size the temporary input vectors, dimension 1
  _Qabs2d.resize (_nsize);
  _Qsca2d.resize (_nsize);
  _g2d.resize (_nsize);

  for (int sz = 0; sz < _nsize; sz++)
    {
      getline (_file, _line); // Blank line
      getline (_file, _line); // This Size definition
      _parts = StringManip::strsplit (_line, "=");
      _size.push_back (atof (_parts[0].c_str ()) * Constant::UM_CM); // Size in cm
      _Const = (Constant::PI)*pow (_size[sz], 2);
      getline (_file, _line); // column headers
      cout.flush ();

      for (int wv = 0; wv < _nwave; wv++)
        {
          getline (_file, _line); // First data line.
          _parts = StringManip::strsplit (_line);
          _offset = _parts.size () - 5;
          if (sz == 0)
            _wave.push_back (atof (_parts[_offset].c_str ()) * Constant::UM_CM);
          _Qabs.push_back (atof (_parts[_offset + 1].c_str ()) * _Const);
          _Qsca.push_back (atof (_parts[_offset + 2].c_str ()) * _Const);
          _g.push_back (atof (_parts[_offset + 4].c_str ()));
        }

      // Interpolate onto MasterWave
      if (interpWave)
        {
          // cout << "INTERPOLATING ON WAVE GRID" << endl;
          NumUtils::interpolr (_Qabs, _wave, MasterWave);
          NumUtils::interpolr (_Qsca, _wave, MasterWave);
          NumUtils::interpolr (_g, _wave, MasterWave);
        }

      // Copy into vector<vector>
      copy (_Qabs.begin (), _Qabs.end (), back_inserter (_Qabs2d[sz]));
      copy (_Qsca.begin (), _Qsca.end (), back_inserter (_Qsca2d[sz]));
      copy (_g.begin (), _g.end (), back_inserter (_g2d[sz]));

      // Clear temporary vectors for push_back.
      _Qabs.erase (_Qabs.begin (), _Qabs.end ());
      _Qsca.erase (_Qsca.begin (), _Qsca.end ());
      _g.erase (_g.begin (), _g.end ());
    }

  // Done with optical "constants" file.
  _file.close ();
  // for (int _sz=0;_sz<MasterSize.size();++_sz) cout << _sz << " " <<
  // MasterSize[_sz] << endl; for (int _sz=0;_sz<_nsize;++_sz) cout << _sz << "
  // " << _size[_sz] << endl;
  //  Assign MasterWave to the object.
  if (interpWave)
    {
      nwave = MasterWave.size ();
      copy (MasterWave.begin (), MasterWave.end (), back_inserter (wave));
    }
  else
    {
      nwave = _nwave;
      copy (_wave.begin (), _wave.end (), back_inserter (wave));
    }

  // We will not extrapolate the size distribution.  If MasterSize
  // comes in outside the bounds of where Q's are defined, MasterSize
  // WILL BE MODIFIED!
  if (interpSize)
    {
      MasterSize.erase (
          remove_if (MasterSize.begin (), MasterSize.end (), bind2nd (less<float> (), _size[0])),
          MasterSize.end ());
      MasterSize.erase (remove_if (MasterSize.begin (), MasterSize.end (),
                                   bind2nd (greater<float> (), _size[_nsize - 1])),
                        MasterSize.end ());
      // trim size of MasterSize if needed.
      vector<float> (MasterSize).swap (MasterSize);

      // Assign correct size grid to class member.
      nsize = MasterSize.size ();
      copy (MasterSize.begin (), MasterSize.end (), back_inserter (size));
      // Trim on amin/max
    }
  else
    {
      // If the user has defined a_min and/or a_max, trim the file tabulated
      // size grid to the new limits. If a_min or a_max are outside the
      // tabulated range, the tabulated range ends up being used... NOTE THAT
      // THIS COULD RESULT IN UNINTENDED CONSEQUENCES FOR EG. A PRE-DEFINED
      // MODEL.
      // Copy input tabulation into class size member.
      copy (_size.begin (), _size.end (), back_inserter (size));

      // If necessary, modify class size member.
      if (a_min != -1)
        { // been set...
          // Stick in the lower limit if it is within defined bounds.
          // Otherwise, keep lower bound as is.
          // if (a_min > size[0])
          //	size.insert(size.begin(),a_min);
          // for (int _sz=0;_sz<size.size();++_sz) cout << _sz << " " << size[_sz]
          // << endl;
          // Remove everything ouside of a_min
          size.erase (remove_if (size.begin (), size.end (), bind2nd (less<float> (), a_min)),
                      size.end ());
          // cout << "post: " << endl;
          // for (int _sz=0;_sz<size.size();++_sz) cout << _sz << " " << size[_sz]
          // << endl;
          //  Since we may have modified the size distribution, do interpolate.
          interpSize = true;
        }
      if (a_max != -1)
        { // been set...
          // Stick in the upper limit is it is within defined bounds
          // Otherwise, keep upper bound as is.
          // if (a_max < size[size.size()-1])
          //	size.push_back(a_max);
          size.erase (remove_if (size.begin (), size.end (), bind2nd (greater<float> (), a_max)),
                      size.end ());
          // Since we may have modified the size distribution, do interpolate.
          interpSize = true;
        }
      // Trim size of size if needed.
      vector<float> (size).swap (size);
      nsize = size.size ();
    }
  // Now, size[] holds the correct size grid.

  // Size the class members to the correct, final size.
  CAbs.resize (nsize);
  CSca.resize (nsize);
  phFunc.resize (nsize);
  mass.resize (nsize);
  for (int sz = 0; sz < nsize; sz++)
    {
      CAbs[sz].resize (nwave);
      CSca[sz].resize (nwave);
      phFunc[sz].resize (nwave);
      mass[sz] = (Constant::VFAC)*size[sz] * size[sz] * size[sz] * density;
    }

  // Loaded up our temporary 2d array, interpolated onto proper
  // wave length scale.  Now interpolate onto proper size scale.
  for (int wv = 0; wv < nwave; wv++)
    {
      // Populate temp vectors with optical constants as f(size)
      // for each wavelength.
      for (int sz = 0; sz < _nsize; sz++)
        {
          _Qabs.push_back (_Qabs2d[sz][wv]);
          _Qsca.push_back (_Qsca2d[sz][wv]);
          _g.push_back (_g2d[sz][wv]);
        }

      // Interpolate onto correct size grid
      if (interpSize)
        {
          NumUtils::interpolr (_Qabs, _size, size);
          NumUtils::interpolr (_Qsca, _size, size);
          NumUtils::interpolr (_g, _size, size);
        }

      // Assign class members the correct values, C in cm^2
      for (int sz = 0; sz < nsize; sz++)
        {
          CAbs[sz][wv] = _Qabs[sz];
          CSca[sz][wv] = _Qsca[sz];
          phFunc[sz][wv] = _g[sz];
        }

      // Erase temporary 1d vectors.
      _Qabs.erase (_Qabs.begin (), _Qabs.end ());
      _Qsca.erase (_Qsca.begin (), _Qsca.end ());
      _g.erase (_g.begin (), _g.end ());
    }

  // ***************************************************************************
  // Get calorimetric properties of grain material.

  // Open calorimetry file and read in data.
  CalorimetryFile = fCalorimetry;
  ifstream _filecal ((Path + fCalorimetry).c_str ());
  if (!_filecal.is_open ())
    {
      cout << "Unable to open calorimetry file " << fCalorimetry << endl;
      exit (8);
    }

  // Get Comment lines
  _StartData = false; // Reset for comment lines
  while (!_StartData)
    {
      getline (_filecal, _line);
      if (_line[0] == '#')
        continue;
      _StartData = true;
    }

  // First non comment line should look like:
  // TMIN,TMAX,NT: 1.0 2500.0  50
  _parts = StringManip::strsplit (_line);
  nTemp = atoi (_parts[3].c_str ());
  for (int nt = 0; nt < nTemp; nt++)
    {
      getline (_filecal, _line);
      _parts = StringManip::strsplit (_line);
      Temperature.push_back (atof (_parts[0].c_str ()));
      SpecificEnthalpy.push_back (atof (_parts[1].c_str ()));
      SpecificHeatCapacity.push_back (atof (_parts[2].c_str ()));
    }
  // Done with calorimetry file
  _filecal.close ();
  // Generate the Enthalpy and specific heat for each grain. Normally let
  // calling code do the conversion for whatever size it wants, but for RT
  // application, this will cut down on the number of multiplies; maybe have
  // a switch for memory conservation.
  Enthalpy.resize (nsize);
  HeatCapacity.resize (nsize);
  for (int sz = 0; sz < nsize; ++sz)
    {
      Enthalpy[sz].resize (nTemp);
      HeatCapacity[sz].resize (nTemp);
      for (int nt = 0; nt < nTemp; ++nt)
        {
          Enthalpy[sz][nt] = SpecificEnthalpy[nt] * mass[sz];
          HeatCapacity[sz][nt] = SpecificHeatCapacity[nt] * mass[sz];
        }
    }
  // ***************************************************************************
}

// ***************************************************************************
// getCAbs(FLOAT): return CAbs for a specified size in cm (rather than at a
// given index.  If a_size is in [size_min,size_max], we will return the CAbs
// If a_size !in [size_min,size_max], it is an error.
// In the range, if there's an exact match, just return CAbs at the correct
// index.  If not, linear interpolation will be done.
vector<float>
Grain::getCAbs (float a_size)

{
  if (a_size < size[0] || a_size > size[nsize - 1])
    {
      cout << "You are attempting to retrieve optical constants outside of the" << endl;
      cout << "size range defined by your input grains.  This is dangerous and" << endl;
      cout << "I VANT TO MAKE AB-SA-LUTE-LY CLEE-AR, STREECTLY VORBIDEN!" << endl;
      cout << "Error in Grain::getCAbs(float)." << endl;
      cout << "Defined size range is " << size[0] << "--" << size[nsize - 1] << endl;
      cout << "Requested size: " << a_size << endl;
      exit (8);
    }

  typedef vector<float>::iterator szIter;
  szIter it_beg = size.begin ();
  szIter it_loc = find (size.begin (), size.end (), a_size);
  int idx = it_loc - it_beg;

  if (it_loc != size.end ())
    { // We have found an exact match
      return CAbs[idx];
    }
  else
    { // do a weighted average.
      it_loc = lower_bound (size.begin (), size.end (), a_size);
      idx = it_loc - it_beg;
      vector<float> Tmp1 (nwave);
      // weight by the square of distance in size space
      float D1 = 1.0 / (pow ((size[idx] - a_size), 2));
      float D2 = 1.0 / (pow ((a_size - size[idx - 1]), 2));
      float DSum = D1 + D2;

      vector<float>::iterator iCAbs1, iCAbs2, iTmp1;
      iCAbs2 = CAbs[idx - 1].begin ();
      iTmp1 = Tmp1.begin ();
      for (iCAbs1 = CAbs[idx].begin (); iCAbs1 != CAbs[idx].end (); iCAbs1++)
        {
          *iTmp1 = ((*iCAbs1) * D1 + (*iCAbs2) * D2) / DSum;
          iTmp1++;
          iCAbs2++;
        }

      return Tmp1;
    }
}

// ***************************************************************************
// getCSca(FLOAT): return CSca for a specified size in cm (rather than at a
// given index.  If a_size is in [size_min,size_max], we will return the CSca
// If a_size !in [size_min,size_max], it is an error.
// In the range, if there's an exact match, just return CSca at the correct
// index.  If not, linear interpolation will be done.
vector<float>
Grain::getCSca (float a_size)

{
  if (a_size < size[0] || a_size > size[nsize - 1])
    {
      cout << "You are attempting to retrieve optical constants outside of the" << endl;
      cout << "size range defined by your input grains.  This is dangerous and" << endl;
      cout << "I VANT TO MAKE AB-SA-LUTE-LY CLEE-AR, STREECTLY VORBIDEN!" << endl;
      cout << "Error in Grain::getCSca(float)." << endl;
      cout << "Defined size range is " << size[0] << "--" << size[nsize - 1] << endl;
      cout << "Requested size: " << a_size << endl;
      exit (8);
    }

  typedef vector<float>::iterator szIter;
  szIter it_beg = size.begin ();
  szIter it_loc = find (size.begin (), size.end (), a_size);
  int idx = it_loc - it_beg;

  if (it_loc != size.end ())
    { // We have found an exact match
      return CSca[idx];
    }
  else
    { // do a weighted average.
      it_loc = lower_bound (size.begin (), size.end (), a_size);
      idx = it_loc - it_beg;
      vector<float> Tmp1 (nwave);

      // weight by the square of distance in size space
      float D1 = 1.0 / (pow ((size[idx] - a_size), 2));
      float D2 = 1.0 / (pow ((a_size - size[idx - 1]), 2));
      float DSum = D1 + D2;

      vector<float>::iterator iCSca1, iCSca2, iTmp1;
      iCSca2 = CSca[idx - 1].begin ();
      iTmp1 = Tmp1.begin ();
      for (iCSca1 = CSca[idx].begin (); iCSca1 != CSca[idx].end (); iCSca1++)
        {
          *iTmp1 = ((*iCSca1) * D1 + (*iCSca2) * D2) / DSum;
          iTmp1++;
          iCSca2++;
        }

      return Tmp1;
    }
}

// ***************************************************************************
// getphFunc(FLOAT): return phFunc for a specified size in cm (rather than at
// a given index.  If a_size is in [size_min,size_max], we will return the
// phFunc. If a_size !in [size_min,size_max], it is an error.
// In the range, if there's an exact match, just return phFunc at the correct
// index.  If not, linear interpolation will be done.
vector<float>
Grain::getphFunc (float a_size)

{
  if (a_size < size[0] || a_size > size[nsize - 1])
    {
      cout << "You are attempting to retrieve optical constants outside of the" << endl;
      cout << "size range defined by your input grains.  This is dangerous and" << endl;
      cout << "I VANT TO MAKE AB-SA-LUTE-LY CLEE-AR, STREECTLY VORBIDEN!" << endl;
      cout << "Error in Grain::getphFunc(float)." << endl;
      cout << "Defined size range is " << size[0] << "--" << size[nsize - 1] << endl;
      cout << "Requested size: " << a_size << endl;
      exit (8);
    }

  typedef vector<float>::iterator szIter;
  szIter it_beg = size.begin ();
  szIter it_loc = find (size.begin (), size.end (), a_size);
  int idx = it_loc - it_beg;

  if (it_loc != size.end ())
    { // We have found an exact match
      return CSca[idx];
    }
  else
    { // do a weighted average.
      it_loc = lower_bound (size.begin (), size.end (), a_size);
      idx = it_loc - it_beg;
      vector<float> Tmp1 (nwave);

      // weight by the square of distance in size space
      float D1 = 1.0 / (pow ((size[idx] - a_size), 2));
      float D2 = 1.0 / (pow ((a_size - size[idx - 1]), 2));
      float DSum = D1 + D2;

      vector<float>::iterator iphFunc1, iphFunc2, iTmp1;
      iphFunc2 = phFunc[idx - 1].begin ();
      iTmp1 = Tmp1.begin ();
      for (iphFunc1 = phFunc[idx].begin (); iphFunc1 != phFunc[idx].end (); iphFunc1++)
        {
          *iTmp1 = ((*iphFunc1) * D1 + (*iphFunc2) * D2) / DSum;
          iTmp1++;
          iphFunc2++;
        }

      return Tmp1;
    }
}
