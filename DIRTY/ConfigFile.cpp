/**************************************************************
 * ConfigFile -- ConfigFile Class, source code.
 *
 * Read in configuration files.  Files should be split into
 * sections with header [SectionName] and entries should be
 * in the form of "key=value".  One can retrieve configuration
 * information by specifying the Section and Key you're after.
 * For example, in a configuration file like the following:
 *
 * # Input data2
 * [Data Files]
 * Enthalpy=Enthalpies.dat
 * Wave=Wave.dat
 *
 * # Grain properties
 * [Grains]
 * Comp1=suvSil_121.dat
 * Comp2=AMC_121.dat
 * Comp3=neutPAH_30.dat
 *
 * # Model parameters
 * [Model]
 * tauv=0.5
 * sfr=10.0
 * geometry=Dusty
 *
 * To retieve the name of the file containing the optical
 * constants of the the second component, the star formation
 * rate and the geometry, one would use calls with the
 * following structure:
 *
 * ConfigFile cf("MyConfigFile");
 * Comp[1]  = cf.SValue("Grains","Comp2");
 * SFR      = cf.DValue("Model","sfr");
 * Geometry = cf.SValue("Model","geometry");
 *
 * Member functions SValue and DValue require string pairs
 * describing "Section" and "Key" for retrieval and return
 * String and Double values, respectively.
 *
 *  TODO: Template-ify so we need only one member function...
 *        ...Possible?
 *
 * History:
 *     Written by:   Putrid, August 2004
 *                   Contact: misselt@as.arizona.edu
 *                   Borrowed liberaly from
 *                       http://www.adp-gmbh.ch/
 *                   Putrid, Mar 2006 - added "Bad" return values
 *                      for not found key/value pairs. Calling code
 *                      should handle those internally.
 *                      - SValue.BAD = NULL
 *                      - IValue.BAD = -99... yeah, yeah, not very robust.
 *                   KDG 30 Jun 2006
 *                      updated all functions to return null values
 *                      added LVALUE function for longs
 *                   Putrid, Late 2007 sometime.  More robust handling
 *                      of NaNs in _BadValues.
 *
 **************************************************************/
#include "ConfigFile.h"

using namespace std;

ConfigFile::ConfigFile(string const &configFile) {

  _BadString = "NULL";
  _BadFloat = strtof("NaN", NULL);
  _BadDouble = strtod("NaN", NULL);
  _BadInt = -99;
  _BadLong = -99;

  ifstream file(configFile.c_str());

  string line;
  string name;
  string value;
  string inSection;
  int posEqual;

  while (getline(file, line)) {

    if (!line.length())
      continue;

    if (line[0] == '#')
      continue;
    if (line[0] == ';')
      continue;

    if (line[0] == '[') {
      inSection = line.substr(1, line.find(']') - 1);

      continue;
    }

    posEqual = line.find('=');
    name = line.substr(0, posEqual);

    value = line.substr(posEqual + 1);

    content_[inSection + '/' + name] = value;
  }
}

bool ConfigFile::BValue(string const &section, string const &entry) const {

  map<string, string>::const_iterator ci = content_.find(section + '/' + entry);

  if (ci == content_.end()) {
    return false;
  }
  string tmp_str = ci->second;
  transform(tmp_str.begin(), tmp_str.end(), tmp_str.begin(),
            (int (*)(int))tolower);

  if (tmp_str == "yes")
    return true;

  return false;
}

string ConfigFile::SValue(string const &section, string const &entry) const {

  map<string, string>::const_iterator ci = content_.find(section + '/' + entry);

  if (ci == content_.end()) {
    // return "NULL";
    return _BadString;
  }
  return ci->second;
}

float ConfigFile::FValue(string const &section, string const &entry) const {

  map<string, string>::const_iterator ci = content_.find(section + '/' + entry);

  if (ci == content_.end()) {
    return _BadFloat;
    // return strtof("NaN",NULL);  // default value
  }
  return atof(ci->second.c_str());
}

double ConfigFile::DValue(string const &section, string const &entry) const {

  map<string, string>::const_iterator ci = content_.find(section + '/' + entry);

  if (ci == content_.end()) {
    return _BadDouble;
    // return strtod("NaN",NULL);  // default value
  }
  return atof(ci->second.c_str());
}

int ConfigFile::IValue(string const &section, string const &entry) const {

  map<string, string>::const_iterator ci = content_.find(section + '/' + entry);

  if (ci == content_.end()) {
    return _BadInt;
    // return -99;
    // return int(strtof("NaN",NULL));  // default value ... have to check
    // this...
  }
  return atoi(ci->second.c_str());
}

long ConfigFile::LValue(string const &section, string const &entry) const {

  map<string, string>::const_iterator ci = content_.find(section + '/' + entry);

  if (ci == content_.end()) {
    return _BadLong;
    // return -99;  // default value
  }
  return atol(ci->second.c_str());
}
