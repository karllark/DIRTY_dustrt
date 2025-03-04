/**************************************************************
 * ConfigFile -- ConfigFile Class, header file.
 *
 * History:
 *     Written by:   Putrid, August 2004
 *                   Contact: misselt@as.arizona.edu
 *
 *
 **************************************************************/
#ifndef __CONFIG_FILE_H__
#define __CONFIG_FILE_H__

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

using namespace std;

class ConfigFile
{
private:
  map<string, string> content_;
  string _BadString;
  double _BadDouble;
  float _BadFloat;
  int _BadInt;
  long _BadLong;

public:
  ConfigFile (string const &configFile);

  bool BValue (string const &section, string const &entry) const;
  string SValue (string const &section, string const &entry) const;
  double DValue (string const &section, string const &entry) const;
  float FValue (string const &section, string const &entry) const;
  int IValue (string const &section, string const &entry) const;
  long LValue (string const &section, string const &entry) const;

  inline string
  BadString (void)
  {
    return _BadString;
  }
  inline bool
  BadFloat (float _val)
  {
    return _val != _val;
  }
  inline bool
  BadDouble (double _val)
  {
    return _val != _val;
  }
  inline int
  BadInt (void)
  {
    return _BadInt;
  }
  inline long
  BadLong (void)
  {
    return _BadLong;
  }

  inline bool
  isBadFloat (float _val)
  {
    if (isnan (_val))
      return true;
    else
      return false;
  }

  inline bool
  isBadDouble (float _val)
  {
    if (isnan (_val))
      return true;
    else
      return false;
  }
};

#endif
