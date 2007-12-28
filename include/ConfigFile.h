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

#include <iostream>
#include <string>
#include <map>
#include <fstream>

using namespace std;

class ConfigFile {

 private: 
  map<string,string> content_;
  
 public:
  ConfigFile(string const& configFile);
  
  bool BValue(string const& section, string const& entry) const;
  string SValue(string const& section, string const& entry) const;
  double DValue(string const& section, string const& entry) const; 
  float  FValue(string const& section, string const& entry) const; 
  int    IValue(string const& section, string const& entry) const;
  long   LValue(string const& section, string const& entry) const;

};

#endif
