/**************************************************************
 * DataFile -- DataFile Class, header file.
 *
 * History:
 *     Written by:   Putrid, August 2004
 *                   Contact: misselt@as.arizona.edu
 *
 * TODO : Make templates?
 **************************************************************/
#ifndef __DATA_FILE_H__
#define __DATA_FILE_H__

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class DataFile {

public:
  DataFile(string const &dataFile, vector<double> &v1);
  DataFile(string const &dataFile, vector<double> &v1, vector<double> &v2);
  DataFile(string const &dataFile, vector<double> &v1, vector<double> &v2,
           vector<double> &v3);
  DataFile(string const &dataFile, vector<double> &v1, vector<double> &v2,
           vector<double> &v3, vector<double> &v4);
  DataFile(string const &dataFile, vector<double> &v1, vector<double> &v2,
           vector<double> &v3, vector<string> &v4);
  DataFile(string const &dataFile, vector<double> &v1, vector<double> &v2,
           vector<double> &v3, vector<double> &v4, vector<double> &v5);
  DataFile(string const &dataFile, vector<double> &v1, vector<double> &v2,
           vector<double> &v3, vector<double> &v4, vector<string> &v5);
};

#endif
