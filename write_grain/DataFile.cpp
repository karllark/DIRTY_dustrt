/**************************************************************
 * DataFile -- DataFile Class, source code.
 *
 * Read in data files 1 - 5 columns of FP data.  DataFile is
 * overloaded; number of vectors passed determines constructor
 * Requires a string holding the filename and then N vector
 * types, 1 for each column in the input file. Lines begining
 * with '#' are considered coments and are ignored.  Blank
 * lines are ignored.  Yeah, yeah, there's probably a more
 * clever elegant way to do this but...
 *
 * History:
 *     Written by:   Putrid, August 2004
 *                   Contact: misselt@as.arizona.edu
 *
 **************************************************************/
#include "DataFile.h"
#include <iostream>

using namespace std;

// One column
DataFile::DataFile(string const &dataFile, vector<double> &v1) {

  ifstream file(dataFile.c_str());

  //  int linecount=0;
  string buf;
  string line;

  while (getline(file, line)) {

    if (!line.length())
      continue;
    if (line[0] == '#')
      continue;

    stringstream ss(line);
    while (ss >> buf) {
      v1.push_back(atof(buf.c_str()));
    }
  }
}

// Two column
DataFile::DataFile(string const &dataFile, vector<double> &v1,
                   vector<double> &v2) {

  ifstream file(dataFile.c_str());

  //  int linecount=0;
  string buf1, buf2;
  string line;

  while (getline(file, line)) {
    if (!line.length())
      continue;
    if (line[0] == '#')
      continue;

    stringstream ss(line);
    while (ss >> buf1 >> buf2) {
      v1.push_back(atof(buf1.c_str()));
      v2.push_back(atof(buf2.c_str()));
    }
  }
}

// Three column
DataFile::DataFile(string const &dataFile, vector<double> &v1,
                   vector<double> &v2, vector<double> &v3) {

  ifstream file(dataFile.c_str());

  //  int linecount=0;
  string buf1, buf2, buf3;
  string line;

  while (getline(file, line)) {

    if (!line.length())
      continue;
    if (line[0] == '#')
      continue;

    stringstream ss(line);
    while (ss >> buf1 >> buf2 >> buf3) {
      v1.push_back(atof(buf1.c_str()));
      v2.push_back(atof(buf2.c_str()));
      v3.push_back(atof(buf3.c_str()));
    }
  }
}

// Four column
DataFile::DataFile(string const &dataFile, vector<double> &v1,
                   vector<double> &v2, vector<double> &v3, vector<double> &v4) {

  ifstream file(dataFile.c_str());

  //  int linecount=0;
  string buf1, buf2, buf3, buf4;
  string line;

  while (getline(file, line)) {

    if (!line.length())
      continue;
    if (line[0] == '#')
      continue;

    stringstream ss(line);
    while (ss >> buf1 >> buf2 >> buf3 >> buf4) {
      v1.push_back(atof(buf1.c_str()));
      v2.push_back(atof(buf2.c_str()));
      v3.push_back(atof(buf3.c_str()));
      v4.push_back(atof(buf4.c_str()));
    }
  }
}

// Four column, last a string
DataFile::DataFile(string const &dataFile, vector<double> &v1,
                   vector<double> &v2, vector<double> &v3, vector<string> &v4) {

  ifstream file(dataFile.c_str());

  //  int linecount=0;
  // int i=0;
  string buf1, buf2, buf3, buf4;
  string line, ThePath;

  while (getline(file, line)) {

    if (!line.length())
      continue;
    if (line[0] == '#')
      continue;

    // First non-zero, non-comment line expected to be a filter trace path.
    ThePath = line.c_str();
    // cout << "the path is " << ThePath << endl;
    break;
  }

  while (getline(file, line)) {

    if (!line.length())
      continue;
    if (line[0] == '#')
      continue;

    stringstream ss(line);
    while (ss >> buf1 >> buf2 >> buf3 >> buf4) {
      v1.push_back(atof(buf1.c_str()));
      v2.push_back(atof(buf2.c_str()));
      v3.push_back(atof(buf3.c_str()));
      if (buf4 == "NULL")
        v4.push_back(buf4.c_str());
      else
        v4.push_back(ThePath + buf4.c_str());
      // cout << i << " v4 is " << v4[i] << endl;
      // i++;
    }
  }
}

// Five column
DataFile::DataFile(string const &dataFile, vector<double> &v1,
                   vector<double> &v2, vector<double> &v3, vector<double> &v4,
                   vector<double> &v5) {

  ifstream file(dataFile.c_str());

  //  int linecount=0;
  string buf1, buf2, buf3, buf4, buf5;
  string line;

  while (getline(file, line)) {

    if (!line.length())
      continue;
    if (line[0] == '#')
      continue;

    stringstream ss(line);
    while (ss >> buf1 >> buf2 >> buf3 >> buf4 >> buf5) {
      v1.push_back(atof(buf1.c_str()));
      v2.push_back(atof(buf2.c_str()));
      v3.push_back(atof(buf3.c_str()));
      v4.push_back(atof(buf4.c_str()));
      v5.push_back(atof(buf5.c_str()));
    }
  }
}

// Five column, Last a string
DataFile::DataFile(string const &dataFile, vector<double> &v1,
                   vector<double> &v2, vector<double> &v3, vector<double> &v4,
                   vector<string> &v5) {

  ifstream file(dataFile.c_str());

  //  int linecount=0;
  string buf1, buf2, buf3, buf4, buf5;
  string line, ThePath;

  while (getline(file, line)) {

    if (!line.length())
      continue;
    if (line[0] == '#')
      continue;

    // First non-zero, non-comment line expected to be a filter trace path.
    ThePath = line.c_str();
    // cout << "the path is " << ThePath << endl;
    break;
  }

  while (getline(file, line)) {

    if (!line.length())
      continue;
    if (line[0] == '#')
      continue;

    stringstream ss(line);
    while (ss >> buf1 >> buf2 >> buf3 >> buf4 >> buf5) {
      v1.push_back(atof(buf1.c_str()));
      v2.push_back(atof(buf2.c_str()));
      v3.push_back(atof(buf3.c_str()));
      v4.push_back(atof(buf4.c_str()));
      if (buf5 == "NULL")
        v5.push_back(buf5.c_str());
      else
        v5.push_back(ThePath + buf5.c_str());
    }
  }
}
