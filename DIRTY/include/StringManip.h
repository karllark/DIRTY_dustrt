// Why does code (vtos()) live here instead of cpp? Because I couldn't
// get an external declaration for the template to work.  This works.
// fucky-you.  Unless you know how to clean it up... then let me know.
#ifndef __STRINGMANIP_H__
#define __STRINGMANIP_H__

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace std;

namespace StringManip
{

vector<string> strsplit (const string &str, const string &delims = " ");

bool FileExists (const string FileName);
inline void
coutStart (const string message)
{
  cout << endl << "---> Start: " << message << endl;
}
inline void
coutStart1 (const string message)
{
  cout << endl << "   ---> Start: " << message << endl;
}
inline void
coutEnd (const string message)
{
  cout << "<--- Done: " << message << endl << endl;
}
inline void
coutEnd1 (const string message)
{
  cout << "   <--- Done: " << message << endl << endl;
}
inline void
coutM1 (const string message)
{
  cout << " * " << message << endl;
}
inline void
coutE1 (const string message)
{
  cout << " * Err: " << message << endl;
}
inline void
coutM2 (const string message)
{
  cout << "   - " << message << endl;
}
inline void
coutE2 (const string message)
{
  cout << "   - Err: " << message << endl;
}
inline void
coutM3 (const string message)
{
  cout << "      " << message << endl;
}
inline void
coutE3 (const string message)
{
  cout << "       Err: " << message << endl;
}

inline string
BoolVal (const bool _abool)
{
  if (_abool)
    return "true";
  else
    return "false";
}
inline string
ToLower (string _s)
{
  transform (_s.begin (), _s.end (), _s.begin (), (int (*) (int))tolower);
  return _s;
}

template <typename T>
string
vtos (const T &value)
{
  stringstream ss;
  ss << value;
  return ss.str ();
}
// Dep check

}
#endif
