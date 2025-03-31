#include "StringManip.h"

using namespace std;

namespace StringManip
{

// Test for file existence.
bool
FileExists (const string FileName)
{
  bool bExists = true;
  struct stat buffer;
  if (stat (FileName.c_str (), &buffer) == -1)
    bExists = false;
  return bExists;
}

vector<string>
strsplit (const string &str, const string &delimiters)

{

  string::size_type lastPos = str.find_first_not_of (delimiters, 0);
  string::size_type pos = str.find_first_of (delimiters, lastPos);
  vector<string> parts;

  while (string::npos != pos || string::npos != lastPos)
    {
      parts.push_back (str.substr (lastPos, pos - lastPos));
      lastPos = str.find_first_not_of (delimiters, pos);
      pos = str.find_first_of (delimiters, lastPos);
    }

  return parts;
}

}
