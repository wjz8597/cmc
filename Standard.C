////////////////////////////////////////////////////////////////////////////
//
// Standard definitions and functions that are used almost everywhere
//
// Burkhard Militzer                                    Urbana 4-9-99
//
///////////////////////////////////////////////////////////////////////////

#include "Standard.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

void Terminate() {
  
#ifdef USE_MPI
  int errorcode=1;
  MPI_Abort(MPI_COMM_WORLD, errorcode);
#endif

   cerr << "Exiting now\n";
   exit(-1);
}

///////////////////////////////////////////////////////////////////////////

ofstream COUT;

// Does this need "ss << ends" ???
string IntToString(const int i) {
  ostringstream ss;
  ss << i;
  return ss.str();
}

// Does this need "ss << ends" ???
string IntToString(const int i, const int d, const char c) {
  ostringstream ss;
  ss << i;
  string s;
  for(int j=0; j<d-int(ss.str().length()); j++) 
    //    s+="0";
    s+=c;
  return s+ss.str();
}

string DoubleToString(const double d) {
  ostringstream ss;
  ss << d;
  return ss.str();
}

bool IsInt(const string & s) {
  const string intChars(" +-0123456789");
  string::size_type j=s.find_first_not_of(intChars);
  return (j == string::npos);
}

int StringToInt(const string & s) {
  const string intChars(" +-0123456789");
  string::size_type j=s.find_first_not_of(intChars);
  if (j != string::npos)
    error("Parser: Not an integer ",s,j+1);
 
  istringstream ss(s);
  int value;
  ss >> value;
  return value;
}

bool IsDouble(const string & s) {
  const string doubleChars(" +-0123456789.edDE");
  string::size_type j=s.find_first_not_of(doubleChars);
  return (j == string::npos);
}

double StringToDouble(const string & s) {
  const string doubleChars(" +-0123456789.edDE");
  string::size_type j=s.find_first_not_of(doubleChars);
  if (j != string::npos)
    error("Parser: Not a double ",s,j+1);
  
  istringstream ss(s);
  double value;
  ss >> value;
  return value;
}

string UpperCase(const string & s) {
  string sl;
  for(string::const_iterator p=s.begin();p!=s.end();++p) {
    sl += toupper(*p);
  }
  return sl;
}
 
string LowerCase(const string & s) {
  string sl;
  for(string::const_iterator p=s.begin();p!=s.end();++p) {
    sl += tolower(*p);
  }
  return sl;
}

