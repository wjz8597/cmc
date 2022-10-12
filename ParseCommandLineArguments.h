////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Parse Command Line Arguments          by B. Militzer, Berkeley, CA 3/04/15 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef _PARSECOMMANDLINEARGUMENTS_
#define _PARSECOMMANDLINEARGUMENTS_

#include "Array.h"

inline Array1 <string> CopyCommandLineArgumentsToArray(const int argc, char **argv) {
  Array1 <string> arg(argc);
  for(int i=0; i<argc; i++) {
    arg[i] = argv[i];
  }
  return arg;
}

inline bool ParseCommandLineFlag(Array1 <string> & arg, const string flag) {
  if (flag.length()==0) error("ParseCommandLineFlag() called with empty string.");
  int ii = arg.Find(flag);
  if (ii==notFound) return false;
  arg.Delete(ii);
  return true;
}

inline bool ParseCommandLineArgument(Array1 <string> & arg, const string name, int & result) {
  if (name.length()==0) error("ParseCommandLineArgument() called with empty string.");
  int ii = arg.Find(name);
  if (ii==notFound) return false;
  if (ii==arg.Size()-1 || !IsInt(arg[ii+1])) error("command line argument must be followed by integer",name);
  result = StringToInt(arg[ii+1]);
  arg.Delete(ii,2);
  if (arg.Find(name)!=notFound) error("Cannot specify the same command line argument twice",name);
  return true;
}

inline bool ParseCommandLineArgument(Array1 <string> & arg, const string name, double & result) {
  if (name.length()==0) error("ParseCommandLineArgument() called with empty string.");
  int ii = arg.Find(name);
  if (ii==notFound) return false;
  if (ii==arg.Size()-1 || !IsDouble(arg[ii+1])) error("command line argument must be followed by floating point number",name);
  result = StringToDouble(arg[ii+1]);
  arg.Delete(ii,2);
  if (arg.Find(name)!=notFound) error("Cannot specify the same command line argument twice",name);
  return true;
}

inline bool ParseCommandLineArgument(Array1 <string> & arg, const string name, string & result) {
  if (name.length()==0) error("ParseCommandLineArgument() called with empty string.");
  int ii = arg.Find(name);
  if (ii==notFound) return false;
  if (ii==arg.Size()-1) error("command line argument must be followed by string",name);
  result = arg[ii+1];
  arg.Delete(ii,2);
  if (arg.Find(name)!=notFound) error("Cannot specify the same command line argument twice",name);
  return true;
}

inline bool ParseCommandLineArguments(Array1 <string> & arg, const string name, int & result1, int & result2) {
  if (name.length()==0) error("ParseCommandLineArgument() called with empty string.");
  int ii = arg.Find(name);
  if (ii==notFound) return false;
  if (ii>=arg.Size()-2 || !IsInt(arg[ii+1]) || !IsInt(arg[ii+2])) error("command line argument must be followed by two integers",name);
  result1 = StringToInt(arg[ii+1]);
  result2 = StringToInt(arg[ii+2]);
  arg.Delete(ii,3);
  if (arg.Find(name)!=notFound) error("Cannot specify the same command line argument twice",name);
  return true;
}

inline bool ParseCommandLineArgumentNoSpace(Array1 <string> & arg, const string name, string & result) {
  if (name.length()==0) error("ParseCommandLineArgument() called with empty string.");
  int ii = arg.Find(name);
  if (ii!=notFound) {

    if (ii==arg.Size()-1) error("command line argument must be followed by string",name);
    result = arg[ii+1];
    arg.Delete(ii,2);
    if (arg.Find(name)!=notFound) error("Cannot specify the same command line argument twice",name);
    return true;

  } else {

    for(int i=0; i<arg.Size(); i++) {
      if (arg[i].substr(0,name.length())==name) {
	result = arg[i].substr(name.length());
	arg.Delete(i);
	return true;
      }
    }
    return false;

  }
}

inline bool ParseCommandLineArgumentNoSpace(Array1 <string> & arg, const string name, int & result) {
  if (name.length()==0) error("ParseCommandLineArgument() called with empty string.");
  int ii = arg.Find(name);
  if (ii!=notFound) {

    if (ii==arg.Size()-1 || !IsInt(arg[ii+1])) error("command line argument must be followed by integer",name);
    result = StringToInt(arg[ii+1]);
    arg.Delete(ii,2);
    if (arg.Find(name)!=notFound) error("Cannot specify the same command line argument twice",name);
    return true;

  } else {

    for(int i=0; i<arg.Size(); i++) {
      if (arg[i].substr(0,name.length())==name) {
	string s = arg[i].substr(name.length());
	if (!IsInt(s)) error("command line argument must be followed by integer",name,s);
	result = StringToInt(s);
	arg.Delete(i);
	return true;
      }
    }
    return false;

  }
}

inline bool ParseCommandLineArgumentNoSpace(Array1 <string> & arg, const string name, double & result) {
  if (name.length()==0) error("ParseCommandLineArgument() called with empty string.");
  int ii = arg.Find(name);
  if (ii!=notFound) {

    if (ii==arg.Size()-1 || !IsDouble(arg[ii+1])) error("command line argument must be followed by floating point number",name);
    result = StringToInt(arg[ii+1]);
    arg.Delete(ii,2);
    if (arg.Find(name)!=notFound) error("Cannot specify the same command line argument twice",name);
    return true;

  } else {

    for(int i=0; i<arg.Size(); i++) {
      if (arg[i].substr(0,name.length())==name) {
	string s = arg[i].substr(name.length());
	if (!IsDouble(s)) error("command line argument must be followed by floating point number",name,s);
	result = StringToDouble(s);
	arg.Delete(i);
	return true;
      }
    }
    return false;

  }
}

#endif // _PARSECOMMANDLINEARGUMENTS_
