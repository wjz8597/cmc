////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Standard definitions and functions that are used almost everywhere         //
//                                                                            //
// Burkhard Militzer                                        Urbana 4-9-99     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef _STANDARD_
#define _STANDARD_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
using namespace std;

/////////////////////////////////////////////////////////////////////////////////

// Define RANGE_CHECKING for testing
// #define RANGE_CHECKING

//For SLOW testing, define STANDARD_DEBUG
// #define STANDARD_DEBUG

////////////////////////////////////////////////////////////////////////////////

string IntToString(const int i);
string IntToString(const int i, const int d, const char c=' '); // was '0' before 9/20/17

inline int DigitsInt(const int i) { return int(ceil(log10(double(i+1)))); }
inline string IntToStringMaxNumber(const int i, const int maxNumber) {
  return IntToString(i,DigitsInt(maxNumber));
}

string DoubleToString(const double d);
bool IsInt(const string & s);
bool IsDouble(const string & s);

inline bool IsNumber(const string & s) {
  return IsInt(s) || IsDouble(s);
}
inline bool IsNotNumber(const string & s) {
  return !(IsNumber(s));
}

int StringToInt(const string & s);
double StringToDouble(const string & s);

string UpperCase(const string & s);
string LowerCase(const string & s);
string Spaces(const int n);

#ifdef __PGI // PG compiler bug
inline bool getline(istream & is, string & s) {
  s="";
  char c;
  while(is.get(c)) {
    if (c=='\n') return true;
    s += c;
  }
  return false;
}
inline bool getline(ifstream & is, string & s) {
  s="";
  char c;
  while(is.get(c)) {
    if (c=='\n') return true;
    s += c;
  }
  return false;
}
#endif

///////////////////////// redefine cout to write to a file /////////////////////
// #ifndef NO_COUT
// #define cout COUT
// extern ofstream COUT;
// #endif // NO_COUT

////////////////////////////////////////////////////////////////////////////////

#ifndef pi
const double pi=3.1415926535897932385;
#endif

// Is required for KCC but not on IBM using KCC3.3
// laters ones require it too.
/* #ifndef __rs6000__  */
/* #ifndef __osf__  */
// #if !defined __GNUC__ || (__GNUC__ == 2 || (__GNUC__ == 3 && __GNUC_MINOR__ < 1))
/*
inline double abs(double x) {  
   return fabs(x);  
} 
*/

// #endif
/* #endif  */
/* #endif  */

inline double sign(double x) {
  return (x>0.0) ? 1.0 : ((x<0.0) ? -1.0 : 0.0);
}

inline int sign(int x) {
  return (x>0) ? 1 : ((x<0) ? -1 : 0);
}

inline double nint(const double x) {
  return int(x+0.5*sign(x));
}

inline double min(double x, double y) {
  return (x<=y) ? x : y;
}

inline int min(int x, int y) {
  return (x<=y) ? x : y;
}

inline double max(double x, double y) {
  return (x>=y) ? x : y;
}

inline int max(int x, int y) {
  return (x>=y) ? x : y;
}

inline double sqr(double x) {
  return (x*x);
}

inline int sqr(int x) {
  return (x*x);
}

/*
inline int Factorial(const int n) { 
  int fac = 1;
  for(int i=2; i<=n; i++) {
    fac *= i;
  }
  return fac;
}
*/

// now returns a double because we had trouble with integers for n=12
inline double Factorial(const int n) { 
  double fac = 1.0;
  for(int i=2; i<=n; i++) {
    fac *= double(i);
  }
  return fac;
}

////////////////////////////////////////////////////////////////////////////////

// Write name of the variable and its value
#define write1(i)                    {cout << " "#i"= " << i; }
#define write2(i,j)                  {write1(i); write1(j);}
#define write3(i,j,k)                {write2(i,j); write1(k); }
#define write4(i,j,k,l)              {write3(i,j,k); write1(l); }
#define write5(i,j,k,l,m)            {write4(i,j,k,l); write1(m); }
#define write6(i,j,k,l,m,n)          {write5(i,j,k,l,m); write1(n); }
#define write7(i,j,k,l,m,n,o)        {write6(i,j,k,l,m,n); write1(o); }
#define write8(i,j,k,l,m,n,o,p)      {write7(i,j,k,l,m,n,o); write1(p); }
#define write9(i,j,k,l,m,n,o,p,q)    {write8(i,j,k,l,m,n,o,p); write1(q); }
#define write10(i,j,k,l,m,n,o,p,q,r) {write9(i,j,k,l,m,n,o,p,q); write1(r); }

#define Write(i)                     {write1(i); cout << endl;}
#define Write2(i,j)                  {write2(i,j); cout << endl;}
#define Write3(i,j,k)                {write3(i,j,k); cout << endl;}
#define Write4(i,j,k,l)              {write4(i,j,k,l); cout << endl;}
#define Write5(i,j,k,l,m)            {write5(i,j,k,l,m); cout << endl;}
#define Write6(i,j,k,l,m,n)          {write6(i,j,k,l,m,n); cout << endl;}
#define Write7(i,j,k,l,m,n,o)        {write7(i,j,k,l,m,n,o); cout << endl;}
#define Write8(i,j,k,l,m,n,o,p)      {write8(i,j,k,l,m,n,o,p); cout << endl;}
#define Write9(i,j,k,l,m,n,o,p,q)    {write9(i,j,k,l,m,n,o,p,q); cout << endl;}
#define Write10(i,j,k,l,m,n,o,p,q,r) {write10(i,j,k,l,m,n,o,p,q,r); cout << endl;}

void Terminate();

inline void Quit(const string & s) {
  cout << endl;
  cout << "Quitting the program early: " << s << endl;
  cout << endl;
  Terminate();
}

inline void WriteError(ostringstream & ss) {
  const string errorString = "Error   ";
  //  ss << ends;
  // cout is redirect into a file which might not yet be opened
  if (cout) {
    cout.precision(16);
    cout << errorString << ss.str() << endl;
  }
  cerr.precision(16);
  cerr << errorString << ss.str() << endl;
  Terminate();
}

inline void error(){
  ostringstream ss;
  ss << " - Program exited by calling error() function.";
  WriteError(ss);
}

inline void error(const string & m){
  ostringstream ss;
  ss << m;
  WriteError(ss);
}

template<class T> inline
void error(const string & m, const T& n){
  ostringstream ss;
  ss << m << " " << n;
  WriteError(ss);
}

template<class T, class U> inline
void error(const string & m, const T& t, const U& u){
  ostringstream ss;
  ss << m << " " << t << " " << u;
  WriteError(ss);
}

template<class T, class U, class V> inline
void error(const string & m, const T& t, const U& u, const V& v){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v;
  WriteError(ss);
}

template<class T, class U, class V, class W> inline
void error(const string & m, const T& t, const U& u, const V& v, const W& w){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w;
  WriteError(ss);
}

template<class T, class U, class V, class W, class X> inline
void error(const string & m, const T& t, const U& u, const V& v, const W& w, const X& x){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x;
  WriteError(ss);
}

template<class T, class U, class V, class W, class X, class Y> inline
void error(const string & m, const T& t, const U& u, const V& v, const W& w, const X& x, const Y& y){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x << " " << y;
  WriteError(ss);
}

template<class T, class U, class V, class W, class X, class Y, class Z> inline
void error(const string & m, const T& t, const U& u, const V& v, const W& w, const X& x, const Y& y, const Z& z){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x << " " << y << " " << z;
  WriteError(ss);
}

inline void WriteWarning(ostringstream & ss) {
  const string warningString = "WARNING   ";
  //  ss << ends;
#ifndef NO_COUT
  cout << warningString << ss.str() << endl;
#endif
  cerr << warningString << ss.str() << endl;
}

inline void warning() {
  ostringstream ss;
  ss << "...";
  WriteWarning(ss);
}

inline void warning(const string & m){
  ostringstream ss;
  ss << m;
  WriteWarning(ss);
}

template<class T> inline
void warning(const string & m, const T& t){
  ostringstream ss;
  ss << m << " " << t;
  WriteWarning(ss);
}

template<class T, class U> inline
void warning(const string & m, const T& t, const U& u){
  ostringstream ss;
  ss << m << " " << t << " " << u;
  WriteWarning(ss);
}

template<class T, class U, class V> inline
void warning(const string & m, const T& t, const U& u, const V& v){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v;
  WriteWarning(ss);
}

template<class T, class U, class V, class W> inline
void warning(const string & m, const T& t, const U& u, const V& v, const W& w){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w;
  WriteWarning(ss);
}

template<class T, class U, class V, class W, class X> inline
void warning(const string & m, const T& t, const U& u, const V& v, const W& w, const X& x){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x;
  WriteWarning(ss);
}

template<class T, class U, class V, class W, class X, class Y> inline
void warning(const string & m, const T& t, const U& u, const V& v, const W& w, const X& x, const Y& y){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x << " " << y;
  WriteWarning(ss);
}

template<class T, class U, class V, class W, class X, class Y, class Z> inline
void warning(const string & m, const T& t, const U& u, const V& v, const W& w, const X& x, const Y& y, const Z& z){
  ostringstream ss;
  ss << m << " " << t << " " << u << " " << v << " " << w << " " << x << " " << y << " " << z;
  WriteWarning(ss);
}

///////////////////////// Functions for array bound checking ///////////////////

// Limits is an inline function that checks indices 
// ok is 0<= n < max
inline void Limits(const int n, const int max) {
#ifdef RANGE_CHECKING
  if ((n<0) || (n>=max)) {
    error("Array index out of range ",n,max);
    cerr << "Array error: Index out of range:  0<= " 
	 << n << " < " << max << "\n" ;
    Terminate();
  }
#endif // RANGE_CHECKING
}

// Limits is an inline function that checks indices 
// ok is 0<= n <= max
inline void LimitsInclusive(const int n, const int max) {
#ifdef RANGE_CHECKING
  if ((n<0) || (n>max)) {
    //    error("Array Error: Index out of range ",n,max);
    cerr << "Array error: Upper limit for index out of range:  0<= " 
	 << n << " <= " << max << "\n" ;
    Terminate();
  }
#endif // RANGE_CHECKING
}

inline void EqualLimits(const int max1, const int max2) {
#ifdef RANGE_CHECKING
  if (max1!=max2) {
    cerr << "Array copy error: array sizes not equal:" 
	 << max1 << "," << max2 << endl;
    Terminate();
  }
#endif // RANGE_CHECKING
}

inline void BiggerLimit(const int lower, const int upper) {
#ifdef RANGE_CHECKING
  if (lower>=upper) {
    cerr << "Sub-array limits error: lower limit not lower " 
	 << lower << "," << upper << endl;
    Terminate();
  }
#endif // RANGE_CHECKING
}


#endif // _STANDARD_
