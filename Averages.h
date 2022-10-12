/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Average scalar or array data over time                              //
//                                                                     //
// Burkhard Militzer                               Washington 12-31-06 //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#ifndef _AVERAGES_
#define _AVERAGES_

#include "PrintErrorBar.h"

class ScalarAverage {
 public:
  double sn,sx,sx2;
  bool shrt;

 ScalarAverage():sn(0.0),sx(0.0),sx2(0.0),shrt(true) {}

  bool Valid() const {
    return sn>0.0;
  }
  void SetToZero() {
    sx = 0.0;
    sx2= 0.0;
    sn = 0.0;
  }

  void Add(const double x, const double f=1.0) {
    sn  += f;
    sx  += x*f;    
    sx2 += x*x*f;
  }

  /*
  void Normalize(const double nn=1.0) {
    if (sn==0.0) error("Norm=0 problem ");
    sn   = nn;
    sx  *= nn/sn;    
    sx2 *= nn/sn;
  }
  */

  double GetAverage() const {
    //    if (sn==0.0) error("Norm=0 problem ");
    if (sn==0.0) return 0.0;
    return sx/sn;
  }

  double GetVariance() const {
    //    if (sn==0.0) error("Norm=0 problem ");
    if (sn==0.0) return 0.0;
    double sxn  = sx/sn;
    double sx2n = sx2/sn;
    double var = sqrt(max(sx2n - sxn*sxn,0.0));
    return var;
  }
 
  double GetErrorBar() const {
    //    if (sn==0.0) error("Norm=0 problem ");
    if (sn==0.0) return 0.0;
    double sxn  = sx/sn;
    double sx2n = sx2/sn;
    double e = sqrt(max(sx2n - sxn*sxn,0.0)/sn);
    return e;
  }

 // pretend x+f was passed instead of x
  ScalarAverage operator+(const double f) const {
    ScalarAverage s=*this;
    s.sx2 += sn*f*f + 2.0*f*sx; // use old 'sx'
    s.sx  += f*sn;
    return s;
  }
 // pretend x-f was passed instead of x
  ScalarAverage operator-(const double f) const {
    return *this+(-f);
  }
 // pretend x*f was passed instead of x
  ScalarAverage operator*(const double f) const {
    ScalarAverage s=*this;
    s.sx  *= f;
    s.sx2 *= f*f;
    return s;
  }
  // pretend x*f was passed instead of x
  ScalarAverage operator/(const double f) const {
    ScalarAverage s=*this;
    s.sx  /= f;
    s.sx2 /= f*f;
    return s;
  }

  void SetShort() {
    shrt=true;
  }
  void SetLong() {
    shrt=false;
  }

  friend ostream& operator<<(ostream &os, const ScalarAverage & a ) {
    if (a.shrt) {
      os << PrintAverageResultShort(a.GetAverage(),a.GetErrorBar()) << " ";
    } else {
      os << a.GetAverage() << " +- " << a.GetErrorBar() << " ";
    }
    return os;
  }

};

class ArrayAverage {
 public:
  ArrayAverage(const int n_=0):n(n_),sn(0.0),sx(n,0.0),sx2(n,0.0),shrt(false) {}
  
  void Init(const string & name_, const int n_) {
    name = name_;
    Init(n_);
  }
  void Init(const int n_) {
    n = n_;
    sn = 0.0;
    sx.Resize(n);
    sx2.Resize(n);
    sx  = 0.0;
    sx2 = 0.0;
  }
  void SetToZero() {
    sx = 0.0;
    sx2= 0.0;
    sn = 0.0;
  }

  void Add(const Array1 <double> & x, const double f=1.0) {
    sn  += f;
    for(int i=0; i<n; i++) {
      sx[i]  += x[i]*f;    
      sx2[i] += x[i]*x[i]*f;
    }
  }

  /*
  void Normalize(const double nn=1.0) {
    if (sn==0.0) error("Norm=0 problem ");
    sn   = nn;
    sx  *= nn/sn;    
    sx2 *= nn/sn;
  }
  */

  double GetAverage(const int i) const {
    //    if (sn==0.0) error("Norm=0 problem ");
    if (sn==0.0) return 0.0;
    return sx[i]/sn;
  }

  Array1 <double> GetAverage() const {
    Array1 <double> a(n);
    for(int i=0; i<n; i++) 
      a[i] = GetAverage(i);
    return a;
  }

  double GetErrorBar(const int i) const {
    //    if (sn==0.0) error("Norm=0 problem ");
    if (sn==0.0) return 0.0;
    double sxn  = sx[i]/sn;
    double sx2n = sx2[i]/sn;
    double e = sqrt(max(sx2n - sxn*sxn,0.0)/sn);
    return e;
  }

  Array1 <double> GetErrorBar() const {
    Array1 <double> e(n);
    for(int i=0; i<n; i++) 
      e[i] = GetErrorBar(i);
    return e;
  }

  // pretend x*f was passed instead of x
  ArrayAverage & operator*(const double f) {
    sx  *= f;
    sx2 *= f*f;
    return *this;
  }
  // pretend x*f was passed instead of x
  ArrayAverage & operator/(const double f) {
    sx  /= f;
    sx2 /= f*f;
    return *this;
  }
  int Size() const {
    return n;
  }
  bool Valid() const {
    return sn>0.0;
  }
  int n;
  double sn;
  Array1 <double> sx,sx2;
  string name;

  bool shrt;
  void SetShort() {
    shrt=true;
  }
  void SetLong() {
    shrt=false;
  }

  friend ostream& operator<<(ostream &os, const ArrayAverage & a ) {
    for(int i=0; i<a.n; i++) {
      os << i << "  ";
      if (a.shrt) {
	os << PrintAverageResultShort(a.GetAverage(i),a.GetErrorBar(i)) << " ";
      } else {
	os << a.GetAverage(i) << " +- " << a.GetErrorBar(i) << " ";
      }
    }
    return os;
  }

};

#endif // _AVERAGES_
