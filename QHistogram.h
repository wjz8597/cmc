//
// Burkhard Militzer                                   Berkeley, 10-02-18
//
// Q histogram
//
#ifndef _QHISTOGRAM_
#define _QHISTOGRAM_

#include "Histogram.h"

class QHistogram {

  // private:
 public:
  Histogram histogram;
  //  int number;

 public:
  string name;

  void Init(const string s, const int n, const double QMin, const double QMax) {
    histogram.Init(n,QMin,QMax);  
    name=s;
  }
  void Add(double y, double norm=1.0) {
    histogram.Add(y,norm);
  }
  void Normalize(const double norm=1.0) {
    histogram.Normalize(norm);
  }
  friend istream& operator>>(istream &is, const QHistogram & qh){
    return is;
  } 
  void SetToZero() {
    histogram.SetToZero();
  }

  double GetPoint(const int i) const { return histogram.GetPoint(i); }
  double GetMidPoint(const int i) const { return histogram.GetMidPoint(i); }
  double GetFirstMidPoint() const { return histogram.GetMidPoint(0); }
  double GetLastMidPoint() const { return histogram.GetMidPoint(histogram.Size()-1); }
  int Size() const { return histogram.Size(); }
  double* Pointer() { return histogram.histogram.v; }
  Array1 <double> & GetArray() { return histogram.histogram; }

  double GetNorm() const {
    return histogram.GetBasicNorm();
  }
  double GetNumber() const {
    return histogram.GetNumber();
  }

  double GetValue(const int i) const {
    return histogram.GetValue(i);
  }
 
};

#endif // _QHISTOGRAM_
