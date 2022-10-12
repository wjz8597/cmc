//
// Burkhard Militzer                                    Urbana 4-9-99
//
// Pair Correlation function g(r)
//
#ifndef _GOFR_
#define _GOFR_

#include "Histogram.h"

class GOfR {
 private:
  Histogram histogram;
  //  int number;

 public:
  string name;

  void Init(const string s, const int n, const double maxSeparation) {
    histogram.Init(n,0.0,maxSeparation);  
    name  =s;
  }
  void Add(double y, double norm=1.0) {
    histogram.Add(y,norm);
  }
  friend istream& operator>>(istream &is, const GOfR &g){
    return is;
  } 
  void SetToZero() {
    histogram.SetToZero();
  }
  double GetPoint(const int i) const {
    return histogram.GetPoint(i);
  }
  double GetMidPoint(const int i) const {
    return histogram.GetMidPoint(i);
  }
  double GetFirstMidPoint() const {
    return histogram.GetMidPoint(0);
  }
  double GetLastMidPoint() const {
    return histogram.GetMidPoint(histogram.Size()-1);
  }
  int Size() const {
    return histogram.Size();
  }
  double* Pointer() {
    return histogram.histogram.v;
  }
  Array1 <double> & GetArray() {
    return histogram.histogram;
  }
  double GetNorm() const {
    return histogram.GetBasicNorm();
  }
  double GetNumber() const {
    return histogram.GetNumber();
  }

  //Normalize as if all entries had been positive
  void Normalize(const double y) {
    // There could be some empty g(r) if there is only one particle of one species
    if (histogram.GetNumber()==0) return; 
    double dR = histogram.GetDeltaR();  
    for(int i=0;i<histogram.Size(); ++i) {
      double r1=i*dR;
      double r2=(i+1)*dR;
      double v=4.0/3.0*pi*(r2*r2*r2-r1*r1*r1);
      double f=y/(v*histogram.GetNumber());
      histogram.histogram[i] *= f;
    }
  }
};

#endif // _GOFR_
