////////////////////////////////////////////////////////////////////////////
//
// Histogram for integer and double data 
//
// Burkhard Militzer                                    Urbana 11-9-99
//
#ifndef _HISTOGRAM_
#define _HISTOGRAM_

#include "Standard.h"        
#include "Array.h"        

class BasicHistogram {
 public:
  Array1 <double> histogram;
 protected:
  int nPts;
  double norm;           // Sum of the weights of the entries made so far
  int number;           // Number of entries made so far

 public:
  BasicHistogram():nPts(0),norm(0.0),number(0){};
  void Multiply(const double f);
  void Normalize();
  void SetToZero();
  int GetNumber() const { return number; }
  int Size() const { return nPts; };
  double* GetPointer() { return histogram.GetPointer(); }
  Array1 <double> & GetArray() { return histogram; }
  double GetValue(const int i) const {
    return histogram[i];
  }
  
};

// Set up bins from n1 to n2 including the boundaries
class IntHistogram : public BasicHistogram {
 public:
  void Init(const int n1, const int n2);
  inline void Add(const int x, const double argNorm=1.0);
  inline void Add(const int x, double argNorm, const int n);
  void Normalize(const double argNorm);
  double GetNorm() const;
  friend ostream& operator<<(ostream &, const IntHistogram &);
  int GetPoint(const int i) const { return i+n1; }

 private:
  int n1,n2;
};

// Set up bins from xx1 to xx2 using dx=(xx2-xx1)/n, extra last bin any more
class Histogram : public BasicHistogram {
 public:
  void Init(const int n, const double xx1, const double xx2);
  inline void Add(const double x, const double norm=1.0);
  void Normalize(const double argNorm=1.0);
  double GetNorm() const;                   // is multiplied by deltaR 
  double GetBasicNorm() const {             // not multiplied by deltaR
    return norm;
  }
  friend ostream& operator<<(ostream &, const Histogram &);
  double GetDeltaR() const { return deltaR; }; 
  double GetPoint(const int i) const { return i*deltaR+x1; }
  double GetMidPoint(const int i) const { return (i+0.5)*deltaR+x1; }

 private:
  double deltaRInv,deltaR,offset,x1,x2;
};

inline void IntHistogram::Add(const int x, const double argNorm) {
  int index=x-n1;
  //  WriteS4(x,n1,n2,index);
  /*   Would prevent array error if index outside
  if (index<0) {
      index=0;
  } else {
    if (index>=nPts) {
	index=nPts-1;
    }
    }*/
  //  cout << "Index = " << index << endl;

  // do a check before adding
  if (index>=0 && index<nPts)
    histogram[index]+=argNorm;
  norm            +=argNorm;
  number++;
}

// make the very same entry n times
inline void IntHistogram::Add(const int x, double argNorm, const int n) {
  int index=x-n1;
  argNorm *= (double)(n);
  // do a check before adding
  if (index>=0 && index<nPts)
    histogram[index]+=argNorm;
  norm            +=argNorm;
  number += n;
}

inline void Histogram::Add(const double x, const double argNorm) {
  double xx=x*deltaRInv+offset;
  int index=int(floor(xx));
  /*
  if (index<0) {
      index=0;
  } else {
    if (index>=nPts) {
	index=nPts-1;
    }
  }
  */

  //  WriteS4(x,xx,index,nPts);
  // Check for overflow of norm since it is a double
  // double norm_old=norm;
  // do a check before adding
  if (index>=0 && index<nPts)
    histogram[index]+=argNorm;
  norm            +=argNorm;
  number++;
}

#endif // _HISTOGRAM_
