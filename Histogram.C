////////////////////////////////////////////////////////////////////////////
//
// Histogram for integers and double data 
//
// Burkhard Militzer                                    Urbana 11-9-99
//

#include "Histogram.h"

void Histogram::Init(const int n, const double xx1, const double xx2) {
  if (xx1>=xx2)
    error("Histogram: R_min must be less than R_max ",xx1,xx2);
  nPts=n;
  x1=xx1;
  x2=xx2;
  deltaRInv = (double)(nPts)/(x2-x1);
  deltaR = 1.0/deltaRInv;
  offset = -deltaRInv*x1;
  histogram.Resize(nPts);
  histogram.Set(0.0);
  norm=0.0;
  number=0;
}

void IntHistogram::Init(const int nn1, const int nn2) {
  if (nn1>nn2)
    error("IntHistogram: I_min must be less than I_max ",nn1,nn2);
  nPts=nn2-nn1+1;
  n1=nn1;
  n2=nn2;
  histogram.Resize(nPts);
  histogram.Set(0.0);
  norm=0.0;
  number=0;
}

void BasicHistogram::SetToZero() {
  histogram.Set(0.0);
  norm  =0.0;
  number=0;
}

void Histogram::Normalize(const double argNorm) {
  if (fabs(norm)<1e-20) {
    //    cout << "Warning Histogram::Normalize: norm is zero\n";
    return;
  }
  double f= argNorm/norm*deltaRInv;
  (*this).Multiply(f);
}

void IntHistogram::Normalize(const double argNorm) {
  if (fabs(norm)<1e-20) {
    //    cout << "Warning Histogram::Normalize: norm is zero\n";
    return;
  }
  double f= argNorm/norm;
  (*this).Multiply(f);
}


//void BasicHistogram::Multiply (const double f) {
void BasicHistogram::Multiply (double f) {
  for(int i=0;i<nPts;++i) {
    histogram[i] *= f;
  }
  norm*=f;
}

double Histogram::GetNorm() const {
  return norm*deltaR;
}

double IntHistogram::GetNorm() const {
  return norm;
}

ostream& operator<<(ostream &os, const Histogram & h) {
  if (h.Size()==0) {
    os << "Histogram undefined\n";
    return os;
  }
  os << "# r1= " << h.x1 << " r2= " << h.x2 << " nPts= " << h.nPts << endl;
  os.setf(ios::scientific);
  for(int i=0;i<h.nPts;i++) {
    double x=h.x1+h.deltaR*(double)(i);
    os << x << " " << h.histogram[i] << endl;
  }
  os << "# ------------------\n";
  return os;
}

ostream& operator<<(ostream &os, const IntHistogram & h) {
  if (h.Size()==0) {
    os << "Histogram undefined\n";
    return os;
  }
  os << "# n1= " << h.n1 << " n2= " << h.n2 << " nPts= " << h.nPts << endl;
  os.setf(ios::scientific);
  for(int i=0;i<h.nPts;i++) {
    os << i+h.n1 << " " << h.histogram[i] << endl;
  }
  os << "# ------------------\n";
  return os;
}
