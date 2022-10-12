/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Modified Spline interpolation with linear extrapolation                 //
//                                                                         //
// Burkhard Militzer                                 Washington 08-03-04   //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _LSPLINE_
#define _LSPLINE_

#include "Array.h"

class LSpline {
 public:
  LSpline():n(0){};
  LSpline(const int m):n(0){ // just pre-allocate memory, do not fill in any elements yet
    x.Allocate(m);
    y.Allocate(m);
    y2.Allocate(m);
  };

  static void CalculateY2(Array1 <double> & x,
			  Array1 <double> & y, 
			  const int n,
			  Array1 <double> & y2,
			  const bool natural=true, 
			  const double ypLow=0.0, 
			  const double ypHigh=0.0);
  // changed the default
  // no longer us natual splines with y''=0
  // but match the first derivative to the linear extrapolation
  // use finite difference to specify dy/dx at the boundaries
  void CalculateY2UsingFiniteDifferences() { 
    if (n<2) error(":Spline::CalculateY2() is called with two few data points",n);
    Sort(x,y,n);
    double ypLow  = (y[1]  -y[0]  )/(x[1]  -x[0]  );
    double ypHigh = (y[n-1]-y[n-2])/(x[n-1]-x[n-2]);
    CalculateY2(x,y,n,y2,false,ypLow,ypHigh);
  }
  void Init() {
    CalculateY2UsingFiniteDifferences();
  }

  void CalculateY2(const bool natural,
		   const double ypLow,
		   const double ypHigh) {
    CalculateY2(x,y,n,y2,natural,ypLow,ypHigh);
  }

  // old natural splines can still be used in the following way:
  void InitNatural() {
    CalculateY2(x,y,n,y2,true,0.0,0.0);
  }
  
  // use finite difference to specify dy/dx at the boundaries
  void InitUsingFiniteDifferences() {
    Init(); // now default
  }

  static void Sort(Array1 <double> & x,
		   Array1 <double> & y, 
		   const int n);

  void Sort() {
    Sort(x,y,n);
  }

  void Reset() {
    n=0;
    x.Resize(0);
    y.Resize(0);
  }

  void PushBack(const double xx, const double yy) {
    n++;
    x.PushBack(xx);
    y.PushBack(yy);
  }

  void PushBack(const Array1 <double> & xx, const Array1 <double> & yy) {
    if (xx.Size()!=yy.Size()) 
      error("Array size error ",xx.Size(),yy.Size());
    for(int i=0; i<xx.Size(); i++) {
      PushBack(xx[i],yy[i]);
    }
  }

  static double Interpolate(const Array1 <double> & x,
			    const Array1 <double> & y, 
			    const int n, 
			    const Array1 <double> & y2,
			    const double xx);
  static double InterpolateLinearly(const Array1 <double> & x,
				    const Array1 <double> & y, 
				    const int n, 
				    const Array1 <double> & y2,
				    const double xx);
  static double Derivative(const Array1 <double> & x,
			   const Array1 <double> & y, 
			   const int n, 
			   const Array1 <double> & y2,
			   const double xx);
  static double SecondDerivative(const Array1 <double> & x,
				 const Array1 <double> & y, 
				 const int n, 
				 const Array1 <double> & y2,
				 const double xx);
  
  double Interpolate(const double xx) const {
    return Interpolate(x,y,n,y2,xx);
  }
  double InterpolateLinearly(const double xx) const {
    return InterpolateLinearly(x,y,n,y2,xx);
  }
  double Derivative(const double xx) const {
    return Derivative(x,y,n,y2,xx);
  }
  double SecondDerivative(const double xx) const {
    return SecondDerivative(x,y,n,y2,xx);
  }

  double f(const double xx) const {
    return Interpolate(xx);
  }

  double operator()(const double xx) const {
    return Interpolate(xx);
  }

  double FirstX() const {
    return x[0];
  }

  double LastX() const {
    return x[x.Size()-1];
  }

  double FirstY() const {
    return y[0];
  }

  double LastY() const {
    return y[y.Size()-1];
  }

  int Size() const {
    return GetNumber();
  }

  int GetNumber() const {
    return x.Size();
  }

  void Print() const {
    for(int i=0; i<GetNumber(); i++) {
      Write3(i,x[i],y[i]);
    }
  }

  // protected:

  static int GetIndex(const Array1 <double> & x, const int n, const double xx);
 
  bool operator==(const LSpline & s) const {
    return ((n==s.n) && (x==s.x) && (y==s.y) && (y2==s.y2));
  }
  bool operator!=(const LSpline & s) const {
    return !(*this==s);
  }

  int n;
  Array1 <double> x;
  Array1 <double> y;
  Array1 <double> y2;
  static Array1 <double> u;

};

#endif // _LSPLINE_
