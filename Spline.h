/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Spline interpolation                                                    // 
//                                                                         //
// Burkhard Militzer                                  Livermore 4-27-01    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _SPLINE_
#define _SPLINE_

#include "Array.h"

class Spline {
 public:
  Spline():n(0){};
  Spline(const int m):n(0){ // just pre-allocate memory, do not fill in any elements yet
    x.Allocate(m);
    y.Allocate(m);
    y2.Allocate(m);
  };
  Spline(const Array1 <double> & xx, const Array1 <double> & yy):n(0) {
    PushBack(xx,yy);
    Init();
  }

  static void CalculateY2(Array1 <double> & x,
			  Array1 <double> & y, 
			  const int n,
			  Array1 <double> & y2,
			  const bool natural=true, 
			  const double ypLow=0.0, 
			  const double ypHigh=0.0);
  void CalculateY2(const bool natural=true, 
		   const double ypLow=0.0, 
		   const double ypHigh=0.0) {
    CalculateY2(x,y,n,y2,natural,ypLow,ypHigh);
  }
  
  void Init(const bool natural=true, 
	    const double ypLow=0.0, 
	    const double ypHigh=0.0) {
    CalculateY2(x,y,n,y2,natural,ypLow,ypHigh);
  }
  
  // use finite difference to specify dy/dx at the boundaries
  void InitUsingFiniteDifferences() {
    if (n<=2) error("Spline:InitUsingFiniteDifferences is called with two few data points",n);
    Sort(x,y,n);
    double ypLow =(y[1]-y[0])    /(x[1]-x[0]);
    double ypHigh=(y[n-1]-y[n-2])/(x[n-1]-x[n-2]);
    Init(false,ypLow,ypHigh);
  }

  static void Sort(Array1 <double> & x,
		   Array1 <double> & y, 
		   const int n);

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

  int GetNumber() const {
    return x.Size();
  }

  void Print() const {
    for(int i=0; i<GetNumber(); i++) {
      Write3(i,x[i],y[i]);
    }
  }

  static double IntegralTerm1(const double x1, const double x2, const double xx);
  static double IntegralTerm2(const double x1, const double x2, const double xx);
  static double IntegralTerm3(const double x1, const double x2, const double xx);
  static double IntegralTerm4(const double x1, const double x2, const double xx);
  static double IntegralInterval1(const double x1, const double x2, const double xx, const double y1, const double y2, const double ys1, const double ys2);
  static double IntegralFullInterval(const double h, const double y1, const double y2, const double ys1, const double ys2);
  static double IntegralInterval2(const double x1, const double x2, const double xx, const double y1, const double y2, const double ys1, const double ys2);

  static double Integrate(const Array1 <double> & x,
			  const Array1 <double> & y, 
			  const int n, 
			  const Array1 <double> & y2,
			  const double xx1, 
			  const double xx2);
  double Integrate(const double xx1, const double xx2) const {
    return Integrate(x,y,n,y2,xx1,xx2);
  }


  // protected:

  static int GetIndex(const Array1 <double> & x, const int n, const double xx);
 
  bool operator==(const Spline & s) const {
    return ((n==s.n) && (x==s.x) && (y==s.y) && (y2==s.y2));
  }
  bool operator!=(const Spline & s) const {
    return !(*this==s);
  }

  int n;
  Array1 <double> x;
  Array1 <double> y;
  Array1 <double> y2;
  static Array1 <double> u;

};

#endif // _SPLINE_
