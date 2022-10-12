/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Some routines for smooth interpolation                                 //
//                                                                         //
//  B. Militzer                                       Oxford, OH, 09-30-17 //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _SWITCH_
#define _SWITCH_

inline double SmoothFastSwitch(const double x, const double x1, const double x2) {
  if (x<=x1) return 0.0;
  if (x>=x2) return 1.0;
  double arg = (x-x1)/(x2-x1);  // cannot be reached in case x1==x2
  double y   = sin(arg*0.5*pi); // sin(0)=0, sin(pi/2)=1
  double yy  = y*y; // sin^2(0) = 0, sin^2(pi/2)=1, derivaitves are zero at arg=0 and arg=1 b/c we square the function
  return yy;
}

inline double SmoothFastSwitch(const double x, const double x1, const double x2, const double y1, const double y2) {
  return y1 + (y2-y1)*SmoothFastSwitch(x,x1,x2);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// return 0.0 for x << mid
// return 1.0 for x >> mid
inline double SmoothSlowSwitch(const double x, const double mid, const double width) {
  double arg = (x-mid)/width;
  double y = (1.0+tanh(arg) )/2.0; // tanh(-inf)=-1, tanh(+inf)=+1
  return y;
}

inline double SmoothSlowSwitch(const double x, const double mid, const double width, const double y1, const double y2) {
  double y = y1 + (y2-y1) * SmoothSlowSwitch(x,mid,width);
  return y;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline double DTanh(const double x) {
  double ex = exp(x);
  return 4.0/sqr(ex+1.0/ex);
}

inline void SmoothSlowSwitchDerivative(const double x, const double mid, const double width, double & y, double & dydx) {
  double arg = (x-mid)/width;
  y    = (1.0+tanh(arg) )/2.0; // tanh(-inf)=-1, tanh(+inf)=+1
  dydx = 0.5/width*DTanh(arg); 
}

inline double SmoothSlowSwitchDerivative(const double x, const double mid, const double width) {
  double arg = (x-mid)/width;
  double dydx = 0.5/width*DTanh(arg); 
  return dydx;
}

#endif // _SWITCH_
