/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// Define a Point in n dimensions without period boundary conditions           //
//                                                                             //
// Burkhard Militzer                                Paris 4-23-99              //
//                                                  Livermore 01-09-03         //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

#ifndef _BASICPOINT_
#define _BASICPOINT_

#include "ShortVector.h"

const int NDIM=3;

template <int DIM>
class BasicPoint : public ShortVector<double,DIM> {
 public:
  static const int nDim = DIM;

  explicit BasicPoint() {};
  explicit BasicPoint(const double p) {
    for(int i=0; i<DIM; i++) 
      this->v[i]=p;
  }
  explicit BasicPoint(const BasicPoint & p) {
    for(int i=0; i<DIM; i++) 
      this->v[i]=p[i];
  }

  BasicPoint & operator=(const double p) {   
    for(int i=0; i<DIM; i++) 
      this->v[i]=p;
    return *this;
  }
  BasicPoint & operator=(const BasicPoint & p) {   
    for(int i=0; i<DIM; i++) 
      this->v[i]=p[i];
    return *this;
  }


  BasicPoint(const ShortVector<double,DIM> & p) {
    for(int i=0; i<DIM; i++) 
      this->v[i]=p[i];
  }
  inline BasicPoint(const double x, const double y);
  inline BasicPoint(const double x, const double y, const double z);
  
  /*
  inline BasicPoint(const double x);
  inline BasicPoint(const double x, const double y);
  inline BasicPoint(const double x, const double y, const double z);
  */

  static double Distance2(const BasicPoint & p1, const BasicPoint & p2) {
    double r2=0.0;
    for(int i=0;i<DIM;++i) {
      register double r=p1.v[i]-p2.v[i];
      r2 += r*r;
    }
    return r2;
  }
  static double Distance(const BasicPoint & p1, const BasicPoint & p2) {
      return sqrt(Distance2(p1,p2));
  }
  static BasicPoint Difference(const BasicPoint & p1, const BasicPoint & p2) {
    return p1-p2;
  }
  static double CosAngle(const BasicPoint & x, const BasicPoint & y) {
    return x*y/(x.Norm()*y.Norm());
  }
  static double Angle(const BasicPoint & x, const BasicPoint & y) {
      return acos(CosAngle(x,y));
  }
  static double AngleSafely(const BasicPoint & x, const BasicPoint & y) {
    double xn=x.Norm();
    double yn=y.Norm();
    if (xn*yn==0.0) error("Cannot calculate angle between vectors if one has length zero.",x,y);
    double ca = x*y/(xn*yn);
    if (ca>+1.0) ca=+1.0;
    if (ca<-1.0) ca=-1.0;
    return acos(ca);  
  }
  static double AngleSafelyDeg(const BasicPoint & x, const BasicPoint & y) {
    return AngleSafely(x,y)*180.0/pi;
  }
  static double AngleSafelyDeg90m90(const BasicPoint & x, const BasicPoint & y) {
    double a = AngleSafelyDeg(x,y);
    return (a>90.0) ? a-180.0 : a;
  }
  inline static BasicPoint CrossProduct(const BasicPoint & a, const BasicPoint & b) {
    BasicPoint p;
    p[0] = a[1]*b[2]-a[2]*b[1];
    p[1] = a[2]*b[0]-a[0]*b[2];
    p[2] = a[0]*b[1]-a[1]*b[0];
    return p;
  }

};

/////////////////////////////////////////////////////////////////////////////////

template <>
inline BasicPoint<1>::BasicPoint(const double x) {
  v[0]=x;
}

template <>
inline BasicPoint<2>::BasicPoint(const double x, const double y) {
  v[0]=x;
  v[1]=y;
}

template <>
inline BasicPoint<3>::BasicPoint(const double x, const double y, const double z) {
  v[0]=x;
  v[1]=y;
  v[2]=z;
}

/////////////////////////////////////////////////////////////////////////////////

#endif  // _BASICPOINT_

