//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Define a Point in n dimensions WITH periodic boundary conditions     //
//                                                                      //
// Burkhard Militzer                                    Paris 4-23-99   //
// Modified class structure                         Livermore 8-08-02   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef _PERIODICPOINT_
#define _PERIODICPOINT_

#include "BasicPoint.h"
#include "Array.h" // just needed for GetLatticeVectors(), which was just added for compatibilty with PeriodicPointGen

// Dirty trick to avoid the slow function floor() (BM 09-20-98)
// wrong if x<-1000.0 and its application fails if
// distances greater than 1000*boxsize occur
// Do not overwrite floor directly since this header could be included
// elsewhere where floor is still needed
inline double Floor(double x) {
  //  return (double(int((x)+1000.0)-1000));
  //   Turn off with
  return floor(x);
}
// the following works everywhere except for Linux
// #define Floor(x) (double(int((x)+1000.0)-1000))


inline double Sign(double x) {
#ifndef __GNUC__
  return 1.0-2.0*(double(int(x))-Floor(x));
#else
  // egcs requires this work around instead of the line above 
  double g = double(int(x));
  return 1.0-2.0*(g-Floor(x));
#endif
}

template<class T>
inline void Swap(T & a, T & b) {
  T c(a);
  a = b;
  b = c;
}

template <int DIM>
class PeriodicPoint:public BasicPoint <DIM> {
 public:
  PeriodicPoint() {};
  PeriodicPoint(const double d) : BasicPoint <DIM> (d) {}
  PeriodicPoint(const BasicPoint <DIM> & p) : BasicPoint <DIM> (p) {}
  PeriodicPoint(const PeriodicPoint & p) : BasicPoint <DIM> (p) {}
  PeriodicPoint(const ShortVector<double,DIM> & p) : BasicPoint <DIM> (p) {}
  PeriodicPoint(const double x, const double y, const double z) : BasicPoint <DIM> (x,y,z) {}

  // Do not apply p.b.c. to operator += -= + - because I forgot about it later
  // this caused serious bugs several times!!!

  // Similarly, disable Distance() and Distance2() and enforce the use of 
  // DistancePBC() and DistancePBC2() to avoid confusion.
  // Using the class name BasicPoint::Distance() and PeriodicPoint::Distance()
  // has proven to be a confusing strategy!!!
  // RE-ENABLE THEM NOW
    /*
  inline double static Distance2(const PeriodicPoint& p1, 
			     const PeriodicPoint& p2);
  inline double static Distance(const PeriodicPoint & p1, 
			    const PeriodicPoint & p2);
    */

  inline double static DistancePBC2(const PeriodicPoint& p1, const PeriodicPoint& p2) {
    double r2=0.0;
    for(int i=0;i<DIM;++i) {
      register double r=p1.v[i]-p2.v[i];
      p1.MapInside(r,i);
      r2 += r*r;
    }
#ifdef STANDARD_DEBUG
    if (r2>p1.GetMaxDistance()*p1.GetMaxDistance())
      error ("Boxsize error",r2,p1.GetMaxDistance()*p1.GetMaxDistance());
#endif
    return r2;
  }
  inline double static DistancePBC(const PeriodicPoint & p1, const PeriodicPoint & p2) {
    return sqrt(DistancePBC2(p1,p2));
  }
  static double AnglePBCSafelyDeg(const PeriodicPoint & r1, const PeriodicPoint & rM, const PeriodicPoint & r2) {
    PeriodicPoint d1 = DifferencePBC(r1,rM);
    PeriodicPoint d2 = DifferencePBC(r2,rM);
    return PeriodicPoint::AngleSafelyDeg(d1,d2);
  }

  inline static void MapInside(double & c, const int i) {
    c -= boxSize[i]*Floor(c*oneOverBoxSize[i]+0.5);
  }
  inline static void MapInside(double & c, const double cell) {
    c -= cell*Floor(c/cell+0.5);
  }
  inline void MapInside(const int i) {
    MapInside(this->v[i],i);
  }
  inline static void MapInside(PeriodicPoint<DIM> & p) {
    for(int i=0; i<DIM; ++i)
      MapInside(p[i],i);
  }
  inline static void MapInside(PeriodicPoint<DIM> & p, const PeriodicPoint<DIM> & cell) {
    for(int i=0; i<DIM; ++i)
      MapInside(p[i],cell[i]);
  }
  inline void MapInside() {
    MapInside(*this);
  }

  static void MapReducedInside(PeriodicPoint & dd) {
    for(int i=0; i<DIM; i++) {
      dd[i] -= floor(dd[i]+0.5);
    }
  }
  static void MapReducedInsideZeroToOne(PeriodicPoint & dd) {
    for(int i=0; i<DIM; i++) {
      dd[i] -= floor(dd[i]);
    }
  }

  static void MapInsideZeroToOne(PeriodicPoint & d) {
     PeriodicPoint dd = Reduced(d);
     MapReducedInsideZeroToOne(dd);
     d = Cartesian(dd);
  }

  void MapInsideZeroToOne() {
    MapInsideZeroToOne(*this);
  }

  inline bool IsInside(const int i, const double eps) const {
    return (fabs(this->v[i]) <= (1.0+eps)*GetHalfBoxSize(i));
  }
  inline bool IsInside(const double eps=1e-8) const {
    for(int i=0;i<DIM;++i) {
      if (!IsInside(i,eps)) return false;
    }
    return true;
  }

  // set this=MapInside(a-b)
  inline void Separation(const PeriodicPoint & p1, const PeriodicPoint & p2) {
    for(int i=0;i<DIM;++i) {
      this->v[i] = p1.v[i]- p2.v[i];
      MapInside(this->v[i],i);
    }
  }

  inline static PeriodicPoint DifferencePBC(const PeriodicPoint & p1, const PeriodicPoint & p2) {
    PeriodicPoint p;
    p.Separation(p1,p2);
    return p;
  }

  // Set this += MapInside(a+b)
  // Do not apply PBC to *this
  inline void AddSeparation(const PeriodicPoint & p1, const PeriodicPoint & p2) {
    register double r;
    for(int i=0; i<DIM; ++i) {
      r = p1.v[i] - p2.v[i];
      MapInside(r,i);
      this->v[i] += r;
    }
  }
  inline void SetMidPoint(const PeriodicPoint& p1, const PeriodicPoint& p2) {
    for(int i=0;i<DIM;++i) {
      // put the difference vector in r
      this->v[i] = p2.v[i] - p1.v[i];
      MapInside(this->v[i],i);
      // add half of the difference vector to p1
      this->v[i] = p1.v[i] + this->v[i]*0.5;
      MapInside(this->v[i],i);
    }
  }

  // Set this=MapInside(a+b)
  inline void SumPBC(const PeriodicPoint& p1, const PeriodicPoint& d) {
    for(int i=0;i<DIM;++i) {
      this->v[i] = p1.v[i] + d.v[i];
      MapInside(this->v[i],i);
    }
  }

  static void SetBoxSize(const PeriodicPoint & size) {
    maxDistance=0.0;
    maxBoxSize =0.0;
    minBoxSize =size[0];
    volume     =1.0;
    for(int i=0;i<DIM;++i) {
      if (size[i]<=0.0) error("Box size cannot be <=0.0",i,size[i]);
      boxSize[i] = size[i];
      halfBoxSize[i] = 0.5*boxSize[i];
      oneOverBoxSize[i] = 1.0/boxSize[i];
      maxDistance += halfBoxSize[i]*halfBoxSize[i];
      if (boxSize[i]>maxBoxSize) maxBoxSize=boxSize[i];
      if (boxSize[i]<minBoxSize) minBoxSize=boxSize[i];
      volume *= boxSize[i];
    }
    maxDistance=sqrt(maxDistance)*(1.0+1e-10);
  }

  static void SetBoxSize(const double size) {
     PeriodicPoint p(size);
     SetBoxSize(p);
  }
  static void PrintCellInfo(const bool shortFlag=false) { // dummy argument
    cout << "|a|= " << boxSize[0] << endl;
    cout << "|b|= " << boxSize[1] << endl;
    cout << "|c|= " << boxSize[2] << endl;
    cout << "alpha= " << 90.0 << endl;
    cout << "beta=  " << 90.0 << endl;
    cout << "gamma= " << 90.0 << endl << endl;
  }

  inline static PeriodicPoint GetBoxSize() {
    PeriodicPoint p;
    for(int i=0;i<DIM;++i)
      p.v[i]=boxSize[i];
    //  cout << "Making a point from box" << p << endl;
#ifdef STANDARD_DEBUG
    if (p.Norm()==0.0)
      error("Box size not set");
#endif
    return p;
  }
  inline static double GetBoxSize(const int iDim) {
    Limits(iDim,DIM);
    return boxSize[iDim];
  }
  //  inline static double* GetBoxSizePointer() { // should be avoided
  //      return boxSize;
  //  }
  inline static double GetMaxBoxSize() {
    return maxBoxSize;
  }
  inline static double GetMinBoxSize() {
    return minBoxSize;
  }
  static double MinimumImageDistance() { // just for compatibility with PeriodicPointGen.h
    return GetMinimumImageDistance();
  }
  static double GetMinimumImageDistance() { // just for compatibility with PeriodicPointGen.h
    return GetMinBoxSize(); 
  }
  inline static double GetIsoBoxSize(const double eps=1e-12) {
    if (fabs(maxBoxSize-minBoxSize)>eps*maxBoxSize) error("Box is not cubic",minBoxSize,maxBoxSize,maxBoxSize-minBoxSize);
    return maxBoxSize;
  }
  inline static PeriodicPoint GetLatticeVector(const int iDim) {
    PeriodicPoint p(0.0);
    p[iDim] = boxSize[iDim];
    return p;
  }
  inline static Array1 <PeriodicPoint> GetLatticeVectors() {
    Array1 <PeriodicPoint> p(DIM);
    for(int d=0; d<DIM; d++) {
      p[d] = GetLatticeVector(d);
    }
    return p;
  }
  inline static PeriodicPoint GetReciprocalLatticeVector(const int iDim) {
    PeriodicPoint p(0.0);
    p[iDim] = 2.0*pi/boxSize[iDim];
    return p;
  }
  inline static PeriodicPoint GetLatticeVectorNorms() {
    return GetBoxSize();
  }
  inline static bool LatticeVectorsMinimal(const PeriodicPoint & a, const PeriodicPoint & b,const PeriodicPoint & c) {
    return true; // do not perform a real check, assuem they are all orthogonal
  }
  // smallest first
  static void Order(PeriodicPoint & a, PeriodicPoint & b) {
    double an=a.Norm();
    double bn=b.Norm();
    if (an>bn) { Swap(a,b); }
  }
  // smallest first
  static void Order(PeriodicPoint & a, PeriodicPoint & b, PeriodicPoint & c) {
    double an=a.Norm();
    double bn=b.Norm();
    double cn=c.Norm();
    if (an>bn) { Swap(a,b); Swap(an,bn); }
    if (an>cn) { Swap(a,c); Swap(an,cn); } // a is now the smallest
    if (bn>bn) { Swap(b,c); }              // c is now the largest
  }
  static double MaximalLatticeVectorRatio(PeriodicPoint & aaa, PeriodicPoint & bbb, PeriodicPoint & ccc) {
    Order(aaa,bbb,ccc); // aa is now the smallest vector
    return ccc.Norm() / aaa.Norm();
  }

  static void SetLatticeVectors(const Array1 <PeriodicPoint> & p, const bool print=false) {
    PeriodicPoint pp;
    double eps=1e-12;
    for(int i=0;i<DIM;++i) {
      if (p[i][i]<=0.0) error("Lattice vector element negative:",p[i]);
      double pn = p[i].Norm();
      if (fabs(pn-p[i][i])>pn*eps) error("Compiled with Point==PeriodicPoint for speed. So lattice vectors cannot have off-diagonal element:",p[i]);
      pp[i]=p[i][i];
    }
    SetBoxSize(pp);
  }
  inline static PeriodicPoint GetHalfBoxSize() {
    PeriodicPoint p;
    for(int i=0;i<DIM;++i)
      p.v[i]=halfBoxSize[i];
    //  cout << "Making a point from box" << p << endl;
#ifdef STANDARD_DEBUG
    if (p.Norm()==0.0)
      error("Box size not set");
#endif
    return p;
  }
  inline static double GetHalfBoxSize(const int iDim) {
    Limits(iDim,DIM);
    return halfBoxSize[iDim];
  }
  //  inline static double* GetHalfBoxSizePointer() { // should be avoided
  //    return halfBoxSize;
  //  }
  inline static double GetUnitCellAngle(const int i) {
    return pi/2;
  }
  inline static double GetUnitCellAngleDegrees(const int i) {
    return GetUnitCellAngle(i)*180.0/pi;
  }
  inline static double GetMaxDistance() {
      return maxDistance;
  }
  inline static double GetVolume() {
      return volume;
  }
  inline static double Reduced(const double x) {
    if (!Isotropic()) error("Cannot scale scalar in non-uniform cell");
    return x/boxSize[0];
  }
  inline static double Cartesian(const double x) {
    if (!Isotropic()) error("Cannot scale scalar in non-uniform cell");
    return x*boxSize[0];
  }
  inline static PeriodicPoint Reduced(const PeriodicPoint & p) {
    PeriodicPoint pp;
    for(int d=0; d<DIM; d++) 
      pp[d] = p[d]/boxSize[d];
    return pp;
  }
  PeriodicPoint Reduced() const {
    return Reduced(*this);
  }
  inline static PeriodicPoint Cartesian(const PeriodicPoint & p) {
    PeriodicPoint pp;
    for(int d=0; d<DIM; d++) 
      pp[d] = p[d]*boxSize[d];
    return pp;
  }
  PeriodicPoint Cartesian() const {
    return Cartesian(*this);
  }
  inline void Reduce() {
    *this = PeriodicPoint::Reduced(*this);
  }
  inline void ToCartesian() {
    *this = PeriodicPoint::Cartesian(*this);
  }

  inline static bool Cubic() {
    return Isotropic() && Orthogonal();
  }
  inline static bool Isotropic() {
    for(int d=1; d<DIM; ++d) 
      if (boxSize[d]!=boxSize[0]) return false;
    return true;
  }
  inline static bool Orthogonal() {
    return true;
  }
 public:
  static double boxSize[DIM];     // L_x, L_y, L_z for periodic box.
  static double oneOverBoxSize[DIM];    
  static double halfBoxSize[DIM]; // L_x/2, L_y/2, L_z/2
  static double maxDistance;            // maximal possible separation
  static double maxBoxSize;             // maximal direction
  static double minBoxSize;             // maximal direction
  static double volume;                 // Lx*Ly*Lz
};

#endif  // _PERIODICPOINT_

