/////////////////////////////////////////////////////////////////////////
//                                                                      //
// Define a Point in n dimensions WITH periodic boundary conditions     //
//                                                                      //
// Burkhard Militzer                                    Paris 4-23-99   //
// Modified class structure                         Livermore 8-08-02   //
// Allow non-orthogonal bases                       Washington 4-10-06  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef _PERIODICPOINTGEN_
#define _PERIODICPOINTGEN_

#include "BasicPoint.h"
#include "Array.h" // only needed for SetBoxSize(Array1 p);
#include "Form.h" // only needed for PrintLatticeVector

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
class PeriodicPointGen:public BasicPoint <DIM> {
 public:
  static PeriodicPointGen <DIM> a,b,c;             // lattice vectors given in carthesian coordinates
  static PeriodicPointGen <DIM> aMin,bMin,cMin;    // (possibly different) lattice vectors point to the nearest images
  static PeriodicPointGen <DIM> aa,bb,cc;          // reciprocal lattice vectors 
  static PeriodicPointGen <DIM> aaMin,bbMin,ccMin; // (possibly different) lattice vectors point to the nearest images
  static double volume;                           // = Lx*Ly*Lz in cubic case only
  static bool isentropic;                         // |a|=|b|=|c|
  static bool orthogonal;                         // all angle are 90 degrees 

  PeriodicPointGen() {};
  PeriodicPointGen(const double d) : BasicPoint <DIM> (d) {}
  PeriodicPointGen(const BasicPoint <DIM> & p) : BasicPoint <DIM> (p) {}
  PeriodicPointGen(const PeriodicPointGen & p) : BasicPoint <DIM> (p) {}
  PeriodicPointGen(const ShortVector<double,DIM> & p) : BasicPoint <DIM> (p) {}
  PeriodicPointGen(const double x, const double y, const double z) : BasicPoint <DIM> (x,y,z) {}

  // from carthesian to reduced
  static PeriodicPointGen Reduced(const PeriodicPointGen & r) {
    PeriodicPointGen rr;
    rr[0] = r*aa;
    rr[1] = r*bb;
    rr[2] = r*cc;
    rr /= 2.0*pi;
    return rr;
  }
  PeriodicPointGen Reduced() const {
    return Reduced(*this);
  }

  // from reduced to cartesian
  static PeriodicPointGen Cartesian(const PeriodicPointGen & rr) {
    return rr[0]*a+rr[1]*b+rr[2]*c;
  }
  PeriodicPointGen Cartesian() const {
    return Cartesian(*this);
  }

  static void MapReducedInside(PeriodicPointGen & dd) {
    for(int i=0; i<DIM; i++) {
      dd[i] -= floor(dd[i]+0.5);
    }
  }
  void MapReducedInside() {
    MapReducedInside(*this);
  }
  static void MapInside(PeriodicPointGen & d) {
    PeriodicPointGen dd = Reduced(d);
    MapReducedInside(dd);
    d = Cartesian(dd);
  }
  static void MapReducedInsideZeroToOne(PeriodicPointGen & dd) {
    for(int i=0; i<DIM; i++) {
      dd[i] -= floor(dd[i]);
    }
  }
  void MapReducedInsideZeroToOne() {
    MapReducedInsideZeroToOne(*this);
  }

  static void MapInsideZeroToOne(PeriodicPointGen & d) {
     PeriodicPointGen dd = Reduced(d);
     MapReducedInsideZeroToOne(dd);
     d = Cartesian(dd);
  }
  void MapInsideZeroToOne() {
    MapInsideZeroToOne(*this);
  }

  static void MapInsideVector(PeriodicPointGen & d, const PeriodicPointGen & v) {
    double q = d*v/v.Norm2();
    //    if (floor(q+0.5)!=0.0) Write4(d,d.Norm(),v.Reduced(),q);
    //    Write4(d,d.Norm(),v.Reduced(),q);
    d -= v*floor(q+0.5);
    //    if (floor(q+0.5)!=0.0) Write4(d,d.Norm(),v.Reduced(),q);
    //    Write4(d,d.Norm(),v.Reduced(),q);
  }
  static void MapInsideMinimumImageTest(PeriodicPointGen & d) {
    //  static void MapInside(PeriodicPointGen & d) {
    MapInsideVector(d,a);
    MapInsideVector(d,b);
    MapInsideVector(d,c);
    MapInsideVector(d,a+b);
    MapInsideVector(d,a-b);
    MapInsideVector(d,b+c);
    MapInsideVector(d,b-c);
    MapInsideVector(d,c+a);
    MapInsideVector(d,c-a);
    MapInsideVector(d,a+b+c);
    MapInsideVector(d,a+b-c);
    MapInsideVector(d,b+c-a);
    MapInsideVector(d,c+a-b);
    /*
    PeriodicPointGen dd = Reduced(d);
    MapReducedInside(dd);
    d = Cartesian(dd);
    */
  }
  static void MapInsideMinimumImageSafely(PeriodicPointGen & d, const int n=1) {
    PeriodicPointGen dd = Reduced(d);
    MapReducedInside(dd);
    d = Cartesian(dd);

    PeriodicPointGen dMin = d;
    double dMinn = dMin.Norm();
    for(int i=-n; i<=+n; i++) {
      for(int j=-n; j<=+n; j++) {
	for(int k=-n; k<=+n; k++) {
	  PeriodicPointGen dijk = d + a*double(i) + b*double(j) + c*double(k);
	  double dijkn=dijk.Norm();
	  if (dijkn<dMinn) {
	    dMinn = dijkn;
	    dMin  = dijk;
	    //	    Write3(i,j,k);
	  }
	}
      }
    }
    d = dMin;
  }

  static double MinimumImageDistanceSafely(const int n=1) {
    double dMin = a.Norm();
    for(int i=-n; i<=+n; i++) {
      for(int j=-n; j<=+n; j++) {
	for(int k=-n; k<=+n; k++) {
	  if (abs(i)+abs(j)+abs(k)>0) {
	    PeriodicPointGen dijk = a*double(i) + b*double(j) + c*double(k);
	    double dijkn=dijk.Norm();
	    if (dijkn<dMin) {
	      dMin = dijkn;
	      //	      Write4(i,j,k,dMin*0.529177);
	    }
	  }
	}
      }
    }
    //    Write(dMin*0.529177);
    return dMin;
  }
  
  // If speed is important, do not call MapInsideMinimumImage() but MapInside() instead
  static void MapInsideMinimumImage(PeriodicPointGen & d, const int n=1) {
    //    error("MapInsideMinimumImage: Is this slow function really needed?");
    //    MapInsideSafely(d);
    //    return;

    PeriodicPointGen d2=d;
    MapInsideMinimumImageSafely(d,n);
    //    MapInsideMinimumImageTest(d2); // take out for speed
    PeriodicPointGen dd = d2-d;
    PeriodicPointGen rr = dd.Reduced();
    for(int iDim=0; iDim<DIM; iDim++) {
      rr[iDim] = rr[iDim]-floor(rr[iDim]+0.5); // set to difference from nearest integer
    }
    if (rr.Norm()>1e-8) error("PBC",d.Reduced(),d2.Reduced(),rr);
  }
  void MapInside() {
    MapInside(*this);
  }
  static void MapInside(double & x, const double cell) { // old function
    x -= cell*floor(x/cell+0.5);
  }
  static void MapInside(double & x, const int d) { // old function
    CheckOrthogonal();
    if (d==0) MapInside(x,a[0]);
    if (d==1) MapInside(x,b[1]);
    if (d==2) MapInside(x,c[2]);
  }
  static bool Inside(const PeriodicPointGen & r, const double eps=1e-08) {
    PeriodicPointGen rr = Reduced(r);
    for(int d=0; d<DIM; d++) {
      if (fabs(rr[d])>0.5+eps) return false;
    }
    return true;
  }
  static bool InsideZeroToOne(const PeriodicPointGen & r, const double eps=1e-08) { // 0 and 1 are ok
    PeriodicPointGen rr = Reduced(r);
    for(int d=0; d<DIM; d++) {
      if (rr[d]<-eps || rr[d]>1.0+eps) return false;
    }
    return true;
  }
  static bool TruelyInsideZeroToOne(const PeriodicPointGen & r, const double eps=1e-08) {// 0 and 1 are not ok
    PeriodicPointGen rr = Reduced(r);
    for(int d=0; d<DIM; d++) {
      if (rr[d]<+eps || rr[d]>1.0-eps) return false;
    }
    return true;
  }
  static bool OneTruelyInsideZeroToOne(const PeriodicPointGen & r, const double eps=1e-08) {// 0 and 1 are not ok
    PeriodicPointGen rr = Reduced(r);
    int nBoundary = 0;
    for(int d=0; d<DIM; d++) {
      if (rr[d]<-eps || rr[d]>1.0+eps) return false;
      if (rr[d]<+eps || rr[d]>1.0-eps) nBoundary++;
    }
    return (nBoundary<DIM);
  }

  static PeriodicPointGen DifferenceMinimumImage(const PeriodicPointGen & r1, const PeriodicPointGen & r2, const int n=1) {
    PeriodicPointGen d = r1-r2;
    MapInsideMinimumImage(d,n);
    return d;
  }
  static double DistanceMinimumImage(const PeriodicPointGen & r1, const PeriodicPointGen & r2, const int n=1) {
    PeriodicPointGen d = DifferenceMinimumImage(r1,r2,n);
    return d.Norm();
  }
  static PeriodicPointGen DifferencePBCReduced(const PeriodicPointGen & rr1, const PeriodicPointGen & rr2) {
    PeriodicPointGen dd = rr1-rr2;
    MapReducedInside(dd);
    return dd;
  }
  // Not always correct for noncubic boxes - use DifferenceMinimumImages
  static PeriodicPointGen DifferencePBC(const PeriodicPointGen & r1, const PeriodicPointGen & r2) {
    PeriodicPointGen rr1 = Reduced(r1);
    PeriodicPointGen rr2 = Reduced(r2);
    PeriodicPointGen dd = DifferencePBCReduced(rr1,rr2);
    PeriodicPointGen d = Cartesian(dd);
    return d;
  }
  static PeriodicPointGen DifferencePBCNormalized(const PeriodicPointGen & r1, const PeriodicPointGen & r2) {
    PeriodicPointGen d = DifferencePBC(r1,r2);
    double dn = d.Norm();
    d /= dn;
    return d;
  }
  void Separation(const PeriodicPointGen & r1, const PeriodicPointGen & r2) {
    *this = DifferencePBC(r1,r2);
  }
  static double DistancePBC2(const PeriodicPointGen & r1, const PeriodicPointGen & r2) {
    return DifferencePBC(r1,r2).Norm2();
  }
  static double DistancePBC2Reduced(const PeriodicPointGen & r1, const PeriodicPointGen & r2) {
    return DifferencePBCReduced(r1,r2).Norm2();
  }
  static double DistancePBC(const PeriodicPointGen & r1, const PeriodicPointGen & r2) {
    return sqrt(DistancePBC2(r1,r2));
  }
  static double DistancePBCReduced(const PeriodicPointGen & r1, const PeriodicPointGen & r2) {
    return sqrt(DistancePBC2Reduced(r1,r2));
  }
  static double AnglePBCSafelyDeg(const PeriodicPointGen & r1, const PeriodicPointGen & rM, const PeriodicPointGen & r2) {
    PeriodicPointGen d1 = DifferencePBC(r1,rM);
    PeriodicPointGen d2 = DifferencePBC(r2,rM);
    return AngleSafelyDeg(d1,d2);
  }
  void SetMidPoint(const PeriodicPointGen & r1, const PeriodicPointGen & r2) {
    PeriodicPointGen rr1 = Reduced(r1);
    PeriodicPointGen rr2 = Reduced(r2);
    PeriodicPointGen dd  = rr2-rr1;
    MapReducedInside(dd);
    dd = rr1 + dd*0.5;
    MapReducedInside(dd);
    *this = Cartesian(dd);
  }
  static void CheckIsotropic() {
    if (!isentropic) error("Basis not isentropic as needed for this function call");
  }
  static void CheckOrthogonal() {
    //    PrintCellInfo();
    if (!orthogonal) error("Basis not orthogonal as needed for this function call");
  }
  static bool RightHanded() {
    return (a*CrossProduct(b,c)>0.0);
  }
  static bool LeftHanded() {
    return !RightHanded();
  }
  static void PrintLatticeVector(const PeriodicPointGen & v, const bool flag=true) {
    Form f=Form(fix,10,6);
    cout << "(" << f(v[0]) << "," << f(v[1]) << "," << f(v[2]) << ")";
    if (flag) cout << endl;
  }
  static void PrintCellInfoForTable() {
    double AUToA=0.529177208607;
    Form fl(fix,6,3);
    cout << fl(a.Norm()*AUToA) << endl;
    cout << fl(b.Norm()*AUToA) << endl;
    cout << fl(c.Norm()*AUToA) << endl;
    Form fa(fix,5,2);
    cout << fa(AngleSafelyDeg(b,c)) << endl;
    cout << fa(AngleSafelyDeg(c,a)) << endl;
    cout << fa(AngleSafelyDeg(a,b)) << endl;
    cout << "Z= " << fix(10,8,sin(AngleSafely(c,a))*c.Norm()*AUToA) << endl;
  }

  static void WritePoint(ostream & of, const Form & f, const PeriodicPointGen & p, const string s) {
    of << s;
    for(int d=0; d<DIM; d++) {
      of << f(p[d]);
      if (d<DIM-1) of << " ";
      else of << endl;
    }
  }
  static void PrintCellInfoVASP() {
    double AUToA=0.529177208607;
    Form f(fix,21,16);
    double s=1.0;
    if (Cubic()) s = GetIsoBoxSize()*AUToA;
    cout << f(s) << endl;
    WritePoint(cout,f,GetLatticeVector(0)*AUToA/s,"  ");
    WritePoint(cout,f,GetLatticeVector(1)*AUToA/s,"  ");
    WritePoint(cout,f,GetLatticeVector(2)*AUToA/s,"  ");
  }

  static void PrintCellInfo(const bool shortFlag=false, const bool semiShortFlag=false) {
    double AUToA=0.529177208607;
    if (shortFlag) {
      cout << "a[A]= " << a.Norm()*AUToA 
	   << " b[A]= " << b.Norm()*AUToA 
	   << " c[A]= " << c.Norm()*AUToA 
	   << " alpha= " << AngleSafelyDeg(b,c)
	   << " beta= "  << AngleSafelyDeg(c,a)
	   << " gamma= " << AngleSafelyDeg(a,b)
	   << " V[A^3]= " << GetVolume()*AUToA*AUToA*AUToA
	   << endl;
    }
    if (semiShortFlag) {
      cout << "a[A]= "; PrintLatticeVector(a*AUToA);
      cout << "b[A]= "; PrintLatticeVector(b*AUToA);
      cout << "c[A]= "; PrintLatticeVector(c*AUToA);
    }
    if (shortFlag || semiShortFlag) return;
    cout << "a   = "; PrintLatticeVector(a);
    cout << "b   = "; PrintLatticeVector(b);
    cout << "c   = "; PrintLatticeVector(c);
    cout << "a[A]= "; PrintLatticeVector(a*AUToA);
    cout << "b[A]= "; PrintLatticeVector(b*AUToA);
    cout << "c[A]= "; PrintLatticeVector(c*AUToA);
    cout << "|a|    = " << a.Norm() << endl;
    cout << "|b|    = " << b.Norm() << endl;
    cout << "|c|    = " << c.Norm() << endl;
    cout << "|a| [A]= " << a.Norm()*AUToA << endl;
    cout << "|b| [A]= " << b.Norm()*AUToA << endl;
    cout << "|c| [A]= " << c.Norm()*AUToA << endl;
    cout << "alpha= " << AngleSafelyDeg(b,c) << endl;
    cout << "beta=  " << AngleSafelyDeg(c,a) << endl;
    cout << "gamma= " << AngleSafelyDeg(a,b) << endl;
    //    cout << "|a|*|b|*|c| = " << a.Norm()*b.Norm()*c.Norm() << endl;
    cout << "V     = " << GetVolume() << endl;
    cout << "V[A^3]= " << GetVolume()*AUToA*AUToA*AUToA << endl;
    cout << "Reciprocal a= "; PrintLatticeVector(aa);
    cout << "Reciprocal b= "; PrintLatticeVector(bb);
    cout << "Reciprocal c= "; PrintLatticeVector(cc);
    cout << "a,b,c form a " << ((RightHanded()) ? "right" : "left") << "handed system." << endl;
    //    Write6((a[1]==0.0),(a[2]==0.0),(b[0]==0.0),(b[2]==0.0),(c[0]==0.0),(c[1]==0.0));
    //    Pause("PrintCellInfo");
    cout << endl;
  }
  static void PrintCellInfoOrdered() {
    // print such that |a| >= |b| >= |c| 
    // print such that gamma < 90 and beta  < 90 
    PeriodicPointGen ap=a;
    PeriodicPointGen bp=b;
    PeriodicPointGen cp=c;
    if (ap.Norm()<bp.Norm()) { PeriodicPointGen p=ap; ap=bp; bp=p; }
    if (ap.Norm()<cp.Norm()) { PeriodicPointGen p=ap; ap=cp; cp=p; }
    if (bp.Norm()<cp.Norm()) { PeriodicPointGen p=bp; bp=cp; cp=p; }
    if (AngleSafelyDeg(ap,bp)<90.0) { bp = -bp; } // gamma > 90
    if (AngleSafelyDeg(ap,cp)<90.0) { cp = -cp; } // beta  > 90 
    cout << "a= " << ap << endl;
    cout << "b= " << bp << endl;
    cout << "c= " << cp << endl;
    cout << "|a|= " << ap.Norm() << endl;
    cout << "|b|= " << bp.Norm() << endl;
    cout << "|c|= " << cp.Norm() << endl;
    cout << "alpha= " << AngleSafelyDeg(bp,cp) << endl;
    cout << "beta=  " << AngleSafelyDeg(cp,ap) << endl;
    cout << "gamma= " << AngleSafelyDeg(ap,bp) << endl << endl;
    //    cout << "alpha-90= " << AngleSafelyDeg(bp,cp)-90.0 << endl;
    //    cout << "beta-90= "  << AngleSafelyDeg(cp,ap)-90.0 << endl;
    //    cout << "gamma-90= " << AngleSafelyDeg(ap,bp)-90.0 << endl << endl;
  }
  // do not switch axes because atoms might be registered already, those must be in reduced coord. however 
  static void RealignCell(const bool print=true) { 
    SetBoxSize(a.Norm(),b.Norm(),c.Norm(),
	       AngleSafelyDeg(b,c),
	       AngleSafelyDeg(c,a),
	       AngleSafelyDeg(a,b),print);
  }
  // do not switch axes because atoms might be registered already, those must be in reduced coord. however 
  static void RealignCellTriclinicRudy() { 
    SetBoxSizeTriclinicRudy(a.Norm(),b.Norm(),c.Norm(),
			    AngleSafelyDeg(b,c),
			    AngleSafelyDeg(c,a),
			    AngleSafelyDeg(a,b));
  }
  // do not switch axes because atoms might be registered already, those must be in reduced coord. however 
  static void RealignCellThirdOrientation() { 
    SetBoxSizeThirdOrientation(a.Norm(),b.Norm(),c.Norm(),
			    AngleSafelyDeg(b,c),
			    AngleSafelyDeg(c,a),
			    AngleSafelyDeg(a,b));
  } 
  // do not switch axes because atoms might be registered already, those must be in reduced coord. however 
  static void RealignCellForthOrientation() { 
    SetBoxSizeForthOrientation(a.Norm(),b.Norm(),c.Norm(),
			       AngleSafelyDeg(b,c),
			       AngleSafelyDeg(c,a),
			       AngleSafelyDeg(a,b));
  }
  static void SetBoxSize(const PeriodicPointGen & a_, const PeriodicPointGen & b_, const PeriodicPointGen & c_, const bool print=false) {
    //    warning("SLOW PBC testing still switched ON");
    if (DIM!=3) error("dimension != 3");
    a = a_;
    b = b_;
    c = c_;

    volume = a*CrossProduct(b,c);
    //    Write2(volume,a.Norm()*b.Norm()*c.Norm());
    aa = 2.0*pi*CrossProduct(b,c)/volume;
    bb = 2.0*pi*CrossProduct(c,a)/volume;
    cc = 2.0*pi*CrossProduct(a,b)/volume;
    volume = fabs(volume);
    double eps = 1e-8;
    // adopt this very strict def of orthogonality - no rotated ortho grids allowed
    //    orthogonal = (a[1]==0.0) && (a[2]==0.0) && (b[0]==0.0) && (b[2]==0.0) && (c[0]==0.0) && (c[1]==0.0);

    double scale = max(max(a.Norm(),b.Norm()),c.Norm());
    orthogonal = true;
    if (fabs(a[1])>scale*eps) orthogonal=false;
    if (fabs(a[2])>scale*eps) orthogonal=false;
    if (fabs(b[0])>scale*eps) orthogonal=false;
    if (fabs(b[2])>scale*eps) orthogonal=false;
    if (fabs(c[0])>scale*eps) orthogonal=false;
    if (fabs(c[1])>scale*eps) orthogonal=false;

    //    orthogonal = (abs(volume-a.Norm()*b.Norm()*c.Norm()) < eps*volume);
    //    Write(a);
    //    Write(b);
    //    Write(c);
    //    Write(orthogonal);

    isentropic = true;
    isentropic &= (fabs(a.Norm()-b.Norm()) < eps * a.Norm());
    isentropic &= (fabs(a.Norm()-c.Norm()) < eps * a.Norm());

    SetMinimumImageLatticeVectors(); // added 01/24/14

    if (print) {
      cout << endl << "Setting the box size:" << endl;
      PrintCellInfo(true);
      //      cout << endl << "Setting the box size [printed such that |a| >= |b| >= |c| and gamma,beta < 90 deg] :" << endl;
      //      PrintCellInfoOrdered();
    }
  }

  static void SetBoxSize(const Array1 <PeriodicPointGen> & p, const bool print=false) {
    SetBoxSize(p[0],p[1],p[2],print);
  }
  static void SetBoxSize(const PeriodicPointGen & p) {
    PeriodicPointGen a(p[0],0.0, 0.0 );
    PeriodicPointGen b(0.0, p[1],0.0 );
    PeriodicPointGen c(0.0, 0.0, p[2]);
    SetBoxSize(a,b,c);
  }
  static void SetBoxSize(const double a, const double b, const double c) {
    PeriodicPointGen p(a,b,c);
    SetBoxSize(p);
  }
  static void SetBoxSize(const double a_) {
    PeriodicPointGen a(a_);
    SetBoxSize(a);
  }
  static bool SetBoxSizeRadians(const double a, const double b, const double c, 
				const double alpha, const double beta, const double gamma, const bool print=true, const bool test=false) {
    double alphaDeg = alpha/pi*180.0;
    double betaDeg  = beta /pi*180.0;
    double gammaDeg = gamma/pi*180.0;
    return SetBoxSize(a,b,c,alphaDeg,betaDeg,gammaDeg,print,test);
  }
  // Dera's recommendation
  // X || a_vector  
  // Z || a cross b --> b must be in (X,Y) plane
  // Y || X cross Z  
  // X = a * (1,0,0) 
  // Y = b * (cos(gamma),sin(gamma),0)
  // Z = c * (...)
  // test = true, just check if parameters are valid, do not set new cell
  static bool SetBoxSize(const double a_, const double b_, const double c_, 
			 const double alphaDeg, const double betaDeg, const double gammaDeg, const bool print=true, const bool test=false) {
    double alpha = alphaDeg*pi/180.0;
    double beta  = betaDeg *pi/180.0;
    double gamma = gammaDeg*pi/180.0;

    if (cos(alpha)*cos(alpha) + cos(beta)*cos(beta) + cos(gamma)*cos(gamma) - 2.0*cos(alpha)*cos(beta)*cos(gamma) >= 1.0) {
      if (test) return false;
      else error("bad angles 1: e.g. alpha+beta>gamma",alphaDeg,betaDeg,gammaDeg);
    }

    PeriodicPointGen a(a_,0.0,0.0);

    PeriodicPointGen b(b_*cos(gamma),b_*sin(gamma),0.0);
    if (gammaDeg==90.0) { b[0]=0.0; b[1]=b_; }

    PeriodicPointGen c;
    c[0] = cos(beta);
    if (sin(gamma)==0.0) {
      if (test) return false;
      else error("Bad gamma angle",gamma);
    }

    c[1] = (cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);

    c[2] = 1.0-c[0]*c[0]-c[1]*c[1];
    if (c[2]<0.0)  {
      if (test) return false;
      else error("bad angles 2: e.g. alpha+beta>gamma",alphaDeg,betaDeg,gammaDeg);
    }
    c[2] = sqrt(max(0.0,c[2]));
		//    if (c[2]==0.0) Write10(alphaDeg,betaDeg,gammaDeg,a,b,c,cos(alpha),cos(beta),cos(gamma),sin(gamma));

    if (alphaDeg==90.0 && betaDeg==90.0 && gammaDeg==90.0) {
      c[0] = 0.0;
      c[1] = 0.0;
      c[2] = 1.0;
    }      

    c *= c_;
    if (!test) SetBoxSize(a,b,c,print);
    return true;
  }

  static void SetBoxSize(const PeriodicPointGen & norm, const PeriodicPointGen & anglesDeg, const bool print=true) { 
    SetBoxSize(norm[0],norm[1],norm[2],anglesDeg[0],anglesDeg[1],anglesDeg[2]);
  }

  // c--> a              b->c   a->b
  // gamma --> alpha     
  static void SetBoxSizeTriclinicRudyOld(const double a_, const double b_, const double c_, 
					 const double alphaDeg, const double betaDeg, const double gammaDeg) {
    double alpha = alphaDeg*pi/180.0;
    double beta  = betaDeg *pi/180.0;
    double gamma = gammaDeg*pi/180.0;
    PeriodicPointGen c(0.0,0.0,c_);
    PeriodicPointGen b(b_*sin(alpha),0.0,b_*cos(alpha));
    PeriodicPointGen a;
    a[2] = cos(beta);
    if (sin(alpha)==0.0) error("Bad alpha angle",alpha);
    a[0] = (cos(gamma)-cos(beta)*cos(alpha))/sin(alpha);
    a[1] = sqrt(max(0.0,1.0-a[2]*a[2]-a[0]*a[0]));
    a *= a_;
    Write3(a,b,c);
    //    warning("Unit cell MAY NOT BE CORRECT, check alpha and beta and gamma",alphaDeg,betaDeg);
    SetBoxSize(a,b,c,false);
    double phi = (bb[0]==0.0) ? ((fabs(bb[1]-pi/2.0)<1e-8) ? -pi/2.0 : +pi/2.0) : atan(bb[1]/bb[0]);
    PeriodicPointGen a2(a);
    PeriodicPointGen b2(b);
    phi = -phi;
    a2[0] =  cos(phi)*a[0]-sin(phi)*a[1];
    a2[1] =  sin(phi)*a[0]+cos(phi)*a[1];
    b2[0] =  cos(phi)*b[0]-sin(phi)*b[1];
    b2[1] =  sin(phi)*b[0]+cos(phi)*b[1];
    SetBoxSize(a2,b2,c,false);
    if (LeftHanded()) {
      a2 *= -1.0;
      b2 *= -1.0;
      c  *= -1.0;
    }
    SetBoxSize(a2,b2,c,true);
  }
  // Z || c_vector  
  // Y || a cross c --> a must be in (X,Z) plane
  // X || Y cross Z  
  // Z = c * (0,0,1) 
  // Y = b * (...)
  // X = a * (sin(beta),0,cos(beta)) 
  static void SetBoxSizeTriclinicRudy(const double a_, const double b_, const double c_, 
				      const double alphaDeg, const double betaDeg, const double gammaDeg) {
    double alpha = alphaDeg*pi/180.0;
    double beta  = betaDeg *pi/180.0;
    double gamma = gammaDeg*pi/180.0;
    PeriodicPointGen c(0.0,0.0,c_);
    PeriodicPointGen b;
    PeriodicPointGen a(a_*sin(beta),0.0,a_*cos(beta));;
    b[2] = cos(alpha); // check this line!!!
    if (sin(beta)==0.0) error("Bad alpha angle",alpha);
    b[0] = (cos(gamma)-cos(alpha)*cos(beta))/sin(beta);
    b[1] = sqrt(max(0.0,1.0-b[2]*b[2]-b[0]*b[0]));
    b *= b_;
    Write3(a,b,c);
    //    warning("Unit cell MAY NOT BE CORRECT, check alpha and beta and gamma",alphaDeg,betaDeg);
    SetBoxSize(a,b,c,true);
  }
  // Y || b_vector  
  // Z || a cross b --> a must be in (X,Y) plane
  // X || Y cross Z  
  // Y = b * (0,1,0) 
  // X = a * (sin(gamma),cos(gamma),0) 
  // Z = c * (...)
  static void SetBoxSizeThirdOrientation(const double a_, const double b_, const double c_, 
						  const double alphaDeg, const double betaDeg, const double gammaDeg) {
    double alpha = alphaDeg*pi/180.0;
    double beta  = betaDeg *pi/180.0;
    double gamma = gammaDeg*pi/180.0;
    PeriodicPointGen b(0.0,b_,0.0);
    PeriodicPointGen a(a_*sin(gamma),a_*cos(gamma),0.0);;
    PeriodicPointGen c;
    c[1] = cos(alpha);
    if (sin(gamma)==0.0) error("Bad gamma angle",gamma);
    c[0] = (cos(beta)-cos(alpha)*cos(gamma))/sin(gamma);
    c[2] = sqrt(max(0.0,1.0-c[0]*c[0]-c[1]*c[1]));
    c *= c_;
    //    Write3(a,b,c);
    //    Quit("W");
    //    warning("Unit cell MAY NOT BE CORRECT, check alpha and beta and gamma",alphaDeg,betaDeg);
    SetBoxSize(a,b,c,false);
  }
  // Dera's recommendation but exchange X and Y
  static void SetBoxSizeForthOrientation(const double a_, const double b_, const double c_, 
					 const double alphaDeg, const double betaDeg, const double gammaDeg) {
    double alpha = alphaDeg*pi/180.0;
    double beta  = betaDeg *pi/180.0;
    double gamma = gammaDeg*pi/180.0;
    PeriodicPointGen a(0.0,a_,0.0);
    PeriodicPointGen b(b_*sin(gamma),b_*cos(gamma),0.0);
    PeriodicPointGen c;
    c[1] = cos(beta);
    if (sin(gamma)==0.0) error("Bad gamma angle",gamma);
    c[0] = (cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);
    c[2] = sqrt(max(0.0,1.0-c[0]*c[0]-c[1]*c[1]));
    c *= c_;
    //    warning("Unit cell MAY NOT BE CORRECT, check alpha and beta and gamma",alphaDeg,betaDeg);
    SetBoxSize(a,b,c,true);
  }

  static double GetMaxBoxSize() {
    CheckOrthogonal();
    return max(a[0],max(b[1],c[2]));
  }
  static double GetMinBoxSize() {
    CheckOrthogonal();
    return min(a[0],min(b[1],c[2]));
    //    return min(a.Norm(),min(b.Norm(),c.Norm()));
  }
  static double GetIsoBoxSize() {
    CheckIsotropic();
    CheckOrthogonal();
    return a.Norm();
  }
  static PeriodicPointGen GetBoxSize(const bool check=true) {
    if (check) CheckOrthogonal();
    return PeriodicPointGen(a[0],b[1],c[2]);
  }
  static double GetBoxSize(const int iDim, const bool check=true) {
    if (check) CheckOrthogonal();
    if (iDim==0) return a[0];
    if (iDim==1) return b[1];
    if (iDim==2) return c[2];
    error("GetBoxSize()");
    return 0.0;
  }
  static PeriodicPointGen GetLatticeVector(const int iDim) {
    if (iDim==0) return a;
    if (iDim==1) return b;
    if (iDim==2) return c;
    return a;
  }
  static Array1 <PeriodicPointGen> GetLatticeVectors() {
    Array1 <PeriodicPointGen> p(DIM);
    for(int d=0; d<DIM; d++) {
      p[d] = GetLatticeVector(d);
    }
    return p;
  }
  static void SetLatticeVector(const PeriodicPointGen & p, const int iDim, const bool print=false) {
    if (iDim==0) SetBoxSize(p,b,c,print);
    if (iDim==1) SetBoxSize(a,p,c,print);
    if (iDim==2) SetBoxSize(a,b,p,print);
  }
  static void SetLatticeVectors(const Array1 <PeriodicPointGen> & p, const bool print=false) {
    SetBoxSize(p[0],p[1],p[2],print);
  }
  static void SetLatticeVectors(const PeriodicPointGen & a_, const PeriodicPointGen & b_, const PeriodicPointGen & c_, const bool print=false) {
    SetBoxSize(a_,b_,c_,print);
  }
  static PeriodicPointGen GetReciprocalLatticeVector(const int iDim) {
    if (iDim==0) return aa;
    if (iDim==1) return bb;
    if (iDim==2) return cc;
    return aa;
  }
  static double GetUnitCellAngle(const int i) {
    if (i==0) return AngleSafely(b,c);
    if (i==1) return AngleSafely(c,a);
    if (i==2) return AngleSafely(a,b);
    error("wrong dimenstion in GetUnitCellAngle()");
    return -1.0;
  }
  static double GetUnitCellAngleDegrees(const int i) {
    return GetUnitCellAngle(i)*180.0/pi;
  }
  static PeriodicPointGen GetUnitCellAngles() {
    PeriodicPointGen a;
    for(int d=0; d<DIM; d++) 
      a[d] = GetUnitCellAngle(d);
    return a;
  }
  static PeriodicPointGen GetUnitCellAnglesDegrees() {
    PeriodicPointGen a;
    for(int d=0; d<DIM; d++) 
      a[d] = GetUnitCellAngleDegrees(d);
    return a;
  }
  static PeriodicPointGen GetLatticeVectorNorms() {
    PeriodicPointGen a;
    for(int d=0; d<DIM; d++) 
      a[d] = GetLatticeVector(d).Norm();
    return a;
  }
  static double GetMaxDistance() { // be careful for non-orthorhmbic cells
    //    PeriodicPointGen pp(0.5);
    //    pp.ToCartesian();
    //    return pp.Norm();
    PeriodicPointGen pp1(0.5, 0.5, 0.5);
    PeriodicPointGen pp2(0.5,-0.5, 0.5);
    PeriodicPointGen pp3(0.5, 0.5,-0.5);
    PeriodicPointGen pp4(0.5,-0.5,-0.5);
    pp1.ToCartesian();
    pp2.ToCartesian();
    pp3.ToCartesian();
    pp4.ToCartesian();
    double n1 = pp1.Norm();
    double n2 = pp2.Norm();
    double n3 = pp3.Norm();
    double n4 = pp4.Norm();
    //    Write4(n1,n2,n3,n4)
    double m=max(n1,n2);
    m=max(m,n3);
    m=max(m,n4);
    return m;
  }
  static double GetMaxDistanceReduced() {
    PeriodicPointGen pp(0.5);
    return pp.Norm();
  }
  static double GetVolume() {
    return volume;
  }
  static double Reduced(const double x) {
    CheckIsotropic();
    CheckOrthogonal();
    return x/a.Norm();
  }
  static double Cartesian(const double x) {
    CheckIsotropic();
    CheckOrthogonal();
    return x*a.Norm();
  }
  void ToCartesian() {
    *this = Cartesian(*this);
  }
  void Reduce() {
    *this = Reduced(*this);
  }

  static bool Cubic() {
    return Isotropic() && Orthogonal();
  }

  /*  static bool IsoOrthogonal() {
    return Isotropic() && Orthogonal();
  }
  */
  static bool Isotropic() {
    return isentropic;
  }
  static bool Orthogonal() {
    return orthogonal;
  }
  // Note: here we do not allow cell that are orthogonal but not aligned with coordinate system
  static bool OrthogonalTolerance(const double eps=1e-8) {
    double d=GetMaxDistance()*eps;
    if (fabs(a[1])>d) return false;
    if (fabs(a[2])>d) return false;
    if (fabs(b[0])>d) return false;
    if (fabs(b[2])>d) return false;
    if (fabs(c[0])>d) return false;
    if (fabs(c[1])>d) return false;
    return true;
  }

  // still not quite sure why I need to added. Understanding this page would help:
  // http://clang.llvm.org/compatibility.html#dep_lookup_bases
  static PeriodicPointGen CrossProduct(const PeriodicPointGen & a, const PeriodicPointGen & b) {
    return BasicPoint<DIM>::CrossProduct(a,b);
  }
  static double AngleSafelyDeg(const PeriodicPointGen & x, const PeriodicPointGen & y) {
    return BasicPoint<DIM>::AngleSafelyDeg(x,y);
  }
  static double AngleSafely(const PeriodicPointGen & x, const PeriodicPointGen & y) {
    return BasicPoint<DIM>::AngleSafely(x,y);
  }

 private:
  // smallest first
  static void Order(PeriodicPointGen & a, PeriodicPointGen & b) {
    double an=a.Norm();
    double bn=b.Norm();
    if (an>bn) { Swap(a,b); }
  }
  // smallest first
  static void Order(PeriodicPointGen & a, PeriodicPointGen & b, PeriodicPointGen & c) {
    double an=a.Norm();
    double bn=b.Norm();
    double cn=c.Norm();
    if (an>bn) { Swap(a,b); Swap(an,bn); }
    if (an>cn) { Swap(a,c); Swap(an,cn); } // a is now the smallest
    if (bn>bn) { Swap(b,c); }              // c is now the largest
  }
  static int RemoveComponent(const PeriodicPointGen & a, PeriodicPointGen & b) {
    //    cout << "1:"; Write(Point::AngleSafelyDeg(a,b));
    double ii = b*a/a.Norm2();
    if (fabs(ii)<=0.5+1e-10) return 0;
    int i = int(round(b*a/a.Norm2()));
    b = b - double(i) * a;
    //    cout << "2:"; Write(Point::AngleSafelyDeg(a,b));
    return i;
  }
  // return true if routine above would remove some i * 'a' from 'b'
  static bool AnyComponentToBeRemoved(const PeriodicPointGen & a, const PeriodicPointGen & b) { 
    double ii = b*a/a.Norm2();
    return (fabs(ii) > 0.5+1e-10);
  }
  // assume aaa,bbb,ccc have been set to some non-minimum lattice vectors
  static void MinimumLatticeVectors(PeriodicPointGen & aaa, PeriodicPointGen & bbb, PeriodicPointGen & ccc) {
    int ib,ic,ibc;
    do {
      //      Write3(aa,bb,cc);
      Order(aaa,bbb,ccc); // aa is now the smallest vector
      //      Write3(aa,bb,cc);
      ib = RemoveComponent(aaa,bbb);
      //      Write4(aa,bb,cc,ib);
      ic = RemoveComponent(aaa,ccc);
      //      Write4(aa,bb,cc,ic);
      Order(bbb,ccc);
      ibc = RemoveComponent(bbb,ccc);
      //      Write4(aa,bb,cc,ibc);
    } while (ib*ib+ic*ic+ibc*ibc>0);
  }
 public:
  // assume aaa,bbb,ccc have been set to some non-minimum lattice vectors
  static bool LatticeVectorsMinimal(PeriodicPointGen & aaa, PeriodicPointGen & bbb, PeriodicPointGen & ccc) {
    Order(aaa,bbb,ccc); // aa is now the smallest vector
    if (AnyComponentToBeRemoved(aaa,bbb)) return false;
    if (AnyComponentToBeRemoved(bbb,ccc)) return false;
    if (AnyComponentToBeRemoved(aaa,ccc)) return false;
    return true;
  }
  static double MaximalLatticeVectorRatio(PeriodicPointGen & aaa, PeriodicPointGen & bbb, PeriodicPointGen & ccc) {
    Order(aaa,bbb,ccc); // aa is now the smallest vector
    return ccc.Norm() / aaa.Norm();
  }
 private:
  static void SetMinimumImageLatticeVectors() {
    aMin=a;
    bMin=b;
    cMin=c;
    MinimumLatticeVectors(aMin,bMin,cMin);
    double volumeMin = aMin*CrossProduct(bMin,cMin);
    aaMin = 2.0*pi*CrossProduct(bMin,cMin)/volumeMin; // needed to map r1-r2 inside (aMin,bMin,cMin) cell
    bbMin = 2.0*pi*CrossProduct(cMin,aMin)/volumeMin;
    ccMin = 2.0*pi*CrossProduct(aMin,bMin)/volumeMin;
  }
  // not sure if we need to keep this as a separate routine
  static double MinimumImageDistanceSub(const PeriodicPointGen & aaa, const PeriodicPointGen & bbb, const PeriodicPointGen & ccc) {
    const int n=1;
    double dMin2=aaa.Norm2();
    for(int ix=-n; ix<=+n; ix++) { // 1...n
      for(int iy=-n; iy<=+n; iy++) {
	for(int iz=-n; iz<=+n; iz++) {
	  PeriodicPointGen r = aaa*double(ix) + bbb*double(iy) + ccc*double(iz);
	  double rn2 = r.Norm2();
	  if ((rn2>1e-10) && (dMin2>rn2)) { // ix=iy=iz=0 excluded
	    dMin2=rn2; 
	    //	    if (print) Write5(ix,iy,iz,dMin,r);
	  }
	}
      }
    }
    double dMin = sqrt(dMin2);
    //    Write(Point::AngleSafelyDeg(aa,bb));
    //    Write(Point::AngleSafelyDeg(bb,cc));
    //    Write(Point::AngleSafelyDeg(cc,aa));
    //    Write(-aa-bb+cc);
    //    Write(dMin);
    return dMin;
  }
 public:
  static double MinimumImageDistance() {
    return MinimumImageDistanceSub(aMin,bMin,cMin);
  }
  // from carthesian to reduced
  static PeriodicPointGen MinimumImageReduced(const PeriodicPointGen & r) {
    PeriodicPointGen rr;
    rr[0] = r*aaMin;
    rr[1] = r*bbMin;
    rr[2] = r*ccMin;
    rr /= 2.0*pi;
    return rr;
  }
  void MinimumImageReduce() {
    *this = MinimumImageReduced(*this);
  }
  static PeriodicPointGen MinimumImageCartesian(const PeriodicPointGen & rr) {
    return rr[0]*aMin+rr[1]*bMin+rr[2]*cMin;
  }
  void MinimumImageToCartesian() {
    *this = MinimumImageCartesian(*this);
  }
  static PeriodicPointGen MinimumImageDifference(const PeriodicPointGen & dOrg) {
    PeriodicPointGen d = dOrg;
    d.MinimumImageReduce();
    d.MapReducedInside();
    d.MinimumImageToCartesian();

    const int n=1; // should be sufficient (after mapping inside)
    //    const int n=2;
    double dMin2=d.Norm2();
    PeriodicPointGen dMin = d;
    for(int ix=-n; ix<=+n; ix++) { // 1...n
      for(int iy=-n; iy<=+n; iy++) {
	for(int iz=-n; iz<=+n; iz++) {
	  PeriodicPointGen dd = aMin*double(ix) + bMin*double(iy) + cMin*double(iz) + d;
	  double ddn2 = dd.Norm2();
	  if (ddn2 < dMin2) { // ix=iy=iz=0 now included
	    dMin = dd;
	    dMin2= ddn2; 
	    //	    if (print) Write5(ix,iy,iz,dMin,r);
	  }
	}
      }
    }
    return dMin;
  }
  static PeriodicPointGen MinimumImageDifference(const PeriodicPointGen & r1, const PeriodicPointGen & r2) {
    return MinimumImageDifference(r1-r2);
  }
  static double MinimumImageDistance(const PeriodicPointGen & d) {
    PeriodicPointGen dd = MinimumImageDifference(d);
    return dd.Norm();
  }
  static double MinimumImageDistance(const PeriodicPointGen & r1, const PeriodicPointGen & r2) {
    return MinimumImageDistance(r1-r2);
  }
};

#endif  // _PERIODICPOINTGEN_

