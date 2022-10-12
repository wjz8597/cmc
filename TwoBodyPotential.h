/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Store the simple pair-potential information                            //
//                                                                         //
//  B. Militzer                                    Washington 07-12-04     //
//  J.Vorberger						01/29/06	//
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _TWOBODYPOTENTIAL_
#define _TWOBODYPOTENTIAL_

#include "Array.h"
#include "Point.h"
#include "Matrix3.h"
#include "FindRoot.h"

class TwoBodyPotential {
 public:
  virtual ~TwoBodyPotential() {} // = 0;
  virtual double OneBodyTerm() const = 0;
  virtual double V(const BPoint & r) const = 0;
  virtual BPoint DV(const BPoint & r) const = 0;
  virtual Matrix3 DDV(const BPoint & r) const = 0;
  
  double V(const Point & r1, const Point & r2) const {
    return V(Point::DifferencePBC(r1,r2));
  }
  BPoint DV(const Point & r1, const Point & r2) const {
    return DV(Point::DifferencePBC(r1,r2));
  }
  Matrix3 DDV(const Point & r1, const Point & r2) const {
    return DDV(Point::DifferencePBC(r1,r2));
  }
};

class IsotropicPotential : public TwoBodyPotential {
 public: 
  virtual double Vm(const double r) const = 0;
  virtual double DVm(const double r) const = 0;
  virtual double DDVm(const double r) const = 0;

  double V(const BPoint & r) const {
    return Vm(r.Norm());
  }
  BPoint DV(const BPoint & r) const {
    double rn=r.Norm();
    return r/rn * DVm(rn);
  }
  Matrix3 DDV(const BPoint & r) const {
    Matrix3 m = 0.0;

    double rn=r.Norm();
    double dv =DVm(rn);
    double ddv=DDVm(rn);

    double term1 = dv/rn;
    double term2 = 1.0/rn/rn*( ddv - dv/rn );
    //    Write5(rn,dv,ddv,term1,term2);
    for(int i=0; i<BPoint::nDim; i++) {
      m(i,i) += term1;
      for(int j=0; j<BPoint::nDim; j++) {
	m(i,j) += r[i]*r[j]*term2;
      }
    }
    return m;
  }
  void TestDVm(const double rn, const double eps=1e-8) {
    double dVm = DVm(rn);
    double d = rn*eps;
    double Vm2 = Vm(rn+d);
    double Vm1 = Vm(rn-d);
    Write4(rn,d,dVm,(Vm2-Vm1)/(2.0*d));
  }

};
  
class LennardJonesPotential : public IsotropicPotential {
 public:
  double A12, C6;
  LennardJonesPotential(const double A12_, const double C6_):A12(A12_),C6(C6_) {}

  double Vm(const double r) const {
    if (r==0.0) error("LJ r=0");
    double r2 = r*r;
    double r6 = 1.0/(r2*r2*r2);
    return A12*r6*r6 - C6*r6;
  }
  double DVm(const double r) const {
    if (r==0.0) error("LJ r=0");
    double r2 = r*r;
    double r6 = 1.0/(r2*r2*r2);
    return -12.0*A12*r6*r6/r + 6.0*C6*r6/r;
  }
  double DDVm(const double r) const {
    if (r==0.0) error("LJ r=0");
    double r2 = r*r;
    double r6 = 1.0/(r2*r2*r2);
    return 12.0*13.0*A12*r6*r6/r/r - 6.0*7.0*C6*r6/r/r;
  }
};

// Not useful:
// V(r) = A(|r|-r0)^2
class HarmonicPotential1D : public IsotropicPotential {
 public:
  double A, r0;
  HarmonicPotential1D(const double A_, const double r0_):A(A_),r0(r0_) {}

  double Vm(const double r) const {
    return A*sqr(r-r0);
  }
  double DVm(const double r) const {
    return 2.0*A*(r-r0);
  }
  double DDVm(const double r) const {
    return 2.0*A;
  }
};

// V(r) = A*(r-r0)^2
class HarmonicPotential {
 public:
  double A;
  BPoint r0;
  HarmonicPotential(const double A_, const BPoint & r0_):A(A_),r0(r0_) {}

  double V(const BPoint & r) const {
    BPoint d=r-r0;
    return A*d.Norm2();
  }
  BPoint DV(const BPoint & r) const {
    return 2.0*A*(r-r0);
  }
  Matrix3 DDV(const BPoint & r) const {
    return Matrix3::Identity()*(2.0*A);
  }
};

class CoulombPotential : public IsotropicPotential {
 public:
  const double C;
  CoulombPotential(const double C_):C(C_) {}

  double Vm(const double r) const {
    if (r==0.0) error("Coulomb r=0 problem");
    return C/r;
  }
  double DVm(const double r) const {
    if (r==0.0) error("Coulomb r=0 problem");
    return -C/r/r;
  }
  double DDVm(const double r) const {
    if (r==0.0) error("Coulomb r=0 problem");
    return 2.0*C/r/r/r;
  }
};

/*
class EwaldPairPotential : public TwoBodyPotential {
 public:
  EwaldPairPotential(Ewald & ewald_, const double q1_, const double q2_):
    ewald(ewald_),q1(q1_),q2(q2_),qq(q1*q2){}
  Ewald & ewald;
  double q1,q2,qq;
  double V(const BPoint & r) const {
    return ewald.FullPotentialInterpolation(r)*qq;
  }
  BPoint DV(const BPoint & r) const {
    BPoint f;
    double v;
    ewald.ForceAndPotentialInterpolation(r,v,f);
    error("check sign of forces: DV needs dV/dr but forces are -dV/dr");
    return f*qq;
  }
  Matrix3 DDV(const BPoint & r) const {
    Matrix3 der;
    BPoint f;
    ewald.FullPotentialDer(r,f,der);
    return der*qq;
  }
};

// oxygen-oxygen potential in TIP4P
class OOPotential {
 public:
  EwaldPairPotential ewaldPP;
  LennardJonesPotential lJP;

  OOPotential(Ewald & ewald_, const double qO_, const double A12_, const double C6_)
    :ewaldPP(ewald_,qO_,qO_),lJP(A12_,C6_){}

  double V(const BPoint & r) const {
    return ewaldPP.V(r)+lJP.V(r);
  }
  BPoint DV(const BPoint & r) const {
    return ewaldPP.DV(r)+lJP.DV(r);
  }
  Matrix3 DDV(const BPoint & r) const {
    return ewaldPP.DDV(r)+lJP.DDV(r);
  }
};
*/

// helps avoiding if statements  
class ZeroPotential : public IsotropicPotential {
  double OneBodyTerm() const {
    return 0.0;
  }
  double Vm(const double r) const {
    return 0.0;
  }
  double DVm(const double r) const {
    return 0.0;
  }
  double DDVm(const double r) const {
    return 0.0;
  }
};

class Exp6PotentialFindInflectionPoint {
 public:
  Exp6PotentialFindInflectionPoint(const double eps_, const double rStar_, const double alpha_)
    :eps(eps_),rStar(rStar_),alpha(alpha_){}
  const double eps, rStar, alpha;
  double operator()(const double r) const {
    double rr = r/rStar;
    double V = 6.0*alpha/rStar*alpha/rStar*exp(alpha*(1.0-rr)) - 42.0*alpha/(rr*rr*rr*rr*rr*rr*r*r);
    V *= eps/(alpha-6.0);
    return V;
  }
};

inline double FindInflectionPoint(const double eps, const double rStar, const double delta) {
  Exp6PotentialFindInflectionPoint fip(eps,rStar,delta);
  return FindRootStepper(rStar*0.1,1e-8,1e-8,fip);
}

class Exp6Potential : public IsotropicPotential {
 public:
  Exp6Potential(const double eps_, const double rStar_, const double alpha_)
    :eps(eps_),rStar(rStar_),alpha(alpha_),W(0.0) {
    //    cout << "Exp-6 potential initialized: eps= " << eps*PC::AUToK << " K  r*= " << rStar*PC::AUToA << " A  alpha= " << alpha << endl;
    double w   = FindInflectionPoint(eps,rStar,alpha);
    double Vw  = Vm(w);
    double DVw = DVm(w);
    B = -DVw/Vw;
    A = Vw/exp(-B*w);
    if (A<0.0) error("Found the wrong inflection point",W,A,B);
    W = w;
    //    Write3(W,A,B);
    //    Write3(W*PC::AUToA,A*PC::AUToK,B/PC::AUToA);
  }
  const double eps, rStar, alpha;
  double A,B,W;

  double OneBodyTerm() const {
    return 0.0;
  }
  double Vm(const double r) const {
    if (r<W) return A*exp(-B*r);
    double rr = r/rStar;
    double V = 6.0*exp(alpha*(1.0-rr)) - alpha/(rr*rr*rr*rr*rr*rr);
    V *= eps/(alpha-6.0); // minus sign added
    return V;
  }
  double DVm(const double r) const {
    if (r<W) return -B*A*exp(-B*r);
    double rr = r/rStar;
    double V = -6.0*alpha/rStar*exp(alpha*(1.0-rr)) + 6.0*alpha/(rr*rr*rr*rr*rr*rr*r);
    V *= eps/(alpha-6.0);
    return V;
  }
  double DDVm(const double r) const {
    if (r<W) return B*B*A*exp(-B*r);
    double rr = r/rStar;
    double V = 6.0*alpha/rStar*alpha/rStar*exp(alpha*(1.0-rr)) - 42.0*alpha/(rr*rr*rr*rr*rr*rr*r*r);
    V *= eps/(alpha-6.0);
    return V;
  }
};

////////////////////////////////////////////////////////////////////////////////////

// Jan's addition: potential acc. Slattery & Hubbard, Icarus 29, 187 (1976)
class SlatteryHubbardPotential : public IsotropicPotential {
 public:
  SlatteryHubbardPotential():C1(8.2),C2(1.74),C3(13.0),C4(116.0),C5(2.7144),C56(C5*C5*C5*C5*C5*C5) {}
  
    const double C1,C2,C3,C4,C5,C56;
  
  double OneBodyTerm() const {
    return 0.0;
  }
  double Vm(const double r) const {
    if (r==0.0) error("Coulomb r=0 problem");
    double r6 = r*r*r*r*r*r;
    double r8 = r6*r*r;
    double V = C1*exp(-C2*r) - (C3/(r6)+C4/(r8))*exp(-C56/r6);
    V *= 2.0;
    return V;
  }
  double DVm(const double r) const {
    if (r==0.0) error("Coulomb r=0 problem");
    double r6 = r*r*r*r*r*r;
    double r7 = r6*r;
    double r8 = r7*r;
    double r9 = r8*r;
    double V = -C1*C2*exp(-C2*r)+(6.0*C3/r7+8.0*C4/r9)*exp(-C56/r6)-6.0*C56/r7*((C3/r6+C4/r8)*exp(-C56/r6));
    V *= 2.0;
    return V;
  }
  double DDVm(const double r) const {
    if (r==0.0) error("Coulomb r=0 problem");
    double r6 = r*r*r*r*r*r;
    double r7 = r6*r;
    double r8 = r7*r;
    double r9 = r8*r;
    double r10 = r9*r;
    double V= -C1*C2*C2*exp(-C2*r)-(42.0*C3/r8+72.0*C4/r10)*exp(-C56/r6)
    +12.0*C56/r7*((6.0*C3/r7+8.0*C4/r9)*exp(-C56/r6))
    +42.0*C56/r8*((C3/r6+C4/r8)*exp(-C56/r6))
    -36.0*C56*C56/(r7*r7)*(C3/r6+C4/r8)*exp(-C56/r6);
    V*= 2.0;
    return V;
  }
};


// Jan's addition: potential acc. Silvera & Goldman, from Ross et al., J. Chem. Phys. 79, 1487 (1983). Eq. (1)
class SilveraGoldmanPotential : public IsotropicPotential {
 public:
  SilveraGoldmanPotential():alpha(1.73),beta(1.5671),
  gamma(0.00993),C6(12.14),C8(215.2),C9(143.1),C10(4813.9),rm(6.444),lim(1.28*rm) {}
  
  const double alpha, beta, gamma;
  const double C6,C8,C9,C10, rm, lim;
  
  double OneBodyTerm() const {
    return 0.0;
  }
  double Vm(const double r) const {
    if (r==0.0) error("Coulomb r=0 problem");
     
    double r6 = r*r*r*r*r*r;
    double r8 = r6*r*r;
    double r9 = r8*r;
    double r10 = r9*r;
    double f;
    if (r<=lim) {
      f = exp(-(1.28*rm/r-1)*(1.28*rm/r-1)); 
    } else {
      f = 1.0; 
    }
    double vp = exp(alpha-beta*r-gamma*r*r)-(C6/r6+C8/r8+C10/r10)*f;
    double V = vp+C9/r9*f;
    V *= 2.0;
    return V;
  }
  double DVm(const double r) const {
    if (r==0.0) error("Coulomb r=0 problem");

    double r6 = r*r*r*r*r*r;
    double r7 = r6*r;
    double r8 = r7*r;
    double r9 = r8*r;
    double r10 = r9*r;
    double r11 = r10*r;
    double f,fs;
    if (r<=lim) {
      f = exp(-(lim/r-1)*(lim/r-1)); 
      fs = 2.56*rm/(r*r)*(lim/r-1)*exp(-(lim/r-1)*(lim/r-1));
    } else {
      f = 1.0; 
      fs= 0.0;
    }
    double vps = (-beta-2.0*gamma*r)*exp(alpha-beta*r-gamma*r*r)-(-6.0*C6/r7-8.0*C8/r9-10.0*C10/r11)*f-(C6/r6+C8/r8+C10/r10)*fs;
    double V = vps-9.0*C9/r10*f+C9/r9*fs;
    V *= 2.0;
    return V;
  }

  double DDVm(const double r) const {
    if (r==0.0) error("Coulomb r=0 problem");

    double r4 = r*r*r*r;
    double r6 = r4*r*r;
    double r7 = r6*r;
    double r8 = r7*r;
    double r9 = r8*r;
    double r10 = r9*r;
    double r11 = r10*r;
    double r12 = r11*r;
    double f,fs,fss;
    if (r <= 1.28*rm) {
      f = exp(-(1.28*rm/r-1)*(1.28*rm/r-1));
      fs = 2.56*rm/(r*r)*(1.28*rm/r-1)*f;
      fss = -3.2768*rm*rm/r4*f-2.0/r*fs+fs*fs/f;
    } else {
      f = 1.0;
      fs = 0.0;
      fss = 0.0;
    }

    double vpss = -2.0*gamma*exp(alpha-beta*r-gamma*r*r)+(-beta-2.0*gamma*r)*(-beta-2.0*gamma*r)*exp(alpha-beta*r-gamma*r*r)
      -2.0*(-6.0*C6/r7-8.0*C8/r9-10.0*C10/r11)*fs
      -(42.0*C6/r8+72.0*C8/r10+110.0*C10/r12)*f
      -(C6/r6+C8/r8+C10/r10)*fss;
    double V = vpss+90.0*C9/r11*f-9.0*C9/r10*fs-9.0*C9/r10*fs+C9/r9*fss;
    V*= 2.0;
    return V;
  }
};

// Jan's addition: improved S&G potential acc. Ross et al., J. Chem. Phys. 79, 1487 (1983), Eq. (5)
class YangRossPotential : public SilveraGoldmanPotential {
 public:
  YangRossPotential():A(0.0018296),B(2.523856),C(0.63134),D(0.14154),E(0.01946),rc(4.8188),r1(2.26767) {}

  const double A,B,C,D,E,rc,r1;

  double OneBodyTerm() const {
    return 0.0;
  }
  double Vm(const double r) const {
    if (r>=rc) { 
      return SilveraGoldmanPotential::Vm(r);
    } else {
      double rc1 = r-rc;
      double rc2 = rc1*rc1;
      double rc3 = rc2*rc1;
      double expx = exp(-B*rc1 -C*rc2 -D*rc3 -E*rc3*(r-r1));
      double V = A*expx;
      return V;
    }
  }
  double DVm(const double r) const {
    if (r>=rc) {
      return SilveraGoldmanPotential::DVm(r);
    } else {
      double rc1 = r-rc;
      double rc2 = rc1*rc1;
      double rc3 = rc2*rc1;
      double expx = exp(-B*rc1 -C*rc2 -D*rc3 -E*rc3*(r-r1));
      double DV = A*(-B -2.0*C*rc1 -3.0*D*rc2 -3.0*E*rc2*(r-r1) -E*rc3) * expx;
      return DV;
    }
  }
  double DDVm(const double r) const {
    if (r>=rc) {
      return SilveraGoldmanPotential::DDVm(r);
    } else {
      double rc1 = r-rc;
      double rc2 = rc1*rc1;
      double rc3 = rc2*rc1;
      double expx = exp(-B*rc1 -C*rc2 -D*rc3 -E*rc3*(r-r1));
      double DDV = A*(-2.0*C -6.0*D*rc1 -6.0*E*rc1*(r-r1) -6.0*E*rc2)      * expx
	        + A*sqr(-B -2.0*C*rc1 -3.0*D*rc2 -3.0*E*rc2*(r-r1) -E*rc3) * expx;
     return DDV;
    }
  }
};

class AzizPotential : public IsotropicPotential {
 public:

  double OneBodyTerm() const {
    return 0.0;
  }

  // r in Angstroem
  // return V in K
  double VmInAngstroemAndKelvin(const double r) const {
    const double eps   = 10.956;
    const double rm    = 2.9683;
    const double a     = 1.86924404e5;
    const double alpha = 10.5717543;
    const double beta  = -2.07758779;
    const double d   = 1.438;
    const double c6  = 1.35186623;
    const double c8  = 0.41495143;
    const double c10 = 0.17151143;
    // R. A. Aziz, A. R. Janzen, M. R. Moldover, Phys. Rev. Letts. 74, 1586 (1995)
 
    const double x  =r/rm;
    const double a2i=1.0/(x*x);
    const double a6i=a2i*a2i*a2i;
    double f  =1.0;
    if (x<d) f=exp(-sqr(d/x-1.0));
    const double pot=eps*(a*exp(-alpha*x+beta*x*x)-f*a6i*(c6+a2i*(c8+a2i*c10)));
    //  Write3("Aziz",r,pot);
    return pot;
  }    
  
  double Vm(const double r) const {
    double v = PC::KToAU * VmInAngstroemAndKelvin(r*PC::ab/PC::A);
    //  Write2(r,v);
    return v;
  }

  double DVm(const double r) const {
    double dr = 1e-8;
    double p1 = Vm(r-dr);
    double p2 = Vm(r+dr);
    return (p2-p1)/(2.0*dr);
  }

  double DDVm(const double r) const {
    error("2nd derivative of Aziz pot not implemented.");
    return 0.0;
  }
};


///////////////////////// Miguel Morales'  Yukawa-type potential //////////////////////////

class MoralesYukawaTypePotential : public IsotropicPotential {
 public:
  MoralesYukawaTypePotential(const double a_, const double b_, const double L_)
    :a(a_),b(b_),L(L_) { 
  }
  const double a,b,L;

  double OneBodyTerm() const {
    return 0.0;
  }
  double Vm(const double r) const {
    if (r>L/2.0) return 0.0;
    double expBr  = exp(-b*r);
    double expBLr = exp(b*(r-L));
    double expBL  = exp(-b*L/2);
    return a*(expBr/r + expBLr/(L-r) - 4.0*expBL/L);
  }
  double DVm(const double r) const {
    if (r>L/2.0) return 0.0;
    double expBr  = exp(-b*r);
    double expBLr = exp(b*(r-L));
    double dVdr = -(b+1.0/r)/r * expBr + (b+1.0/(L-r))/(L-r) * expBLr;
    return a*dVdr;
  }
  double DDVm(const double r) const {
    error("2nd derivative was not implemented");
    return 0.0;
  }
};



#endif // _TWOBODYPOTENTIAL_
