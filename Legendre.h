/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Associated Legendre Polynomials                                         //
//                                                                         //
// adopted from Numerical Recipes                                          //
//                                                                         //
// Burkhard Militzer                               Oxford, OH, 01-03-18    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _LEGENDRE_
#define _LEGENDRE_

#include "Standard.h"

inline double AssociatedLegendrePolynomial(const int l, const int m, const double x) {
  //	float fact,pll,pmm,pmmp1,somx2;
  //	int i,ll;

  if (m<0) {
    double f = double(Factorial(l+m)) / double(Factorial(l-m)); // avoid int division
    if (m%2==0) return f*AssociatedLegendrePolynomial(l,-m,x);
    else  return -f*AssociatedLegendrePolynomial(l,-m,x);
  }

  //  if (m < 0 || m > l || fabs(x) > 1.0) 
  if ( m>l || fabs(x)>1.0 ) 
    error("AssociatedLegendrePolynomial: bad arguments",l,m,x);

  double pmm=1.0;
  if (m > 0) {
    const double somx2 = sqrt((1.0-x)*(1.0+x));
    double fact = 1.0;
    for (int i=0; i<m; i++) {
      pmm *= -fact*somx2;
      fact += 2.0;
      //      Write4(x,somx2,pmm,fact);
    }
  }
  //  Write2(m,pmm);

  if (l == m) {

    return pmm;

  } else {

    double pmmp1 = x*(2*m+1)*pmm;
    if (l == (m+1)) {
      return pmmp1;
    } else {
      //      Write3(m,pmm,pmmp1);
      double pll;
      for (int ll=m+2; ll<=l; ll++) {
	pll = ( x*(2*ll-1)*pmmp1-(ll+m-1)*pmm ) / (ll-m);
	pmm=pmmp1;
	pmmp1=pll;
      }
      return pll;
    }

  }
}

// computes all Plm[m] for m=0...l
inline void AssociatedLegendrePolynomials(const int l, const double x, Array1 <double> & Plm) {
  if ( fabs(x)>1.0 ) error("AssociatedLegendrePolynomials: bad argument x=",x);

  Plm[0] = 1.0;
  if (l==0) return; // otherwise trobule below in the line Plm[l-1]=...
  
  double pmm=1.0;
  const double somx2 = sqrt((1.0-x)*(1.0+x));
  double fact = 1.0;
  for (int i=0; i<l; i++) {
    pmm *= -fact*somx2;
    fact += 2.0;
    Plm[i+1] = pmm;
    //      Write4(x,somx2,pmm,fact);
  }
  //  Write(Plm);

  Plm[l-1] *= x*(2*l-1);

  for(int m=0; m<l-1; m++) { // maybe there is smart way to combine the 'm' and 'll' loops but I do not see it yet.
    double pmm_   = Plm[m];
    double pmmp1_ = x*(2*m+1)*pmm_;
    //    Write3(m,pmm_,pmmp1_ );

    double pll;
    for (int ll=m+2; ll<=l; ll++) {
      pll    = ( x*(2*ll-1)*pmmp1_ - (ll+m-1)*pmm_ ) / (ll-m);
      pmm_   = pmmp1_;
      pmmp1_ = pll;
    }
    Plm[m] = pll;
  }
}

#endif // _LEGENDRE_
