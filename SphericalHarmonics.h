/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Spherical Harmonics Routine                                             //
//                                                                         //
// Jizhou Wu                                                               //
// Felipe Gonzalez                                                         //
// Burkhard Militzer                               Oxford, OH, 01-03-18    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _BONDORDER_
#define _BONDORDER_

#include <complex>
#include "Standard.h"
#include "Legendre.h"

typedef std::complex<double> Complex;

/*
double Legendre(int l,int m,double x){

 if(l<abs(m)){std::cout  << "Error"<< std::endl;
 return 0;
 }

 if(m==0){
  if(l==0)return 1.0;
  else if (l==1) return x;
  return ((2*l-1)*x*Legendre(l-1,m,x)-(l+m-1)*Legendre(l-2,m,x))/(l-m);
 }

 else if(m<0){
  // std::cout  << "m<0" << std::endl;
   return cos(M_PI*m)*Factorial(l+m)/Factorial(l-m)*Legendre(l,-m,x);
 }

 else {
  return ((l-m+1)*x*Legendre(l,m-1,x)-(l+m-1)*Legendre(l-1,m-1,x))/sqrt(1-x*x);
}

}

inline Complex SphericalHarmonicsOrg(const int l, const int m, const double theta, const double phi) {
  const Complex I(0.0,1.0);
  return sqrt((2*l+1)*Factorial(l-m)/(4*M_PI*Factorial(l+m))) * Legendre(l,m,cos(theta)) * exp(m*phi*I) ;
}
*/

inline Complex SphericalHarmonics(const int l, const int m, const double theta, const double phi){
  double factor = sqrt( (2*l+1) * Factorial(l-m) / (4.0*pi*Factorial(l+m)) );
  double legendre = AssociatedLegendrePolynomial(l,m,cos(theta));
  const Complex I(0.0,1.0);
  return factor * legendre * exp(m*phi*I) ;
}

inline Complex SphericalHarmonicsCosTheta(const int l, const int m, const double cosTheta, const double phi) {
  /*
  Array1 <double> Plm(l+1);
  AssociatedLegendrePolynomials(l,cosTheta,Plm);
  for(int m=0; m<=l; m++) {
    double legendre = AssociatedLegendrePolynomial(l,m,cosTheta);
    Write4(l,m,legendre,Plm[m]);
  }
  */

  double factor = sqrt( (2*l+1) * Factorial(l-m) / (4.0*pi*Factorial(l+m)) );
  double legendre = AssociatedLegendrePolynomial(l,m,cosTheta);
  const Complex I(0.0,1.0);
  //  if (l==-m) Write5(l,m,factor,legendre,exp(m*phi*I));
  return factor * legendre * exp(m*phi*I) ;
}

// Array Ylm has size 2*l+1
// store m=-l result in Ylm[0]
// store m= 0 result in Ylm[l]
// store m=+l result in Ylm[2*l]
//
// For l=6, Ylm has the size 13
// store m=-6 result in Ylm[0]
// store m= 0 result in Ylm[6]
// store m=+6 result in Ylm[12]
inline void SphericalHarmonicsCosTheta(const int l, const double cosTheta, const double phi, Array1 <Complex> & Ylm) {
  Array1 <double> Plm(l+1);
  AssociatedLegendrePolynomials(l,cosTheta,Plm);

  double factorM0 = sqrt( (2*l+1) / (4.0*pi) );
  Ylm[l] = Complex(factorM0*Plm[0],0.0); // case m=0

  for(int m=1; m<=l; m++) {
    double factor = factorM0 * sqrt( Factorial(l-m) / Factorial(l+m) ); // remember Factorial() now returns a double
    const Complex I(0.0,1.0);
    Ylm[l+m] = factor  * Plm[m] * exp( m*phi*I);
    if (m%2!=0) factor *= -1.0; // flip sign only for negative and odd 'm'
    Ylm[l-m] = factor  * Plm[m] * exp(-m*phi*I);
  }
}

inline void AddSphericalHarmonicsCosTheta(const int l, const double cosTheta, const double phi, 
					  Array1 <Complex> & Ylm, const double f=1.0) {
  Array1 <double> Plm(l+1);
  AssociatedLegendrePolynomials(l,cosTheta,Plm);

  double factorM0 = sqrt( (2*l+1) / (4.0*pi) );
  Ylm[l] += f*Complex(factorM0*Plm[0],0.0); // case m=0

  for(int m=1; m<=l; m++) {
    double factor = f*factorM0 * sqrt( Factorial(l-m) / Factorial(l+m) ); // remember Factorial() now returns a double
    const Complex I(0.0,1.0);
    Ylm[l+m] += factor  * Plm[m] * exp( m*phi*I);
    if (m%2!=0) factor *= -1.0; // flip sign only for negative and odd 'm'
    Ylm[l-m] += factor  * Plm[m] * exp(-m*phi*I);
  }
}

inline void AddSphericalHarmonicsCosTheta(const int l, const double cosTheta, const double phi, 
					  Array2 <Complex> & Ylm, const int i, const int j, const double f=1.0) {//ylm[#particle][#m] as the element of Array1<Array2<Complex>>Q_(i)lm
  Array1 <double> Plm(l+1);
  AssociatedLegendrePolynomials(l,cosTheta,Plm);

  double factorM0 = sqrt( (2*l+1) / (4.0*pi) );
  Complex term = f*Complex(factorM0*Plm[0],0.0); // case m=0
  Ylm(i,l) += term;
  Ylm(j,l) += term;

  for(int m=1; m<=l; m++) {
    double factor = f * factorM0 * sqrt( Factorial(l-m) / Factorial(l+m) ); // remember Factorial() now returns a double
    const Complex I(0.0,1.0);
    const Complex term1 = factor  * Plm[m] * exp( m*phi*I);
    Ylm(i,l+m) += term1;
    Ylm(j,l+m) += term1;
    if (m%2!=0) factor *= -1.0; // flip sign only for negative and odd 'm'
    const Complex term2 = factor  * Plm[m] * exp(-m*phi*I);
    Ylm(i,l-m) += term2;
    Ylm(j,l-m) += term2;
  }
}

#endif // _BONDORDER_
