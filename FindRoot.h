/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  1D search routines for roots of equation                               //
//  uses templates                                                         //
//                                                                         //
//  B. Militzer                                     Livermore 09-03-01     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _FINDROOT_
#define _FINDROOT_

#include "Standard.h"

template <class T>
double FindRootBisection(const double x1, 
			 const double x2, 
			 const double tolX,
			 const T & f) {
  if (tolX<0.0)
    error("The accuracy parameter must be positive",tolX);
  
  double f1=sign(f(x1));
  double f2=sign(f(x2));
  
  if (f1==f2)
    error("FindRootBisection: same sign at both interval boundaries",x1,x2,f1,f2);
  
  const double acc=fabs(x2-x1)*tolX;
  double xx1=x1;
  double xx2=x2;
  double xx;
  do {

    xx=(xx2+xx1)*0.5;
    double ff=sign(f(xx));

    if (ff==f1) xx1=xx;
    else xx2=xx;

    //    Write3(xx1,xx2,fabs(xx2-xx1)/acc);

  } while (fabs(xx2-xx1)>acc);
  return xx;
}

template <class T>
double FindRootRegulaFalsi(const double x1, 
			   const double x2, 
			   const double tolX,
			   const T & f) {
  if (tolX<0.0)
    error("The accuracy parameter must be positive",tolX);

  double f1=f(x1);
  double f2=f(x2);
  
  const double acc=fabs(x2-x1)*tolX;

  double xx1=x1;
  double xx2=x2;
  double xx;
  do {
    xx=(f1*xx2-f2*xx1)/(f1-f2);
    double ff=f(xx);
    xx1=xx2; f1=f2;
    xx2=xx;  f2=ff;
  } while (fabs(xx2-xx1)>acc);
  return xx;
}

template <class T>
double FindRootStepper(const double x1, 
		       const double x2, 
		       const double tolX,
		       const double tolF,
		       const T & f) {
  if (tolX<0.0)
    error("The accuracy parameter must be positive",tolX);
  if (tolF<0.0)
    error("The accuracy parameter must be positive",tolF);

  const double r=0.2;
  double f1=f(x1);
  double f2=f(x2);

  if (sign(f1)!=sign(f2)) 
    return FindRootBisection(x1,x2,tolX,f);

  double xx,dx,ff;
  if (fabs(f1)<fabs(f2)) {
    xx=x1;
    dx=(x2-x1)*r;
    ff=f1;
  } else {
    xx=x2;
    dx=(x1-x2)*r;
    ff=f2;
  }

  for(;;) {
    const double xxn=xx+dx;
    const double ffn=f(xxn);
    //    Write5(xx,xxn,dx,ff,ffn);

    if (sign(ffn)!=sign(ff))
      return FindRootBisection(xx,xxn,tolX,f);

    if (sign(ffn)!=sign(ff) || fabs(ffn)>fabs(ff)) {
      //      Write2(int(sign(ffn)!=(ff)),int(fabs(ffn)>fabs(ff)));
      dx *= -r;
      if (fabs(dx)<tolX && fabs(ffn)<tolF) {
	return (xxn+xx)*0.5;
      }

      const double overrun=1e-8;
     if (fabs(dx)<tolX*overrun) 
	error("FindRootStepper: Could not find root (dx got too small)",dx,ffn);
    }
    xx=xxn;
    ff=ffn;
  }
  return 0.0;
}

template <class T>
double FindRootStepper(const double x, 
		       const double tolX,
		       const double tolF,
		       const T & f) {
  const double r=0.1;
  return FindRootStepper( (1.0-r)*x, (1.0+r)*x, tolX, tolF, f );
}

template <class T>
double FindRootBrancher(const double x1, 
			const double x2, 
			const double tolX,
			const double tolF,
			const T & f) {
  if (tolX<0.0)
    error("The accuracy parameter must be positive",tolX);
  if (tolF<0.0)
    error("The accuracy parameter must be positive",tolF);

  double f1=f(x1);
  double f2=f(x2);

  if (sign(f1)!=sign(f2)) 
    return FindRootBisection(x1,x2,tolX,f);

  // are there any solution inside the interval?
  const int nn=10;
  for(int i=1; i<nn; i++) {
    double xx = x1+(x2-x1)*i/nn;
    double ff = f(xx);
    if (sign(ff)!=sign(f1)) {
      if (i<=nn/2) return FindRootBisection(x1,xx,tolX,f);
      else return FindRootBisection(xx,x2,tolX,f);
    }
  }

  // now look outside
  const int nTrials=100;
  double xx1=x1;
  double xx2=x2;
  for(int i=0; i<nTrials; i++) {
    xx1 /= (1.0+1.0/nn);
    xx2 *= (1.0+1.0/nn);
    double ff1 = f(xx1);
    double ff2 = f(xx2);
    if (sign(f1)!=sign(ff1)) return FindRootBisection(x1,xx1,tolX,f);
    if (sign(f2)!=sign(ff2)) return FindRootBisection(x2,xx2,tolX,f);
  }

  error("Could not find any root in and outside of interval",x1,x2,xx1,xx2);
  return 0.0;
}

template <class T>
double FindRootBrancher(const double x, 
			const double tolX,
			const double tolF,
			const T & f) {
  const double r=0.1;
  return FindRootBrancher( (1.0-r)*x, (1.0+r)*x, tolX, tolF, f );
}

#endif // _FINDROOT_
