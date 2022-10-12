/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Spline interpolation                                                    // 
//                                                                         //
// Burkhard Militzer                                  Livermore 4-27-01    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include "Spline.h"
#include "Vector.h"

class XY {
public:
  double x,y;
  bool operator< (const XY & xy) const {
    return x<xy.x;
  }
};

void Spline::Sort(Array1 <double> & x,
		  Array1 <double> & y, 
		  const int n) {
  bool sortFlag=false;
  for (int i=0; i<n-1; i++) {
    if (x[i]>x[i+1]) {
      sortFlag=true;
      break;
    }
  }
  if (!sortFlag) return;

  Vector <XY> v(n);
  for(int i=0; i<n; i++) {
    v[i].x = x[i];
    v[i].y = y[i];
  }
  sort(v.begin(),v.end());
  for(int i=0; i<n; i++) {
    x[i] = v[i].x;
    y[i] = v[i].y;
  }
  
}

void Spline::CalculateY2(Array1 <double> & x,
			 Array1 <double> & y, 
			 const int n,
			 Array1 <double> & y2,
			 const bool natural, 
			 const double ypLow, 
			 const double ypHigh) {
  if (n<2) 
    error("Spline::CalculateY2 call for empty spline",n);
  if (n!=x.Size() || n!=y.Size()) 
    error("Array size problem",n,x.Size(),y.Size());
  
  Sort(x,y,n);
  
  Array1 <double> u(n-1);
  
  y2.Resize(n);

  if (natural) {
    u[0] =0.0;
    y2[0]=0.0;
  } else {
    y2[0] = -0.5;
    u[0]  = 3.0/(x[1]-x[0]) * ((y[1]-y[0])/(x[1]-x[0])-ypLow);
  }
  
  for (int i=1; i<n-1; i++) {
    if (x[i-1]>=x[i] || x[i-1]>=x[i+1]) 
      error("Spline:x values are not monotonic",i,x[i-1],x[i],x[i+1]);

    const double sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    const double p  =sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  
  double qn,un;
  if (natural) {
    qn= 0.0;
    un= 0.0;
  } else {
    qn= 0.5;
    un= 3.0/(x[n-1]-x[n-2]) * (ypHigh-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }

  y2[n-1] = (un-qn*u[n-2]) / (qn*y2[n-2]+1.0);
  for (int k=n-2; k>=0; k--)
    y2[k]=y2[k]*y2[k+1]+u[k];

  //  y2=0.0;
  //  warning("0.0=y2");
}

inline int Spline::GetIndex(const Array1 <double> & x, const int n, const double xx) {
  int k1=0;
  int k2=n-1;
  
  while (k2-k1 > 1) {
    int k=(k2+k1) >> 1;
    if (x[k] > xx) k2=k;
    else k1=k;
  }

  if (x[k2] == x[k1]) error("SplineInterpolation: Bad x-vector input",k1,k2,n,x[k1],x[k2],xx);
  return k1;
}

double Spline::Interpolate(const Array1 <double> & x,
			   const Array1 <double> & y, 
			   const int n, 
			   const Array1 <double> & y2,
			   const double xx) {
  int k1=GetIndex(x,n,xx);
  int k2=k1+1;
  double h=x[k2]-x[k1];
  double a=(x[k2]-xx)/h;
  double b=1.0-a;
  return a*y[k1] + b*y[k2] + 
    ( (a*a*a-a)*y2[k1]+(b*b*b-b)*y2[k2] ) * (h*h)/6.0;
}

double Spline::Derivative(const Array1 <double> & x,
			  const Array1 <double> & y, 
			  const int n, 
			  const Array1 <double> & y2,
			  const double xx) {
  int k1=GetIndex(x,n,xx);
  int k2=k1+1;
  double h=x[k2]-x[k1];
  double a=(x[k2]-xx)/h;
  double b=1.0-a;
  
  return (y[k2]-y[k1])/h + h/6.0* 
    ( (1.0-3.0*a*a)*y2[k1] + (3.0*b*b-1.0)*y2[k2] ); // bug fixed 6/18/06: y2[k2] instead of y2[k1]
}

double Spline::SecondDerivative(const Array1 <double> & x,
				const Array1 <double> & y, 
				const int n, 
				const Array1 <double> & y2,
				const double xx) {
  int k1=GetIndex(x,n,xx);
  int k2=k1+1;
  double h=x[k2]-x[k1];
  double a=(x[k2]-xx)/h;
  double b=(xx-x[k1])/h;
  
  return a*y2[k1]+b*y2[k2];
}

double Spline::IntegralTerm1(const double x1, const double x2, const double xx) {
  //  Simplify[Integrate[a*y1,{x,x1,xx}]]
  //         -((x1 - xx) (x1 - 2 x2 + xx) y1)
  //Out[37]= --------------------------------
  //                   2 (x1 - x2)
  return 0.5*(xx-x1)/(x1-x2) * (x1 - 2.0*x2 + xx);
}

double Spline::IntegralTerm2(const double x1, const double x2, const double xx) {
  //Simplify[Integrate[b*y2,{x,x1,xx}]]
  //                    2
  //         -((x1 - xx)  y2)
  //Out[38]= ----------------
  //           2 (x1 - x2)
  return sqr(xx-x1)/(2.0*(x2-x1));
}

double Spline::IntegralTerm3(const double x1, const double x2, const double xx) {
  //In[32]:= Integrate[(a^3-a)*ys1*h*h/6,{x,x1,xx}]
  //                  2                 2
  //         (x1 - xx)  (x1 - 2 x2 + xx)  ys1
  //Out[32]= --------------------------------
  //                   24 (x1 - x2)
  return sqr(xx-x1)*sqr(x1 - 2.0*x2 + xx) / ( 24.0 * (x1-x2) );
}

double Spline::IntegralTerm4(const double x1, const double x2, const double xx) {
  //In[33]:= Integrate[(b^3-b)*ys2*h*h/6,{x,x1,xx}]
  //                  2    2                 2               2
  //         (x1 - xx)  (x1  - 4 x1 x2 + 2 x2  + 2 x1 xx - xx ) ys2
  //Out[33]= ------------------------------------------------------
  //                              24 (x1 - x2)
  return sqr(xx-x1)*(x1*x1 - 4.0*x1*x2 + 2.0*x2*x2 + 2.0*x1*xx - xx*xx) / ( 24.0 * (x1-x2) );
}

// integral dx_{x1...xx}
double Spline::IntegralInterval1(const double x1, const double x2, const double xx, const double y1, const double y2, const double ys1, const double ys2) {
  return 
    y1 *IntegralTerm1(x1,x2,xx)+y2 *IntegralTerm2(x1,x2,xx)+
    ys1*IntegralTerm3(x1,x2,xx)+ys2*IntegralTerm4(x1,x2,xx);
}

// integral dx_{x1...x2}
double Spline::IntegralFullInterval(const double h, const double y1, const double y2, const double ys1, const double ys2) {
  //  In[39]:= Integrate[f,{x,x1,x2}]
  //                                              2
  //         (x1 - x2) (-12 y1 - 12 y2 + (x1 - x2)  (ys1 + ys2))
  //Out[39]= ---------------------------------------------------
  //                                 24
  return h * ( 12.0*(y1+y2) - h*h*(ys1+ys2) ) / 24.0;
}

// integral dx_{xx...x2}
double Spline::IntegralInterval2(const double x1, const double x2, const double xx, const double y1, const double y2, const double ys1, const double ys2) {
  return 
    IntegralFullInterval(x2-x1,y1,y2,ys1,ys2) - IntegralInterval1(x1,x2,xx,y1,y2,ys1,ys2);
}

double Spline::Integrate(const Array1 <double> & x,
			 const Array1 <double> & y, 
			 const int n, 
			 const Array1 <double> & y2,
			 const double xx1, 
			 const double xx2) {
  if (xx2<xx1) return -Integrate(x,y,n,y2,xx2,xx1);

  int k1=GetIndex(x,n,xx1);
  int k2=GetIndex(x,n,xx2); // definition differs from above!
  int kk1 = k1+1;
  int kk2 = k2+1;

  //  Write2(k1,k2);

  if (k1==k2) {
    /*
    int nn=100000;
    double ss=0.0;
    for(int i=0; i<nn; i++) {
      ss+=Interpolate(x,y,n,y2,xx1+(xx2-xx1)*(0.5+i)/nn)*(xx2-xx1)/nn;
    }
    Write6(xx1,xx2,k1,k2,n,ss);
    */

    //    double f1 = IntegralInterval1(x[k1],x[kk1],xx2,y[k1],y[kk1],y2[k1],y2[kk1])
    //      -IntegralInterval1(x[k1],x[kk1],xx1,y[k1],y[kk1],y2[k1],y2[kk1]);
    //    double f2 = IntegralFullInterval(x[kk1]-x[k1],y[k1],y[kk1],y2[k1],y2[kk1]);
    //    Write2(f1,f2);

    return 
       IntegralInterval1(x[k1],x[kk1],xx2,y[k1],y[kk1],y2[k1],y2[kk1])
      -IntegralInterval1(x[k1],x[kk1],xx1,y[k1],y[kk1],y2[k1],y2[kk1]);
  }
  double sInt = IntegralInterval2(x[k1],x[kk1],xx1,y[k1],y[kk1],y2[k1],y2[kk1]);
  sInt += IntegralInterval1(x[k2],x[kk2],xx2,y[k2],y[kk2],y2[k2],y2[kk2]);

  for(int k=k1+1; k<k2; k++) {
    int kk = k+1;
    sInt += IntegralFullInterval(x[kk]-x[k],y[k],y[kk],y2[k],y2[kk]);
  }

  /*
  int nn=100000;
  double ss=0.0;
  for(int i=0; i<nn; i++) {
    ss+=Interpolate(x,y,n,y2,xx1+(xx2-xx1)*(0.5+i)/nn)*(xx2-xx1)/nn;
  }
  Write7(xx1,xx2,k1,k2,n,sInt,ss);
  */

  return sInt;
}

/*
int main() {
  const int n=10;
  Array1 <double> x(n);
  Array1 <double> y(n);

  Spline s;
  for(int i=0; i<n; i++) {
    x[i] = i*1.5;
    //    x[i] = i*1.0;
    y[i] = sin(i)+cos(10*i);
    Write2(x[i],y[i]);
    s.PushBack(x[i],y[i]);
  }


  Array1 <double> y2;
  Spline ss;
  s.CalculateY2(x,y,n,y2);
  //  ss.CalculateY2(x,y,n,y2,false,5,5);
  s.Init();
  
  //  Write(s.Integrate(x[1],x[2]*1.5));
  //  Write(s.Integrate(x[1]*0.1,1.5));
  Write(s.Integrate(0.5,13.0+1e-15));
  error("Q");

  const int m=10;
  for(int i=0; i<n*m; i++) {
    double xx = double(i)/double(m);
    double yy   = ss.Interpolate(x,y,n,y2,xx);
    double yyd  = ss.Derivative(x,y,n,y2,xx);
    double yydd = ss.SecondDerivative(x,y,n,y2,xx);
    double y   = s.Interpolate(xx);
    double yd  = s.Derivative(xx);
    double ydd = s.SecondDerivative(xx);
    double yII = s.Integrate(x[1],x[n/2]);
    Write4(xx,yy,yyd,yydd);
    Write5(xx,y,yd,ydd,yII);
  }
}
*/
