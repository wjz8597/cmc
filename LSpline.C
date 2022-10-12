/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Modified Spline interpolation with linear extrapolation                 //
//                                                                         //
// Burkhard Militzer                                 Washington 08-03-04   //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include "LSpline.h"
#include "Vector.h"

class XY {
public:
  double x,y;
  bool operator< (const XY & xy) const {
    return x<xy.x;
  }
};

void LSpline::Sort(Array1 <double> & x,
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

void LSpline::CalculateY2(Array1 <double> & x,
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
}

inline int LSpline::GetIndex(const Array1 <double> & x, const int n, const double xx) {
  if (xx<x[0]) warning("Spline: x too low",xx,x[0]);
  if (xx>x.Last()) warning("Spline x too high",xx,x.Last());
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

inline double LinearInterpolation(const double x1, const double x2, 
				  const double y1, const double y2, 
				  const double xx) {
  double q2 = (xx-x1)/(x2-x1);
  double q1 = (x2-xx)/(x2-x1);
  return q1*y1+q2*y2;
}

double LSpline::Interpolate(const Array1 <double> & x,
			   const Array1 <double> & y, 
			   const int n, 
			   const Array1 <double> & y2,
			   const double xx) {
  if (xx<x[0])   return LinearInterpolation(x[0],  x[1],  y[0],  y[1],  xx);
  if (xx>x[n-1]) return LinearInterpolation(x[n-2],x[n-1],y[n-2],y[n-1],xx);
  int k1=GetIndex(x,n,xx);
  int k2=k1+1;
  double h=x[k2]-x[k1];
  double a=(x[k2]-xx)/h;
  double b=1.0-a;
  return a*y[k1] + b*y[k2] + 
    ( (a*a*a-a)*y2[k1]+(b*b*b-b)*y2[k2] ) * (h*h)/6.0;
}

double LSpline::InterpolateLinearly(const Array1 <double> & x,
				    const Array1 <double> & y, 
				    const int n, 
				    const Array1 <double> & y2,
				    const double xx) {
  if (xx<x[0])   return LinearInterpolation(x[0],  x[1],  y[0],  y[1],  xx);
  if (xx>x[n-1]) return LinearInterpolation(x[n-2],x[n-1],y[n-2],y[n-1],xx);
  int k1=GetIndex(x,n,xx);
  int k2=k1+1;
  return LinearInterpolation(x[k1],  x[k2],  y[k1],  y[k2],  xx);
}

double LSpline::Derivative(const Array1 <double> & x,
			  const Array1 <double> & y, 
			  const int n, 
			  const Array1 <double> & y2,
			  const double xx) {
  if (xx<x[0])   return (y[1]  -y[0]  )/(x[1]  -x[0]  );
  if (xx>x[n-1]) return (y[n-1]-y[n-2])/(x[n-1]-x[n-2]);
  int k1=GetIndex(x,n,xx);
  int k2=k1+1;
  double h=x[k2]-x[k1];
  double a=(x[k2]-xx)/h;
  double b=1.0-a;
  
  return (y[k2]-y[k1])/h + h/6.0* 
    ( (1.0-3.0*a*a)*y2[k1] + (3.0*b*b-1.0)*y2[k2] ); // bug fixed 6/18/06: y2[k2] instead of y2[k1]
}

double LSpline::SecondDerivative(const Array1 <double> & x,
				const Array1 <double> & y, 
				const int n, 
				const Array1 <double> & y2,
				const double xx) {
  if (xx<x[0])   return 0.0;
  if (xx>x[n-1]) return 0.0;

  int k1=GetIndex(x,n,xx);
  int k2=k1+1;
  double h=x[k2]-x[k1];
  double a=(x[k2]-xx)/h;
  double b=(xx-x[k1])/h;
  
  return a*y2[k1]+b*y2[k2];
}

/*
int main() {
  const int n=10;
  Array1 <double> x(n);
  Array1 <double> y(n);

  for(int i=0; i<n; i++) {
    x[i] = i;
    y[i] = sin(i)+cos(10*i);
    Write2(x[i],y[i]);
  }


  Array1 <double> y2;
  Spline s;
  //  s.CalculateY2(x,y,n,y2);
  s.CalculateY2(x,y,n,y2,false,5,5);

  const int m=10;
  for(int i=0; i<n*m; i++) {
    double xx = double(i)/double(m);
    double yy = s.Interpolate(x,y,n,y2,xx);
    double yyd = s.Derivative(x,y,n,y2,xx);
    double yydd = s.SecondDerivative(x,y,n,y2,xx);
    Write4(xx,yy,yyd,yydd);
  }
}
*/
