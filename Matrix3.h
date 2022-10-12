//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Matrix operations for (3x3) matrices                                 //
//                                                                      //
// Burkhard Militzer                              Washington, 10-26-05  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef _MATRIX3_
#define _MATRIX3_

#include "Array.h"
#include "MatrixAlgebra.h"
#include "Point.h"

class Matrix3 : public Array2 <double> {
 public:
  Matrix3():Array2 <double>(Point::nDim,Point::nDim) {};
  Matrix3(const double x):Array2 <double>(Point::nDim,Point::nDim,x) {};

  BPoint operator[] (const int i) const {
    BPoint v;
    for(int d=0; d<Point::nDim; d++) {
      v[d]=(*this)(i,d);
    }						
    return v;
  }
  
  static Matrix3 Identity() {
    Matrix3 id(0.0);
    for(int d=0; d<Point::nDim; d++) {
      id(d,d)=1.0;
    }
    return id;
  }

  Matrix3 Inverse() const {
    Matrix3 A1;
    ::InvertMatrix((*this),A1,Point::nDim);
    Matrix3 one = (*this)*A1;
    double detOne = one.Determinant();
    if (fabs(detOne-1.0)>1e-10) warning("Inverse problem");
    return A1;
  }

  Matrix3 Transpose() const {
    Matrix3 AT(*this);
    TransposeMatrix(AT,Point::nDim);
    return AT;
  }

  double Determinant() const {
    return ::Determinant(*this,Point::nDim);
  }

  Matrix3 operator*(const double x) {
    for(int d=0; d<Point::nDim; d++) {
      for(int dd=0; dd<Point::nDim; dd++) {
	(*this)(dd,d)*=x;
      }
    }
    return *this;
  }

  BPoint operator*(const BPoint & x) {
    BPoint y(0.0);
    for(int d=0; d<Point::nDim; d++) {
      for(int dd=0; dd<Point::nDim; dd++) {
	y[dd] += (*this)(dd,d)*x[d];
      }
    }
    return y;
  }

  Matrix3 operator*(const Matrix3 & A) const {
    Matrix3 C;
    MultiplyMatrices(*this,A,C,Point::nDim);
    return C;
  }

  // "this *= A" equals "this=this*A" but not "this=A*this" !!!
  Matrix3 & operator*=(const Matrix3 & A) {
    Matrix3 C;
    MultiplyMatrices(*this,A,C,Point::nDim);
    *this = C;
    return *this;
  }

  Matrix3 operator+(const Matrix3 & A) const {
    Matrix3 C;
    for(int i=0; i<Point::nDim; i++) {
      for(int j=0; j<Point::nDim; j++) {
	C(i,j) = (*this)(i,j) + A(i,j);
      }
    }
    return C;
  }

  friend ostream& operator<<(ostream & os, const Matrix3 & A) {
    os << endl;
    for(int i=0; i<Point::nDim; ++i) {
      os << "( ";
      for(int j=0; j<Point::nDim; ++j) {
	if (j>0) 
	  os << " , ";
	os << A(i,j);
      }
      os << " )" << endl;
    }
    return os;
  }
};

#endif // _MATRIX3_
