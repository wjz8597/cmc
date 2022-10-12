/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Burkhard Militzer                                  Dresden, 9-4-04  //
//                                                                     //
// Initial position on a lattice                                       //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#include "Standard.h"
#include "InitialPos.h"

enum LatticeType {SC, BCC, FCC};

int NDivisions(const int np, const double q) {
  return int(floor(pow(np/q,1.0/3.0)+1.0-1e-10));
}

// return number of empty lattice site
int NEmpty(const int np, const int q) {
  int m = NDivisions(np,double(q));
  int nE = q*m*m*m-np;
  //  Write4(np,q,m,nE);
  return nE;
}

inline void Periodic(double &c, const double cell) {
  c -= cell*floor(c/cell+0.5);
}

void Set(const double ix, const double iy, const double iz, const double cell, const int nD, 
	 int & ii, const int np, const double rRandom, Array1 <Point> & r) {
  if (ii>=np) return;
  
  Point rr;
  rr[0] = (ix+0.5)*cell/nD -cell/2.0;
  rr[1] = (iy+0.5)*cell/nD -cell/2.0;
  rr[2] = (iz+0.5)*cell/nD -cell/2.0;

  rr[0] += (2.0*drand48()-1.0)*rRandom;
  rr[1] += (2.0*drand48()-1.0)*rRandom;
  rr[2] += (2.0*drand48()-1.0)*rRandom;
  Periodic(rr[0],cell);
  Periodic(rr[1],cell);
  Periodic(rr[2],cell);

//    cout << rr[0] << "  " 
//         << rr[1] << "  " 
//         << rr[2] << endl;
//    cout << rr[0]/cell << "  " 
//         << rr[1]/cell << "  " 
//         << rr[2]/cell << endl;
  r.PushBack(rr);

  ii++;
}

void Set(const int np, const double q, const double cell, 
	 const LatticeType lt, const double rRandom, Array1 <Point> & r) {
  int m = NDivisions(np,q);
  cout << "Number of divisions = " << m << endl;
  //  Write2(q,m);
  
  int ii=0;
  for(int ix=m-1; ix>=0; ix--) {
    for(int iy=m-1; iy>=0; iy--) {
      for(int iz=m-1; iz>=0; iz--) {
	if (lt==SC) {
	  Set(ix,iy,iz,cell,m,ii,np,rRandom,r);
	}
	if (lt==BCC) {
	  Set(ix    ,iy    ,iz    ,cell,m,ii,np,rRandom,r);
	  Set(ix+0.5,iy+0.5,iz+0.5,cell,m,ii,np,rRandom,r);
	}
	if (lt==FCC) {
	  Set(ix    ,iy    ,iz    ,cell,m,ii,np,rRandom,r);
	  Set(ix+0.5,iy+0.5,iz    ,cell,m,ii,np,rRandom,r);
	  Set(ix+0.5,iy    ,iz+0.5,cell,m,ii,np,rRandom,r);
	  Set(ix    ,iy+0.5,iz+0.5,cell,m,ii,np,rRandom,r);
	}
      }
    }
  }
}

void InitialPos(const int np, const double cell, 
		const double rRandom, Array1 <Point> & r) {
  if (cell<=0.0) error("cell",cell);
  if (np<=0) error("NP",np);
  if (rRandom<0.0) error("random L",rRandom);

  int nSC = NEmpty(np,1);
  int nBCC= NEmpty(np,2);
  int nFCC= NEmpty(np,4);

  r.Resize(0); // delete any old coord since we use PushBack later on

  if (nSC<=nBCC && nSC<=nFCC) {
    cout << "Setting to s.c. lattice" << endl;
    Set(np,1.0,cell,SC,rRandom,r);
  } else if (nBCC<=nSC && nBCC<=nFCC) {
    cout << "Setting to b.c.c. lattice" << endl;
    Set(np,2.0,cell,BCC,rRandom,r);
  } else {
    cout << "Setting to f.c.c. lattice" << endl;
    Set(np,4.0,cell,FCC,rRandom,r);
  }
}
      
