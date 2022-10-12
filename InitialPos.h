/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Burkhard Militzer                                  Dresden, 9-4-04  //
//                                                                     //
// Initial position on a lattice                                       //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#ifndef _INITIALPOS_
#define _INITIALPOS_

#include "Array.h"
#include "Point.h"

void InitialPos(const int np, const double cell, 
		const double rRandom, Array1 <Point> & r);

#endif // _INITIALPOS_
