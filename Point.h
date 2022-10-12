//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Provide TYPEDEFs for the actual Point class being used everywhere    //
//                                                                      //
// Burkhard Militzer                                    Paris 4-23-99   //
// Modified class structure                         Livermore 8-08-02   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef _POINT_
#define _POINT_

//#define POINTGEN
#ifdef POINTGEN

#include "PeriodicPointGen.h"
typedef PeriodicPointGen<NDIM> Point;

#else

#include "PeriodicPoint.h"
typedef PeriodicPoint<NDIM> Point;

#endif

typedef BasicPoint<NDIM>    BPoint;

#endif  // _POINT_

