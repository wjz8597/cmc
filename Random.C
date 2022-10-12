/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// Interface to different RANDOM NUMBER GENERATORS                             //
// Make generators exchangeable                                                //                 
//                                                                             //
// Burkhard Militzer                                    Urbana 1998            //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

#include "Random.h"

#ifndef SIMPLE_SPRNG 
int* sprngStream;
#endif
