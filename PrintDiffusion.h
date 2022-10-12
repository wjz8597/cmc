/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Print diffusion distance for each type                                  //
//                                                                         //
// Burkhard Militzer                                 Berkeley, 02-19-09    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _PRINTDIFFUSION_
#define _PRINTDIFFUSION_

#include "State.h"

void PrintDiffusion(const State & s1, const State & s2, const int t, const double dHansen, const double time, const bool final=true);
void PrintDiffusion(const State & s1, const State & s2, const int t, const double time, const bool final=true);
void PrintDiffusionAllTypes(const State & s1, const State & s2, const double time, const bool final=true);

#endif // _PRINTDIFFUSION_
