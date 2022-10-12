/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Classical MC+MD program                                             // 
//                                                                     //
// Burkhard Militzer                              Berkeley, 07-01-2009 //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#ifndef _WRITEVASPPOSCARFILE_
#define _WRITEVASPPOSCARFILE_

#include "Standard.h"
#include "State.h"

void WriteVASPPOSCARFile(ostream & of, const State & s, const State & sp, const double dt, const Array1 <int> & typeOrder);
void WriteVASPPOSCARFile(const string & fn, const State & s, const State & sp, const double dt, const Array1 <int> & typeOrder);
void WriteVASPPOSCARFile(const string & fn, const State & s, const State & sp, const double dt);
void WriteVASPPOSCARFile(const string & fn, const State & s);

#endif // _WRITEVASPPOSCARFILE_
