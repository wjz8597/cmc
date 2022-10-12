////////////////////////////////////////////////////////////////////
//
// Burkhard Militzer                                  Urbana 10-1-99
//
// Simple class to keep track of the acceptance ratios
//
#include <iomanip>
#include "Standard.h"
#include "AccRatio.h"

int AccRatio::nameLengthMax=0;

ostream& operator<<(ostream &os, const AccRatio & ar ) {
  if (!ar.Flag()) return os;

  const ios::fmtflags current = os.flags();
  const int p = os.precision();

  if (!ar.infoFlag) error("Cannot output accRatio without labels set\n");
  os << " " << setw(AccRatio::nameLengthMax) << ar.name;
  if (ar.nMovers>=0) 
    os << " movers=" << setw(2) << ar.nMovers+1;
  //  else 
  //    os << " moves=all";
  os << " " << setw(4) << ar.label;
  os << " tri=" << setw(7) << ar.nTrial;
  os << " acc=" << setw(7) << ar.nAccept;
  if (ar.nTrial>0) {
    os.setf(ios::fixed);
    os << " ratio= " << setprecision(4) << double(ar.nAccept)/double(ar.nTrial);
  }
  os << endl;

  os.precision(p);
  os.setf(current);

  return os;
}
