/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Print Error Bars in Short Chemistry Notation                        //
//                                                                     //
// Burkhard Militzer                                Washington 8-17-04 //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#include "PrintErrorBar.h"

string PrintAverageResultShort(const double av, const double e) {
  if (e<=0.0) {
    ostringstream os;
    os.setf(ios::fixed,ios::floatfield);
    os.precision(6);
    os << av;
    return os.str();
  }

  // how many digit 1 --> 1020(20) 1.02(2) 
  //             or 2 --> 1018(23) 1.018(23)
  const int errorDigits=2;    
  //  const int errorDigits=1;    

  // add an extra digit if error is 116 --> 120 rather than 100
  const double errorLimit=2.0; // default
  //  const double errorLimit=1.0; 

  int pp=0;
  if (errorDigits == 1) {
    if (e/pow(10.0,floor(log10(e)))<errorLimit) pp=1; // add extra digit for error between 1.0 and 2.0
  }

  int p = int ( floor(log10(e))-(errorDigits-1) ); // e>10^p
  p -= pp;
  double f = pow(10.0,p);
  
  double ee = f*nint(e/f);  // round the error
  //  double aa = f*nint(av/f);
  
  //  printf("ef= %f p= %d pp=%d\n",ef,p,pp);
  
  int pe,pa;
  if (p<0) {
    pa = -p; // print too many decimal digits
    pe = -p;

    if (ee<1.0) { // this condition means we print 102.0(1.1) instead of 102.0(11)
      pe = 0;
      ee /= f;
    }

  } else {
    pa = 0;   // no decimal digit 
    pe = 0;   // no decimal digit 
  }

  ostringstream os;
  os.setf(ios::fixed,ios::floatfield);
  os.precision(pa);
  os << av;
  os << "(";
  os.precision(pe);
  os << ee;
  os << ")";

  return os.str();
}

