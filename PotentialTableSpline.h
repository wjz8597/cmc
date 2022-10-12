/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Store the simple pair-potential derived from V(r) table                //
//  Now use LSpline (slower) instead of linear interpolation to get        //
//  second derivative for lattic dynamics calculations                     //
//                                                                         //
//  V(r) has units Ha/a_B                                                   //
//  B. Militzer                                    Washington 05-11-07     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _POTENTIALTABLESPLINE_
#define _POTENTIALTABLESPLINE_

// Do not inherite this of TwoBodyPotential (as is should be) because TwoBodyPotential uses Point/BasicPoint/PeriodicPoint
// and both conflict with the Ewald class
// do nto waste time on fixing this!
// --> This is not na issue any longer since we elimiated Ewald from the list of potentials

#include "TwoBodyPotential.h"
#include "LSpline.h"
#include "Parser.h"

class PotentialTableSpline : public IsotropicPotential {
  // class PotentialTableSpline {
 public:
  LSpline s;

  PotentialTableSpline() {} // can now be called with a constructor without arguments
  PotentialTableSpline(const string & filename, const int col1=1, const int col2=2, const double f=1.0) {
    Init(filename,col1,col2,f);
  }
  void Init(const string & filename, const int col1=1, const int col2=2, const double f=1.0) {
    Parser p(filename);
    while (p.ReadLine()) {
      double r = p.GetDouble(col1-1);
      double v = p.GetDouble(col2-1);
      if (s.n>0 && s.x.Last()>=r) error("r points are in increasing",s.x.Last(),r);
      s.PushBack(r,v*f);
    }
    
    //    s.Init();
    s.InitUsingFiniteDifferences(); // should call this routine so that dV/dr is continous at least (L0 but not L1)

    cout << "Read in potential file with " << s.n
	 << " points successfully" << endl;

    cout << "Setting up V(r) spline." << endl;
    if (f!=1.0) cout << "Using " << f << "*V(r) instead of V(r)." << endl;
  };

  double OneBodyTerm() const {
    // enter any one-body term here
    return 0.0;
  }

  double GetRMin() const {
    return s.FirstX();
  }
  double GetRMax() const {
    return s.LastX();
  }

  double Vm(const double r) const {
    return s.Interpolate(r);
  }

  double DVm(const double r) const {
    return s.Derivative(r);
  }

  double DDVm(const double r) const {
    return s.SecondDerivative(r);
  }


};

#endif // _POTENTIALTABLESPLINE_
