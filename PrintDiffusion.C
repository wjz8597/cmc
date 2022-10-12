/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Print diffusion distance for each type                                  //
//                                                                         //
// Burkhard Militzer                                 Berkeley, 02-19-09    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "PrintDiffusion.h"

void PrintDiffusion(const State & s1, const State & s2, const int t, const double dHansen, const double time, const bool final) {
  if (final) cout << "Diffusion analysis by only comparing the initial and final snapshots:" << endl;
  const Point rCM1 = s1.CenterOfMassType(t);
  const Point rCM2 = s2.CenterOfMassType(t);
  //  Write2(t,rCM1);
  //  Write2(t,rCM2);
  //  Write2(t,rCM2-rCM1);
  //  Write(s1.CenterOfMassEinsteinAtoms());
  //  Write(s2.CenterOfMassEinsteinAtoms());
 Form f(gen,2,2);
  f.left();
  
  double d2=0.0;
  int n=0;
  for(int i=0; i<s1.atomsByType[t].Size(); i++) {
    int ii = s1.atomsByType[t][i];
    if (s1.names[ii]!=s2.names[ii]) error("Type mismatch for particle",ii);
    Point di = s2.r[ii]-s1.r[ii]-(rCM2-rCM1); // would yield infinity at some large block - unclear why
    d2 += di.Norm2();
    //      Write3(i,di,d2);
    n++;
  }
  //  Write2(name,d2);
  if (n==0) return;
  d2 /= (double) n;

  string name = s1.listOfTypes[t];
  if (final) {
    cout << f(name)  << ": Final diff sqrt(<r^2>)/d_min_image    = " << sqrt(d2)/Point::MinimumImageDistance() << endl;
    cout << f(name)  << ": Final diff <r^2>                      = " << d2 << endl;
    if (d2>0.0) 
      cout << f(name)  << ": Estimated Einstein k = 3 kb T / <r^2> = " << 3.0*s1.GetTemperature()/d2 << endl;
    
    if (time>0.0) {
      double dCoeff   = d2/(6.0*time);
      double dCoeffSI = dCoeff*sqr(PC::ab)/PC::AUToS;
      cout << " D[a.u.]= " << dCoeff;
      //    cout << " D[SI]= "   << dCoeffSI;
      cout << " D[cgs]= "  << dCoeffSI*1e4;
      if (dHansen>0.0) cout << " D/DHansen= "  << dCoeff/dHansen;
    }
    cout << endl;
  } else {
    cout << " " << f(name)  << ": diffusion so far <r^2> = " << d2 << endl;
  }
  //  Write3(dCoeff,dCoeffSI,dCoeffSI*1e4);
}

void PrintDiffusion(const State & s1, const State & s2, const int t, const double time, const bool final) {
  PrintDiffusion(s1,s2,t,0.0,time,final);
}

void PrintDiffusionAllTypes(const State & s1, const State & s2, const double time, const bool final) {
  for(int t=0; t<s1.GetNTypes(); t++)  
    PrintDiffusion(s1,s2,t,time,final);
}

