/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Physical Constants                                                      // 
//                                                                         //
// Burkhard Militzer                                  Livermore 4-27-01    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _PHYSICS_
#define _PHYSICS_

#include <math.h>
#include "Standard.h"

#define HELIUM

class PC {
 public:
  static const double eConst;
  static const double mu0;
  static const double c;
  static const double eps0;
  static const double h;
  static const double hBar;
  static const double kb;
  static const double fConst;
  static const double u;
  static const double me;
  static const double mp;
  static const double md;
  static const double ab;
  static const double Ha;
  static const double Ry;
  static const double NA;
  static const double A;
  
  /*
  static const double e0;
  static const double e0Standard;
  static const double rho0D2gcc;
  static const double rho0D2gccStandard;
  static const double rho0D2kgm3;
  static const double nn0;
  static const double p0;
  */

  static const double rho0D2gccStandard;
  static const double rho0D2gcc;
  static const double e0Standard;
  static const double e0;
  static const double e0DFT;
  static const double V0Hecm3molStandard;
  static const double V0Hecm3mol;
  static const double V0Hem3mol;
  static const double V0Hem3atom;
  static const double V0HeAUatom;
  static const double V0HeAUelectron;
  static const double rho0Hegcc;
  static const double nn0;
  static const double p0;

  static const double GPaToAU;
  static const double AUToGPa;
  static const double MbarToAU;
  static const double AUToMbar;
  static const double AUToK;
  static const double KToAU;
  static const double UToAU;
  static const double AUToU;
  static const double K4ToAU;
  static const double AUToeV;
  static const double eVToAU;
  static const double timeAU;
  static const double AUToPS;
  static const double AUToFS;
  static const double sToAU;
  static const double AUToS;
  static const double AToAU;
  static const double AUToA;
  static const double A3ToAU;
  static const double AUToA3;
  static const double e6cmPerSToAU;
  static const double kmPerSToAU;
  static const double AUTokmPerS;
  static const double AUToe6cmPerS;
  static const double mdAU;
  static const double AUToD2gcc;
  static const double D2gccToAU;
  static const double AUTogcc;
  static const double gccToAU;
  static const double rho0Algcc;
  static const double rho0AlAU;
  static const double rho0Quartzgcc;
  static const double rho0QuartzAU;
  static const double kcalMolToAU;
  static const double AUTokcalMol;
};

inline double DensityFromRs(double rs) {
  double v=4.0/3.0*pi*rs*rs*rs;
  return 1.0/v;
}

inline double RsFromDensity(double nn) {
  return pow(3.0/(4.0*pi*nn),1.0/3.0);
}

// volume for 1 mol of H2 molecules or He atoms per cm^3, 
// same function for H and He
inline double VolumeCCMolFromDensity(double nn) {
  return 2e+6*PC::NA*PC::ab*PC::ab*PC::ab/nn;
}

// number density of electrons given the volume of 1 mol of H2 molecules or He atoms per cm^3
// same function for H and He
inline double DensityFromVolumeCCMol(double vol) {
  return 2e+6*PC::NA*PC::ab*PC::ab*PC::ab/vol;
}

// same function for DEUTERIUM (not hydrogen) and HELIUM
inline double MassDensityFromElectronDensity(double nn) {
  return nn*PC::AUToD2gcc;
}

inline double HugoniotFunction(const double e2,
			       const double p2,
			       const double nn2,
			       const double e1,
			       const double p1,
			       const double nn1) {
  return (e2-e1)+0.5*(p2+p1)*(1.0/nn2-1.0/nn1);
}

inline double HugoniotFunction(const double e,
			       const double p,
			       const double nn) {
  
  return HugoniotFunction(e,p,nn,PC::e0,PC::p0,PC::nn0);
}

#endif // PHYSICS
