/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Physical Constants                                                      // 
//                                                                         //
// Burkhard Militzer                                  Livermore 4-27-01    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "Standard.h"
#include "Physics.h"

const double PC::eConst    = 1.602176462e-19;
const double PC::mu0       = 4.0*pi*1e-7;
const double PC::c         = 299792458.0;
const double PC::eps0      = 1.0/(PC::mu0*PC::c*PC::c);
const double PC::h         = 6.62606876e-34;
const double PC::hBar      = PC::h/(2.0*pi);
const double PC::kb        = 1.3806503e-23;
const double PC::fConst    = PC::eConst*PC::eConst/(4.0*pi*PC::eps0);
const double PC::u         = 1.66053873e-27; // atomic mass unit in kg
const double PC::me        = 9.10938188e-31;
const double PC::mp        = 1.67262158e-27;
const double PC::md        = 1.99900750083*PC::mp;
const double PC::ab        = PC::hBar*PC::hBar/(PC::fConst*PC::me);
const double PC::Ha        = PC::fConst*PC::fConst*PC::me/PC::hBar/PC::hBar;
const double PC::Ry        = PC::Ha/2.0;
const double PC::NA        = 6.0221367e23;
const double PC::A         = 1e-10; // Angstroem

const double PC::GPaToAU  = 1e9/(PC::Ha/PC::ab/PC::ab/PC::ab);
const double PC::AUToGPa  = 1.0/PC::GPaToAU;
const double PC::MbarToAU = 100.0*PC::GPaToAU;
const double PC::AUToMbar = 1.0/PC::MbarToAU;
const double PC::AUToK    = PC::Ha/PC::kb;
const double PC::KToAU    = PC::kb/PC::Ha;
const double PC::UToAU    = PC::u/PC::me;
const double PC::K4ToAU   = 10000.0*PC::kb/PC::Ha;
const double PC::AUToeV   = PC::Ha/PC::eConst;
const double PC::eVToAU   = 1.0/PC::AUToeV;
const double PC::timeAU   = PC::hBar/PC::Ha;
const double PC::sToAU    = PC::Ha/PC::hBar;               // == 1/tAU(sec)
const double PC::AUToS    = PC::timeAU;
const double PC::AUToPS   = PC::timeAU*1e12; 
const double PC::AUToFS   = PC::timeAU*1e15; 
const double PC::AToAU    = PC::A/PC::ab;
const double PC::AUToA    = PC::ab/PC::A;
const double PC::A3ToAU   = PC::AToAU*PC::AToAU*PC::AToAU;
const double PC::AUToA3   = PC::AUToA*PC::AUToA*PC::AUToA;
const double PC::e6cmPerSToAU = 1e6*0.01/PC::ab/PC::sToAU; // 10^4*tAU/lAU
const double PC::kmPerSToAU = 1e3/PC::ab/PC::sToAU;        // 10^3*tAU/lAU
const double PC::AUTokmPerS = 1.0/PC::kmPerSToAU;
const double PC::AUToe6cmPerS = 1.0/PC::e6cmPerSToAU;
const double PC::AUToD2gcc  = PC::md/(PC::ab*PC::ab*PC::ab)/1000.0;
const double PC::D2gccToAU  = 1.0/PC::AUToD2gcc;
const double PC::AUTogcc  = PC::me/(PC::ab*PC::ab*PC::ab)/1000.0;
const double PC::gccToAU  = 1.0/PC::AUTogcc;
const double PC::kcalMolToAU = 4.184*1e3/PC::NA/PC::Ha;
const double PC::AUTokcalMol = 1.0/PC::kcalMolToAU;

/////////////// materials constants + initial D2 conditions /////////////////////

const double PC::mdAU     = PC::md/PC::me;
#ifndef HELIUM
//const double PC::rho0Algcc = 2.710;  // Nellis paper Tab. III
const double PC::rho0Algcc = 2.744;  // Nellis paper Tab. III
#else
const double PC::rho0Algcc = 2.733;  // Nellis 19834 PRL Tab. II
#endif
const double PC::rho0AlAU = PC::rho0Algcc * PC::gccToAU;

const double PC::rho0Quartzgcc = 2.650;  // Hicks' viewgraphs
const double PC::rho0QuartzAU = PC::rho0Quartzgcc * PC::gccToAU;

#ifndef HELIUM
/////////////////////////////////// initial D2 conditions /////////////////////

const double PC::e0Standard= -0.5*1.1676;      // Ha per atom
//const double PC::e0      = -0.5*1.1676;      // Ha per atom
const double PC::e0        = -15.87*PC::eVToAU;// Stanimir Bonve's E0 (Ha per atom)

const double PC::rho0D2gccStandard =0.171;     // rho0 D2 g/cc
const double PC::rho0D2gcc =0.171;             // rho0 D2 g/cc
//const double PC::rho0D2gcc =0.37;            // rho0 D2 g/cc
const double PC::rho0D2kgm3=PC::rho0D2gcc * 1000.0;        // rho0 D2 kg/m^3

const double PC::nn0       =PC::rho0D2gcc * PC::D2gccToAU; // number of D atoms in ab^3
const double PC::p0        =0.0;

#else
////////////////////////////////He conditions ///////////////////////////////

const double PC::rho0D2gccStandard =0.171;     // rho0 D2 g/cc
const double PC::rho0D2gcc =0.0;

const double PC::e0Standard= -0.5*2.90337;      // Ha per electron
const double PC::e0DFT     = -1.4336452664335;  // Ha per electron
const double PC::e0        = PC::e0DFT;

const double PC::V0Hecm3molStandard =32.4;     // rho0 He cm3/mol
const double PC::V0Hecm3mol         =32.1;     // for double shock state Nellis 1984 PRL table II
// const double PC::V0Hecm3mol         =32.4;     // rho0 He cm3/mol
// const double PC::V0Hecm3mol         =18.0;
// const double PC::V0Hecm3mol         =11.75; // V at 1.2 kBar = 0.12 GPa

const double PC::V0Hem3mol          =PC::V0Hecm3mol * 1e-6;
const double PC::V0Hem3atom         =PC::V0Hem3mol / PC::NA;
const double PC::V0HeAUatom         =PC::V0Hem3atom / (PC::ab*PC::ab*PC::ab);
const double PC::V0HeAUelectron     =0.5*PC::V0HeAUatom;
const double PC::rho0Hegcc          =4.0*PC::mp/(PC::V0Hecm3mol/PC::NA)*1e3;
const double PC::nn0       =1.0 / PC::V0HeAUelectron;
const double PC::p0        =0.0;

#endif

