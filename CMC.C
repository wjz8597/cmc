/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Classical MC+MD program                                             // 
//                                                                     //
// Burkhard Militzer                              Berkeley, 07-01-2009 //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#include "Physics.h"
#include "Sampling.h"
#include "Timer.h"
#include "Form.h"
#include "ReadInTable.h"
#include "State.h"
#include "PotentialTableSpline.h"
#include "PrintDiffusion.h"
#include "ParseCommandLineArguments.h"
#include "WriteVASPPOSCARFile.h"

//////////////////////////////////////////////////////////////////////////////////////////////////

// make a true copy reference of state, do not modify the original
void EstimateMCStepSize(const State & state, const double PTarget, 
			Array1 <double> & ddType, double & dVBest, double & dCBest) {
  double rs    = state.GetRs();
  double d1    = rs*rs*rs/1000.0;
  double d2    = Point::MinimumImageDistance(); // rs >= 4.0
  double dV1   = 0.0001;
  double dV2   = 0.20;
  double dC1   = 0.001;
  double dC2   = 1.0;
  int n = 6;
  //  int n = 20;

  Array1 <double> d; // all step sizes (same for all types)
  Array1 < Array1 <double> > acc; // acc[distances][types]
  Array1 <double> dV; 
  Array1 <double> accVol;
  Array1 <double> dC; 
  Array1 <double> accCell; // acc[distances] - average over all spacial directions

  //  State state1=state; // do not restart from same conf. in each test
  for(int i=0; i<=n; i++) {
    cout << endl;
    double di  = d1 *exp(log(d2/d1)  *double(i)/double(n));
    double dVi = dV1*exp(log(dV2/dV1)*double(i)/double(n));
    double dCi = dC1*exp(log(dC2/dC1)*double(i)/double(n));
    d.PushBack(di);
    dV.PushBack(dVi);
    dC.PushBack(dCi);
    State state1=state; // restart from same conf each time
    Sampling sampling(state1,"","","","",1000);   // nothing written --> just compute TD averages
    sampling.SetPressure(PTarget);
    //    sampling.nSteps =20; warning("reduced nSteps");
    sampling.nSteps =200; // production
    sampling.nBlocks=1;
    sampling.SetDistance(di);
    sampling.SetRelativeVolumeStepSize(dVi);
    sampling.SetLatticeVectorStepSize(dCi);
    Write2(di,di/rs);
    sampling.RunAllBlocks(false,true); // false b/c no step size correction, true b/c this is a test block
    acc.PushBack(    sampling.AverageAcceptanceRatios());
    accVol.PushBack( sampling.AverageVolumeAcceptanceRatio());
    accCell.PushBack(sampling.AverageCellVectorAcceptanceRatio());
    //    Write2(di,acci);
  }
  
  for(int i=0; i<n; i++) {
    Write6(i,d[i],acc[i],dV[i],accVol[i],accCell[i]);
  }
  //  Write(d);
  //  Write(acc);
  
  ddType.Clear(); // step size estimate for each type separately
  for(int it=0; it<state.nTypes; it++) {
    if (acc[0][it]<0.1) error("Acceptance ratio at smallest step size is too small - check step size and state of system.",it,d[0],acc[0][it]);
    
    double dd=0.0;
    if (acc[0][it]<0.5) {
      cout << "Acceptance ratio at smallest step size is already less than 50%, a bit strange but could be ok." << endl;
      dd = d[0];
    }
    if (acc[n][it]>0.5) {
      cout << "Acceptance ratio at largest step size is higher than 50%. This is expected for very high T." << endl;
      dd = d[n];
    }
    for(int i=0; i<n; i++) {
      if (acc[i+1][it]<0.5) {
	double y2 = acc[i+1][it]-0.5;
	double y1 = acc[i][it]  -0.5;
	double x2 = log(d[i+1]);
	double x1 = log(d[i]  );
	double xx = (y2*x1-y1*x2)/(y2-y1);
	dd = exp(xx);
	//      Write6(x1,x2,y1,y2,xx,dd);
	break;
      }
    }
    Write4(it,dd,rs,dd/rs);
    ddType.PushBack(dd);
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  dVBest=0.0;
  if (accVol.Max()>0.0) {
 
    if (accVol[0]<0.1) error("Volume acceptance ratio at smallest step size is too small - check step size and state of system.",dV[0],accVol[0]);
    
    if (accVol[0]<0.5) {
      cout << "Volume acceptance ratio at smallest step size is already less than 50%, a bit strange but could be ok." << endl;
      dVBest = dV[0];
    }
    
    if (accVol[n]>0.5) {
      cout << "Volume acceptance ratio at largest step size is higher than 50%. This is expected for very high T." << endl;
      dVBest = dV[n];
    }
    
    for(int i=0; i<n; i++) {
      if (accVol[i+1]<0.5) {
	double y2 = accVol[i+1]-0.5;
	double y1 = accVol[i]  -0.5;
	double x2 = log(dV[i+1]);
	double x1 = log(dV[i]  );
	double xx = (y2*x1-y1*x2)/(y2-y1);
	dVBest = exp(xx);
	//      Write6(x1,x2,y1,y2,xx,dd);
	break;
      }
    }
  }
  Write(dVBest);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  dCBest=0.0;
  if (accCell.Max()>0.0) {
  
    if (accCell[0]<0.1) error("Cell vector changing move acceptance ratio at smallest step size is too small - check step size and state of system.",dC[0],accCell[0]);
    
    if (accCell[0]<0.5) {
      cout << "Cell move acceptance ratio at smallest step size is already less than 50%, a bit strange but could be ok." << endl;
      dCBest = dC[0];
    }
    
    if (accCell[n]>0.5) {
      cout << "Cell move acceptance ratio at largest step size is higher than 50%. This is expected for very high T." << endl;
      dCBest = dC[n];
    }
    
    for(int i=0; i<n; i++) {
      if (accCell[i+1]<0.5) {
	double y2 = accCell[i+1]-0.5;
	double y1 = accCell[i]  -0.5;
	double x2 = log(dC[i+1]);
	double x1 = log(dC[i]  );
	double xx = (y2*x1-y1*x2)/(y2-y1);
	dCBest = exp(xx);
	//      Write6(x1,x2,y1,y2,xx,dd);
	break;
      }
    }
  }
  Write(dCBest);
  
}

void PrintInfo(const State & state, const Sampling & sampling, const int blockMin) {
  cout << "cmc ";
  //  for(int t=0; t<state.GetNTypes(); t++)  cout << fix(2,state.listOfTypes[t]) << " " << fix(3,state.atomsByType[t].Size())<< " ";
  cout << " rho= " << fix(8,3,state.MassDensityGCC());
  cout << " T[K]= " << fix(8,0,state.GetTemperature()*PC::AUToK);
  //  cout << " m= " << fix(6,0,sampling.nSteps*sampling.nBlocks);
  //  cout << " " << fix(6,0,sampling.nSteps*(sampling.nBlocks-blockMin));
}

// use a reference of state, do not make copy, modify and use later in MD 
double RunMC(State & state, const Array1 <double> & d, const double dV, const double dC, 
	     const int nBlocks, const double PTarget, const bool stepSizeCorrectionFlag,
	     const string POSCARFileName="", const string XBSFileName="") {
  State state0 = state;
  state.SetBoxSize(); // for step size estimation procedure to work, we now copy the box size inside 'state' to Point
  state.Init();  // this bos zie change, require to recalculate other variables
  
  //  state.PrintPotentials();
  //  error("quit");

  Timer timerTotal(Timer::total,true);
  Timer timerUser(Timer::total,true);

  //int blockMin=5;
    int blockMin=50;
 //   int blockMin=100;
   // Sampling sampling(state,"mc","","","",blockMin); // write mc.xyz file contain a snapshot from the end of each block
    Sampling sampling(state,"","POSCAR","","",blockMin); // one POSCAR files per block written --> run DFT on classical snapshots
  //  Sampling sampling(state,"","","",blockMin);   // nothing written --> just compute TD averages
  //Sampling sampling(state,"","","",XBSFileName,blockMin);   // nothing written --> just compute TD averages
  //  Sampling sampling(state,"","","TRAJECTORYn1",blockMin); // save pos+forces fo FM_CODE --> test FM code
  sampling.SetPressure(PTarget);
  sampling.SetRelativeVolumeStepSize(dV);
  sampling.SetLatticeVectorStepSize(dC);
  sampling.Check();

  cout << endl;
  cout << "Ideal kinetic pressure       (GPa): " << sampling.P0*PC::AUToGPa << endl;
  cout << "Ideal kinetic energy per cell (Ha): " << sampling.E0             << endl; // now per cell, no longer per particle

  //  sampling.nSteps = 1;
  //  sampling.nSteps = 50;
  //  sampling.nSteps = 200; // production
  //  sampling.nSteps = 500; 
  sampling.nSteps = 1000; // warning("extra long blocks"); // extra long 

  //  sampling.nBlocks=10*runLongerFactor; warning("Low #blocks in MC");
  //  sampling.nBlocks=20*runLongerFactor; warning("Low #blocks in MC");
  //  sampling.nBlocks=50*runLongerFactor; 
  //  sampling.nBlocks=50;
  sampling.nBlocks=nBlocks;
  //  sampling.nBlocks=200*runLongerFactor;  // standard
  //  sampling.nBlocks=1000*runLongerFactor; 
  //  sampling.nBlocks=5000*runLongerFactor; 

  sampling.SetDistances(d);
  cout << endl;
  Write3(sampling.nBlocks,sampling.nSteps,blockMin);
  Write2(d,dV);

  sampling.RunAllBlocks(stepSizeCorrectionFlag,false); // MC
  cout << endl;
  State state1 = state;
  PrintDiffusionAllTypes(state0,state1,0.0,true); // true because it is after the last block
  cout << endl;

  double tMC = timerUser.ReadAndRestart();
  cout << "User time:  " << fix(2,2,tMC) << " sec = " << fix(2,2,tMC/60.0) << " minutes = " << fix(2,2,tMC/3600.0) << " hours." << endl; // not yet tested
  tMC = timerTotal.ReadAndRestart();
  cout << "Total time: " << fix(2,2,tMC) << " sec = " << fix(2,2,tMC/60.0) << " minutes = " << fix(2,2,tMC/3600.0) << " hours." << endl; // not yet tested
  cout << endl;
  cout << " --------------------------------------------------------------- " << endl;
  cout << endl;
  double ePotAvParticle = sampling.analysis.totalEnergy.GetAverage(); // is per particle (H2 or He) 
  //  Write(ePotAvParticle);

  double P  = sampling.analysis.totalPressure.GetAverage();
  double PE = sampling.analysis.totalPressure.GetErrorBar();
  PrintInfo(state,sampling,blockMin);
  if (P>=0.01 && PE>=0.0003) cout << "   P[AU]=  " << fix(11,4,P) << " " << fix(7,4,PE) << endl;
  else cout << "   P[AU]=  " << sci(11,4,P) << " " << sci(7,1,PE) << endl;

  double PGPa  = P *PC::AUToGPa;
  double PEGPa = PE*PC::AUToGPa;
  PrintInfo(state,sampling,blockMin);
  if (PGPa>=0.01 && PEGPa>=0.0003) cout << "   P[GPa]= " << fix(11,4,PGPa) << " " << fix(7,4,PEGPa) << endl;
  else cout << "   P[GPa]= " << sci(11,4,PGPa) << " " << sci(7,1,PEGPa) << endl;

  double PV  = P *Point::GetVolume();
  double PVE = PE*Point::GetVolume();
  PrintInfo(state,sampling,blockMin);
  if (PVE>=0.000003) cout << "   PV[Ha/cell]=     " << fix(10,6,PV) << " " << fix(9,6,PVE) << endl;
  else cout << "   PV[Ha/cell]=     " << sci(10,6,PV) << " " << sci(9,3,PVE) << endl;

  PrintInfo(state,sampling,blockMin);
  // remember sampling.analysis.totalEnergy.GetAverage() is the POTENTIAL ENERGY PER CELL (H or He atom) 
  // ['total' means it is not block energy]
  double E = (sampling.analysis.totalEnergy.GetAverage()+sampling.analysis.K); // +E0Atom;
  double EE = sampling.analysis.totalEnergy.GetErrorBar(); 
  if (EE>=0.000003) cout << "  (E=V+K)[Ha/cell]= " << fix(10,6,E) << " " << fix( 9,6,EE) << endl;
  else cout << "  (E=V+K)[Ha/cell]= " << fix(10,6,E) << " " << sci( 9,3,EE) << endl;

  cout << endl;
  //if (POSCARFileName.length()>0) WriteVASPPOSCARFile("POSCAR_MC_"+POSCARFileName,state);
  if (POSCARFileName.length()>0) WriteVASPPOSCARFile(POSCARFileName,state);
  return ePotAvParticle;
}

int main(int argc, char *argv[]) {
  cout << endl;
  cout << "Burkhard Militzer's classical MC+MD program" << endl;
  cout << "It uses pair potentials obtained with force matching along DFT-MD trajectories." << endl;
  cout << endl;
  cout << "Oxford, 08-18-2018" << endl;
  cout << endl;

  Array1 <string> arg = CopyCommandLineArgumentsToArray(argc,argv);
  bool quitFlag       = ParseCommandLineFlag(arg,"-q");
  bool leakErrorFlag  = !(ParseCommandLineFlag(arg,"-noLeakError"));
  bool doNotMoveFlag  = ParseCommandLineFlag(arg,"-doNotMove");
  double boxScale = 1.0;
  bool boxScaleFlag   = ParseCommandLineArgumentNoSpace(arg,"boxScale=",boxScale);
  double potScale = 1.0;
  ParseCommandLineArgumentNoSpace(arg,"pairPotentialFactor=",potScale);

  string WriteFinalPOSCARFileName;
  ParseCommandLineArgumentNoSpace(arg,"POSCAR=",WriteFinalPOSCARFileName);
  string XBSFileName;
  ParseCommandLineArgumentNoSpace(arg,"XBS=",XBSFileName);

  if (arg.Size()<4) { 
    cout << "Usage: " << argv[0] << " POSCAR type1 [type2 ...] T=temperature_in_K P=pressure_in_GPa QFactor=f QStar=f FP=f [-q] [boxScale=f] [POSCAR=POSCARFileName] [pairPotentialFactor=f] [-noLeakError] [-doNotMove] [XBS=XBSFileName" << endl;
    cout << "Usage: " << argv[0] << " POSCAR Fe      T=20000 P=3974.11138518 " << endl;
    cout << "Usage: " << argv[0] << " POSCAR Fe Mg O T=20000 P=3974.11138518 -q" << endl;
    Quit("Please correct the command line parameters.");
  }

  double QFactor=0.0;
  ParseCommandLineArgumentNoSpace(arg,"QFactor=",QFactor);
  Write(QFactor);

  double QStar=0.0;
  ParseCommandLineArgumentNoSpace(arg,"QStar=",QStar);
  Write(QStar);

  double fourierPotentialFactor=0.0;
  ParseCommandLineArgumentNoSpace(arg,"FP=",fourierPotentialFactor);
  Write(fourierPotentialFactor);

  int nBlocks = 50;
  ParseCommandLineArgumentNoSpace(arg,"nb=",nBlocks);
  Write(nBlocks);
  int nBlocksRelax = 0;
  ParseCommandLineArgumentNoSpace(arg,"nbr=",nBlocksRelax);
  Write(nBlocksRelax);

  double PGPa;
  bool PressureFlag = ParseCommandLineArgumentNoSpace(arg,"P=",PGPa);
  if (PressureFlag==false) error("Must specify pressure in GPa in command line.");
  double P = PGPa * PC::GPaToAU;
  Write2(PGPa,P);

  double TinK;
  bool temperatureFlag = ParseCommandLineArgumentNoSpace(arg,"T=",TinK);
  if (temperatureFlag==false) error("Must specify temperature in Kelvin in command line with T=.");
  Write(TinK);

  //  cout.precision(12);
  //  Write(PC::AUToA);
  //  Write(PC::AUToeV);
  //  Write(PC::AUToK);
  //  Write(PC::AUToGPa);
  //  Write(PC::hBar);
  //  Write(PC::fConst);
  //  cout << endl;
  //  Quit("Q");
  cout.precision(8);

  string ReadPOSCARfile   = arg[1];
  Array1 <string> types;
  for(int i=2; i<arg.Size();) {
    types.PushBack(arg[i++]); // is element name
  }

  //  InitRandom(true);
  cout << "Running with fixed default random number seed." << endl;
  Write(Random());

  State state;
  //  state.ParseVaspPOSCARFile(ReadPOSCARfile,types,true); // sets the box/cell size in Point and inside State !!!
  state.ParseVaspPOSCARFile(ReadPOSCARfile,types,false); // sets the box/cell size in Point and inside State !!!
  
  if (boxScaleFlag) {
    state.ScaleBox(boxScale);
    state.SetBoxSize();
  }
  //  state.CopyCurrentBoxSizeFromPoint(); // unclear what purpose this once served
  if (doNotMoveFlag) state.doNotMove = true;
  
    //  state.PrintNeighbors();
  state.PrintNeighbors(0.75);

  //  int t1=0;
  //  int t2=1;
  //  int t3=2;
  int nT = state.GetNTypes();

  double T     = TinK*PC::KToAU;
  double beta  = 1.0/T;
  double V     = Point::GetVolume();
  double N     = state.nParticles;
  //  double nn    = double(N)/V;                 // number density N/V
  double V1    = V/double(N);                 // volume for 1 H or He atom in a.u.
  //  double VSImol= V1*PC::NA*PC::ab*PC::ab*PC::ab;
  //  double Vccmol= VSImol*1000000.0;
  double rs    = pow(V1*3.0/(4.0*pi),1.0/3.0); // not the true 'rs' excpet for H

  state.SetBeta(beta);

  cout << endl;
  Write2(rs,V1);
  Write2(N,V);
  cout << endl;
  Write2(TinK,T);
  cout << endl;

  state.TotalFreeEnergy(true);
  if (quitFlag) Quit("Stop after printing F and before reading any potential files");
  //  Quit("Stop after printing F and before reading any potential files");

  Array1 <PotentialTableSpline> pTable(nT*(nT+1)/2);
  int ip=0;
  for(int it=0; it<nT; it++) {
    for(int jt=0; jt<=it; jt++) {
      string file = "Poten_"+IntToString(jt+1)+IntToString(it+1)+".dat";
      pTable[ip].Init(file,1,2,potScale);                      // read in pot file, store a real copy, rather than using a pointer
      state.RegisterPotential(it,jt,pTable[ip]);  
      if (pTable[ip].GetRMax()*(1.0-1e-08)>Point::MinimumImageDistance()*0.5) {
	if (leakErrorFlag) error("Potential leaks out:",it,jt,pTable[ip].GetRMax(),Point::MinimumImageDistance()*0.5);
	else warning("Potential leaks out:",it,jt,pTable[ip].GetRMax(),Point::MinimumImageDistance()*0.5);
      }
      ip++;
    }      
  }

  //  state.PrintPotentials(); Quit("Print pot");

  /*
  if (state.nParticles==384) {
    // good for N = 384 = 4*6x4x4 cell
    state.AddKPoint(12, 0, 0, 1.0);
    state.AddKPoint( 6,12, 0, 1.0);
    state.AddKPoint( 0, 0, 8, 1.0);
  } else if (state.nParticles==144) {
    // good for N = 144 = 4*4x3x3 cell
    state.AddKPoint(8,0,0, 1.0);
    //    state.AddKPoint(0,6,0, 1.0);
    state.AddKPoint(4,9,0, 1.0);
    state.AddKPoint(0,0,6, 1.0);
  } else {
    if (fourierPotentialFactor!=0.0) error("k-point grid only specified for N=144 and 384.");
  }
  */

  state.SetFourierPotentialFactor(fourierPotentialFactor);

  state.SetQFactorOverall(QFactor);
  //  state.RegisterQTerm(2,-1.5233869);
  //  state.RegisterQTerm(4,2.2231621);
  //  state.RegisterQTerm(8,-2.4415873);
  //  state.RegisterQTerm(12,2.0442374);
  //  state.RegisterQTerm(12,1.0);
  //  for(int L=2; L<=24; L+=2) {
  //    state.RegisterQTerm(L,1.0);
  //  }
  // L=  2 QTerm/(2L+1)=   0.0000 QTerm/(2L+1)/N= 0.00000000
  // L=  4 QTerm/(2L+1)=   2.6329 QTerm/(2L+1)/N= 0.00914196
  // L=  6 QTerm/(2L+1)=  10.9231 QTerm/(2L+1)/N= 0.03792728 **4298**
  // L=  8 QTerm/(2L+1)=   6.2461 QTerm/(2L+1)/N= 0.02168800
  // L= 10 QTerm/(2L+1)=   0.1803 QTerm/(2L+1)/N= 0.00062598
  // L= 12 QTerm/(2L+1)=   9.1802 QTerm/(2L+1)/N= 0.03187553
  // L, QTerm target value per N and divided by (2L+1), weight
  //  state.RegisterQTerm( 2, 3.8, 4.2, 0.00000000, 1.0);
  //  state.RegisterQTerm( 4, 3.8, 4.2, 0.00914196, 1.0);
  //  state.RegisterQTerm( 6, 3.8, 4.2, 0.03792728, 1.0);
  //  state.RegisterQTerm( 8, 3.8, 4.2, 0.02168800, 1.0);
  //  state.RegisterQTerm(10, 3.8, 4.2, 0.00062598, 1.0);
  //  state.RegisterQTerm(12, 3.8, 4.2, 0.03187553, 1.0);
  //  state.RegisterQTerm(10, 3.8, 4.2, 0.000, 1.0);

  state.RegisterQTerm(12, 0.0, 0.0, 3.8, 4.2, QStar, 1.0);
  //  state.RegisterQTerm( 8, 4.0, 4.2, 4.8, 5.0, 0.038556446, 1.0);
  //  state.RegisterQTerm( 8, 0.0, 0.0, 0.0, 5.0, QStar, 1.0);
      
  state.Init();
  cout << endl;
  cout << "Initial pair    potential energy: V[Ha/cell]= " << fix(14,8,state.GetTotalEnergy())            << endl;
  //  cout << "Initial pair potential energy: V[eV/cell]= " << fix(14,8,state.GetTotalEnergy()*PC::AUToeV) << endl; 
  cout << "Initial Fourier potential energy: V[Ha/cell]= " << fix(14,8,state.GetFourierPotentialEnergy())            << endl;
  cout << endl;

  //  state.CalculateManySOfK(12);
  //  Quit("sOfK");
  
  /*
  cout << " ----------- " << endl;
  Write(state.SumOfAllSOfKTermsSlow());
  cout << " ----------- " << endl;
  state.SetAllRhoK();
  Write(state.SumOfAllSOfKTermsFromRhoK());
  Quit("S(k)");
  */
  
  Array1 <double> d(state.nTypes,0.3);
  double dV = 0.03;
  double dC = 0.01;
  //  cout << "A"; state.Check();
  if (doNotMoveFlag==false) EstimateMCStepSize(state,P,d,dV,dC);
  // note: this state.boxSize may not be conistent with Point::boxSize b/c EstimateMCStepSize used its own temporary copy of state
  // state.Check() would fail right here.
  
  //////////////////////////////////////////////////////////////////////////////////////////////////

  //  dV = 0.0001; // overwrite estimate
  //  cout << "B"; state.Check();
  //  cout << "B*" << endl;
  RunMC(state,d,dV,dC,nBlocks,P,false,WriteFinalPOSCARFileName,XBSFileName);
  state.CopyCurrentBoxSizeFromPoint(); // copy current boxSize so that the relaxation can continue w/o problems

  //////////////////////////////////////////////////////////////////////////////////////////////////

  if (nBlocksRelax>0) {
  
    cout << endl << "Relaxing current structure" << endl << endl;
    
    double f = 0.1;
    d  *= f;
    dV *= f;
    dC *= f;
    
    state.SetTemperature(state.GetTemperature()*f*f*f);
    //  nBlocks = 50;
    RunMC(state,d,dV,dC,nBlocksRelax,P,true,WriteFinalPOSCARFileName,XBSFileName+"_relax"); // careful: now with step size correction
    
    state.CalculateManySOfK(12);
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////

  cout << endl << "CMC simulation finished properly." << endl << endl;
  return 0;
}
  
