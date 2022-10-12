/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Analysis code of classical MC and MD simulations                    //
//                                                                     //
// Burkhard Militzer                               Washington 12-31-06 //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#ifndef _ANALYSIS_
#define _ANALYSIS_

#include "Averages.h"
#include "GOfR.h"
#include "QHistogram.h"

class Analysis {
 public:
  ScalarAverage blockEnergy;
  ScalarAverage totalEnergy;
  ScalarAverage blockPairEnergy;
  ScalarAverage totalPairEnergy;
  ScalarAverage blockQEnergy;
  ScalarAverage totalQEnergy;
  const int nq;
  ArrayAverage  blockQ;
  ArrayAverage  totalQ;
  ScalarAverage blockFPEnergy;
  ScalarAverage totalFPEnergy;
  ScalarAverage blockFP;
  ScalarAverage totalFP;
  ScalarAverage blockDeltaV;
  ScalarAverage totalDeltaV;
  ScalarAverage blockPressure;
  ScalarAverage totalPressure;
  ScalarAverage blockDisp;
  ScalarAverage totalDisp;
  ScalarAverage blockVolume;
  ScalarAverage totalVolume;
  ScalarAverage blockVolumeChanges;
  ScalarAverage totalVolumeChanges;
  ArrayAverage  blockCellVectors;
  ArrayAverage  totalCellVectors;
  ArrayAverage  blockCellAngles;
  ArrayAverage  totalCellAngles;
  ScalarAverage blockCellVectorChanges;
  ScalarAverage totalCellVectorChanges;
  ArrayAverage  blockSOfK;
  ArrayAverage  totalSOfK;
  
  //  static const int nBinsGOfR = 100;
  static const int nBinsGOfR = 50;
  int nTypes,nGOfR;
  Array1 <GOfR> blockGOfR;
  Array1 <ArrayAverage> totalGOfR;

  //  static const int nBinsQHistogram = 100;
  static const int nBinsQHistogram = 1000;
  Array1 <QHistogram> blockQHistogram;
  Array1 <ArrayAverage> totalQHistogram;

  double K;
  double P0;
  const int blockMin;


  // let the main program decide to pass in K_total or K/N
  // P0 is the pressure for noninteracting system
 Analysis(const int dim, const int nk, const int nq_, const int nTypes_, const int ib, const double dMax, const State & state, const double K_, const double P0_):
  nq(nq_),
  blockQ(nq),
  totalQ(nq),
  blockCellVectors(dim),
  totalCellVectors(dim),
  blockCellAngles(dim),
  totalCellAngles(dim),
  blockSOfK(nk),
  totalSOfK(nk),
  nTypes(nTypes_),
  nGOfR(nTypes*(nTypes+1)/2),
  blockGOfR(nGOfR),
  totalGOfR(nGOfR),
  blockQHistogram(nq),
  totalQHistogram(nq),
  K(K_),
  P0(P0_),
  blockMin(ib) {
    //    Write(dMax);
    int t1=0;
    int t2=0;
    for(int i=0; i<nGOfR; i++) {
      //      Write6(nTypes,nGOfR,i,t1,t2,GOfRIndex(t1,t2));
      if (GOfRIndex(t1,t2)!=i) error("g(r) index mismatch",i,t1,t2,GOfRIndex(t1,t2));
      string name = "gr_" + state.listOfTypes[t1] + "-" + state.listOfTypes[t2];
      blockGOfR[i].Init(name,nBinsGOfR,dMax);
      totalGOfR[i].Init(name,nBinsGOfR);
      t2++;
      if (t2==nTypes) {
	t1++;
	t2=t1;
      }
    }
    for(int iq=0; iq<nq; iq++) {
      double QMin = 0.00; // divided by particle and by (2*L+1)
      double QMax = 0.04; // divided by particle and by (2*L+1)
      //       double QMin = state.QTargets[iq]-state.QInterval/2.0;
      //       double QMax = state.QTargets[iq]+state.QInterval/2.0;
      blockQHistogram[iq].Init("block_hist_Q"+DoubleToString(state.L[iq]),nBinsQHistogram,QMin,QMax); 
      totalQHistogram[iq].Init("total_hist_Q"+DoubleToString(state.L[iq]),nBinsQHistogram);
    }
    blockEnergy.SetLong();
    totalEnergy.SetLong();
    blockPairEnergy.SetLong();
    totalPairEnergy.SetLong();
    blockQEnergy.SetLong();
    totalQEnergy.SetLong();
    blockQ.SetLong();
    totalQ.SetLong();
    blockFPEnergy.SetLong();
    totalFPEnergy.SetLong();
    blockFP.SetLong();
    totalFP.SetLong();
    blockDeltaV.SetLong();
    totalDeltaV.SetLong();
    blockPressure.SetLong();
    totalPressure.SetLong();
    blockDisp.SetLong();
    totalDisp.SetLong();
  }

  void StartBlock(const int ib) {
    blockEnergy.SetToZero();
    blockPairEnergy.SetToZero();
    blockQEnergy.SetToZero();
    blockQ.SetToZero();
    blockFPEnergy.SetToZero();
    blockFP.SetToZero();
    blockDeltaV.SetToZero();
    blockPressure.SetToZero();
    blockDisp.SetToZero();
    blockVolume.SetToZero();
    blockVolumeChanges.SetToZero();
    blockCellVectorChanges.SetToZero();
    blockCellVectors.SetToZero();
    blockCellAngles.SetToZero();
    for(int i=0; i<nGOfR; i++) {
      blockGOfR[i].SetToZero();
    }
    for(int iq=0; iq<nq; iq++)
      blockQHistogram[iq].SetToZero();
    blockSOfK.SetToZero();
    
  }
  void EndBlock(const int ib) {
    for(int iq=0; iq<nq; iq++) {
      blockQHistogram[iq].Normalize();
    }
    if (ib>=blockMin) {
      totalEnergy.Add(blockEnergy       .GetAverage());
      totalPairEnergy.Add(blockPairEnergy   .GetAverage());
      totalQEnergy.Add(blockQEnergy     .GetAverage());
      totalQ.Add(blockQ           .GetAverage());
      totalFPEnergy.Add(blockFPEnergy     .GetAverage());
      totalFP.Add(blockFP           .GetAverage());
      totalDeltaV.Add(blockDeltaV       .GetAverage());
      totalPressure.Add(blockPressure     .GetAverage());
      totalDisp.Add(blockDisp         .GetAverage());
      totalVolume.Add(blockVolume       .GetAverage());
      totalVolumeChanges.Add(blockVolumeChanges.GetAverage());
      totalCellVectorChanges.Add(blockCellVectorChanges.GetAverage());
      totalCellVectors.Add(blockCellVectors.GetAverage());
      totalCellAngles.Add(blockCellAngles.GetAverage());
      for(int i=0; i<nGOfR; i++) {
	blockGOfR[i].Normalize(Point::GetVolume());
	totalGOfR[i].Add(blockGOfR[i].GetArray());
      }
      for(int iq=0; iq<nq; iq++) {
	totalQHistogram[iq].Add(blockQHistogram[iq].GetArray());
      }
      totalSOfK.Add(blockSOfK.GetAverage());
    }
    cout << *this;
    PrintQHistograms();
  }
  void EndSimulation() {
    PrintGOfRs();
    blockEnergy.SetLong();
    totalEnergy.SetLong();
    blockPairEnergy.SetLong();
    totalPairEnergy.SetLong();
    blockQEnergy.SetLong();
    totalQEnergy.SetLong();
    blockQ.SetLong();
    totalQ.SetLong();
    blockFPEnergy.SetLong();
    totalFPEnergy.SetLong();
    blockFP.SetLong();
    totalFP.SetLong();
    blockDeltaV.SetLong();
    totalDeltaV.SetLong();
    blockPressure.SetLong();
    totalPressure.SetLong();
    blockDisp.SetLong();
    totalDisp.SetLong();
    blockVolume.SetLong();
    totalVolume.SetLong();
    blockVolumeChanges.SetLong();
    totalVolumeChanges.SetLong();
    blockCellVectorChanges.SetLong();
    totalCellVectorChanges.SetLong();
    blockCellVectors.SetLong();
    totalCellVectors.SetLong();
    blockCellAngles.SetLong();
    totalCellAngles.SetLong();
    blockSOfK.SetLong();
    totalSOfK.SetLong();
    cout << *this;
  }

  void AddTotalEnergy(const double e) {
    blockEnergy.Add(e);
  }
  void AddPairEnergy(const double e) {
    blockPairEnergy.Add(e);
  }
  void AddQEnergy(const double e) { // QEnergy=QTerm*QFactor
    blockQEnergy.Add(e);
  }
  void AddQTerms(const Array1 <double> & Q, const Array1 <double> & QForHistogram) { // QTerm=QEnergy/QFactor
    blockQ.Add(Q);
    for(int iq=0; iq<nq; iq++) {
      blockQHistogram[iq].Add(QForHistogram[iq]);
    }
  }
  void AddFPEnergy(const double e) { // QEnergy=QTerm*QFactor
    blockFPEnergy.Add(e);
  }
  void AddFPTerm(const double t) { // QTerm=QEnergy/QFactor
    blockFP.Add(t);
  }
  void AddPressure(const double p) {
    blockPressure.Add(p);
  }
  void AddDeltaV(const double dV) {
    blockDeltaV.Add(dV);
  }
  void AddDisplacement(const double d) {
    blockDisp.Add(d);
  }
  void AddVolume(const double V) {
    blockVolume.Add(V);
  }
  void AddVolumeChange(const double dV) {
    blockVolumeChanges.Add(dV);
  }
  void AddCellVectorChange(const double dc) {
    blockCellVectorChanges.Add(dc);
  }
  void AddCellVectors(const Array1 <double> & v) {
    blockCellVectors.Add(v);
  }
  void AddCellAngles(const Array1 <double> & v) {
    blockCellAngles.Add(v);
  }
  int GOfRIndex(int t1, int t2) const {
    // needs t2 >= t1 
    if (t2<t1) { int t=t1; t1=t2; t2=t; }
    return nTypes*(nTypes+1)/2-(nTypes-t1)*(nTypes-t1+1)/2 + t2-t1;
  }
  void AddDistanceToGOfR(const int t1, const int t2, const double r) {
    int i = GOfRIndex(t1,t2);
    blockGOfR[i].Add(r);
  }
  void AddSOfKTerm(const Array1 <double> & sOfK) {
    blockSOfK.Add(sOfK);
  }

  void PrintGOfRs() const {
    //    Form f(sci,12,4);
    Form f(fix,12,6);
    f.center();
    /*
    cout << "g(r) block averages" << endl;
    cout << f("r=");
    for(int g=0; g<nGOfR; g++) 
      cout << f(blockGOfR[g].name);
    cout << endl;

    for(int i=0; i<blockGOfR[0].Size(); i++) {
      cout << f(blockGOfR[0].GetMidPoint(i));
      for(int g=0; g<nGOfR; g++) 
	cout << blockGOfR[g].GetAverage(i);
      cout << endl;
    }
    */
    cout << endl;
    cout << "# total averages g(r)" << endl;
    cout << "#" << f("r=");
    f.right();
    for(int g=0; g<nGOfR; g++) 
      cout << f(totalGOfR[g].name);
    cout << endl;

    int nBinMax = blockGOfR[0].Size();
    //    int nBinMax = int(double(blockGOfR[0].Size())/sqrt(3.0));
    for(int i=0; i<nBinMax; i++) {
      cout << "gr ";
      cout << f(blockGOfR[0].GetMidPoint(i));
      for(int g=0; g<nGOfR; g++) 
	cout << f(totalGOfR[g].GetAverage(i));
      cout << endl;
    }

  }

  void PrintQHistograms() const {
    //    Form f(sci,12,4);
    Form f(fix,14,6);

    cout << endl;
    cout << "# Block and total Q histograms" << endl;

    int iBinMin = 0;
    int iBinMax = nBinsQHistogram-1;
    if (nBinsQHistogram>100) {
      for(int i=0; i<nBinsQHistogram; i++) {
	bool flag = false;
	for(int iq=0; iq<nq; iq++)  
	  flag |= blockQHistogram[iq].GetValue(i)>0.0 || totalQHistogram[iq].GetAverage(i)>0.0;
	if (flag) {
	  iBinMin=i;
	  break;
	}
      }
      for(int i=nBinsQHistogram-1; i>=0; i--) {
	bool flag = false;
	for(int iq=0; iq<nq; iq++)  
	  flag |= blockQHistogram[iq].GetValue(i)>0.0 || totalQHistogram[iq].GetAverage(i)>0.0;
	if (flag) {
	  iBinMax=i;
	  break;
	}
      }
    }
    cout << "# Printing " << iBinMax-iBinMin+1 << " out of " << nBinsQHistogram << " bins." << endl;
    cout << "# First bin: " << iBinMin+1 << "    Last bin: " << iBinMax+1 << endl;

    f.right();
    cout << "#    i=  " << f("Q=   ");
    for(int iq=0; iq<nq; iq++) {
      cout << " " << f(blockQHistogram[iq].name);
      cout << " " << f(totalQHistogram[iq].name);
    }
    cout << endl;

    for(int i=iBinMin; i<=iBinMax; i++) {
      cout << "PQ: ";
      cout << IntToStringMaxNumber(i+1,nBinsQHistogram);
      cout << " " << f(blockQHistogram[0].GetMidPoint(i));
      for(int iq=0; iq<nq; iq++) { 
	cout << " " << f(blockQHistogram[iq].GetValue(i)); // assume the histogram has been normalized beforehand
      	cout << " " << f(totalQHistogram[iq].GetAverage(i));
      }
      cout << endl;
    }
    cout << endl;

  }

  static string PFS(const ScalarAverage & s, const Form & f, const Form & fe) {
    ostringstream os;
    os << f(s.GetAverage()) << " +- " << fe(s.GetErrorBar()) << " ";
    return os.str();
  }

  static string PFS(const ScalarAverage & s, const Form & f) {
    return PFS(s,f,f);
  }

  static string PFS(const ArrayAverage & a, const int i, const Form & f, const Form & fe) {
    ostringstream os;
    os << f(a.GetAverage(i)) << " +- " << fe(a.GetErrorBar(i)) << " ";
    return os.str();
  }

  /*
  friend ostream& operator<<(ostream &os, const Analysis & a ) {
    Form ff(fix,14,8);
    Form fe(fix,11,8);
    os << endl;
    os << "DR   [AU] : block: " << PFS(a.blockDisp      ,ff,fe) << " total: " << PFS(a.totalDisp      ,ff,fe) << endl;
    os << "V/N  [AU] : block: " << PFS(a.blockEnergy    ,ff,fe) << " total: " << PFS(a.totalEnergy    ,ff,fe) << endl;
    os << "E/N  [AU] : block: " << PFS(a.blockEnergy+a.K,ff,fe) << " total: " << PFS(a.totalEnergy+a.K,ff,fe) << endl;
    os << "P    [AU] : block: " << PFS(a.blockPressure  ,ff,fe) << " total: " << PFS(a.totalPressure  ,ff,fe) << endl;
    os << endl;
    os << "V/N  [K]  : block: " << PFS(a.blockEnergy      *PC::AUToK,ff,fe) << " total: " << PFS(a.totalEnergy      *PC::AUToK,ff,fe) << endl;
    os << "E/N  [K]  : block: " << PFS((a.blockEnergy+a.K)*PC::AUToK,ff,fe) << " total: " << PFS((a.totalEnergy+a.K)*PC::AUToK,ff,fe) << endl;
    os << "P    [GPa]: block: " << PFS(a.blockPressure*PC::AUToGPa  ,ff,fe) << " total: " << PFS(a.totalPressure*PC::AUToGPa  ,ff,fe) << endl;
    os << endl;
    os << "(E-E0)/E0 : block: " << PFS( a.blockEnergy        /a.K ,ff,fe) << " total: " << PFS( a.totalEnergy         /a.K ,ff,fe) << endl;
    os << "(P-P0)/P0 : block: " << PFS((a.blockPressure-a.P0)/a.P0,ff,fe) << " total: " << PFS((a.totalPressure-a.P0) /a.P0,ff,fe) << endl;
    //    os << a.P0 << endl;
    //    os << a.P0*PC::AUToGPa << endl;

    return os;
  }
  */

  friend ostream& operator<<(ostream &os, const Analysis & a ) {
    Form ff(fix,14,8);
    Form fe(fix,11,8);
    Form ff2(sci,14,4);
    Form fe2(sci,11,4);
    os << endl;
    os << "DR             [AU] : block: " << PFS(a.blockDisp      ,ff,fe) << " total: " << PFS(a.totalDisp      ,ff,fe) << endl;
    //    os << "dV/cell        [AU] : block: " << PFS(a.blockDeltaV    ,ff,fe) << " total: " << PFS(a.totalDeltaV    ,ff,fe) << endl;
    os << "Pot_pair/cell  [AU] : block: " << PFS(a.blockPairEnergy,ff,fe) << " total: " << PFS(a.totalPairEnergy,ff,fe) << endl;
    os << "Q*factor/cell  [AU] : block: " << PFS(a.blockQEnergy  ,ff,fe) << " total: " << PFS(a.totalQEnergy  ,ff,fe) << endl;
    for(int iq=0; iq<a.nq; iq++) {
      os << "QTerm/cell     [AU] : block: " << PFS(a.blockQ,iq   ,ff,fe) << " total: " << PFS(a.totalQ,iq     ,ff,fe) << endl;
    }
    os << "FP*factor/cell [AU] : block: " << PFS(a.blockFPEnergy  ,ff,fe) << " total: " << PFS(a.totalFPEnergy  ,ff,fe) << endl;
    os << "FP/cell        [AU] : block: " << PFS(a.blockFP        ,ff2,fe2) << " total: " << PFS(a.totalFP        ,ff2,fe2) << endl;
    os << "Pot_total/cell [AU] : block: " << PFS(a.blockEnergy    ,ff,fe) << " total: " << PFS(a.totalEnergy    ,ff,fe) << endl;
    os << "Pot+kin/cell   [AU] : block: " << PFS(a.blockEnergy+a.K,ff,fe) << " total: " << PFS(a.totalEnergy+a.K,ff,fe) << endl;
    os << "Pressure       [AU] : block: " << PFS(a.blockPressure  ,ff,fe) << " total: " << PFS(a.totalPressure  ,ff,fe) << endl;
    //    os << endl;
    //    os << "V/cell [K]  : block: " << PFS(a.blockEnergy      *PC::AUToK,ff,fe) << " total: " << PFS(a.totalEnergy      *PC::AUToK,ff,fe) << endl;
    //    os << "E/cell [K]  : block: " << PFS((a.blockEnergy+a.K)*PC::AUToK,ff,fe) << " total: " << PFS((a.totalEnergy+a.K)*PC::AUToK,ff,fe) << endl;
    os << "Pressure      [GPa] : block: " << PFS(a.blockPressure*PC::AUToGPa ,ff,fe) << " total: " << PFS(a.totalPressure*PC::AUToGPa,ff,fe) << endl;
    os << "Volume/cell    [AU] : block: " << PFS(a.blockVolume     ,ff,fe) << " total: " << PFS(a.totalVolume  ,ff,fe) << endl;

    os << "V_changes/cell [AU] : block: " << PFS(a.blockVolumeChanges,ff,fe) << " total: " << PFS(a.totalVolumeChanges,ff,fe) << endl;
    for(int d=0; d<Point::nDim; d++) {
      os << "cell[" << d+1 << "]        [AU] : block: " << PFS(a.blockCellVectors,d,ff,fe) << " total: " << PFS(a.totalCellVectors,d,ff,fe) << endl;
    }
    for(int d=0; d<Point::nDim; d++) {
      if (a.blockCellAngles.Valid()) {
	os << "angle[" << d+1 << "]       [AU] : block: " << PFS(a.blockCellAngles,d,ff,fe) << " total: " << PFS(a.totalCellAngles,d,ff,fe) << endl;
      }
    }
    if (a.blockCellVectorChanges.Valid()) {
      os << "cell_changes   [AU] : block: " << PFS(a.blockCellVectorChanges,ff,fe) << " total: " << PFS(a.totalCellVectorChanges,ff,fe) << endl;
    }

    for(int ik=0; ik<a.blockSOfK.Size(); ik++) {
      os << "SOfK[" << ik+1 << "]             : block: " << PFS(a.blockSOfK,ik,ff2,fe2) << " total: " << PFS(a.totalSOfK,ik,ff2,fe2) << endl;
    }

    //    os << endl;
    //    os << "(E-E0)/E0    : block: " << PFS( a.blockEnergy        /a.K ,ff,fe) << " total: " << PFS( a.totalEnergy         /a.K ,ff,fe) << endl;
    //    os << "(P-P0)/P0    : block: " << PFS((a.blockPressure-a.P0)/a.P0,ff,fe) << " total: " << PFS((a.totalPressure-a.P0) /a.P0,ff,fe) << endl;
    //    os << a.P0 << endl;
    //    os << a.P0*PC::AUToGPa << endl;

    return os;
  }


};

class AnalysisMD {
 public:
  AnalysisMD(const State & state_, const double dt_, const int nSteps_, const int lHist_,const int nStepsAvMin_):
    state(state_),dt(dt_),nSteps(nSteps_),nStepsAvMin(nStepsAvMin_),
    lHist(lHist_),nt(state.GetNTypes()),nStore(0),iStore(0) {
    v.Resize(lHist,state.nParticles);
    vtvt.Resize(state.GetNTypes());
    for(int t=0; t<state.GetNTypes(); t++) {
      vtvt[t].Init(lHist);
    }
    vtvtParticle.Resize(lHist);
  }
    
  bool VelocitiesNeeded(const int i) {
    return (i%nSteps==0 && i>nStepsAvMin);
  }
  void ReportVelocities(const Array1 <BPoint> vt) {
    for(int i=0; i<state.nParticles; i++) {
      v(iStore,i) = vt[i];
    }
    nStore++;
    if (nStore>=lHist) RecordVVCorrelation();
    iStore = (iStore+1) % lHist;
  }
  void RecordVVCorrelation() {
    for(int i=0; i<state.nParticles; i++) {
      for(int j=0; j<lHist; j++) {
	int jStore = (iStore-j+lHist) % lHist;
	vtvtParticle[j] = v(iStore,i)*v(jStore,i); // / v(iStore,i).Norm2();
      }
      int t = state.typeIndex[i];
      vtvt[t].Add(vtvtParticle);
    }
  }

  void PrintVVCorrelation() const {
    cout << "#Vel-vel autocorr";
    for(int t=0; t<nt; t++) {
      cout << "         " << state.GetTypeName(t) << "     ";
    }
    cout << endl;
    for(int i=0; i<lHist; i++) {
      double tt = dt*i*nSteps;
      cout << tt << "       ";
      for(int t=0; t<nt; t++) {
	cout << vtvt[t].GetAverage(i)/vtvt[t].GetAverage(0);
	cout << " " << vtvt[t].GetErrorBar(i)/vtvt[t].GetAverage(0);
	if (t<nt-1) cout << "       ";
      }
      cout << endl;
    }
  }

  void PrintVVCorrelationIntegral(Array1 <double> & sInt, const int lMeasure) const {
    sInt.Resize(nt);
    cout << "#Diffusion coeff";
    for(int t=0; t<nt; t++) {
      cout << "         " << state.GetTypeName(t) << "     ";
    }
    cout << endl;
    Array1 <Spline> s(nt);
    for(int i=0; i<lHist; i++) {
      double tt = dt*i*nSteps;
      for(int t=0; t<nt; t++) {
	s[t].PushBack(tt,vtvt[t].GetAverage(i)/vtvt[t].GetAverage(0));
      }
    }
    for(int t=0; t<nt; t++) {
      s[t].Init();
    }
    for(int i=0; i<lHist; i++) {
      double tt = dt*i*nSteps;
      cout << tt << "       ";
      for(int t=0; t<nt; t++) {
	double sIT = s[t].Integrate(0.0,tt);
	if (i==lMeasure) sInt[t]=sIT;
	cout << sIT;
	if (t<nt-1) cout << "       ";
      }
      cout << endl;
    }
  }

  const State & state;
  double dt;
  const int nSteps;      // record vtvt every so mamy steps
  const int nStepsAvMin; // start recording after so many steps
  const int lHist;       // #records in histogram, max(t-t') = dt*nSteps*(lHist-1)
  const int nt;
  int nStore,iStore;
  Array2 <BPoint> v;
  Array1 <double> vtvtParticle; // [iHist]
  Array1 <ArrayAverage> vtvt; // vtvt[Type].GetAverage(dt-index)
};

#endif // _ANALYSIS_
