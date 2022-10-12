/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Burkhard Militzer                                Washington 8-17-04 //
//                                                                     //
// Perform the classical MC sampling                                   //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#ifndef _SAMPLING_
#define _SAMPLING_

#include "Standard.h"
#include "Random.h"
#include "Point.h"
#include "State.h"
#include "AccRatio.h"
#include "Analysis.h"
#include "SphericalHarmonics.h" // includes typedef of Complex
#include "Switch.h"

class Sampling {
 public:
  State & state;
  State state0; // copy initial state for printing diffusion info
  const int nParticles;
  const Array1 <Point> & r;
  //  const Array1 <double> & charge;
  //  Array2 <const TwoBodyPotential*> & p; // instead go through state.p 
  double P0,E0;

  int nBlocks;
  int nBlocksMin;
  int nSteps;

  Array1 <double> distances; // displacement [iType] // earlier, was a Point rather than double
  Array1 <AccRatio> acc;
  AccRatio accVol;
  AccRatio accCellVector; // average over all spatial directions

  Analysis analysis;
  bool writeXYZFileFlag;
  ofstream ofXYZ;
  bool writeXBSFileFlag;
  ofstream ofXBS;
  bool writeFM_CODEFlag;
  ofstream ofFM_CODE;
  string fnBlock;

  const int nq;
  Array1 <int> L, M;
   // not a neighbor for r<R1 it r>R4. neighbor with full weight R2<=r<=R3
  Array1 <double> R_cutoff1, R_cutoff2, R_cutoff3, R_cutoff4;
  Array1 < Array1 <double> > NN,NN_new;
  Array1 < Array2 <Complex> > Q_ilm,Q_ilm_new; // size [#QTerms][#particles, M]
  Array1 < Array1 <Complex> > Q_lm, Q_lm_new;  // size [M]
  Array1 <double> Q, Q_new;
  double QTerm, QTerm_new;
  Array1 <double> QWeights;
  Array1 <double> QTargets;
  double QFactor;
  double P,logVolumeStepSize,latticeVectorStepSize;
  const bool doNotMove;
  
  Sampling(State & state_, const string & fn, const string & fnBlock_, const string & fnFM_CODE, const string & fnXBS, const int nBlocksMin_): 
    state(state_),
    nParticles(state.nParticles),
    r(state.r),
    //    charge(state.charge),
    //    p(state.p),
    P0(nParticles*state.GetTemperature()/Point::GetVolume()),
    //    E0(3.0/2.0*state.GetTemperature()), // ideal kinetic per particle
    E0(3.0/2.0*nParticles*state.GetTemperature()), // now total ideal kinetic energy
    nBlocksMin(nBlocksMin_),
    distances(state.nTypes),
    acc(state.nTypes),
    analysis(Point::nDim,state.kPoints.Size(),state.nq,state.nTypes,nBlocksMin_,Point::MinimumImageDistance()/2.0,state,E0,P0), // pass in kinetic energy (now per cell not per particle)
    fnBlock(fnBlock_),
    nq(state.nq),
    L(state.L),
    M(state.M),
    R_cutoff1(state.R_cutoff1),
    R_cutoff2(state.R_cutoff2),
    R_cutoff3(state.R_cutoff3),
    R_cutoff4(state.R_cutoff4),
    QWeights(state.QWeights),
    QTargets(state.QTargets),
    QFactor(state.QFactor),
    P(-1.0),
    logVolumeStepSize(0.0),    // disable volume changing moves until step size has been set
    latticeVectorStepSize(0.0), // disable cell-vector changing moves until step size has been set
    doNotMove(state.doNotMove)  
  {
    //    Write5(nq,L,M,QWeights,QFactor); Quit("Q");
    //    Write(Point::GetIsoBoxSize());

    writeXYZFileFlag = (fn.length()>0);
    if (writeXYZFileFlag) {
      string fnXYZ = fn + ".xyz";
      ofXYZ.open(fnXYZ.c_str());
      if (!ofXYZ) error("Could not open file",fnXYZ);
    }

    writeFM_CODEFlag = (fnFM_CODE.length()>0);
    if (writeFM_CODEFlag) {
      string fnFM_CODE = "TRAJECTORYn1";
      ofFM_CODE.open(fnFM_CODE.c_str());
      if (!ofFM_CODE) error("Could not open file",fnFM_CODE);
    }

    writeXBSFileFlag = (fnXBS.length()>0);
    if (writeXBSFileFlag) {
      state.WriteCompleteXBSFile(fnXBS+".bs");
      string XBSFileName = fnXBS+".mv";
      ofXBS.open(XBSFileName.c_str());
      if (!ofXBS) error("Could not open XBS file ",XBSFileName);
    }

    // do not use the box size inside 'state' for now
    state.SetBoxSize();  // state now include is own Point::boxSize(), which may need to be activated
    state.Init();        // since the boxSize may have changed, recalculate the pair potential and Fourier potential energy
    AllocateAllQVariables();
    InitializeAllQVariables();
    Check();
  }

  void RunAllBlocks(const bool stepSizeCorrectionFlag, const bool testBlockFlag);  
  void RunOneBlock(const int iBlock, const bool print=true);  
  void MakeOneMove(const int ip, const double & stepSize, const int tp);
  void MakeOneMoveFixedCMS(const int ip, const double & stepSize, const int tp);
  void MakeOneVolumeMove();
  void MakeOneCellVectorMove();
  void CalcEnergyDifference(const int ip, const Point & rOld, const Point & rNew, double & dEPair, double & dEQn, double & dEFP) {
    dEPair = CalcEnergyDifferencePairPotentials(ip,rOld,rNew);
    dEQn   = CalcEnergyDifferenceQPotentials(ip,rOld,rNew);
    dEFP   = CalcEnergyDifferenceFourierPotential(ip,rOld,rNew);
   }
  double CalcEnergyDifferencePairPotentials(const int ip, const Point & rOld, const Point & rNew);
  void CorrectAllStepSizes();
  void CorrectOneStepSize(const AccRatio & ar, double & d) const;

  void SetToZeroAccRations();
  void AddGOfR();
  void AddGOfR(const int j);
  void AddPressure();
  void AddDeltaV();
  void AddQTerms(const Array1<double> & Q) {
    Array1 <double> QNormalized(nq);
    for(int iq=0; iq<nq; iq++) {
      QNormalized[iq] = Q[iq] / (M[iq]*nParticles);
    }
    analysis.AddQTerms(Q,QNormalized);
  }
  //  void AddQ6();
  //  void AddQ6Org();
  void AllocateAllQVariables();
  void Add_One_Q_ilm(const int M, const Array2 <Complex> & Q6_ilm, const Array1 <double> & NN, Array1 <Complex> & Q6_lm, double & Q6) const;
  void CalculateAllQVariables(Array1 <Array1 <double> > & NN, Array1 < Array2 <Complex> > & Q_ilm, Array1 < Array1 <Complex> > & Q_lm, Array1 <double> & Q, double & QTerm) const;
  void InitializeAllQVariables();
  void CopyOneSetOfQVariables(const int M,
			      const Array1 <double> & NN1, Array1 <double> & NN2, 
			      const Array2 <Complex> Q_ilm1, Array2 <Complex> & Q_ilm2, 
			      const Array1 <Complex> Q_lm1,  Array1 <Complex> & Q_lm2, 
			      const double Q1, double & Q2) const;
  void CopyAllQVariables(const Array1 < Array1 <double> > & NN1, Array1 < Array1 <double> > & NN2, 
			 const Array1 < Array2 <Complex> > & Q_ilm1, Array1 < Array2 <Complex> > & Q_ilm2,
			 const Array1 < Array1 <Complex> > & Q_lm1,  Array1 < Array1 <Complex> > & Q_lm2,
			 const Array1 <double> & Q1, Array1 <double> & Q2,
			 const double QTerm1, double & QTerm2) const;

  double CalcEnergyDifferenceQPotentials(const int ip, const Point & rOld, const Point & rNew);
  double CalculateOneQTerm(const int iq, const int L) const;
  double CalculateOneQTerm(const int iq) const {
    return CalculateOneQTerm(iq,L[iq]);
  }
  void CheckAllQTerms() const;
  void CalculateAndPrintManyQTerms() const;
  double CalcEnergyDifferenceFourierPotential(const int ip, const Point & rOld, const Point & rNew);
  void Check() const {
    state.Check();
    CheckAllQTerms();
  }
  double OneQTerm(const int iq, const double QPerCell) const {
    double f = double(M[iq]*nParticles);
    // QPerCell is assumed to be the total Q term normalized "per cell" or "per N particles"
    // currently QTargets[iq] are normalized "per particle" and have been divived by (2*L+1)
    // Add 1/f factor because we an energy per normalized per cell (or per N partciles), not per N^2 particles
    //    Write6(iq,QPerCell,QTargets[iq],f,QWeights[iq],QFactor);
    return sqr(QPerCell - QTargets[iq]*f) * QWeights[iq] * QFactor / f;
  }
  double SumOfAllQTerms(const Array1 <double> & Q) const {
    double QTerm = 0.0;
    for(int iq=0; iq<nq; iq++) {
      QTerm += OneQTerm(iq,Q[iq]);
    }
    return QTerm;
  }
    
  double NeighborWeight(const int iq, const double r) const {
    if (r<=R_cutoff2[iq]) 
      return SmoothFastSwitch(r,R_cutoff1[iq],R_cutoff2[iq],0.0,1.0); // 0.0 for r<=Rcutoff1 and 1.0 for r>=R_cutoff2
    else
      return SmoothFastSwitch(r,R_cutoff3[iq],R_cutoff4[iq],1.0,0.0); // 1.0 for r<=Rcutoff3 and 0.0 for r>=R_cutoff4
  }

  void PrintXYZ(const int iBlock) const {
    state.WriteXYZGeometryFile(cout,iBlock);
  }
  void WriteXYZGeometryFile(const int iBlock) {
    if (writeXYZFileFlag) state.WriteXYZGeometryFile(ofXYZ,iBlock);
  }
  void WriteGeometryToXBSFile() {
    if (writeXBSFileFlag) state.WriteGeometryToXBSFile(ofXBS);
  }

  void SetDistance(const double d) {
    distances = d;
  }
  void SetDistances(const Array1 <double> & d) {
    if (d.Size()!=state.nTypes) error("SetDistances array size",d.Size(),state.nTypes);
    distances = d;
  }

  double AverageAcceptanceRatio() const {
    int nTrial=0;
    int nAccept=0;
    for(int t=0; t<state.nTypes; t++) {
      nTrial += acc[t].nTrial;
      nAccept += acc[t].nAccept;
    }
    return (nTrial==0) ? 0.0 : double(nAccept)/double(nTrial);
  }
  Array1 <double> AverageAcceptanceRatios() const {
    Array1 <double> ar(state.nTypes);
    for(int t=0; t<state.nTypes; t++) {
      ar[t] = acc[t].Ratio();
    }
    return ar;
  }
  double AverageVolumeAcceptanceRatio() const {
    return accVol.Ratio();
  }
  double AverageCellVectorAcceptanceRatio() const {
    return accCellVector.Ratio();
  }

  void SetPressure(const double P_) {
    P = P_;
  }

  void SetRelativeVolumeStepSize(const double relDV) {
    SetVolumeStepSize(Point::GetVolume()*relDV);
  }
  void SetVolumeStepSize(const double dV) { // pass in absolute delta V, derive logV step size for logV sampling
    //    volumeStepSize = dV;
    // 'logVolumeStepSize' is the step size in logV space, which is all the matters. As long as this is fixed, logV sampling is ok.
    // derive 'logVolumeStepSize' from an absolute 'dV'. This conversion may be approximate but the MC samplig is exact nevertheless
    logVolumeStepSize = dV / Point::GetVolume();
    Write2(dV,logVolumeStepSize);
  }

  double GetAverageCellVector() const {
    return pow(Point::GetVolume(),1.0/double(Point::nDim));
  }
  void SetRelativeCellVectorStepSize(const double relDC) {
    SetVolumeStepSize(GetAverageCellVector()*relDC);
  }
  void SetLatticeVectorStepSize(const double dc) { // pass in absolute delta c, derive log(c) step size for log(c) sampling
    latticeVectorStepSize = dc;
  }
};

#endif
