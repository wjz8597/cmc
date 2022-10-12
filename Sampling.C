/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Burkhard Militzer                                Washington 8-17-04 //
//                                                                     //
// Perform the MC sampling                                             //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#include "Sampling.h"
#include "PrintDiffusion.h"

void Sampling::CorrectOneStepSize(const AccRatio & ar, double & d) const {
  if (ar.Flag()==false) return; // no correct if nTrial==0

  double dNew = d;
  const double factor = 1.5;
  if (ar.Ratio()<0.40) dNew /= factor; // accRatio too low  --> make smaller steps
  if (ar.Ratio()>0.60) dNew *= factor; // accRatio too high --> make larger steps

  if (dNew!=d) {
    cout << "Step size correction for " << ar.label << " from " << d << " to " << dNew << endl;
  }
  d = dNew;
}

void Sampling::CorrectAllStepSizes() {
  for(int t=0; t<state.nTypes; t++) {
    CorrectOneStepSize(acc[t],distances[t]);
  }
  CorrectOneStepSize(accVol,logVolumeStepSize);
  CorrectOneStepSize(accCellVector,latticeVectorStepSize);
}

void Sampling::RunAllBlocks(const bool stepSizeCorrectionFlag, const bool testBlockFlag) {
  state0 = state; // make a copy to monitor diffusion
  for(int t=0; t<state.nTypes; t++) {
    acc[t].SetNameLabelMovers(state.listOfTypes[t]);
  }
  accVol.SetNameLabelMovers("vol");
  accCellVector.SetNameLabelMovers("cell");
  for(int i=0; i<nBlocks; i++) {
    RunOneBlock(i,testBlockFlag);
    if (stepSizeCorrectionFlag) CorrectAllStepSizes();
  }
  if (!testBlockFlag) analysis.EndSimulation();
  if (writeXYZFileFlag) ofXYZ.close();
  if (writeFM_CODEFlag) ofFM_CODE.close();
  if (writeXBSFileFlag) ofXBS.close();
}

/*
void Sampling::RunOneBlock(const int iBlock) {
  //  cout << "-------- Starting block " << iBlock+1 << " --------------" << endl;

  analysis.StartBlock(iBlock);
  for(int t=0; t<state.nTypes; t++) 
    acc[t].Reset();
  for(int i=0; i<nSteps; i++) {
    for(int ip=0; ip<state.nParticles; ip++) {
      MakeOneMove(ip);
    }
    AddGOfR();
    AddPressure();
  }

  state.Check();
  for(int t=0; t<state.nTypes; t++) 
    cout << acc[t];
  //  PrintXYZ(iBlock);
  WriteXYZGeometryFile(iBlock+1);
  analysis.EndBlock(iBlock);

  cout << "\\_________________________________ Finished block " << iBlock+1 << " _________________________________/" << endl;
  cout << endl;
}
*/

void WriteVASPPOSCARFile(const string & fn, const State & s);
void WriteFM_CODEFile(ostream & of, const State & s);

void Sampling::RunOneBlock(const int iBlock, const bool testBlockFlag) {
  //  cout << "-------- Starting block " << iBlock+1 << " --------------" << endl;

  //  PrintDiffusionAllTypes(state0,state,0.0,false); // print intermediate diffusion information 3/19/13
  analysis.StartBlock(iBlock);
  for(int t=0; t<state.nTypes; t++) 
    acc[t].Reset();
  accVol.Reset();
  accCellVector.Reset();

  for(int i=0; i<nSteps; i++) {
    //    for(int ip=0; ip<state.nParticles; ip++) {
    //      MakeOneMove(ip);
    //    }
    for(int t=0; t<state.nTypes; t++) {
      int nT = state.atomsByType[t].Size();
      int n  = state.nParticles;
      // move particles of each type n/nT times
      // type 1: N=99 --> move 50 of them --> 1 sweep with 50% selection prob.
      // type 2: N= 1 --> move 50 times   --> 50 sweeps with 100% selection prob.
      double f = double(n)/double(state.nTypes)/double(nT);
      int nSweeps = int(ceil(f)); // nSweeps could be a little too high
      for(int is=0; is<nSweeps; is++) {
	for(int i=0; i<nT; i++) { // loop over all particles of type 't'
	  int ii = state.atomsByType[t][i];
	  if (Random()<f/double(nSweeps)) {
	    MakeOneMove(ii,distances[t],t);
	    //	    Check(); cout << "Checking after single particle move." << endl;
	    AddGOfR(state.atomsByType[t][i]);
	  }
	}
	//	AddGOfR();
	//        AddPressure();
      }
      if (P>0.0) MakeOneVolumeMove(); // for now, use P==0 as flag for constant volume simulations
      //      Check(); cout << "Checking after volume changing move move." << endl;
             MakeOneCellVectorMove(); // fixed step size for now
      //      Check(); cout << "Checking after cell vector changing move move." << endl;
    }
    //    AddGOfR();
    AddPressure();
    //    AddQ6(); // now done inside MakeOneMove
  }

  Check();
  for(int t=0; t<state.nTypes; t++) 
    cout << acc[t];
  cout << accVol;
  cout << accCellVector;
  PrintDiffusionAllTypes(state0,state,0.0,false); // print intermediate diffusion information 3/19/13
  
  if (!testBlockFlag) {

    analysis.EndBlock(iBlock);
  
    //  PrintXYZ(iBlock);
    WriteXYZGeometryFile(iBlock+1);//True if Sampling.writeXYZFileFlag==true so that the fn.length()>0
    WriteGeometryToXBSFile();//True if writeXBSFileFlag==True
    if (fnBlock!="" && iBlock>=nBlocksMin) {
      string fn = fnBlock+"_"+IntToString(iBlock-nBlocksMin+1,NumberOfDigits(nBlocks));
      cout << endl;
      //    cout << "Block " << iBlock+1 << " : writing VASP POSCAR: " << fn << endl;
      double E = state.CalcTotalEnergy();
      cout << "Snapshot EPot/cell[Ha]= " << E            << " EPot/N[Ha]= " << E/state.nParticles            << endl;
      cout << "Snapshot EPot/cell[eV]= " << E*PC::AUToeV << " EPot/N[eV]= " << E/state.nParticles*PC::AUToeV << endl;
      cout << "Block " << iBlock+1 << " : ";
      WriteVASPPOSCARFile(fn,state);
      cout << endl;
    }
    cout << "\\_________________________________ Finished block " << iBlock+1 << " _________________________________/" << endl;
    cout << endl;
  } else {
    cout << "\\___________________ Finished test block _______________________/" << endl;
  }
}
 
void Sampling::AddGOfR() {
  for(int i=0; i<state.nParticles-1; i++) {
    for(int j=i+1; j<state.nParticles; j++) {
      double rijn = Point::DistancePBC(r[i],r[j]);
      analysis.AddDistanceToGOfR(state.typeIndex[i],state.typeIndex[j],rijn);
    }
  }
}

void Sampling::AddGOfR(const int j) {
  for(int i=0; i<state.nParticles; i++) {
    if (i!=j) {
      double rijn = Point::DistancePBC(r[i],r[j]);
      analysis.AddDistanceToGOfR(state.typeIndex[i],state.typeIndex[j],rijn);
    }
  }
}

void Sampling::AddPressure() {
  double virial = 0.0;
  for(int i=0; i<state.nParticles-1; i++) {
    int ti = state.typeIndex[i];
    for(int j=i+1; j<state.nParticles; j++) {
      int tj = state.typeIndex[j];
      Point rij = Point::DifferencePBC(r[i],r[j]);
      virial += rij * state.p(ti,tj)->DV(rij); // scalar product
    }
  }
  double kin      = 1.5*state.GetTemperature()*double(nParticles); // 3/2 N kb T
  double pressure = 2.0*kin - virial;
  pressure /= 3.0*Point::GetVolume();
  analysis.AddPressure(pressure);
}

/*
void Sampling::AddQ6Org() {

  double Q6=0.0;
  const int M = 6;
  for(int m=-M; m<=M; m++) {
    Complex Q6_lm = Complex(0.0,0.0); //initialize the Q6_lm

    for(int i=0; i<state.nParticles; i++) {
      Complex Q6_ilm = Complex(0.0,0.0); //q6 contribution of i'th particle
      int NN=0; //number of nearest neighbors' of the i'th particle
      for(int j=0; j<state.nParticles; j++) {
	if (i!=j) {
	  Point r_ij = Point::DifferencePBC(r[i],r[j]); 
	  double r_ijn = r_ij.Norm();
	  const double R_cutoff=4.0; 
	  if (r_ijn<R_cutoff) {
	    NN++;
	    double theta = Point::Angle(r_ij,Point(0,0,1));
	    double phi      = atan2(r_ij[1],r_ij[0]); 
	    Q6_ilm += SphericalHarmonicsOrg(M,m,theta,phi);
	    //	    Q6_ilm += SphericalHarmonicsCosTheta(M,m,cosTheta,phi);
	  }
	}
      }
      Q6_ilm /= double(NN); // divide by number of neighbors (double)
      Q6_lm += Q6_ilm;
      //   cout << "Q6_lm=" <<Q6_lm <<endl;
    }
    //    Q6=Q6+((Q6_lm)*conj(Q6_lm)).real(); 
    Q6 += norm(Q6_lm);
  }
  Q6 = sqrt(Q6);
  Q6 /= state.nParticles;
  //  cout <<"Q6=" <<Q6 <<endl; 
  cout << "org:"; Write(Q6);
  analysis.AddQ6(Q6);
}
*/

// this calculations the Q term from scratch for the current configuration
// iq specifies the set of radii 1,2,3,4
double Sampling::CalculateOneQTerm(const int iq, const int L) const {
  const int M = 2*L+1;
  
  Array1 <Complex> Q_lm(M,Complex(0.0,0.0)); // allocate memory and initialize the Q6_lm for all m values from -M...+M 
  Array1 <Complex> Q_ilm(M);                 // allocate memory only for Q6 contribution of i'th particle

  for(int i=0; i<state.nParticles; i++) {
    for(int m=0; m<M; m++) Q_ilm[m] = Complex(0.0,0.0); // initialize Q6 contribution of i'th particle

    double NNi=0.0; //number of nearest neighbors' of the i'th particle -  now a double because we have weighted neighbors
    for(int j=0; j<state.nParticles; j++) {
      if (i!=j) {
	Point r_ij = Point::DifferencePBC(r[i],r[j]); 
	double r_ijn = r_ij.Norm();
	if (r_ijn>R_cutoff1[iq] && r_ijn<R_cutoff4[iq]) {
	  double weight = NeighborWeight(iq,r_ijn); // 1.0 for r<=Rcutoff1 and 0.0 for r>=R_cutoff2
	  NNi += weight;
	  //	  if (i==0) { cout << "Re: "; Write3(j,weight,NN); }
	  double cosTheta = r_ij[2] / r_ijn;
	  double phi      = atan2(r_ij[1],r_ij[0]); // atan2(y,x) returns arctan(y/x)
	  AddSphericalHarmonicsCosTheta(L,cosTheta,phi,Q_ilm,weight);
	}
      }
    }
    if (NNi>1e-13) { // when we remove the last neighbor NN[i] will not beome exactly 0 again.
      for(int m=0; m<M; m++) {
	Q_ilm[m] /= double(NNi); // divide by number of neighbors (double)
	Q_lm[m]  += Q_ilm[m];
      }
    }
  }

  double Q=0.0;
  for(int m=0; m<M; m++) {
    Q += norm(Q_lm[m]);
  }

  Q = sqrt(Q);
  //  Q /= state.nParticles; // this line would have converted Q6 from "per cell" to "per particle", which we do not want 7/3/18

  //  cout << "CalculateOneQTerm()"; Write(Q);
  //  analysis.AddQ6(Q6);
  return Q;
}

void Sampling::CheckAllQTerms() const {
  for(int iq=0; iq<nq; iq++) {
    double q = CalculateOneQTerm(iq);
    if (fabs(q-Q[iq])>1e-6) error("Error in Q term wrong",iq, q, Q[iq], q-Q[iq]);
  }
}

void Sampling::CalculateAndPrintManyQTerms() const {
  for(int iq=0; iq<nq; iq++) {
    for(int L=0; L<=24; L+=2) {
      double Q = CalculateOneQTerm(iq,L);
      Write5(L,iq,Q,Q/(2*L+1),Q/(2*L+1)/nParticles);
    }
  }
}

/*
void Sampling::AddQ6() {
  analysis.AddQ6(Q6*Q6Factor);
}
*/

void Sampling::AllocateAllQVariables() {
  NN.Resize(nq);
  Q_ilm.Resize(nq);
  Q_lm.Resize(nq);
  Q.Resize(nq);
  for(int iq=0; iq<nq; iq++) {
    NN[iq].Resize(state.nParticles);
    Q_ilm[iq].Resize(state.nParticles,M[iq]);
    Q_lm[iq].Resize(M[iq]);
  }

  NN_new.Resize(nq);
  Q_ilm_new.Resize(nq);
  Q_lm_new.Resize(nq);
  for(int iq=0; iq<nq; iq++) {
    NN_new[iq].Resize(state.nParticles);
    Q_ilm_new[iq].Resize(state.nParticles,M[iq]);
    Q_lm_new[iq].Resize(M[iq]);
  }
  Q_new.Resize(nq);

  //  CalculateAndPrintManyQTerms();  // Quit("after printing many Q terms");
}

void Sampling::Add_One_Q_ilm(const int M, const Array2 <Complex> & Q_ilm, const Array1 <double> & NN, Array1 <Complex> & Q_lm, double & Q) const {
  Q_lm = Complex(0.0,0.0); // initialize the Q6_lm for all m values from -M...+M 

  for(int i=0; i<state.nParticles; i++) {
    if (NN[i]>1e-13) { // when we remove the last neighbor NN[i] will not beome exactly 0 again.
      for(int m=0; m<M; m++) {
	Q_lm[m]  += Q_ilm(i,m) / NN[i];
      }
    }
  }

  Q=0.0;
  for(int m=0; m<M; m++) {
    Q += norm(Q_lm[m]);
  }

  Q = sqrt(Q);
  //  Q6 /= state.nParticles; // this line would have converted Q6 from "per cell" to "per particle", which we do not want 7/3/18

  //  Write(Q6);
  //  analysis.AddQ6(Q6);
}

void Sampling::CalculateAllQVariables(Array1 < Array1 <double> > & NN, Array1 < Array2 <Complex> > & Q_ilm, Array1 < Array1 <Complex> > & Q_lm, Array1 <double> & Q, double & QTerm) const {

  for(int iq=0; iq<nq; iq++) {
    for(int i=0; i<state.nParticles; i++) {
      NN[iq][i]=0.0;                         // no neighbors counted yet
      for(int m=0; m<M[iq]; m++) {
	Q_ilm[iq](i,m) = Complex(0.0,0.0); // initialize Q6 contribution of i'th particle
      }
    }
  }
  
  for(int i=0; i<state.nParticles-1; i++) {
    for(int j=i+1; j<state.nParticles; j++) {
      Point r_ij = Point::DifferencePBC(r[i],r[j]); 
      double r_ijn = r_ij.Norm();
      for(int iq=0; iq<nq; iq++) {
	if (r_ijn>R_cutoff1[iq] && r_ijn<R_cutoff4[iq]) {
	  double weight = NeighborWeight(iq,r_ijn); // 1.0 for r<=Rcutoff1 and 0.0 for r>=R_cutoff2
	  NN[iq][i] += weight;
	  NN[iq][j] += weight;
	  double cosTheta = r_ij[2] / r_ijn;
	  double phi      = atan2(r_ij[1],r_ij[0]); // atan2(y,x) returns arctan(y/x)
	  AddSphericalHarmonicsCosTheta(L[iq],cosTheta,phi,Q_ilm[iq],i,j,weight);
	}
      }
    }
  }
  
  for(int iq=0; iq<nq; iq++) {
    Add_One_Q_ilm(M[iq], Q_ilm[iq], NN[iq], Q_lm[iq], Q[iq]);
  }
  QTerm = SumOfAllQTerms(Q);
}

void Sampling::CopyOneSetOfQVariables(const int MM,
				      const Array1 <double> & NN1, Array1 <double> & NN2, 
				      const Array2 <Complex> Q_ilm1, Array2 <Complex> & Q_ilm2, 
				      const Array1 <Complex> Q_lm1,  Array1 <Complex> & Q_lm2, 
				      const double Q1, double & Q2) const {
  for(int i=0; i<state.nParticles; i++) {
    NN2[i] = NN1[i];
  }
  
  for(int i=0; i<state.nParticles; i++) {
    for(int mm=0; mm<MM; mm++) {
      Q_ilm2(i,mm) = Q_ilm1(i,mm);
    }
  }

  for(int mm=0; mm<MM; mm++) {
    Q_lm2[mm] = Q_lm1[mm];
  }

  Q2 = Q1;
}

void Sampling::CopyAllQVariables(const Array1 < Array1 <double> > & NN1, Array1 < Array1 <double> > & NN2, 
				 const Array1 < Array2 <Complex> > & Q_ilm1, Array1 < Array2 <Complex> > & Q_ilm2,
				 const Array1 < Array1 <Complex> > & Q_lm1,  Array1 < Array1 <Complex> > & Q_lm2,
				 const Array1 <double> & Q1, Array1 <double> & Q2,
				 const double QTerm1, double & QTerm2) const {
  for(int iq=0; iq<nq; iq++) {
    CopyOneSetOfQVariables(M[iq], NN1[iq], NN2[iq], Q_ilm1[iq], Q_ilm2[iq], Q_lm1[iq], Q_lm2[iq], Q1[iq], Q2[iq]);
  }
  
  QTerm2 = QTerm1;
}
  
void Sampling::InitializeAllQVariables() {
  CalculateAllQVariables(NN,Q_ilm,Q_lm,Q,QTerm);
  CopyAllQVariables(NN, NN_new, Q_ilm, Q_ilm_new, Q_lm, Q_lm_new, Q, Q_new, QTerm, QTerm_new); // copy contents to standard variables to 'new' variables
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ip -- index of moving particle
// return ENew - EOld
inline double Sampling::CalcEnergyDifferencePairPotentials(const int ip, const Point & rOld, const Point & rNew) {
  const int tp = state.typeIndex[ip]; // type of moving particle
  double de=0.0;
  for(int t=0; t<state.GetNTypes(); t++) { // type of other particle
    for(int it=0; it<state.atomsByType[t].Size(); it++) {
      int i = state.atomsByType[t][it]; // index of other particle
      if (i!=ip) {
	de += state.p(tp,t)->V(rNew,r[i]);
	de -= state.p(tp,t)->V(rOld,r[i]);
      }
    }
  }
  return de;
}

// ip -- index of moving particle
// assumes all 'new' variable are identical to the standard variables
inline double Sampling::CalcEnergyDifferenceQPotentials(const int ip, const Point & rOld, const Point & rNew) {
  //  if (QFactor==0.0) return 0.0;

  for(int i=0; i<state.nParticles; i++) {
    if (i!=ip) {

      Point dOld = Point::DifferencePBC(r[i],rOld);
      double dnOld = dOld.Norm();
      for(int iq=0; iq<nq; iq++) {
	if (dnOld>R_cutoff1[iq] && dnOld<R_cutoff4[iq]) {
	  double weight = -NeighborWeight(iq,dnOld); // negative weight to remove old contributions
	  NN_new[iq][i]  += weight;
	  NN_new[iq][ip] += weight;
	  //	Write3(i,weight,NN[ip]);
	  double cosTheta = dOld[2] / dnOld;
	  double phi      = atan2(dOld[1],dOld[0]); // atan2(y,x) returns arctan(y/x)
	  AddSphericalHarmonicsCosTheta(L[iq],cosTheta,phi,Q_ilm_new[iq],i,ip,weight); // subtract out old contributions
	}
      }

      Point dNew = Point::DifferencePBC(r[i],rNew);
      double dnNew = dNew.Norm();
      for(int iq=0; iq<nq; iq++) {
	if (dnNew>R_cutoff1[iq] && dnNew<R_cutoff4[iq]) {
	  double weight = +NeighborWeight(iq,dnNew);
	  NN_new[iq][i]  += weight;
	  NN_new[iq][ip] += weight;
	  //	cout << "Diff Q6: "; Write3(i,weight,NN[ip]);
	  double cosTheta = dNew[2] / dnNew;
	  double phi      = atan2(dNew[1],dNew[0]); // atan2(y,x) returns arctan(y/x)
	  AddSphericalHarmonicsCosTheta(L[iq],cosTheta,phi,Q_ilm_new[iq],i,ip,weight); // add new contributions
	}
      }
    }
  }
  
  for(int iq=0; iq<nq; iq++) {
    Add_One_Q_ilm(M[iq], Q_ilm_new[iq], NN_new[iq], Q_lm_new[iq], Q_new[iq]);
    //  Write2(iq,Q[iq],Q_new[iq]);
  }
  return SumOfAllQTerms(Q_new) - SumOfAllQTerms(Q);

}

// ip -- index of moving particle
// assumes all 'new' variable are identical to the standard variables
inline double Sampling::CalcEnergyDifferenceFourierPotential(const int ip, const Point & rOld, const Point & rNew) {
  if (state.fourierPotentialFactor==0.0) return 0.0;
  
  for(int ik=0; ik<state.kPoints.Size(); ik++) {
    state.rhoKNew[ik] = state.rhoK[ik] + state.CalculateRhoKi(rNew,state.kPoints[ik]) - state.CalculateRhoKi(rOld,state.kPoints[ik]);
  }
  return state.fourierPotentialFactor * (state.SumOfAllSOfKTermsFromRhoKNew() - state.SumOfAllSOfKTermsFromRhoK());
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void Sampling::MakeOneMove(const int ip, const double & stepSize, const int tp) {
  Point d=state.Uniform(stepSize);
  Point rOld = r[ip];  
  Point rNew = r[ip] + d;

  double dEPair,dEQ,dEFP;
  CalcEnergyDifference(ip,rOld,rNew,dEPair,dEQ,dEFP);
  double action = state.beta*(dEPair+dEQ+dEFP);
  QTerm_new = QTerm + dEQ;
  //  Write3(QTerm,dEQ,QTerm_new);
  
  acc[tp].Try();
  if (Random() <= exp(-action) && doNotMove==false) {
    acc[tp].Accept();
    state.UpdateTotalEnergy(dEPair);
    state.UpdateCoordinate(ip,rNew);
    //    cout << "Accepted:"; Write4(ip,tp,Q6,Q6_new);
    analysis.AddDisplacement(d.Norm());
    CopyAllQVariables(NN_new, NN, Q_ilm_new, Q_ilm, Q_lm_new, Q_lm, Q_new, Q, QTerm_new, QTerm); // copy contents of 'new' variables to standard variables
    state.UpdateFourierPotentialEnergyAndRhoK(dEFP);
  } else {
    analysis.AddDisplacement(0.0);
    CopyAllQVariables(NN, NN_new, Q_ilm, Q_ilm_new, Q_lm, Q_lm_new, Q, Q_new, QTerm, QTerm_new); // copy contents of standard variables to 'new' variables
    //    cout << "Rejected:"; Write4(ip,tp,Q6,Q6_new);
  }

  //  RecalculateQ6Term();
  //  cout << "MakeOneMove()"; Write(Q6);
  /*
  double QTerm = 0.0;
  for(int iq=0; iq<nq; iq++) {
    QTerm += Q[iq]*QWeights[iq];
  }
  QTerm *= Q6Factor;
  */
  
  analysis.AddTotalEnergy( state.GetTotalEnergy() + QTerm + state.GetFourierPotentialEnergy());  
  analysis.AddPairEnergy(  state.GetTotalEnergy()              );  
  analysis.AddQEnergy(     QTerm                               );  
  analysis.AddFPEnergy(    state.GetFourierPotentialEnergy()   );
  analysis.AddFPTerm(      state.GetFourierPotential()         );
  analysis.AddSOfKTerm(    state.ArrayOfSOfKTerms()            );
  AddQTerms(Q);
  //  analysis.PrintQHistograms(); Quit("PQ");
  //  if (ip==0) { Write(Q6); RecalculateQ6Term(); }
  //  Check();
}

void Sampling::MakeOneVolumeMove() {
  if (logVolumeStepSize<=0.0) return; // disable volume-changing moves until the step size has been set
  accVol.Try();

  double V = Point::GetVolume();
  //  double dV=state.Uniform1D(volumeStepSize);
  //  double dV=0.0;
  double dV=( exp(state.Uniform1D(logVolumeStepSize)) - 1.0 ) * V; // sample log(V) - see notes from 1/06/18, better **4243**
  double VNew = V+dV;
  if (VNew<0.0) return; // reject - should never happen b/c we sample log(V)
    
  double f = exp(log(VNew/V)/3.0);
  //  Point boxSize    = Point::GetBoxSize();
  //  Point boxSizeNew = boxSize*f;
  //  Point::SetBoxSize(boxSizeNew);
  Array1 <Point> l = Point::GetLatticeVectors();
  Array1 <Point> lNew(Point::nDim);
  for(int d=0; d<Point::nDim; d++) lNew[d] = l[d]*f;
  Point::SetLatticeVectors(lNew);
  
  for(int i=0; i<state.nParticles; i++) {
    state.r[i] *= f;
  }

  CalculateAllQVariables(NN_new,Q_ilm_new,Q_lm_new,Q_new,QTerm_new);
  double EPair    = state.GetTotalEnergy();
  double EPairNew = state.CalcTotalEnergy();
  double E    = EPair    + QTerm;
  double ENew = EPairNew + QTerm_new;

  //  Write6(E,ENew,ENew-E,Q6,Q6_new,Q6_new-Q6);

  double G    = E    + P*V;
  double GNew = ENew + P*VNew;
  double dG   = GNew - G;
  double action = state.beta*dG;

  action += state.nParticles * log(VNew/V); // extra factor only b/c we sample log(V), see Schultz & Kofke PRE 84 (2011) 046712, Eq. (12)

  if (Random() <= exp(-action) && doNotMove==false) {
    //    cout << "Accepted V:"; Write10(f,V,VNew,dV,G,GNew,dG,Q6,Q6_new,action);
    accVol.Accept();
    state.UpdateTotalEnergy(EPairNew-EPair);
    state.UpdateKPointsForNewCell();
    CopyAllQVariables(NN_new, NN, Q_ilm_new, Q_ilm, Q_lm_new, Q_lm, Q_new, Q, QTerm_new, QTerm); // copy contents of 'new' variables to standard variables
    analysis.AddVolumeChange(fabs(dV));
    V = VNew;
    state.CopyCurrentBoxSizeFromPoint(); // this is not really needed but BETTER keep the boxSize inside state up to date
  } else {
    //    cout << "Rejected V:"; Write10(f,V,VNew,dV,G,GNew,dG,Q6,Q6_new,action);
    //    Point::SetBoxSize(boxSize);
    Point::SetLatticeVectors(l);
    f = 1.0/f;
    for(int i=0; i<state.nParticles; i++) state.r[i] *= f;
    CopyAllQVariables(NN, NN_new, Q_ilm, Q_ilm_new, Q_lm, Q_lm_new, Q, Q_new, QTerm, QTerm_new); // copy contents of standard variables to 'new' variables
    analysis.AddVolumeChange(0.0);
  }

  //  RecalculateQ6Term();
  //  cout << "MakeOneMove()"; Write(Q6);

  analysis.AddTotalEnergy( state.GetTotalEnergy() + QTerm);  
  analysis.AddPairEnergy(  state.GetTotalEnergy()        );  
  analysis.AddQEnergy(     QTerm                         );  
  analysis.AddVolume(      V                             );
  AddQTerms(Q);
  // do not add any of those since they do not change during volume moves
  //  analysis.AddFPEnergy(    state.GetFourierPotentialEnergy()   );
  //  analysis.AddFPTerm(      state.GetFourierPotential()         );
  //  analysis.AddSOfKTerm(    state.ArrayOfSOfKTerms()            );

  Array1 <double> cell(Point::nDim);
  //  for(int d=0; d<Point::nDim; d++) 
  //    cell[d] = Point::GetBoxSize(d);
  Point cellVectorNorms = Point::GetLatticeVectorNorms();
  for(int d=0; d<Point::nDim; d++)
    cell[d] = cellVectorNorms[d];
  analysis.AddCellVectors(cell); // volume changing moves also affect the cell vectors

  //  if (ip==0) { Write(Q6); RecalculateQ6Term(); }
  //  Check();
}

void Sampling::MakeOneCellVectorMove() {
  if (latticeVectorStepSize<=0.0) return; // disable cell-vector-changing moves until the step size has been set
  accCellVector.Try();

  double V    = Point::GetVolume();
  Array1 <Point> l    = Point::GetLatticeVectors();
  Array1 <Point> lNew = l;
  
  int v = int(Point::nDim * Random()); // for every move, pick a new lattice vector at random
  int d = int(Point::nDim * Random()); // for every move, pick a new spatial direction at random
  double dl = state.Uniform1D(latticeVectorStepSize);
  lNew[v][d] += dl;

  bool ok = Point::LatticeVectorsMinimal(lNew[0],lNew[1],lNew[2]);
  // so far, cells with a=tiny and beta=gamma=90 are not yet excluded, therefore add another criteria
  if (Point::MaximalLatticeVectorRatio(lNew[0],lNew[1],lNew[2])>2.0) ok=false;
  if (ok) {

    state.ConvertAllAtomsToReducedCoordinates();
    Point::SetLatticeVectors(lNew);
    state.ConvertAllAtomsToCartesianCoordinates();
    double VNew = Point::GetVolume();

    CalculateAllQVariables(NN_new,Q_ilm_new,Q_lm_new,Q_new,QTerm_new);
    double EPair    = state.GetTotalEnergy();
    double EPairNew = state.CalcTotalEnergy();
    double E    = EPair    + QTerm;
    double ENew = EPairNew + QTerm_new;

    //  Write6(E,ENew,ENew-E,Q6,Q6_new,Q6_new-Q6);

    double G    = E    + P*V;
    double GNew = ENew + P*VNew;
    double dG   = GNew - G;
    double action = state.beta*dG;
    
    if (Random() <= exp(-action) && doNotMove==false) {
      //    cout << "Accepted V:"; Write9(V,VNew,dV,G,GNew,dG,Q6,Q6_new,action);
      accCellVector.Accept();
      state.UpdateTotalEnergy(EPairNew-EPair);
      CopyAllQVariables(NN_new, NN, Q_ilm_new, Q_ilm, Q_lm_new, Q_lm, Q_new, Q, QTerm_new, QTerm); // copy contents of 'new' variables to standard variables
      analysis.AddCellVectorChange(fabs(dl));
      V = VNew;
      //    state.CopyCurrentBoxSize(); // this is not really needed but keep the boxSize inside state up to date

    } else {
      
      //    cout << "Rejected V:"; Write9(V,VNew,dV,G,GNew,dG,Q6,Q6_new,action);
      state.ConvertAllAtomsToReducedCoordinates();
      Point::SetLatticeVectors(l);
      state.ConvertAllAtomsToCartesianCoordinates();
      CopyAllQVariables(NN, NN_new, Q_ilm, Q_ilm_new, Q_lm, Q_lm_new, Q, Q_new, QTerm, QTerm_new); // copy contents of standard variables to 'new' variables
      analysis.AddCellVectorChange(0.0);
    }
  } else {
    analysis.AddCellVectorChange(0.0);
  }
    
  //  RecalculateQ6Term();
  //  cout << "MakeOneMove()"; Write(Q6);
  
  analysis.AddTotalEnergy( state.GetTotalEnergy() + QTerm);  
  analysis.AddPairEnergy(  state.GetTotalEnergy()              );  
  analysis.AddQEnergy(     QTerm                               );  
  analysis.AddVolume(      V                                   );  
  AddQTerms(Q);
  //  if (ip==0) { Write(Q6); RecalculateQ6Term(); }
  //  Check();

  Array1 <double> cell(Point::nDim);
  //  for(int d=0; d<Point::nDim; d++) 
  //    cell[d] = Point::GetBoxSize(d);
  Point cellVectorNorms = Point::GetLatticeVectorNorms();
  for(int d=0; d<Point::nDim; d++)
    cell[d] = cellVectorNorms[d];
  analysis.AddCellVectors(cell);

  Array1 <double> angles(Point::nDim);
  for(int d=0; d<Point::nDim; d++) 
    angles[d] = Point::GetUnitCellAngleDegrees(d);
  analysis.AddCellAngles(angles);
  
}

