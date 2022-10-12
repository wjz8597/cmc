/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Burkhard Militzer                                    Urbana 4-9-99  //
//                                                                     //
// State of a classical MC simulation                                  //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#ifndef _STATE_
#define _STATE_

#include "Point.h"
#include "InitialPos.h"
#include "Physics.h"
#include "Spline.h"
#include "TwoBodyPotential.h"
#include "Form.h"
#include "Parser.h"
#include "Random.h"
#include "SphericalHarmonics.h"

/*
class Potential {
 public:
  Potential(const double cell):ewald(cell) {
    //    warning("Quadratic potential and force interpolation");
  };

  double OneBodyPot() const {
    return ewald.MadelungConstant();
  }

  double operator() (const Point & r) const {
    //    return 1.0/r.Norm();
    return ewald.FullPotentialInterpolation(r);
    //    return ewald.FullPotentialInterpolationQuadratic(r);
    //    return ewald.FullPotentialCalculation(r);
  }
  
  double operator() (const Point & r1, const Point & r2) const {
    Point r;
    r.Separation(r1,r2);
    return (*this)(r);
  }

  BPoint DV(const BPoint & r) const {
    BPoint f;
    ewald.FullForceInterpolation(r,f);
    //    ewald.FullForceInterpolationQuadratic(r,f);
    //    ewald.FullForceCalculation(r,f);
    return f;
  }

  BPoint DV(const Point & r1, const Point & r2) const {
    Point r;
    r.Separation(r1,r2);
    return DV(r);
  }

  Ewald ewald;
};
*/

inline double AtomicMassNumber(const string & typeName) {
  double mm = 0.0;
  //  if (typeName == "H")  { mm =  2.0; warning("Doing a DEUTERIUM calculation instead HYDROGEN!!!"); }
  if (typeName == "H")  mm =  1.00794;
  if (typeName == "He") mm =  4.00260;
  if (typeName == "Li") mm =  6.94100; // 7
  if (typeName == "Be") mm =  9.01218;
  if (typeName == "B")  mm = 10.81100; // 11
  if (typeName == "C")  mm = 12.0107;
  if (typeName == "N")  mm = 14.00674; // 14.0;
  if (typeName == "O")  mm = 15.9994;  // 16.0;
  if (typeName == "F")  mm = 18.9984;  // 19.0;
  if (typeName == "Mg") mm = 24.3050;  
  if (typeName == "Ne") mm = 20.1797;  // 22.0; 
  if (typeName == "Na") mm = 22.98977; // 23.0;
  if (typeName == "Al") mm = 26.9815386;
  if (typeName == "Si") mm = 28.0855;  // 28
  if (typeName == "P")  mm = 30.973762;
  if (typeName == "S")  mm = 32.065; // ± 0.005 u
  if (typeName == "Ar") mm = 39.948;  // heavier than "K"
  if (typeName == "K")  mm = 39.0983;
  if (typeName == "Ca") mm = 40.078;
  if (typeName == "Fe") mm = 55.845;
  if (typeName == "Co") mm = 58.9332;
  if (typeName == "Ga") mm = 69.723;
  if (typeName == "Kr") mm = 83.798;
  if (typeName == "I")  mm = 126.90447;
  if (typeName == "Xe") mm = 131.29;
  if (typeName == "Gd") mm = 157.25;
  if (typeName == "Pb") mm = 207.2;
  if (typeName == "Po") mm = 209.0;
  if (typeName == "H2O") mm = 18.01528; // wiki
  if (typeName == "D2O") mm = 18.01528+2.0; // need to look up correct mass!!!
  if (mm==0.0) warning("Do not know the mass for the this type",typeName);
  return mm;
}

inline double AtomicMass(const string & s) {
  return AtomicMassNumber(s)*PC::u/PC::me;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// input: lambda = hb^2/(2*m) in A.U.
//        beta   = 1/(kb T)   in A.U.
inline double LambdaDeBroglieLambdaBeta(const double lambda, const double beta) {
  double laDB = sqrt(4.0*pi*lambda*beta);
  return laDB;
}


// input: mass in m_e and T in [Ha/kb]
inline double LambdaDeBroglieMassTemperature(const double m, const double T) {
  //  return PC::h/sqrt(2.0*pi*PC::kb*T*m) / PC::ab;
  //  return 2.0*pi/sqrt(2.0*pi*T*m); // in atmomic units
  return LambdaDeBroglieLambdaBeta(0.5/m,1.0/T);
}

// input: lambda = hb^2/(2*m) in A.U.
// input: nn = N/V with V=ab^3 
// returns  f = F/N 
inline double FreeEnergyBoltzmannGas(const double lambda, const double beta, const double nn, const double g=1.0, const double e=0.0) {
  if (nn<0.0) error("nn",nn);
  if (nn==0.0) return 0.0;
  if (beta<0.0) error("beta",beta);
  double laDB = LambdaDeBroglieLambdaBeta(lambda,beta);
  double f = ( log ( laDB*laDB*laDB /g*nn ) - 1.0 ) / beta + e; // F/N
  //  f *= nn;  // F/N -> F/V
  return f;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline double LogFactorial(const int n) {
  double lf = 0.0;
  for(int i=2; i<=n; i++) 
    lf += log(double(i));
  return lf;
}

class State {
  //  Array1 <TwoBodyPotential> pList; // array of abstract class
  //  "TwoBodyPotential" is not allowed: In C++ a class with at
  //  least one pure virtual function is called Abstract class and you
  //  cannot create objects of that class, You can only have pointers
  //  or references to it.
 public:
  Array2 <const TwoBodyPotential*> p;
  ZeroPotential zeroPotential;

  int nParticles,nTypes;
  Array1 <Point> r;
  Array1 <string> names;
  Array1 <double> mass;
  double totalMass;
  Array1 <int> typeIndex;              // atom belongs to type j [#atoms]
  Array1 <string> listOfTypes;         // names of types [#types]
  Array1 <Array1 <int> > atomsByType;  // atom indices of given type [#types][#atom of this type] 

  Array1 < ShortVector<int,Point::nDim> > kPointIndices;
  Array1 <BPoint> kPoints;
  Array1 <double> kPointWeights;
  Array1 <Complex> rhoK;
  Array1 <Complex> rhoKNew;

  double totalEnergy;
  double fourierPotentialEnergy,fourierPotentialFactor;
  double beta;
  Point boxSize;

  int nq;
  Array1 <int> L, M;
   // not a neighbor for r<R1 it r>R4. neighbor with full weight R2<=r<=R3
  Array1 <double> R_cutoff1, R_cutoff2, R_cutoff3, R_cutoff4;
  Array1 <double> QTargets;
  Array1 <double> QWeights;
  double QFactor;
  bool doNotMove;
  
 State():nParticles(0),nTypes(0),totalMass(0.0),fourierPotentialFactor(0.0),nq(0),doNotMove(false) {
  }

 State(const double cell):nParticles(0),nTypes(0),totalMass(0.0),fourierPotentialFactor(0.0),nq(0),doNotMove(false) {
    Point::SetBoxSize(cell);
  }
    
  // put them in as they come - not in alphabetical order 
  int RegisterType(const string & typeName) {
    int i=listOfTypes.Find(typeName);
    if (i!=notFound) return i; // type has been entered before, return this index

    atomsByType.PushBack();
    listOfTypes.PushBack(typeName);
    nTypes++;
    p.CopyResize(nTypes,nTypes);
    for(int t=0; t<nTypes; t++) {
      p(t,nTypes-1) = &zeroPotential;
      p(nTypes-1,t) = &zeroPotential;
    }
    return listOfTypes.Size()-1;
  }

  int FindType(const string & typeName, const bool flag) const {
    int i=0;
    while (i<listOfTypes.Size() && typeName!=listOfTypes[i]) 
      i++;

    // type has been enter before, return this index
    if (i<listOfTypes.Size()) 
      return i;
    
    if (flag) error("Could not find type in list",typeName,listOfTypes);
    return notFound;
  }

  void RegisterPotential(const int t1, const int t2, const TwoBodyPotential & pp) {
    if (t1<0 || t1>=nTypes) error("Type index problem 1",t1,nTypes);
    if (t2<0 || t2>=nTypes) error("Type index problem 2",t2,nTypes);
    //#ifdef POT_COPY // << array of abstract class "TwoBodyPotential" is not allowed
    //    pList.PushBack(pp);
    //    p(t1,t2) = & pList.Last();
    //    p(t2,t1) = & pList.Last();
    //#else
    p(t1,t2) = & pp;
    p(t2,t1) = & pp;
    //#endif
  }
  int GetNTypes() const {
    return nTypes; // == listOfTypes.Size();
  }
  int GetType(const int i) const {
    return typeIndex[i];
  }
  string GetTypeName(const int t) const {
    return listOfTypes[t];
  }
  string GetAtomName(const int i) const {
    return listOfTypes[typeIndex[i]];
  }
  double GetTypeMass(const int t) const {
    if (GetNAtomsPerType(t)==0) error("Cannot give mass for type without atoms",t,GetTypeName(t));
    return mass[atomsByType[t][0]];
  }
  int GetNAtomsPerType(const int t) const {
    if (t==notFound) return 0;
    return atomsByType[t].Size();
  }

  //  void RegisterNAtoms(const string & name, const double c, const double m, const int n) {
  int RegisterNAtoms(const string & name, const double m, const int n) {
    //    if (n==0) return notFound;  // runs without He still need helium type!
    nParticles += n;
    int iType = RegisterType(name);
    for(int i=0; i<n; i++) {
      r.PushBack(Point(0.0));
      mass.PushBack(m);
      totalMass+=m;
      names.PushBack(name);
      //      charge.PushBack(c);
      int index = mass.Size()-1;
      typeIndex.PushBack(iType);
      atomsByType[iType].PushBack(index);
    }
    return iType;
  }

  Array1 <int> GetDefaultTypeArray() const { // need for VASP POSCAR files
    Array1 <int> typeOrder(GetNTypes());
    for(int t=0; t<GetNTypes(); t++)
      typeOrder[t]=t;
    return typeOrder;
  }

  double CalcTotalEnergy() const { 
    double e=0.0;
    for(int i=0; i<nParticles; i++) {
      int t = typeIndex[i]; 
      e += (p(t,t))->OneBodyTerm();
      //      Write3(i,t,e);
    }
    for(int i=0; i<nParticles-1; i++) {
      int ti = typeIndex[i];
      for(int j=i+1; j<nParticles; j++) {
	int tj = typeIndex[j];
	e += (p(ti,tj))->V(r[i],r[j]);
	//	Write5(i,j,ti,tj,e);
      }
    }
    //    Write2(nParticles,e);
    return e;
  }

  void InitTotalEnergy() {
    totalEnergy = CalcTotalEnergy();
  }
  void InitFourierPotentialEnergy() {
    SetAllRhoK();
    fourierPotentialEnergy = CalcFourierPotentialEnergy();
    //    Write2(fourierPotentialFactor,fourierPotentialEnergy);
  }
  void SetFourierPotentialFactor(const double fourierPotentialFactor_) {
    fourierPotentialFactor = fourierPotentialFactor_;
  }
  
  double GetTotalEnergy() const {
    return totalEnergy;
  }
  double GetFourierPotentialEnergy() const {
    return fourierPotentialEnergy;
  }
  double GetFourierPotential() const {
    return SumOfAllSOfKTermsFromRhoK();
  }
  double CalcFourierPotentialEnergy() {
    return fourierPotentialFactor * SumOfAllSOfKTermsFromRhoK();
  }
  
  void UpdateTotalEnergy(const double de) {
    totalEnergy += de;
  }
  void UpdateCoordinate(const int ip, const Point & rNew) {
    r[ip] = rNew;
  }

  // disable this safety feature for now because 'de' include the Q6 contribution while CalcTotalEnergy() does not
  void Check() { 
    double e=CalcTotalEnergy();
    double de = fabs(e-totalEnergy);
    if (de>1e-8*fabs(e) && de>1e-12) error("Total energy problem",e,totalEnergy,de);
    totalEnergy = e;

    for(int ik=0; ik<kPoints.Size(); ik++) {
      Complex rhoKi = CalculateOneRhoK(kPoints[ik]);
      Complex dRhoKi = rhoKi - rhoK[ik];
      if (norm(dRhoKi)>1e-8) error("Rho_k[i] incorrect",ik,rhoK[ik],rhoKi,dRhoKi);
    }
    double eFP = CalcFourierPotentialEnergy();
    double deFP = fabs(eFP-fourierPotentialEnergy);
    if (deFP>1e-8*fabs(eFP) && deFP>1e-12) error("Fourier potential energy problem",eFP,fourierPotentialEnergy,deFP);
    fourierPotentialEnergy = eFP;
  }
  
  void CopyCurrentBoxSizeFromPoint() {
    boxSize = Point::GetBoxSize();
  }
  void SetBoxSize(const double cell) const {
    Point::SetBoxSize(cell);
  }
  void SetBoxSize() const {
    Point::SetBoxSize(boxSize);
  }

  void SetCoordinates() {
    //    for(int i=0; i<nParticles; i++) {
    //      r[i] = Uniform(Point::GetBoxSize());
    //    }
    InitialPos(nParticles,Point::GetMinBoxSize(),0.0,r);
    for(int i=0; i<nParticles; i++) Write2(i,r[i]);
  }
  void Init() {
    //    SetCoordinates();
    InitTotalEnergy();
    InitFourierPotentialEnergy();
    SetAllRhoK();
  }
  
  void SetTemperature(const double t) {
    //    T=t;
    beta=1.0/t;
  }

  double GetTemperature() const {
    return 1.0/beta;
  }
  double GetTotalMass() const {
    //    return mass.Sum();
    return totalMass;
  }
  double GetRs() const {
    double V     = Point::GetVolume();
    double V1    = V/double(nParticles);
    double rs    = pow(V1*3.0/(4.0*pi),1.0/3.0); // not the true 'rs' except for H simulations
    return rs;
  }
  double MassDensity() const {
    return GetTotalMass()/Point::GetVolume();
  }
  double MassDensitySI() const {
    return MassDensity()*PC::me/PC::ab/PC::ab/PC::ab;
  }
  double MassDensityGCC() const {
    return MassDensitySI()/1000.0;
  }

  void SetBeta(const double b) {
    beta = b;
  }

  double Uniform1D(const double d) const { // uniform random numbers from from -d/2 to d/2.
    return d*(Random()-0.5);
  }
  Point Uniform(const double d) const {
    Point p(d);
    return Uniform(p);
  }
  Point Uniform(const Point & d) const {
    Point r;
    for(int i=0; i<Point::nDim; i++) {
      r[i] = Uniform1D(d[i]);
    }
    return r;
  }

  Point CenterOfMass() const {
    Point rCM(0.0);
    for(int i=0; i<nParticles; i++) {
      rCM += r[i]*mass[i];
    }
    rCM /= GetTotalMass();
    return rCM;
  }
  Point CenterOfMassType(const int t) const { // assume all atoms of this type have the same mass
    Point rCM(0.0);
    for(int i=0; i<atomsByType[t].Size(); i++) {
      int ii = atomsByType[t][i];
      rCM += r[ii];
    }
    rCM /= double(atomsByType[t].Size());
    return rCM;
  }

  // Write in Ang instead of AU
  void WriteXYZGeometryFile(ostream & of, const int nc) const {
    of.precision(10);    
    of << nParticles << endl;
    of << nc << endl;
    for(int i=0; i<nParticles; i++) {
      of << r[i][0]*PC::AUToA << " " << r[i][1]*PC::AUToA << " " << r[i][2]*PC::AUToA << endl;
    }
  }

  void WriteXBSAtom(ofstream & of, const string & name, const Point & r, const bool flag) const {
    Point rr = r*PC::AUToA;
    if (flag) {
      of << "atom " << name << " ";
    }
    //    of << fix(16,10,rr[0]) << " " << fix(16,10,rr[1])<< " " << fix(16,10,rr[2]) << endl;
    of << fix(9,6,rr[0]) << " " << fix(9,6,rr[1])<< " " << fix(9,6,rr[2]) << endl;
  }
  void WriteXBSCellInfo(ofstream & of, const bool flag) const { // write coordinates of all 8 cell corners
    for(int i=0; i<=1; i++) {
      for(int j=0; j<=1; j++) {
	for(int k=0; k<=1; k++) {
	  string name = "c"+IntToString(i)+IntToString(j)+IntToString(k);
	  //	  double shift = 0.5;
	  double shift = 0.0;
	  WriteXBSAtom(of, name,
		       Point::GetLatticeVector(0)*(i-shift)+
		       Point::GetLatticeVector(1)*(j-shift)+
		       Point::GetLatticeVector(2)*(k-shift),flag);
	  if (flag) of << "spec " << name << " 0.01  blue" << endl;
	}
      }
    }
  }
  void WriteGeometryToXBSFile(ofstream & of) const {
    of << "frame" << endl;
    WriteXBSCellInfo(of,false);
    for(int i=0; i<nParticles; i++) {
      //      Point ri=r[i];
      //      ri.MapInsideZeroToOne();
      WriteXBSAtom(of,names[i],r[i],false);
    }
  }
  void WriteCompleteXBSFile(const string & filename) const {
    ofstream of;
    //  string fn = "in.bs";
    of.open(filename.c_str());
    if (!of) error("Could not open XBS file",filename);
    WriteCompleteXBSFile(of);
    of.close();
  }
  void WriteCompleteXBSFile(ofstream & of) const;
  
  void PrintPotentials() const {
    Form sci4(sci,12,4);
    cout << "              ";
    for(int t1=0; t1<GetNTypes(); t1++) {
      for(int t2=t1; t2<GetNTypes(); t2++) {
	cout << "       V_" << GetTypeName(t1) << "-" << GetTypeName(t2);
      }
    }
    cout << endl;
    double r2 = Point::GetMaxBoxSize()/2.0;
    double r1 = r2*1e-4;
    int nn=1000;
    for(int i=0; i<=nn; i++) {
      double r = r1*exp(log(r2/r1)*double(i)/double(nn));
      cout << "rr= " << sci4(r);
      for(int t1=0; t1<GetNTypes(); t1++) {
	for(int t2=t1; t2<GetNTypes(); t2++) {
	  BPoint rp(r/sqrt(Point::nDim)); // must be a point since it should also work for Ewald
	  double V    = (p(t1,t2))->V(rp);
	  cout << "  " << sci4(V);
	  //	  cout << "  " << sci4(V*PC::AUToeV);
	  //	  double dVdr = (p(t1,t2))->DV(rp).Norm();
	  //	  cout << "  " << sci4(dVdr); 
	  //	  double V = (p(t1,t2))->DDV(rp)(0,0);
	}
      }
      cout << endl;
    }
    cout << endl;
  }

  void PrintPotentialsMiguel() const {
    cout.precision(12);
    //    Form sci4(sci,12,4);
    Form sci4(sci,18,12);
    cout << "              ";
    for(int t1=0; t1<GetNTypes(); t1++) {
      for(int t2=t1; t2<GetNTypes(); t2++) {
	cout << "       V_" << GetTypeName(t1) << "-" << GetTypeName(t2);
      }
    }
    cout << endl;
    double r2 = Point::GetMaxBoxSize()/2.0;
    double r1 = 0.0;
    int iMin= 1;
    int nn=29;
    for(int i=iMin; i<=nn; i++) {
      double r = r1+(r2-r1)*double(i)/double(nn);
      r = 8.153273448275862E-002/PC::AUToA;
      cout << "rr= " << sci4(r*PC::AUToA);
      for(int t1=0; t1<GetNTypes(); t1++) {
	for(int t2=t1; t2<GetNTypes(); t2++) {
	  BPoint rp(r/sqrt(double(Point::nDim))); // must be a point since it should also work for Ewald
	  double V    = (p(t1,t2))->V(rp);
	  //	  cout << "  " << sci4(V);
	  cout << "  " << sci4(V*PC::AUToeV);
	  double dVdr = (p(t1,t2))->DV(rp).Norm();
	  cout << "  " << sci4(dVdr*PC::AUToeV/PC::AUToA/(r*PC::AUToA)); // prints [du/dr * 1/r] in eV/A^2
	  //	  double V = (p(t1,t2))->DDV(rp)(0,0);
	}
      }
      cout << endl;
    }
    cout << endl;
  }

  void PrintPotentialsMiguelOld() const {
    Form sci4(sci,12,8);
    cout << "              ";
    for(int t1=0; t1<GetNTypes(); t1++) {
      for(int t2=t1; t2<GetNTypes(); t2++) {
	cout << "       V_" << GetTypeName(t1) << "-" << GetTypeName(t2);
      }
    }
    cout << endl;
    double r2 = 2.3644493*PC::AToAU;
    double r1 = 0.0;
    int nn=29;
    for(int i=1; i<=nn; i++) {
      double r = r1+(r2-r1)*double(i)/double(nn);
      //      cout << "rr= " << sci4(r);
      cout << "rr= " << sci4(r*PC::AUToA);
      for(int t1=0; t1<GetNTypes(); t1++) {
	for(int t2=t1; t2<GetNTypes(); t2++) {
	  BPoint rp(r/sqrt(Point::nDim)); // must be a point since it should also work for Ewald
	  double V = (p(t1,t2))->V(rp);
	  //	  cout << "  " << sci4(V);
	  cout << "  " << sci4(V*PC::AUToeV);
	}
      }
      cout << endl;
    }
    cout << endl;
  }
  int GetNAtoms() const {
    return nParticles;
  }

  void ParseVaspPOSCARFile(Parser & p, const Array1 <string> types, const bool shift=false);
  void ParseVaspPOSCARFile(const string & file, const Array1 <string> types, const bool shift=false) {
    Parser p(file);
    p.UnSetIgnoreEmptyLines();
    ParseVaspPOSCARFile(p,types,shift);
  }
  void ParseVaspREFCARFile(Parser & p, const Array1 <string> types, const bool shift=false);
  void ParseVaspREFCARFile(const string & file, const Array1 <string> types, const bool shift=false) {
    Parser p(file);
    p.UnSetIgnoreEmptyLines();
    ParseVaspREFCARFile(p,types,shift);
  }
  void ParseVaspPOSCARFileSkipOneType(Parser & p, const Array1 <string> types, const string & typeNameSkip, const bool shift=false);
  void ParseVaspPOSCARFileSkipOneType(const string & file, const Array1 <string> types, const string & typeNameSkip, const bool shift=false) {
    Parser p(file);
    p.UnSetIgnoreEmptyLines();
    ParseVaspPOSCARFileSkipOneType(p,types,typeNameSkip,shift);
  }

  void PrintFreeEnergyTerm(const string & s, const int t, const double F1, const bool eVFlag) const {
    double F = double(GetNAtomsPerType(t)) * F1;
    cout << s << gen(2,GetTypeName(t)) << ": F[Ha/particle]= " << fix(14,8,F1)            << " F[Ha/cell]= " << fix(14,8,F) << endl;
    if (eVFlag) {
      cout << s << gen(2,GetTypeName(t)) << ": F[eV/particle]= " << fix(14,8,F1*PC::AUToeV) << " F[eV/cell]= " << fix(14,8,F*PC::AUToeV) << endl;
    }
  }

  double TotalFreeEnergyType(const int t, double & FTypeDeprecated, const bool print=false) const { 
    double F1;
    double F1Deprecated;
    int NT= GetNAtomsPerType(t);

    double V = Point::GetVolume();
    double nn = double(NT)/V;
    //      double m  = 1.660538921e-27/PC::me; // test for Shuai
    double m  = GetTypeMass(t);
    //      double laDB = LambdaDeBroglieMassTemperature(m,1.0/beta);
    double la = 0.5/m;
    double g  = 1.0; // no spin degree of freedom
    double E0 = 0.0;
    F1Deprecated = FreeEnergyBoltzmannGas(la,beta,nn,g,E0); // free energy per particle of this type
    PrintFreeEnergyTerm("Boltzmann free energy of type (Stirling's approx.) ",t,F1Deprecated,false);
    
    double Z1 = V*pow(4.0*pi*la*beta,-3.0/2.0); // partition function of 1 particle in 3D (-3.0/2.0 for 3D)
    //      double F1 = -log(Z1          )/beta; // distinguishable particles
    //      double F1 = -log(Z1*exp(1)/NT)/beta; // same as F1 = FreeEnergyBoltzmannGas(la,beta,nn,g,E0); (used Stirling's formula)
    F1 = -log(Z1          )/beta + LogFactorial(NT)/double(NT)/beta; // avoid Stirling's formula
    PrintFreeEnergyTerm("Boltzmann free energy of type (exact)              ",t,F1,false);
    
    Write5(V*PC::AUToA3,beta,m,m*PC::me,NT);
    
    FTypeDeprecated = F1Deprecated*double(NT);
    double FType = F1*double(NT);
    return FType;
  }

  double TotalFreeEnergy(const bool print=false) const { 
    double F           = 0.0; // in Ha/cell in 3D
    double FDeprecated = 0.0; // in Ha/cell in 3D

    for(int t=0; t<GetNTypes(); t++) {
      double F1Deprecated;
      F += TotalFreeEnergyType(t,F1Deprecated,print);
      FDeprecated += F1Deprecated;
      //      if (print)  cout << endl;
    }

    if (print) cout << endl;
    if (print) cout << "Total free energy, F:" << endl;
    if (print) cout << "F without fixed-particle correction with Stirling's approx.: F[Ha/cell]= " << fix(14,8,FDeprecated)            << endl;
    if (print) cout << "F without fixed-particle correction without Stirling:        F[Ha/cell]= " << fix(14,8,F)            << endl;
    if (print) cout << "Final F                                                      F[Ha/cell]= " << fix(14,8,F)            << endl;
    if (print) cout << endl;

    return F;
  }

  void ConvertAllAtomsToReducedCoordinates() {
    for(int i=0; i<nParticles; i++) {
      r[i].Reduce();
    }
  }
  
  void ConvertAllAtomsToCartesianCoordinates() {
    for(int i=0; i<nParticles; i++) {
      r[i].ToCartesian();
    }
  }

  void ScaleBox(const double f) {
    ConvertAllAtomsToReducedCoordinates();
    Array1 <Point> l = Point::GetLatticeVectors();
    for (int d=0; d<Point::nDim; d++) {
      l[d] *= f;
    }
    Point::SetLatticeVectors(l);
    ConvertAllAtomsToCartesianCoordinates();
  }

  int FindClosestNeighbor(const int ii) const;
  void PrintNeighbors(const double x=0.1, const bool relAbs=true, const bool printPositions=false, const int nMin=0, const int nAtoms=-1) const;

  BPoint DeriveKVector(const ShortVector<int,Point::nDim> & N) { // kx = 2*pi/L * Nx 
    BPoint k(0.0);
    for(int d=0; d<Point::nDim; d++) {
      k += N[d]*Point::GetReciprocalLatticeVector(d);
    }
    return k;
  }
  void AddKPoint(const int Nx, const int Ny, const int Nz, const double weight) {
    ShortVector<int,Point::nDim> N(Nx,Ny,Nz);
    kPointIndices.PushBack(N);
    kPoints.PushBack(DeriveKVector(N));
    kPointWeights.PushBack(weight);
    rhoK.PushBack();
    rhoKNew.PushBack();
  }
  void UpdateKPointsForNewCell() {
    for(int ik=0; ik<kPoints.Size(); ik++) {
      kPoints[ik] = DeriveKVector(kPointIndices[ik]);
    }
  }
  double SOfKTermSlow(const BPoint & k) const {
    double c=0.0;
    double s=0.0;
    //    for(int i=0; i<nParticles-1; i++) {
      //      for(int j=i+1; j<nParticles; j++) {
    for(int i=0; i<nParticles; i++) {
      for(int j=0; j<nParticles; j++){ 
	Point rij = Point::DifferencePBC(r[i],r[j]);
	double rk = rij*k;
	c += cos(rk);
	s += sin(rk);
	//	Write4(i,j,r[i],r[j]);
	//	Write7(i,j,rij,j,rk,c,s);
      }
    }
    //    return sqrt(c*c+s*s);
    //    return 0.25*(c*c+s*s);
    return c*c+s*s;
  }
  double SumOfAllSOfKTermsSlow() const {
    double sum = 0.0;
    for(int ik=0; ik<kPoints.Size(); ik++) {
      sum += kPointWeights[ik] * SOfKTermSlow(kPoints[ik]);
    }
    //    Write(fourierPotentialFactor*sum);
    return sum;
  }
  Complex CalculateRhoKi(const Point & r, const BPoint & k) const {
    double rk = r*k;
    double c = cos(rk);
    double s = sin(rk);
    Complex rhoKi(c,s);
    return rhoKi;
  }
  Complex CalculateOneRhoK(const BPoint & k) const {
    Complex rhoK(0.0,0.0);
    for(int i=0; i<nParticles; i++) {
      rhoK += CalculateRhoKi(r[i],k);
      //      Write6(i,r[i],k,rk,rhoKi,rhoK);
    }
    return rhoK;
  }
  void SetAllRhoK() {
    rhoK.Resize(kPoints.Size());
    rhoKNew.Resize(kPoints.Size());
    for(int ik=0; ik<kPoints.Size(); ik++) {
      rhoK[ik] = CalculateOneRhoK(kPoints[ik]);
    }
  }
  double SOfKTermFromRhoKi(const Complex & rhoKi) const {
    return sqr(norm(rhoKi));
  }
  double SumOfAllSOfKTermsFromRhoK(const Array1 <Complex> & myRhoK) const {
    double sum = 0.0;
    for(int ik=0; ik<kPoints.Size(); ik++) {
      sum += kPointWeights[ik] * SOfKTermFromRhoKi(myRhoK[ik]);
    }
    return sum;
  }
  Array1 <double> ArrayOfSOfKTerms() const {
    Array1 <double> sOfK2(kPoints.Size());
    for(int ik=0; ik<kPoints.Size(); ik++) {
      sOfK2[ik] = SOfKTermFromRhoKi(rhoK[ik]);
    }
    return sOfK2;
  }
  double SumOfAllSOfKTermsFromRhoK() const {
    return SumOfAllSOfKTermsFromRhoK(rhoK);
  }
  double SumOfAllSOfKTermsFromRhoKNew() const {
    return SumOfAllSOfKTermsFromRhoK(rhoKNew);
  }
  void UpdateFourierPotentialEnergyAndRhoK(const double dEFP) {
    for(int ik=0; ik<kPoints.Size(); ik++) {
      rhoK[ik] = rhoKNew[ik];
    }
    fourierPotentialEnergy += dEFP;
  }
  void CalculateManySOfK(const int nMax) {
    for(int ix=0; ix<=nMax; ix++) {
      for(int iy=0; iy<=nMax; iy++) {
	for(int iz=0; iz<=nMax; iz++) {
	  ShortVector<int,Point::nDim> iN(ix,iy,iz);
	  BPoint k = DeriveKVector(iN);
	  double sOfK = SOfKTermSlow(k);
	  if (sOfK>1e-04) {
	    cout << " ix= " << IntToStringMaxNumber(ix,nMax);
	    cout << " iy= " << IntToStringMaxNumber(iy,nMax);
	    cout << " iz= " << IntToStringMaxNumber(iz,nMax);
	    cout << " |k|= " << fix(8,6,k.Norm());
	    cout << " S(k)= " << sOfK;
	    cout << endl;
	  }
	}
      }
    }
  }

  void SetQFactorOverall(const double f) {
    QFactor=f;
  }
  void RegisterQTerm(const int l, const double r_cutoff1, const double r_cutoff2, const double r_cutoff3, const double r_cutoff4,
		     const double target, const double weight) {
    nq++;
    L.PushBack(l);
    M.PushBack(2*l+1);
    R_cutoff1.PushBack(r_cutoff1);
    R_cutoff2.PushBack(r_cutoff2);
    R_cutoff3.PushBack(r_cutoff3);
    R_cutoff4.PushBack(r_cutoff4);
    QTargets.PushBack(target);
    QWeights.PushBack(weight);
  }
};

#endif
