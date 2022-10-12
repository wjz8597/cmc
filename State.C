/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Burkhard Militzer                                    Urbana 4-9-99  //
//                                                                     //
// State of a classical MC simulation                                  //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#include "State.h"
#include "Sort.h"

// sets the box/cell size in Point and here in State !!!
void State::ParseVaspPOSCARFile(Parser & p, const Array1 <string> typeNames, const bool shift) {
  p.ReadLineSafely();
  cout << "Reading VASP POSCAR/CONTCAR file with headline: " << p.GetLineString() << endl;

  //  p.ReadLineSafely(1);
  p.ReadLineSafely();
  double f=p.GetDouble(0);
  f *= PC::AToAU;

  p.ReadLineSafely(3);
  Point a(p.GetDouble(0),p.GetDouble(1),p.GetDouble(2));
  p.ReadLineSafely(3);
  Point b(p.GetDouble(0),p.GetDouble(1),p.GetDouble(2));
  p.ReadLineSafely(3);
  Point c(p.GetDouble(0),p.GetDouble(1),p.GetDouble(2));
#ifdef POINTGEN
  Point::SetBoxSize(a*f,b*f,c*f,true);
#else
  if (a[1]!=0.0 || a[2]!=0.0) error("Compiled without #define POINTGEN for orthorhombic cells only",a);
  if (b[0]!=0.0 || b[2]!=0.0) error("Compiled without #define POINTGEN for orthorhombic cells only",b);
  if (c[0]!=0.0 || c[1]!=0.0) error("Compiled without #define POINTGEN for orthorhombic cells only",c);
  Point cell(a[0],b[1],c[2]);
  cell *= f;
  Point::SetBoxSize(cell);
  //  error("Cannot call Point::SetBoxSize() with 3 vectors without #define POINTGEN");
#endif
  CopyCurrentBoxSizeFromPoint();

  p.ReadLineSafely();
  int nTypes=p.GetNWords();
  if (nTypes>typeNames.Size()) error("Not enough type string available",nTypes,typeNames.Size());
  if (nTypes<typeNames.Size()) warning("More type strings given than needed",nTypes,typeNames.Size());

  for(int i=0; i<nTypes; i++) {
    string name = typeNames[i];
    double mass = AtomicMass(name);
    int n = p.GetInt(i);
    RegisterNAtoms(name,mass,n);
  }

  p.ReadLineSafely(1);
  if (p.GetStringLowerCase(0)!="direct") error("Use reduced coordinates and 'direct' in POSCAR file");

  double s=(shift) ? -0.5 : 0.0;

  for(int i=0; i<nParticles; i++) {
    p.ReadLineSafely(3);
    Point ri(p.GetDouble(0)+s,p.GetDouble(1)+s,p.GetDouble(2)+s);
    ri.ToCartesian();
    r[i] = ri;
  }
  cout << "POSCAR file with " << GetNAtoms() << " atoms parsed successfully." << endl;
}

void State::ParseVaspPOSCARFileSkipOneType(Parser & p, const Array1 <string> typeNames, const string & typeNameSkip, const bool shift) {
  p.ReadLineSafely();
  cout << "Reading VASP POSCAR/CONTCAR file with headline: " << p.GetLineString() << endl;
  cout << "Ignoring type " << typeNameSkip << endl;

  int typeIndexSkip = typeNames.Find(typeNameSkip);
  if (typeIndexSkip==notFound) error("Tyep name to be skipped not in list of types", typeNameSkip, typeNames);

  //  p.ReadLineSafely(1);
  p.ReadLineSafely();
  double f=p.GetDouble(0);
  f *= PC::AToAU;

  p.ReadLineSafely(3);
  Point a(p.GetDouble(0),p.GetDouble(1),p.GetDouble(2));
  p.ReadLineSafely(3);
  Point b(p.GetDouble(0),p.GetDouble(1),p.GetDouble(2));
  p.ReadLineSafely(3);
  Point c(p.GetDouble(0),p.GetDouble(1),p.GetDouble(2));
#ifdef POINTGEN
  Point::SetBoxSize(a*f,b*f,c*f,true);
#else
  if (a[1]!=0.0 || a[2]!=0.0) error("Compiled without #define POINTGEN for orthorhombic cells only",a);
  if (b[0]!=0.0 || b[2]!=0.0) error("Compiled without #define POINTGEN for orthorhombic cells only",b);
  if (c[0]!=0.0 || c[1]!=0.0) error("Compiled without #define POINTGEN for orthorhombic cells only",c);
  Point cell(a[0],b[1],c[2]);
  cell *= f;
  Point::SetBoxSize(cell);
  //  error("Cannot call Point::SetBoxSize() with 3 vectors without #define POINTGEN");
#endif

  p.ReadLineSafely();
  int nTypes=p.GetNWords();
  if (nTypes>typeNames.Size()) error("Not enough type string available",nTypes,typeNames.Size());
  //  if (nTypes<typeNames.Size()) warning("More type strings given than needed",nTypes,typeNames.Size());
  if (nTypes<typeNames.Size()) error("More type strings given than needed",nTypes,typeNames.Size());

  Array1 <int> nAtomsPerType(nTypes);
  for(int t=0; t<nTypes; t++) {
    string name = typeNames[t];
    double mass = AtomicMass(name);
    nAtomsPerType[t]=p.GetInt(t);
    if (name!=typeNameSkip) RegisterNAtoms(name,mass,nAtomsPerType[t]);
  }

  p.ReadLineSafely(1);
  if (p.GetStringLowerCase(0)!="direct") error("Use reduced coordinates and 'direct' in POSCAR file");

  double s=(shift) ? -0.5 : 0.0;

  int ii=0;
  for(int t=0; t<nTypes; t++) {
    for(int i=0; i<nAtomsPerType[t]; i++) {
      p.ReadLineSafely(3);
      Point rr(p.GetDouble(0)+s,p.GetDouble(1)+s,p.GetDouble(2)+s);
      rr.ToCartesian();
      if (t!=typeIndexSkip) {
	r[ii] = rr;
	Write4(t,i,ii,r[ii]);
	ii++;
      }
    }
  }
  cout << "POSCAR file with " << nAtomsPerType.Sum() << " atoms parsed successfully. Only " << GetNAtoms() << " used."<< endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// probe all atoms, ignore types
int State::FindClosestNeighbor(const int ii) const {
  int iMin = notFound;
  double dMin = Point::GetMaxDistance();
  //  double dMin = Point::MinimumImageDistance();
  for(int i=0; i<nParticles; i++) {
    if (i!=ii) {
      double d=Point::DistancePBC(r[i],r[ii]);
      //      double d=Point::MinimumImageDistance(r[i],r[ii]);
      if (d<dMin) {
	dMin=d;
	iMin=i;
      }
    }
  }

  //  if (number>=2 && iMin==notFound) {
  //    WriteVASPPOSCARFile("xxx");
  //    error("iMin==notFound");
  //  }
  // Careful, a particle can be close to one of its images than to any other particle
  // in this case return 'notFound'
  return iMin;
}


// x = tolerance
// relAbs = true  --> tolerance=relative
// relAbs = false --> tolerance=absolute
// printPositions = true --> print fraction coordincates of each atom before the NN distances
// nMin>0    --> print at least 'nMin' neighbors
// nAtoms>=0 --> do not print all atoms, stop after nAtoms
void State::PrintNeighbors(const double x, const bool relAbs, const bool printPositions, const int nMin, const int nAtoms) const {
  cout << " ___                 ___ " << endl;
  cout << "/   Neighbor analysis   \\" << endl;

  int n=(nAtoms<0) ? nParticles : nAtoms;
  if (n>nParticles) error("PrintNeighbors:: nAtoms too large",nAtoms,nParticles);

  if (n<=1) {
    cout << "Only " << n << " particle(s). No neighbor info to be printed." << endl;
    return;
  }

  int ii=10;
  int digits = 1;
  while (ii<n) { digits++; ii*=10;} // fix number of digits to print based on particles in cell
  Form f(fix,digits);

  for(int i=0; i<n; i++) { // main loop over all particles to be printed

    cout << f(i+1) << " "; // could be turned off
    //    double dj = Point::DistancePBC(r[i],r[j]);
    int j = FindClosestNeighbor(i);
    //    double dj = (j!=notFound) ? Point::MinimumImageDistance(r[i],r[j]) : Point::MinimumImageDistance();
    double dj = (j!=notFound) ? Point::DistancePBC(r[i],r[j]) : Point::MinimumImageDistance();
    double dMax = (relAbs) ? dj*(1.0+x) : x;
    //    Write4(i,dMax,relAbs,x);

    if (nMin>0) { // adjust dMax in case nMin>0 is provided (everything else stays the same)
      Array1 <double> d;
      for(int j=0; j<nParticles; j++) {
	if (i!=j) {
	  d.PushBack(Point::DistancePBC(r[i],r[j]));
	  //	  d.PushBack(Point::MinimumImageDistance(r[i],r[j]));
	}
      }
      SortArray(d);
      dMax = max(dMax,(1.0+1e-10)*d[nMin-1]);
    }

    Array1 <int> neighbors;
    Array1 <double> d;
    for(int j=0; j<nParticles; j++) {
      if (i!=j) {
	double dd = Point::DistancePBC(r[i],r[j]);
	//	double dd = Point::MinimumImageDistance(r[i],r[j]);
	//	Write3(i,j,dd*PC::AUToA);
	if (dd<dMax) {
	  neighbors.PushBack(j);
	  d.PushBack(dd);
	}
      }
    }

    Array1 <int> index = SortArrayIntoIndex(d);
    Array1 <int> sortedNeighbors(neighbors.Size());
    for(int j=0; j<neighbors.Size(); j++) sortedNeighbors[j] = neighbors[index[j]];

    cout << GetAtomName(i); 
    cout << ": ";
    cout << sortedNeighbors.Size() << " ";
    //    cout << " <"; for(int j=0; j<neighbors.Size(); j++) { cout << GetAtomName(neighbors(j)); } cout << "> ";
    //    cout << " <"; for(int j=0; j<neighbors.Size(); j++) { cout << neighbors(j) << "_"; } cout << "> ";
    int k=0;
    int nt=0;
    while (k<sortedNeighbors.Size()) {
      k++;
      nt++;
      if (k==sortedNeighbors.Size() || GetType(sortedNeighbors[k-1])!=GetType(sortedNeighbors[k])) {
	cout << GetAtomName(sortedNeighbors[k-1]);
	if (nt>1) {
	  cout << "_" << nt;
	}
	if (k!=sortedNeighbors.Size()) cout << " ";
	nt = 0;
      };
    }

    if (printPositions) { // prints fractional coordinates of each atom
      cout << "   ";
      Point rr = r[i].Reduced();
      Form f(fix,8,4);
      for(int d=0; d<Point::nDim; d++) {
	cout << f(rr[d]) << " ";
      }
      cout << "   ";
    }

    /*
    double dMin=dMax;
    dMax=0.0;
    for(int j=0; j<neighbors.Size(); j++) {
      double d = Point::DistancePBC(r[i],r[neighbors[j]]);
      if (d<dMin) dMin=d;
      if (d>dMax) dMax=d;
    }

    cout << " dist[A]= " << dMin*PC::AUToA;
    if (neighbors.Size()>1) 
      cout << " ... " << dMax*PC::AUToA;
    */

    /* is now simpler
    Array1 <double> dn(neighbors.Size());
    for(int j=0; j<neighbors.Size(); j++) {
      dn[j]=Point::DistancePBC(r[i],r[neighbors[j]]);
    }
    SortArray(dn);
    for(int j=0; j<dn.Size(); j++) {
      cout << " " << dn[j]*PC::AUToA;
    }
    */

    cout << " ";
    Form fn(fix,4,2); // fixed format, added 12-30-12
    for(int j=0; j<sortedNeighbors.Size(); j++) {
      int jj = sortedNeighbors[j];
      //      cout << " " << fn(Point::DistancePBC(r[i],r[jj])*PC::AUToA); // in Ang
      cout << " " << fn(Point::DistancePBC(r[i],r[jj])); // in aBohr
      //      cout << " " << fn(Point::MinimumImageDistance(r[i],r[jj])*PC::AUToA);
    }

    //    bool printAngles=true;
    bool printAngles=false;
    if (printAngles) {
      cout << "  :  ";
      Form fa(fix,3,0); 
      for(int j=0; j<sortedNeighbors.Size()-1; j++) {
	int jj = sortedNeighbors[j];
	for(int k=j+1; k<sortedNeighbors.Size(); k++) {
	  int kk = sortedNeighbors[k];
	  double angle = Point::AnglePBCSafelyDeg(r[jj],r[i],r[kk]);
	  cout << " " << fa(angle);
	}
      }
    }

    cout << endl;

  } // loop over particle index 'i'
  cout << "\\___Neighbor analysis___/" << endl; // << endl;
}

void State::WriteCompleteXBSFile(ofstream & of) const {
  WriteXBSCellInfo(of,true);
  for(int i=0; i<nParticles; i++) {
    //    Point ri=r[i];
    //    ri.MapInsideZeroToOne();
    WriteXBSAtom(of,names[i],r[i],true);
  }

  if (FindType("H",false)>=0) of << " spec H 0.0 1.0 1.0 1.0 " << endl;
  if (FindType("He",false)>=0) of << " spec He 0.15 0.0 0.0 1.0 " << endl;
  if (FindType("Ne",false)>=0) of << " spec Ne 0.2 1.0 0.2 0.5 " << endl;
  if (FindType("Ar",false)>=0) of << " spec Ar 0.3 green" << endl;
  if (FindType("N",false)>=0) of << " spec N 0.27 darkgreen" << endl;
  if (FindType("O",false)>=0) of << " spec O 0.23 1.0 0.0 0.0 " << endl;
  if (FindType("C",false)>=0) of << " spec C 0.23 0.2 0.2 1.0" << endl;
  if (FindType("K",false)>=0) of << " spec K 0.5 blue" << endl;
  if (FindType("Mg",false)>=0) of << " spec Mg 0.4 0.5 0.5 1.0" << endl;
  if (FindType("Si",false)>=0)  of << " spec Si 0.4 0.3 0.3 0.3" << endl;
  if (FindType("Fe",false)>=0)  of << " spec Fe 0.2 0.3 0.3 0.3" << endl;

  if (FindType("Mg",false)>=0 && FindType("O",false)>=0) of << "bonds Mg O 0.0 2.0 0.03 white" << endl;
  if (FindType("Si",false)>=0 && FindType("O",false)>=0) of << "bonds Si O 0.0 2.3 0.03 white" << endl;
  if (FindType("Mg",false)>=0 && FindType("Si",false)>=0) of << "bonds Mg Si 0.0 2.0 0.03 orange" << endl;
  if (FindType("C",false)>=0 && FindType("C",false)>=0) of << "bonds C C 0.0 1.7 0.05 red " << endl;

  if (FindType("O",false)>=0 && FindType("O",false)>=0) of << "bonds O O 0.0 2.0 0.03 green" << endl;
  if (FindType("Fe",false)>=0 && FindType("Fe",false)>=0) of << "bonds Fe Fe 0.0 2.0 0.03 green" << endl;
  //  if (FindType("O",false)>=0 && FindType("O",false)>=0) of << "bonds O O 0.0 1.70 0.03 green" << endl;
  //  if (FindType("O",false)>=0 && FindType("O",false)>=0) of << "bonds O O 0.0 1.60 0.03 green" << endl;

  if (FindType("N",false)>=0 && FindType("N",false)>=0) of << "bonds N N 0.0 2.0 0.03 yellow" << endl;
  if (FindType("H",false)>=0 && FindType("N",false)>=0) of << "bonds H N 0.0 1.2 0.03 white " << endl;
  if (FindType("K",false)>=0 && FindType("K",false)>=0) of << "bonds K K 0.0 3.0 0.03 orange " << endl;

  double L = Point::GetMaxDistance()*PC::AUToA*5.0;
  of << "bonds c000 c001 0.0 " << L << " 0.0 blue" << endl;
  of << "bonds c000 c010 0.0 " << L << " 0.0 blue" << endl;
  of << "bonds c000 c100 0.0 " << L << " 0.0 blue" << endl;
  of << "bonds c111 c110 0.0 " << L << " 0.0 blue" << endl;
  of << "bonds c111 c101 0.0 " << L << " 0.0 blue" << endl;
  of << "bonds c111 c011 0.0 " << L << " 0.0 blue" << endl;
  of << "bonds c001 c011 0.0 " << L << " 0.0 blue" << endl;
  of << "bonds c001 c101 0.0 " << L << " 0.0 blue" << endl;
  of << "bonds c010 c011 0.0 " << L << " 0.0 blue" << endl;
  of << "bonds c010 c110 0.0 " << L << " 0.0 blue" << endl;
  of << "bonds c100 c101 0.0 " << L << " 0.0 blue" << endl;
  of << "bonds c100 c110 0.0 " << L << " 0.0 blue" << endl;

  of <<  "inc 10" << endl;
  of <<  "scale 50" << endl;
}
