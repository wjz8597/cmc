/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Classical MC+MD program                                             // 
//                                                                     //
// Burkhard Militzer                              Berkeley, 07-01-2009 //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#include "State.h"
#include "Sort.h"

void writePoint(ostream & of, const Form & f, const Point & p, const string s="") {
  of << s;
  for(int d=0; d<Point::nDim; d++) {
    of << f(p[d]);
    if (d<Point::nDim-1) of << " ";
  }
}

void WritePoint(ostream & of, const Form & f, const Point & p, const string s="") {
  writePoint(of,f,p,s);
  of << endl;
}

void ComputeVelocities(Array1 <BPoint> & v, const State & s1, const State & s2, const double dt) {
  for(int i=0; i<s1.nParticles; i++) {
    v[i] = (s2.r[i]-s1.r[i])/dt;
  }
}

// Remember VASP wants to center-of-mass motion subtracted out!!!
void WriteVASPPOSCARFile(ostream & of, const State & s, const State & sp, const double dt, const Array1 <int> & typeOrder) {
  int nt = s.GetNTypes();
  if (typeOrder.Size()!=nt) error("Type list error I");
  Array1 <int> tmp=typeOrder;
  SortArray(tmp);
  for(int t=0; t<nt; t++) 
    if (tmp[t]!=t) error("Type list error II",t,tmp[t]);
  //  Write(typeOrder);
  //  Print();

  Form f(fix,21,16);
  Form ff(sci,8);
  ff.SetShowPos();

  // My old way of using the first comment line
  //  of << " VASP POSCAR/CONTCAR file created by \'ta\' program written by B. Militzer " << endl;
  // Stephen's way recording the type there
  for(int t=0; t<nt; t++) {
    int tt = typeOrder[t];
    of << " " << s.listOfTypes[tt];
  }
  of << endl;

  double scale=1.0;
  if (Point::Cubic()) scale = Point::GetIsoBoxSize()*PC::AUToA;
  of << "  " << f(scale) << endl;
  WritePoint(of,f,Point::GetLatticeVector(0)*PC::AUToA/scale,"   ");
  WritePoint(of,f,Point::GetLatticeVector(1)*PC::AUToA/scale,"   ");
  WritePoint(of,f,Point::GetLatticeVector(2)*PC::AUToA/scale,"   ");
  for(int t=0; t<s.GetNTypes(); t++) {
    int tt = typeOrder[t];
    of << " " << s.GetNAtomsPerType(tt) << " ";
  }
  of << endl;

  of << "Direct" << endl; // reduced coordinates
  for(int t=0; t<s.listOfTypes.Size(); t++) {
    int tt = typeOrder[t];
    for(int it=0; it<s.atomsByType[tt].Size(); it++) {
      int i = s.atomsByType[tt][it];
      Point r = s.r[i];
      r.MapInside(); // do not forget this because MD takes atoms outside (unrolles traj)
      r.Reduce(); // VASP wants reduced coord!
      WritePoint(of,f,r);
    }
  }

  int n = s.nParticles;
  Array1 <BPoint> v(n);
  ComputeVelocities(v,s,sp,dt);

  Point vCMS = 0.0;
  double v2Sum = 0.0;
  double m = 0.0;
  for(int i=0; i<n; i++) {
    m += s.mass[i];
    vCMS  += v[i]*s.mass[i];
    v2Sum += v[i].Norm();
  }
  vCMS /= m;

  if (v2Sum==0.0) return; // skip writing velocities if case all are zero!

  //  warning("Do not subtract out v_CMS");
  //  vCMS = 0.0;

  //  of << endl; // VASP version 4.6: an empty lines mean cartesian velocities with CMS motion in Ang/fs - could change in new version
//  of << "Cartesian" << endl; // just to make sure

//  for(int t=0; t<s.listOfTypes.Size(); t++) {
//    int tt = typeOrder[t];
//    for(int it=0; it<s.atomsByType[tt].Size(); it++) {
//      int i = s.atomsByType[tt][it];
      //    WritePoint(of,ff,v[i].Reduced()*dt," ");
//      WritePoint(of,ff,(v[i]-vCMS)*PC::AUToA/PC::AUToFS," ");
//    }
//  }
}

void WriteVASPPOSCARFile(const string & fn, const State & s, const State & sp, const double dt, const Array1 <int> & typeOrder) {
  ofstream of;
  of.open(fn.c_str());
  if (!of) error("Could not open file",fn);
  WriteVASPPOSCARFile(of,s,sp,dt,typeOrder);
  cout << "VASP POSCAR file written with name: " << fn << endl;
  of.close();
}

void WriteVASPPOSCARFile(const string & fn, const State & s, const State & sp, const double dt) {
  WriteVASPPOSCARFile(fn, s, sp, dt, s.GetDefaultTypeArray());
}

void WriteVASPPOSCARFile(const string & fn, const State & s) {
  WriteVASPPOSCARFile(fn, s, s, 1.0);
}

