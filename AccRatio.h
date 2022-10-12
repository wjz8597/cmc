////////////////////////////////////////////////////////////////////
//
// Burkhard Militzer                                  Urbana 10-1-99
//
// Simple class to keep track of the acceptance ratios
//
#ifndef _ACCRATIO_
#define _ACCRATIO_

class AccRatio {
 public:
  int nTrial,nAccept; 
  string name,label;
 private:
  int nMovers;
  bool infoFlag;
  static int nameLengthMax;

 public:
  AccRatio():nTrial(0),nAccept(0),infoFlag(false){};
  AccRatio(const int i):nTrial(i),nAccept(i){};
  AccRatio& operator=(const int x) { nTrial=x; nAccept=x; return *this; }
  AccRatio& operator=(const AccRatio& a) { nTrial=a.nTrial; nAccept=a.nAccept; return *this; }

  void Try() { nTrial++; }
  void Accept() { nAccept++; } 

  void Reset() { nTrial=0; nAccept=0; }
  double Ratio() const { return (nTrial>0) ? double(nAccept)/double(nTrial) : 0.0; }
  bool Flag() const { return nTrial>0; }

  void SetNameLabelMovers(const string argName, const string argLabel, const int argN=-1) {
    if (infoFlag) error("Cannot set accRatio info twice");
    infoFlag= true;
    name    = argName;
    nMovers = argN;
    label   = argLabel;
    nameLengthMax = max(name.length(),nameLengthMax);
  }
  void SetNameLabelMovers(const string argLabel, const int argN=-1) {
    SetNameLabelMovers("Acceptances type",argLabel,argN);
  }
  friend ostream& operator<<(ostream &os, const AccRatio & ar );
};


#endif // _ACCRATIO_
