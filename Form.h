////////////////////////////////////////////////////////////////////////////////////
//                                                                                //
// Manipulate the floating point output temporarily                               //
//                                                                                //
// B. Militzer                                        Livermore 5-30-01           //
//                                                                                //
////////////////////////////////////////////////////////////////////////////////////

#ifndef _FORM_
#define _FORM_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
using namespace std;

#include "Standard.h"

class Form;

enum StringAdjustment {LEFT,RIGHT,CENTER};

class BoundForm {
 public:
  inline BoundForm(const Form & f, const double v);
  inline BoundForm(const Form & f, const int v);
  inline BoundForm(const BoundForm &);
  inline BoundForm(const Form & f, const char *ss);
  inline BoundForm(const Form & f, const string & ss);

  friend ostream& operator<<(ostream& os, const BoundForm& bf) {
    os << bf.s.str(); // s must be terminated by ends otherwise str() cannot be called!!!
    return os;
  }
  int GetLength() {
    return s.str().length();
  }
  string GetString() {
    return s.str();
  }

 private:
  ostringstream s;

};


class Form {
  friend class BoundForm;

  int prc;              // precision
  int wdt;              // width = number of fields
  ios::fmtflags fmt;   // extra flags with "+"
  StringAdjustment sa;      // LEFT,RIGHT,CENTER
  
 public:
  // different constructors to specify a certain output format
  Form(const Form & f)                            :prc(f.prc),wdt(f.wdt),fmt(f.fmt),sa(f.sa) {}
  Form(const Form & f, const int p)               :prc(p),    wdt(f.wdt),fmt(f.fmt),sa(f.sa) {}
  Form(const Form & f, const int w, const int p)  :prc(p),    wdt(w),    fmt(f.fmt),sa(f.sa) {}   // w and p were swapped before
  Form(const Form & f, const StringAdjustment sa1):prc(f.prc),wdt(f.wdt),fmt(f.fmt),sa(sa1)  {}
  explicit Form(const ios::fmtflags f=ios::fmtflags(0), 
  		const int p=6, 
		const int w=0):
    prc(p),
    wdt(w),
    fmt(f),
    sa(RIGHT)
    {}

  // different functions with format arguments 
  // all return a BoundForm which can be output to a stream
  BoundForm operator() (const double d) const {
    return BoundForm(*this,d);
  }
  BoundForm operator() (const int p, const double d) const {
    Form f(*this,p,(*this).wdt);
    return BoundForm(f,d);
  }
  BoundForm operator() (const int w, const int p, const double d) const {
    Form f(*this,w,p);
    return BoundForm(f,d);
  }
  BoundForm operator() (const int i) const {
    return BoundForm(*this,i);
  }
  BoundForm operator() (const int w, const int i) const {
    Form f(*this);
    f.width(w);
    return BoundForm(f,i);
  }
  BoundForm operator() (const StringAdjustment sa, const string & ss) const {
    Form f(*this,sa);
    return BoundForm(f,ss);
  }
  BoundForm operator() (const string & ss) const {
    return BoundForm(*this,ss);
  }
  BoundForm operator() (const int w, const string & ss) const {
    Form f(*this);
    f.width(w);
    return BoundForm(f,ss);
  }

  // All the following function modify the current output format (even though they return *this)

  void SetWidthToScientificWidth() {
    BoundForm bf(*this,(double)(1.0));
    wdt = bf.GetLength();
  }

  void SetShowPos() {
    fmt |= ios::showpos;
  }

  // function to specify the format
  Form& scientific() {
    fmt = ios::scientific;
    return *this;
  }

  Form& fixed() {
    fmt = ios::fixed;
    return *this;
  }

  Form& general() {
    fmt = ios::fmtflags(0);
    return *this;
  }

  Form& width(int w) {
    wdt = w;
    return *this;
  }

  Form& left() {
    sa=LEFT;
    return *this;
  }

  Form& right() {
    sa=RIGHT;
    return *this;
  }

  Form& center() {
    sa=CENTER;
    return *this;
  }

  Form& precision(int p) {
    prc = p;
    return *this;
  }
};

inline BoundForm::BoundForm(const Form & f, const double v) {
  s.precision(f.prc);
  s.width(f.wdt);
  //  s.setf(f.fmt,ios::floatfield);
  s.setf(f.fmt);
  s << v;
#if !defined __GNUC__ && !defined __PGI
  s << ends; // must be terminated by ends other s.str() cannot be called
#endif
}

inline BoundForm::BoundForm(const Form & f, const int v) {
  s.width(f.wdt);
  s.setf(f.fmt);
  s << v;
#if !defined __GNUC__ && !defined __PGI
  s << ends; // must be terminated by ends other s.str() cannot be called
#endif
}

inline BoundForm::BoundForm(const Form & f, const char *ss) {
  BoundForm bf(f,string(ss));
  s << bf;
#if !defined __GNUC__ && !defined __PGI
  s << ends;
#endif
}

// all string formatting is done HERE
inline BoundForm::BoundForm(const Form & f, const string & ss) {
  if (f.sa!=CENTER) s.width(f.wdt);
  s.setf(f.fmt);
  //  Write3(f.sa,f.wdt,f.fmt);
  if (f.sa==RIGHT) {
    s.setf(ios::right,ios::adjustfield);
    s << ss;
  } else if (f.sa==LEFT) {
    s.setf(ios::left,ios::adjustfield);
    s << ss;
  } else {
    int i = f.wdt-ss.length();
    int il = i/2;
    int ir = i-il;
    for (int j=0; j<il; j++) 
      s << " ";
    s << ss;
    for (int j=0; j<ir; j++) 
      s << " ";
  }
#if !defined __GNUC__ && !defined __PGI
  s << ends; // must be terminated by ends other s.str() cannot be called
#endif
}

inline BoundForm::BoundForm(const BoundForm & b) {
  s << b.s;
#if !defined __GNUC__ && !defined __PGI
  s << ends;
#endif
}

inline int NumberOfDigits(int n) {
  int i=1;
  if (n<0) {
    i++;
    n=-n;
  }
  while (n>=10) {
    n/=10;
    i++;
  }
  return i;
}

extern Form gen,fix,sci;

#endif // _FORM_
