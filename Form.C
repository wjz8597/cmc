////////////////////////////////////////////////////////////////////////////////////
//                                                                                //
// Manipulate the floating point output temporarily                               //
//                                                                                //
// B. Militzer                                        Livermore 5-30-01           //
//                                                                                //
////////////////////////////////////////////////////////////////////////////////////

#include "Form.h"

Form gen(ios::fmtflags(0));
Form fix(ios::fixed);
Form sci(ios::scientific);

/*
int main() {
  //  double a=12345.678912345678912345;
  double a=20.3277111;
  int i=102;

  cout << a << endl;
  // simple way: use predefined formats gen, fix and sci
  // gen - precision = total nubmer of digits
  // fix - precision = total nubmer of digits after decimal point
  // sci - precision = total nubmer of digits after decimal point
  cout << gen(a) << "#" << endl;    // standand way with flexible nubmer of digits
  cout << fix(a) << "#" << endl;    // fixed size with default precision = 6
  cout << sci(a) << "#" << endl;    // scientific format with default precision = 6

  cout << gen(3,a) << "#" << endl;  // precision=3, total digits=3
  cout << fix(3,a) << "#" << endl;  // precision=3, digits after decimal point=3
  cout << sci(3,a) << "#" << endl;  // precision=3, digits after decimal point=3
  cout << gen(15,3,a) << "#" << endl; // as above but right bound with width=15
  cout << fix(15,3,a) << "#" << endl; // as above but right bound with width=15
  cout << sci(15,3,a) << "#" << endl; // as above but right bound with width=15

  cout << "------ AA ------" << endl;

  cout << i << endl;
  cout << gen(i) << "#" << endl;    // just print the integer
  cout << fix(i) << "#" << endl;    // same as above
  cout << sci(i) << "#" << endl;    // same as above
  cout << gen(10,i) << "#" << endl; // as above but right bound with width=10
  cout << fix(10,i) << "#" << endl; // as above but right bound with width=10
  cout << sci(10,i) << "#" << endl; // as above but right bound with width=10

  cout << "------ BB ------" << endl;

  // Define the formats to be used repeatedly.
  // print 2 digits (total, not just after the floating point)
  Form gen2(gen,2);

  // print 2 total digits, width=10
  Form gen10_2(gen,2,10);

  // print with precision=2
  Form fix2(fix,2);
  // print with precision=2, width=10
  Form fix10_2(fix,2,10);

  // print in scientic format with precision=4
  Form sci4(sci,4);
  // print in scientic format with precision=4
  Form sci8(sci,8);

  cout << a << endl;
  cout << gen2(a) << endl;
  cout << gen10_2(a) << endl;
  cout << fix2(a) << endl;
  cout << fix10_2(a) << endl;
  cout << sci4(a) << endl;
  cout << sci8(a) << endl;

  cout << endl;

  cout << "----- CC -------" << endl;

  cout << i << "#" << endl;
  cout << gen2(i) << "#" << endl;
  cout << gen10_2(i) << "#" << endl;
  cout << fix2(i)  << "#" << endl;
  cout << fix10_2(i) << "#" << endl;
  cout << sci4(i) << "#" << endl;
  cout << sci8(i) << "#" << endl;

  i = 2;
  cout << i  << "#" << endl;
  cout << gen2(i) << "#" << endl;
  cout << gen10_2(i) << "#" << endl;
  cout << fix2(i)  << "#" << endl;
  cout << fix10_2(i) << "#" << endl;
  cout << sci4(i) << "#" << endl;
  cout << sci8(i) << "#" << endl;

  cout << "----- DD -------" << endl;

  string s = "abc";
  Form sl(gen,5,5);
  sl.left();
  Form sr(gen,5,5);
  sr.right();
  cout << sl(s) << "#" << endl;
  cout << sr(s) << "#" << endl;
  cout << sl(5) << "# <<< left adjustment of numbers does not work -- not implemented?" << endl;
  cout << sr(5) << "#" << endl;

  cout << sr(LEFT,s) << "#" << endl;
  cout << sr(RIGHT,s) << "#" << endl;
  cout << sr(CENTER,s) << "#" << endl;
}
*/
