////////////////////////////////////////////////////////////////////////////////
//
// Measure the time like using a stop watch
//
// Burkhard Militzer                                     Livermore, 5-10-01
//       
////////////////////////////////////////////////////////////////////////////////

#include <sys/times.h>
#include <sys/time.h>
#include <iostream>

#ifndef _TIMER_
#define _TIMER_

class Timer {
 public:
  enum mode {user, system, total, off};
  // one might want to change this back later so that Array1 get initialized with Timer::user
  //  Timer(const mode & ma=user, const bool flag=true);
  Timer(const mode & ma=off, const bool flag=true);
  void Start();
  double Read();
  double ReadAndRestart();
  friend ostream& operator<<(ostream &os, Timer & t );
  double Stop();
  void Continue();
  void Reset();
  mode GetMode() const {
    return m;
  }
  void SetMode(const mode ma) {
    m = ma;
  }
  double Fraction(const double t);
  double Fraction(Timer & t);
 private:
  long GetHertz() const;
  double GetUserTime();
  double GetSystemTime();
  double GetTotalTime();
  double GetTime();

  struct tms buf;
  mode m;
  const double f;
  double time0,runTime;
  bool running;
  
  struct timeval val;
  struct timezone zone;
};

#endif // _TIMER_
