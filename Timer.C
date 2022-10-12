////////////////////////////////////////////////////////////////////////////////
//
// Measure the time like using a stop watch
//
// Burkhard Militzer                                     Livermore, 5-10-01
//       
////////////////////////////////////////////////////////////////////////////////

#include <unistd.h> 
#include <sys/times.h>
#include <math.h>

#include <iomanip>
using namespace std;

#include "Standard.h"
#include "Timer.h"

Timer::Timer(const mode & ma, const bool flag):
  m(ma),
  f(1.0/(double)(GetHertz())),
  runTime(0.0),
  running(flag)
{
  zone.tz_minuteswest=0;
  if (flag) Start();
}

void Timer::Start() {
  time0   = GetTime();
  running = true;
}

void Timer::Reset() {
  running = false;
  runTime = 0.0;
}

long Timer::GetHertz() const {
  return sysconf (_SC_CLK_TCK) ;
}

double Timer::GetUserTime() {
  times(&buf);
  return
    (double) (buf.tms_utime)*f;
}

double Timer::GetSystemTime() {
  times(&buf);
  return
    (double) (buf.tms_stime)*f;
}

double Timer::GetTotalTime() {
  gettimeofday(&val,&zone);
  return (double)val.tv_sec + 1e-6*(double) val.tv_usec;

  //  times(&buf);
  //  return (double) (buf.tms_utime+buf.tms_stime)*f;
}

double Timer::GetTime() {
  switch (m) {
  case user:
    return GetUserTime();
  case system:
    return GetSystemTime();
  case total:
    return GetTotalTime();
  default: // case "off"
    return 0.0;
  }
}

double Timer::Read() {
  return (running) ? (GetTime()-time0) : runTime;
}

double Timer::ReadAndRestart() {
  double t =Read();
  Start();
  return t;
}

double Timer::Stop() {
  if (!running)
    error("Timer::Stop: Timer already stopped.");
  double t =GetTime();
  runTime = t-time0;
  running = false;
  return runTime;
}

void Timer::Continue() {
  if (running) 
    error("Cannot continue a timer before having stopped it.");
  double t = GetTime();
  time0    = t-runTime;
  running  = true;
  runTime  = 0.0;
}

double Timer::Fraction(const double t) {
  if (t==0.0) return 0.0;
  return Read()/t;
}

double Timer::Fraction(Timer & t) {
  return Fraction(t.Read());
}

ostream& operator<<(ostream &os, Timer & t ) {
  if (t.m!=Timer::off) {
    int p=int(log10( (double)(t.GetHertz()) ))+1;
    double tt= max(0.0,t.Read());
    const ios::fmtflags current = os.flags();
    const streamsize ss = os.precision();
    os.precision(p);
    os.setf(ios::fixed,ios::floatfield);
    os << tt;
    os.setf(current);
    os.precision(ss);
  }
  return os;
}
