/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Function that contains function f(x), f'(x) and flag if df/dx>0         //
//                                                                         //
// Burkhard Militzer                                  Berkeley 01-25-08    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _FUNCTION_
#define _FUNCTION_

class Function {
 public:
  double operator()(const double x) const {
    return f(x);
  }
  double f(const double x) const {
    return x;
  }
  double Inverse(const double x) const {
    return x;
  }
  double Derivative(const double x) const {
    return 1.0;
  }
  double SecondDerivative(const double x) const {
    return 0.0;
  }
  bool Rising() const {
    return true;
  }
  string Name() const {
    return "f(x) = x";
  }
};

class LogFunction : public Function {
 public:
  double operator()(const double x) const {
    return f(x);
  }
  double f(const double x) const {
    return log(x);
  }
  double Inverse(const double x) const {
    return exp(x);
  }
  double Derivative(const double x) const {
    return 1.0/x;
  }
  double SecondDerivative(const double x) const {
    return -1.0/(x*x);
  }
  bool Rising() const {
    return true;
  }
  string Name() const {
    return "f(x) = log(x)";
  }
};

class PowerFunction : public Function {
 public:
  const double alpha;
  PowerFunction(const double alpha_):alpha(alpha_) {}

  double operator()(const double x) const { // must be declared again, otherwise wrong f() called
    return f(x);
  }
  double f(const double x) const {
    return pow(x,alpha);
  }
  double Inverse(const double x) const { // ignore any negative roots
    return pow(x,1.0/alpha);
  }
  double Derivative(const double x) const {
    return alpha*pow(x,alpha-1.0);
  }
  double SecondDerivative(const double x) const {
    return alpha*(alpha-1)*pow(x,alpha-2.0);
  }
  bool Rising() const {
    return alpha>0.0;
  }
  string Name() const {
    return "f(x) = x^"+DoubleToString(alpha);
  }
};

class Power3Function : public PowerFunction {
 public:
  Power3Function():PowerFunction(3.0){}
};

class Power2Function : public PowerFunction {
 public:
  Power2Function():PowerFunction(2.0){}
};

class Power1Function : public PowerFunction {
 public:
  Power1Function():PowerFunction(1.0){}
};

#endif _FUNCTION_
