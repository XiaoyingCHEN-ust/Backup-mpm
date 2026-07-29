#include "sin_function.h"
// Constructor
mpm::SinFunction::SinFunction(unsigned id, const Json& properties)
    : mpm::FunctionBase(id, properties) {
  x0_ = properties.at("x0");
  a_ = properties.at("a");
  b_ = properties.at("b");
  c_ = properties.at("c");  
}

// Return f(x) for a given x
double mpm::SinFunction::value(double x_input) const {

  double f = 0;
  if (x_input <= b_) f = std::sin(a_ * (x_input - x0_));
  if (x_input > b_) f = c_;
  return f;  
  // return std::sin(a_ * (x_input - x0_));
}
