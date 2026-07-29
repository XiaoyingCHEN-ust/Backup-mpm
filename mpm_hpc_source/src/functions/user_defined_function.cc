#include "user_defined_function.h"
// Constructor
mpm::UserDefinedFunction::UserDefinedFunction(unsigned id, const Json& properties)
    : mpm::FunctionBase(id, properties) {
  a_ = properties.at("a");
  b_ = properties.at("b");
  c_ = properties.at("c");
  d_ = properties.at("d");  
}

// Return f(x) for a given x
double mpm::UserDefinedFunction::value(double x_input) const {
  double f = 0;
  if (x_input <= a_) f = b_;
  if (x_input > a_) f = std::sin(c_ * (x_input - d_));
  return f;
}
