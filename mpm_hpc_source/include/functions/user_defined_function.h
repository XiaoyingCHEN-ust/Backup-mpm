#ifndef MPM_USER_DEFINED_FUNCTION_H_
#define MPM_USER_DEFINED_FUNCTION_H_

#include <cmath>

#include "function_base.h"

namespace mpm {

//! UserDefinedFunction class
//! \brief class that computes a mathematical UserDefinedFunction
//! \details UserDefinedFunction computes the value of a UserDefinedFunction
class UserDefinedFunction : public FunctionBase {
 public:
  // Construct a UserDefinedFunction with a unique id
  //! \param[in] id Global id
  //! \param[in] json object of function properties
  UserDefinedFunction(unsigned id, const Json& properties);

  //! Default destructor
  ~UserDefinedFunction() override = default;

  //! Return the value of the UserDefinedFunction at given input
  double value(double x) const override;

 private:
  //! function id
  using FunctionBase::id_;
  double a_{0.};
  double b_{0.};
  double c_{0.};
  double d_{0.};  
};  // UserDefinedFunction class
}  // namespace mpm

#endif  //MPM_UDER_DEFINED_FUNCTION_H_
