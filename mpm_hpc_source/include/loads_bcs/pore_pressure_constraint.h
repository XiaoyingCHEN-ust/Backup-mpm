#ifndef MPM_PORE_PRESSURE_CONSTRAINT_H_
#define MPM_PORE_PRESSURE_CONSTRAINT_H_

#include "function_base.h"

namespace mpm {

//! PorepressureConstraint class to store pore pressure constraint on a set
//! \brief PorepressureConstraint class to store a constraint on a set
//! \details PorepressureConstraint stores the constraint as a static value
class PorepressureConstraint {
 public:
  // Constructor
  //! \param[in] setid  set id
  //! \param[in] pore_pressure Constraint  pore_pressure
  PorepressureConstraint(int setid, double pore_pressure,
                     const std::shared_ptr<FunctionBase>& Pfunction)
      : setid_{setid}, pore_pressure_{pore_pressure}, Pfunction_{Pfunction} {};

  // Set id
  int setid() const { return setid_; }
  // Return pore_pressure
  double pore_pressure(double current_time = -1) const {
    const double scalar = ((Pfunction_ != nullptr) & (current_time >= 0))
                              ? Pfunction_->value(current_time)
                              : 1.0;

    return pore_pressure_ * scalar;
  }

 private:
  // ID
  int setid_;
  // Pore pressure
  double pore_pressure_;
  // Mathematic function for pore pressure;
  std::shared_ptr<FunctionBase> Pfunction_;
};  // namespace mpm
}  // namespace mpm
#endif  // MPM_PORE_PRESSURE_CONSTRAINT_H_
