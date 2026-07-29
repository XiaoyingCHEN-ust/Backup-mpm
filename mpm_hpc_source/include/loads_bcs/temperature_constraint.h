#ifndef MPM_TEMPERATURE_CONSTRAINT_H_
#define MPM_TEMPERATURE_CONSTRAINT_H_

#include "function_base.h"

namespace mpm {

//! TemperatureConstraint class to store temperature constraint on a set
//! \brief TemperatureConstraint class to store a constraint on a set
//! \details TemperatureConstraint stores the constraint as a static value
class TemperatureConstraint {
 public:
  // Constructor
  //! \param[in] setid  set id
  //! \param[in] temperature Constraint  temperature
  TemperatureConstraint(int setid, double temperature,
                     const std::shared_ptr<FunctionBase>& Tfunction)
      : setid_{setid}, temperature_{temperature}, Tfunction_{Tfunction} {};

  // Set id
  int setid() const { return setid_; }
  // Return temperature
  double temperature(double current_time = -1) const {
    const double scalar = ((Tfunction_ != nullptr) & (current_time >= 0))
                              ? Tfunction_->value(current_time)
                              : 1.0;

    return temperature_ * scalar;
  }

 private:
  // ID
  int setid_;
  // Temperature
  double temperature_;
  // Mathematic function for temperature;
  std::shared_ptr<FunctionBase> Tfunction_;
};  // namespace mpm
}  // namespace mpm
#endif  // MPM_TEMPERATURE_CONSTRAINT_H_
