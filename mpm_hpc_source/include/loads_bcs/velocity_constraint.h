#ifndef MPM_VELOCITY_CONSTRAINT_H_
#define MPM_VELOCITY_CONSTRAINT_H_

#include "function_base.h"

namespace mpm {

//! VelocityConstraint class to store velocity constraint on a set
//! \brief VelocityConstraint class to store a constraint on a set
//! \details VelocityConstraint stores the constraint as a static value
class VelocityConstraint {
 public:
  // Constructor
  //! \param[in] setid  set id
  //! \param[in] dir Direction of constraint load
  //! \param[in] velocity Constraint  velocity
  VelocityConstraint(int setid, unsigned dir, double velocity,
                     std::shared_ptr<FunctionBase>& vfunction)
      : setid_{setid}, dir_{dir}, velocity_{velocity}, vfunction_{vfunction} {};

  // Set id
  int setid() const { return setid_; }
  // Direction
  unsigned dir() const { return dir_; }
  // Return velocity
  double velocity(double current_time = -1) const {
    const double scalar = ((vfunction_ != nullptr) & (current_time >= 0))
                              ? vfunction_->value(current_time)
                              : 1.0;

    return velocity_ * scalar;
  }

 private:
  // ID
  int setid_;
  // Direction
  unsigned dir_;
  // Velocity
  double velocity_;
  // Mathematic function for velocity;
  std::shared_ptr<FunctionBase> vfunction_;
};  // namespace mpm
}  // namespace mpm
#endif  // MPM_VELOCITY_CONSTRAINT_H_
