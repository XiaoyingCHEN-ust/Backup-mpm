#ifndef MPM_HEAT_SOURCE_H_
#define MPM_HEAT_SOURCE_H_

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "function_base.h"

namespace mpm {

//! Heat source class to store the load on a set
//! \brief Heat_source class to store the load on a set
//! \details Heat_source stores the load on a set using mathematical functions, so
//! the load can vary dynamically with time
class Heat_source {
 public:
  // Constructor
  //! \param[setid] setid  set id
  //! \param[in] mfunction Math function if defined
  //! Heat source
  Heat_source(int setid, const std::shared_ptr<mpm::FunctionBase>& heat_source_fn, double heat_source)
      : setid_{setid},
        heat_source_fn_{heat_source_fn},
        heat_source_{heat_source} {};

  // Set id
  int setid() const { return setid_; }

  // Return heat source
  double heat_source(double current_time) const {
    // Static load when no math function is defined
    double scalar = (this->heat_source_fn_ != nullptr)
                        ? (this->heat_source_fn_)->value(current_time)
                        : 1.0;
    return heat_source_ * scalar;
  }

 private:
  // ID
  int setid_;
  // Math function
  std::shared_ptr<mpm::FunctionBase> heat_source_fn_;
  // Heat source
  double heat_source_;
};
}  // namespace mpm
#endif  // MPM_HEAT_SOURCE_H_
