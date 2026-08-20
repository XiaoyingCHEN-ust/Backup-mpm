#ifndef MPM_PRESCRIBED_PRESSURE_VALIDATION_H_
#define MPM_PRESCRIBED_PRESSURE_VALIDATION_H_

#include <stdexcept>
#include <string>

#include "data_types.h"

namespace mpm::prescribed_pressure::detail {

inline void validate_write_max_step(mpm::Index max_step,
                                    mpm::Index analysis_steps) {
  if (max_step > analysis_steps)
    throw std::invalid_argument(
        "prescribed_phase_pressures.max_step (" + std::to_string(max_step) +
        ") exceeds analysis.nsteps (" + std::to_string(analysis_steps) +
        ") while writing the pressure database");
}

}  // namespace mpm::prescribed_pressure::detail

#endif  // MPM_PRESCRIBED_PRESSURE_VALIDATION_H_
