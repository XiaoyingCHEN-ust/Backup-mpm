#ifndef MPM_LOGGER_H_
#define MPM_LOGGER_H_

#include <memory>

// Speed log
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

//! MPM namespace
namespace mpm {

// Create an stdout colour sink
const std::shared_ptr<spdlog::sinks::stdout_color_sink_mt> stdout_sink =
    std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

struct Logger {
  // Create a logger for IO
  static const std::shared_ptr<spdlog::logger> io_logger;

  // Create a logger for reading mesh
  static const std::shared_ptr<spdlog::logger> io_mesh_logger;

  // Create a logger for reading ascii mesh
  static const std::shared_ptr<spdlog::logger> io_mesh_ascii_logger;

  // Create a logger for point generator
  static const std::shared_ptr<spdlog::logger> point_generator_logger;

  // Create a logger for MPM
  static const std::shared_ptr<spdlog::logger> mpm_logger;

  // Create a logger for MPM Base
  static const std::shared_ptr<spdlog::logger> mpm_base_logger;

  // Create a logger for MPM Explicit
  static const std::shared_ptr<spdlog::logger> mpm_explicit_logger;

  // Create a logger for Thermo-mechanical MPM Explicit
  static const std::shared_ptr<spdlog::logger> thermo_mpm_explicit_logger;

  // Create a logger for Thermo-mechanical MPM Explicit
  static const std::shared_ptr<spdlog::logger> thermo_mpm_implicit_logger;  

  // Create a logger for MPM Explicit Two Phase
  static const std::shared_ptr<spdlog::logger> mpm_explicit_two_phase_logger;

  // Create a logger for Thermal MPM Explicit Two Phase
  static const std::shared_ptr<spdlog::logger> thermo_mpm_explicit_two_phase_logger;

  // Create a logger for Thermal MPM Explicit Three Phase
  static const std::shared_ptr<spdlog::logger> thermo_mpm_explicit_three_phase_logger;

  // Create a logger for Thermal MPM Explicit Three Phase New
  static const std::shared_ptr<spdlog::logger> thermo_mpm_explicit_three_phase_new_logger;

  // Create a logger for Thermal MPM Explicit Three Phase New
  static const std::shared_ptr<spdlog::logger> thermo_mpm_explicit_three_phase_lag_logger;

  // Create a logger for THM-MPM Explicit for saturated frozen soil
  static const std::shared_ptr<spdlog::logger> thm_mpm_explicit_sat_frozen_logger;

  // Create a logger for MPM Semi-implicit Two Phase
  static const std::shared_ptr<spdlog::logger> mpm_semi_implicit_two_phase_logger;

  // Create a logger for thermal MPM Semi-implicit Two Phase
  static const std::shared_ptr<spdlog::logger> thermo_mpm_semi_implicit_two_phase_logger;

  // Create a logger for THM-MPM Semi-implicit for saturated frozen soil
  static const std::shared_ptr<spdlog::logger> thm_mpm_semi_implicit_sat_frozen_logger;

  // Create a logger for THM-MPM Semi-implicit for unsaturated frozen soil
  static const std::shared_ptr<spdlog::logger> thm_mpm_semi_implicit_unsat_frozen_logger;

  // Create a logger for THM-MPM Semi-implicit for saturated frozen soil
  static const std::shared_ptr<spdlog::logger> thermo_mpm_phase_change_logger;

  // Create a logger for THM-MPM Semi-implicit for saturated frozen soil
  static const std::shared_ptr<spdlog::logger> thm_mpm_explicit_MHBS_logger;

  // Create a logger for 
  static const std::shared_ptr<spdlog::logger> hydro_mpm_explicit_logger; 

  // Create a logger for 
  static const std::shared_ptr<spdlog::logger> hydro_mpm_explicit_MHBS_logger;   

  // Create a logger for MPM Explicit USF
  static const std::shared_ptr<spdlog::logger> mpm_explicit_usf_logger;

  // Create a logger for MPM Explicit USL
  static const std::shared_ptr<spdlog::logger> mpm_explicit_usl_logger;
};

}  // namespace mpm

#endif  // MPM_LOGGER_H_
