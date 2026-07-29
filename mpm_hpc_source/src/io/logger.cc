#include "logger.h"

// Create a logger for IO
const std::shared_ptr<spdlog::logger> mpm::Logger::io_logger =
    spdlog::stdout_color_st("IO");

// Create a logger for reading mesh
const std::shared_ptr<spdlog::logger> mpm::Logger::io_mesh_logger =
    spdlog::stdout_color_st("IOMesh");

// Create a logger for reading ascii mesh
const std::shared_ptr<spdlog::logger> mpm::Logger::io_mesh_ascii_logger =
    spdlog::stdout_color_st("IOMeshAscii");

// Create a logger for point generator
const std::shared_ptr<spdlog::logger> mpm::Logger::point_generator_logger =
    spdlog::stdout_color_st("PointGenerator");

// Create a logger for MPM
const std::shared_ptr<spdlog::logger> mpm::Logger::mpm_logger =
    spdlog::stdout_color_st("MPM");

// Create a logger for MPM Base
const std::shared_ptr<spdlog::logger> mpm::Logger::mpm_base_logger =
    spdlog::stdout_color_st("MPMBase");

// Create a logger for MPM Explicit
const std::shared_ptr<spdlog::logger> mpm::Logger::mpm_explicit_logger =
    spdlog::stdout_color_st("MPMExplicit");

// Create a logger for Thermo-mechanical MPM Explicit
const std::shared_ptr<spdlog::logger> mpm::Logger::thermo_mpm_explicit_logger =
    spdlog::stdout_color_st("ThermoMPMExplicit");

// Create a logger for Thermal MPM Explicit
const std::shared_ptr<spdlog::logger> mpm::Logger::thermo_mpm_implicit_logger =
    spdlog::stdout_color_st("ThermoMPMImplicit");

// Create a logger for MPM Explicit Two Phase
const std::shared_ptr<spdlog::logger>
    mpm::Logger::mpm_explicit_two_phase_logger =
        spdlog::stdout_color_st("MPMExplicitTwoPhase");

// Create a logger for Thermal MPM Explicit Two Phase
const std::shared_ptr<spdlog::logger>
    mpm::Logger::thermo_mpm_explicit_two_phase_logger =
        spdlog::stdout_color_st("ThermoMPMExplicitTwoPhase");

// Create a logger for Thermal MPM Explicit Three Phase
const std::shared_ptr<spdlog::logger>
    mpm::Logger::thermo_mpm_explicit_three_phase_logger =
        spdlog::stdout_color_st("ThermoMPMExplicitThreePhase");

// Create a logger for Thermal MPM Explicit Three Phase New
const std::shared_ptr<spdlog::logger>
    mpm::Logger::thermo_mpm_explicit_three_phase_new_logger =
        spdlog::stdout_color_st("ThermoMPMExplicitThreePhaseNew");

// Create a logger for Thermal MPM Explicit Three Phase Lag
const std::shared_ptr<spdlog::logger>
    mpm::Logger::thermo_mpm_explicit_three_phase_lag_logger =
        spdlog::stdout_color_st("ThermoMPMExplicitThreePhaseLag");

// Create a logger for Explicit THM-MPM For Saturated Frozen Soil
const std::shared_ptr<spdlog::logger>
    mpm::Logger::thm_mpm_explicit_sat_frozen_logger =
        spdlog::stdout_color_st("THMMPMExplicitSatFrozen"); 

// Create a logger for MPM Semi Implicit Two Phase
const std::shared_ptr<spdlog::logger>
    mpm::Logger::mpm_semi_implicit_two_phase_logger =
        spdlog::stdout_color_st("MPMSemiImplicitTwoPhase");

// Create a logger for Thermal MPM Semi Implicit Two Phase
const std::shared_ptr<spdlog::logger>
    mpm::Logger::thermo_mpm_semi_implicit_two_phase_logger =
        spdlog::stdout_color_st("ThermoMPMSemiImplicitTwoPhase");

// Create a logger for Semi Implicit THM-MPM For Saturated Frozen Soil
const std::shared_ptr<spdlog::logger>
    mpm::Logger::thm_mpm_semi_implicit_sat_frozen_logger =
        spdlog::stdout_color_st("THMMPMSemiImplicitSatFrozen"); 

// Create a logger for Semi Implicit THM-MPM For Unsaturated Frozen Soil
const std::shared_ptr<spdlog::logger>
    mpm::Logger::thm_mpm_semi_implicit_unsat_frozen_logger =
        spdlog::stdout_color_st("THMMPMSemiImplicitUnsatFrozen"); 

// Create a logger for 
const std::shared_ptr<spdlog::logger>
    mpm::Logger::thermo_mpm_phase_change_logger =
        spdlog::stdout_color_st("ThermoMPMPhaseChange");

// Create a logger for 
const std::shared_ptr<spdlog::logger>
    mpm::Logger::thm_mpm_explicit_MHBS_logger =
        spdlog::stdout_color_st("THMMPMExplicitMHBS");

// Create a logger for 
const std::shared_ptr<spdlog::logger>
    mpm::Logger::hydro_mpm_explicit_logger =
        spdlog::stdout_color_st("HydroMPMExplicit");

// Create a logger for 
const std::shared_ptr<spdlog::logger>
    mpm::Logger::hydro_mpm_explicit_MHBS_logger =
        spdlog::stdout_color_st("HydroMPMExplicitMHBS"); 

// Create a logger for MPM Explicit USF
const std::shared_ptr<spdlog::logger> mpm::Logger::mpm_explicit_usf_logger =
    spdlog::stdout_color_st("MPMExplicitUSF");

// Create a logger for MPM Explicit USL
const std::shared_ptr<spdlog::logger> mpm::Logger::mpm_explicit_usl_logger =
    spdlog::stdout_color_st("MPMExplicitUSL");
