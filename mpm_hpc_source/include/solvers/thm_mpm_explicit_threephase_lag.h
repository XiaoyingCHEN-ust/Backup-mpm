#ifndef MPM_THERMO_MPM_EXPLICIT_THREEPHASE_LAG_H_
#define MPM_THERMO_MPM_EXPLICIT_THREEPHASE_LAG_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "solvers/mpm_base.h"
#include "solvers/prescribed_pressure_validation.h"
#include "pipeline/rigid_pipeline2d.h"
#ifdef USE_VTK
#include "vtk_writer.h"
#endif

namespace mpm {

//! ThermoMPMExplicitThreePhaseLag class
//! \brief A class that implements the fully explicit THREE phase mpm
//! \details A THREE-phase explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ThermoMPMExplicitThreePhaseLag : public MPMBase<Tdim> {
 public:
  //! Default constructor
  ThermoMPMExplicitThreePhaseLag(const std::shared_ptr<IO>& io);

  //! Solve
  bool solve() override;

  //! Compute stress strain
  void compute_stress_strain();

  //! Pressure smoothing
  //! \param[in] phase Phase to smooth pressure
  void pressure_smoothing(unsigned phase) override;

  // Compute time step size
  void compute_critical_timestep_size(double dt); 

  //! Initialise externally prescribed phase pressures
  bool initialise_prescribed_phase_pressures();

  //! Assign prescribed phase pressures for current time
  void assign_prescribed_phase_pressures();

#ifdef USE_VTK
  //! Prepare lag-only seepage data for VTK output
  void prepare_seepage_velocity_vtk();

  //! Re-write particle VTK output with lag-only seepage data
  void write_particle_vtk_with_seepage_velocity(mpm::Index step,
                                                mpm::Index max_steps);
#endif

 protected:
  // Generate a unique id for the analysis
  using mpm::MPMBase<Tdim>::uuid_;
  //! Time step size
  using mpm::MPMBase<Tdim>::dt_;
  //! Current step
  using mpm::MPMBase<Tdim>::step_;
  //! Number of steps
  using mpm::MPMBase<Tdim>::nsteps_;
  //! Output steps
  using mpm::MPMBase<Tdim>::output_steps_;
  //! A unique ptr to IO object
  using mpm::MPMBase<Tdim>::io_;
  //! JSON analysis object
  using mpm::MPMBase<Tdim>::analysis_;
  //! JSON post-process object
  using mpm::MPMBase<Tdim>::post_process_;
  //! Logger
  using mpm::MPMBase<Tdim>::console_;
  //! Stress update
  using mpm::MPMBase<Tdim>::stress_update_;
  //! pic value
  using mpm::MPMBase<Tdim>::pic_;
  //! PIC value
  using mpm::MPMBase<Tdim>::pic_t_; 
  //! Gravity
  using mpm::MPMBase<Tdim>::gravity_;
  //! Mesh object
  using mpm::MPMBase<Tdim>::mesh_;
  //! Materials
  using mpm::MPMBase<Tdim>::materials_;
  //! VTK attributes
  using mpm::MPMBase<Tdim>::vtk_attributes_;
  //! Write VTK
  using mpm::MPMBase<Tdim>::write_vtk_;
  //! Write hdf5
  using mpm::MPMBase<Tdim>::write_hdf5_;   
  //! Damping factor
  using mpm::MPMBase<Tdim>::damping_factor_;

  // Free surface detection
  std::string free_surface_particle_{"detect"};
  //! Volume tolerance for free surface
  double volume_tolerance_{0.25};

  // Time step matrix
  using mpm::MPMBase<Tdim>::dt_matrix_;
  // Steps
  using mpm::MPMBase<Tdim>::nsteps_matrix_;
  // Time step matrix size
  using mpm::MPMBase<Tdim>::dt_matrix_size;
  //! Interface
  bool interface_{false};
  // Current time
  double current_time_{0};
    // Output number
  int No_output{0};

 private:
  //! Pressure smoothing
  bool pressure_smoothing_{false};
  //! Apply pressure smoothing inside each time step and write it back to state
  bool pressure_smoothing_in_loop_{false};
  //! Number of PIC pressure smoothing passes for VTK diagnostics
  unsigned pressure_smoothing_iterations_{1};
    //! Variable timestep
  bool variable_timestep_{false};
  //! Log output steps
  bool log_output_steps_{false};
  // DEBUG
  bool debug_{false};

  struct PressureSample {
    mpm::Index id{std::numeric_limits<mpm::Index>::max()};
    Eigen::Matrix<double, Tdim, 1> coordinates =
        Eigen::Matrix<double, Tdim, 1>::Zero();
    double liquid_pressure{0.0};
    double gas_pressure{0.0};
    bool has_coordinates{false};
  };

  struct PressureFrame {
    mpm::Index step{0};
    std::vector<double> liquid_pressures;
    std::vector<double> gas_pressures;
  };

  bool prescribed_phase_pressures_{false};
  bool write_prescribed_phase_pressures_{false};
  std::string prescribed_pressure_path_;
  std::string prescribed_pressure_prefix_{"pressure"};
  std::string prescribed_pressure_mapping_{"id"};
  double prescribed_pressure_source_dt_{0.0};
  double prescribed_pressure_time_interval_{0.0};
  mpm::Index prescribed_pressure_step_interval_{0};
  mpm::Index prescribed_pressure_max_step_{0};
  double prescribed_pressure_coordinate_bin_size_{0.02};
  std::size_t prescribed_pressure_particle_count_{0};
  //! V1 databases stored source frame zero after the first source update.
  bool prescribed_pressure_legacy_step_offset_{false};
  std::vector<PressureSample> prescribed_pressure_samples_;
  std::unordered_map<mpm::Index, std::size_t> prescribed_pressure_id_to_sample_;
  //! Stable dynamic-particle id to pressure-database sample binding. This is
  //! built once after checkpoint resume so large deformation cannot change the
  //! pressure history assigned to a material point.
  std::unordered_map<mpm::Index, std::size_t>
      prescribed_pressure_particle_to_sample_;
  std::unordered_map<long long, std::vector<std::size_t>>
      prescribed_pressure_coordinate_bins_;
  Eigen::Matrix<double, Tdim, 1> prescribed_pressure_min_coordinates_ =
      Eigen::Matrix<double, Tdim, 1>::Zero();
  bool prescribed_pressure_has_spatial_index_{false};
  //! Two-frame LRU cache used by temporal interpolation.
  std::array<PressureFrame, 2> prescribed_pressure_frame_cache_;
  std::array<bool, 2> prescribed_pressure_frame_cache_valid_{{false, false}};
  std::array<std::uint64_t, 2> prescribed_pressure_frame_cache_use_{{0, 0}};
  std::uint64_t prescribed_pressure_frame_cache_counter_{0};

  std::chrono::time_point<std::chrono::steady_clock> solver_begin;

  std::string prescribed_pressure_points_filename() const;
  std::string prescribed_pressure_values_filename() const;
  std::streamoff prescribed_pressure_frame_offset(mpm::Index source_step) const;
  void write_prescribed_pressure_points();
  void write_prescribed_pressure_binary_header();
  void read_prescribed_pressure_points();
  void validate_prescribed_pressure_binary_header();
  void write_prescribed_pressure_frame(mpm::Index source_step);
  const PressureFrame& prescribed_pressure_frame(mpm::Index source_step);
  PressureFrame read_prescribed_pressure_frame(mpm::Index source_step);
  void build_prescribed_pressure_spatial_index();
  void bind_prescribed_pressure_samples();
  bool initial_pressure_sample_index(
      mpm::Index particle_id,
      const Eigen::Matrix<double, Tdim, 1>& particle_coordinates,
      std::size_t* sample_index) const;
  bool pressure_sample_index(
      const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle,
      std::size_t* sample_index) const;
  bool pressure_from_frame(
      const PressureFrame& frame,
      const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle,
      double* liquid_pressure, double* gas_pressure) const;
  long long pressure_bin_key(const Eigen::Matrix<double, Tdim, 1>& coordinates,
                             const Eigen::Matrix<double, Tdim, 1>& min_coordinates,
                             double bin_size) const;

  //! Configure the optional circular rigid pipeline from the analysis JSON.
  void initialise_rigid_pipeline(int mpi_rank);

  //! Map pipe contact to soil nodes and return the equal-and-opposite reaction.
  RigidCircleContactResult<Tdim> map_rigid_pipeline_contact();

  //! Write pipe motion/reaction history and its current circular VTP geometry.
  void write_rigid_pipeline_output(
      const RigidCircleContactResult<Tdim>& contact,
      const Eigen::Vector2d& total_force, double total_moment, int mpi_rank);

  bool rigid_pipeline_enabled_{false};
  bool rigid_pipeline_fixed_{true};
  double rigid_pipeline_release_time_{0.};
  bool rigid_pipeline_allow_rotation_{true};
  double rigid_pipeline_particle_radius_{0.};
  double rigid_pipeline_normal_penalty_{0.};
  double rigid_pipeline_normal_damping_{0.};
  double rigid_pipeline_tangential_damping_{0.};
  double rigid_pipeline_friction_{0.};
  double rigid_pipeline_fluid_density_{0.};
  double rigid_pipeline_translational_damping_{0.};
  double rigid_pipeline_rotational_damping_{0.};
  mpm::Index rigid_pipeline_history_interval_{1};
  std::unique_ptr<mpm::pipeline::RigidPipeline2D> rigid_pipeline_;
  std::ofstream rigid_pipeline_history_;
};
}  // namespace mpm

#include "thm_mpm_explicit_threephase_lag.tcc"

#endif  // MPM_THERMO_MPM_EXPLICIT_THREEPHASE_LAG_H_
