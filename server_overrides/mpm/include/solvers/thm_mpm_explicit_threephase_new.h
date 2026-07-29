#ifndef MPM_THERMO_MPM_EXPLICIT_THREEPHASE_NEW_H_
#define MPM_THERMO_MPM_EXPLICIT_THREEPHASE_NEW_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include <Eigen/SparseCholesky>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <vector>

#include "particles/threephase_pressure_state.h"
#include "solvers/mpm_base.h"

namespace mpm {

//! ThermoMPMExplicitThreePhaseNew class
//! \brief A class that implements the fully explicit TWO phase mpm
//! \details A TWO-phase explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ThermoMPMExplicitThreePhaseNew : public MPMBase<Tdim> {
 public:
  //! Default constructor
  ThermoMPMExplicitThreePhaseNew(const std::shared_ptr<IO>& io);

  //! Solve
  bool solve() override;

  //! Compute stress strain
  void compute_stress_strain();

  //! Pressure smoothing
  //! \param[in] phase Phase to smooth pressure
  void pressure_smoothing(unsigned phase) override;

  // Compute time step size
  void compute_critical_timestep_size(double dt);

  //! Solve the coupled liquid/gas pressure equations by backward Euler.
  bool solve_semi_implicit_pressure(double dt);

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
    //! Variable timestep
  bool variable_timestep_{false};
  //! Log output steps
  bool log_output_steps_{false};
  //! Use the global semi-implicit liquid/gas pressure solve.
  bool semi_implicit_pressure_{false};
  //! Weak particle-pressure boundary penalty used by the global solve.
  double pressure_boundary_penalty_{1.0e6};
  //! Log the implicit pressure residual at each output step.
  bool log_pressure_solver_{false};
  //! Reconstruct gradients from the current nodal pressure solution.
  bool reconstruct_pressure_gradient_{false};
  //! Reconstruct only the pressure used by the phase momentum internal force.
  bool reconstruct_pressure_force_{false};
  //! Reconstruct phase velocities from the Darcy relation instead of FLIP.
  bool reconstruct_darcy_velocity_{false};
  //! Apply the diagnostic local nodal bound to transferred pressures.
  bool bounded_pressure_transfer_{false};
  // DEBUG
  bool debug_{false};

  std::chrono::time_point<std::chrono::steady_clock> solver_begin;
};
}  // namespace mpm

#include "thm_mpm_explicit_threephase_new.tcc"

#endif  // MPM_THERMO_MPM_EXPLICIT_THREEPHASE_NEW_H_
