#ifndef MPM_HYDRO_MPM_EXPLICIT_MHBS_H_
#define MPM_HYDRO_MPM_EXPLICIT_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "solvers/mpm_base.h"
#include <Eigen/Dense>
namespace mpm {

//! HydroMPMExplicitMHBS class
//! \brief A class that implements the fully explicit one phase thermal-mechanical mpm
//! \details A single-phase thermal-mechanical explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class HydroMPMExplicitMHBS : public MPMBase<Tdim> {
 public:
  //! Default constructor
  HydroMPMExplicitMHBS(const std::shared_ptr<IO>& io);

  //! Solve
  bool solve() override;

  // Compute time step size
  void compute_critical_timestep_size(double dt);  


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
  //! write hdf5
  using mpm::MPMBase<Tdim>::write_hdf5_;
  //! write vtk
  using mpm::MPMBase<Tdim>::write_vtk_;
  //! A unique ptr to IO object
  using mpm::MPMBase<Tdim>::io_;
  //! JSON analysis object
  using mpm::MPMBase<Tdim>::analysis_;
  //! JSON post-process object
  using mpm::MPMBase<Tdim>::post_process_;
  //! Logger
  using mpm::MPMBase<Tdim>::console_;

#ifdef USE_GRAPH_PARTITIONING
  //! Graph
  using mpm::MPMBase<Tdim>::graph_;
#endif

  //! PIC value
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
  //! Node concentrated force
  using mpm::MPMBase<Tdim>::set_node_concentrated_force_;
  //! Damping type
  using mpm::MPMBase<Tdim>::damping_type_;
  //! Damping factor
  using mpm::MPMBase<Tdim>::damping_factor_;

 private:
  //! Pressure smoothing
  bool pressure_smoothing_{false};
  //! Variable timestep
  bool variable_timestep_{false};
  //! Log output steps
  bool log_output_steps_{false};
  //! Interface
  bool interface_{false};
  //! Solver begin time
  std::chrono::time_point<std::chrono::steady_clock> solver_begin;

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
  // Current time
  double current_time_{0};
  // Output number
  int No_output{0};

};  // HydroMPMExplicitMHBS class
}  // namespace mpm

#include "hydro_mpm_explicit_MHBS.tcc"

#endif  // MPM_HYDRO_MPM_EXPLICIT_MHBS_H_
