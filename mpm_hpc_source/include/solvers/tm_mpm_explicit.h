#ifndef MPM_THERMO_MPM_EXPLICIT_H_
#define MPM_THERMO_MPM_EXPLICIT_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "solvers/mpm_base.h"
#include <Eigen/Dense>
namespace mpm {

//! ThermoMPMExplicit class
//! \brief A class that implements the fully explicit one phase thermal-mechanical mpm
//! \details A single-phase thermal-mechanical explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ThermoMPMExplicit : public MPMBase<Tdim> {
 public:
  //! Default constructor
  ThermoMPMExplicit(const std::shared_ptr<IO>& io);

  //! Domain decomposition
  void mpi_domain_decompose();

  //! Solve
  bool solve() override;

  //! Pressure smoothing
  //! \param[in] phase Phase to smooth pressure
  void pressure_smoothing(unsigned phase) override;

  //! Compute stress strain
  //! \param[in] phase Phase to smooth pressure
  void compute_stress_strain(unsigned phase);
  
  //! Apply temperature constraints
  void apply_nodal_temperature_constraints(unsigned phase);

  // Compute time step size
  void compute_critical_timestep_size(double dt);  

  //! Pre process for MPM-DEM
  bool pre_process() override;

  //! Get deformation gradient for MPM-DEM
  bool get_deformation_task() override;

  //! Get analysis information
  void get_info(unsigned& dim, bool& resume,
                unsigned& checkpoint_step) override;

  //! get step, time
  void get_status(double& dt, unsigned& step, unsigned& nsteps,
                  unsigned& output_steps) override;

  //! send deformation task
  bool send_deformation_task(
      std::vector<unsigned>& id,
      std::vector<Eigen::MatrixXd>& displacement_gradients) override;

  //! send temperature task
  bool send_temperature_task(
      std::vector<unsigned>& id,
      std::vector<double>& particle_temperature) override;    

  //! Set particle stress
  bool set_stress_task(const Eigen::MatrixXd& stresses,
                       bool increment) override;

  //! Set particle porosity
  bool set_porosity_task(const Eigen::MatrixXd& porosities) override;

  // Set particle fabric
  bool set_fabric_task(std::string fabric_type,
                       const Eigen::MatrixXd& fabrics) override;

  // Set particle rotation
  bool set_rotation_task(const Eigen::MatrixXd& rotations) override;

  // Update particle state eg position, velocity
  bool update_state_task() override;

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
  // Free surface detection
  std::string free_surface_particle_{"detect"};
  //! Volume tolerance for free surface
  double volume_tolerance_{0.25}; 
  //! Virtual flux
  bool virtual_flux_{false};
  // Flux type
  std::string flux_type_{"convective"};
  // Heat transfer coeff
  double heat_transfer_coeff_{0};
  // Ambient temperature
  double ambient_temperature_{0};
  // COnductive flux
  double flux_{0};

};  // ThermoMPMExplicit class
}  // namespace mpm

#include "tm_mpm_explicit.tcc"

#endif  // MPM_THERMO_MPM_EXPLICIT_H_
