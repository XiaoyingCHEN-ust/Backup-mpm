#ifndef MPM_THM_MPM_EXPLICIT_MHBS_H_
#define MPM_THM_MPM_EXPLICIT_MHBS_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "solvers/mpm_base.h"

namespace mpm {

//! THMMPMExplicitMHBS class
//! \brief A class that implements the fully explicit TWO phase mpm
//! \details A TWO-phase explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class THMMPMExplicitMHBS : public MPMBase<Tdim> {
 public:
  //! Default constructor
  THMMPMExplicitMHBS(const std::shared_ptr<IO>& io);

  //! Solve
  bool solve() override;

  //! Compute stress strain
  void compute_stress_strain();

  //! Pressure smoothing
  //! \param[in] phase Phase to smooth pressure
  void pressure_smoothing(unsigned phase) override;

  // Compute time step size
  void compute_critical_timestep_size(double dt); 

  bool pre_process() override;

  bool get_deformation_task() override;

  //! Get analysis information
  void get_info(unsigned& dim, bool& resume,
                unsigned& checkpoint_step) override;

  //! get step, time
  void get_status(double& dt, unsigned& step, unsigned& nsteps,
                  unsigned& output_steps) override;

  bool send_deformation_task(
      std::vector<unsigned>& id,
      std::vector<Eigen::MatrixXd>& displacement_gradients) override;

  //! send temperature task
  bool send_saturation_task(
      std::vector<unsigned>& id,
      std::vector<double>& hydrate_saturation) override;

  //! Set particle stess
  bool set_stress_task(const Eigen::MatrixXd& stresses,
                        bool increment) override;

  //! Set particle porosity
  bool set_porosity_task(const Eigen::MatrixXd& porosities) override;

  // Set particle fabric
  bool set_fabric_task(std::string fabric_type,
                        const Eigen::MatrixXd& fabrics) override;

  // Set particle rotation
  bool set_rotation_task(const Eigen::MatrixXd& rotations) override;

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
  // DEBUG
  bool debug_{false};

  std::chrono::time_point<std::chrono::steady_clock> solver_begin;
}; 
}  // namespace mpm

#include "thm_mpm_explicit_MHBS.tcc"

#endif  // MPM_THERMO_MPM_EXPLICIT_TWOPHASE_H_
