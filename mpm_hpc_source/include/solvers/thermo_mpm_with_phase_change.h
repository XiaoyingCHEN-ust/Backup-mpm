#ifndef MPM_THERMO_MPM_PHASE_CHANGE_H_
#define MPM_THERMO_MPM_PHASE_CHANGE_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include <iostream>
#include "solvers/mpm_base.h"
#include <cmath>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "matrix/solver_base.h"

namespace mpm {

//! ThermoMPMPhaseChange class
//! \brief A class that implements the semi-implicit one phase mpm
//! \details A two-phase semi-implicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ThermoMPMPhaseChange : public MPMBase<Tdim> {
 public:
  //! Default constructor
  ThermoMPMPhaseChange(const std::shared_ptr<IO>& io);

  //! Domain decomposition
  void mpi_domain_decompose();

  //! Solve
  bool solve() override;

  // Compute time step size
  void compute_critical_timestep_size(double dt);    

  //! Class private variables
 private:
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
  //! Variable timestep
  bool variable_timestep_{false};
  //! Log output steps
  bool log_output_steps_{false};

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
  double volume_tolerance_{0};

  std::chrono::time_point<std::chrono::steady_clock> solver_begin;
};  // 
}  // namespace mpm

#include "thermo_mpm_with_phase_change.tcc"

#endif  // MPM_THERMO_MPM_PHASE_CHANGE_H_
