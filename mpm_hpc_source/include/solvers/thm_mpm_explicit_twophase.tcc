//! Constructor
template <unsigned Tdim>
mpm::ThermoMPMExplicitTwoPhase<Tdim>::ThermoMPMExplicitTwoPhase(
    const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("ThermoMPMExplicitTwoPhase");
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////                                                                ////////
////////                  HM-MPM Explicit Functions                     ////////
////////                                                                ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//! Thermo-hydro-mechanical MPM Explicit compute stress strain
template <unsigned Tdim>
void mpm::ThermoMPMExplicitTwoPhase<Tdim>::compute_stress_strain() {
  // Iterate over each particle to calculate mechancial strain
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::update_particle_strain, std::placeholders::_1, dt_));

  // Iterate over each particle to calculate thermal strain
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::update_particle_thermal_strain, std::placeholders::_1));

  // // Iterate over each particle to update particle volume
  // mesh_->iterate_over_particles(std::bind(
  //     &mpm::ParticleBase<Tdim>::update_particle_volume, std::placeholders::_1));

  // // Iterate over each particle to update material density of particle
  // mesh_->iterate_over_particles(std::bind(
  //     &mpm::ParticleBase<Tdim>::update_particle_density, std::placeholders::_1, dt_));

  // // Iterate over each particle to calculate particle porosity
  // mesh_->iterate_over_particles(std::bind(
  //     &mpm::ParticleBase<Tdim>::update_particle_porosity, std::placeholders::_1, dt_));

  // Iterate over each particle to compute stress
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::update_particle_stress, std::placeholders::_1));

  // Iterate over each particle to compute pore pressure
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_thermal_pore_pressure,
                std::placeholders::_1, dt_));  
                                            
}

//! MPM Explicit two-phase pressure smoothing
template <unsigned Tdim>
void mpm::ThermoMPMExplicitTwoPhase<Tdim>::pressure_smoothing(unsigned phase) {

  if (phase == mpm::ParticlePhase::Solid) {
    // Assign pressure to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_pressure_to_nodes,
                  std::placeholders::_1, current_time_));
  } else if (phase == mpm::ParticlePhase::Liquid) {
    // Assign pore pressure to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_pore_pressure_to_nodes,
                  std::placeholders::_1, current_time_));
  }

#ifdef USE_MPI
  int mpi_size = 1;

  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Run if there is more than a single MPI task
  if (mpi_size > 1) {
    // MPI all reduce nodal pressure
    mesh_->template nodal_halo_exchange<double, 1>(
        std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1, phase),
        std::bind(&mpm::NodeBase<Tdim>::assign_pressure, std::placeholders::_1,
                  phase, std::placeholders::_2));
  }
#endif

  if (phase == mpm::ParticlePhase::Solid) {
    // Smooth pressure over particles
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_smoothing,
                  std::placeholders::_1));
  } else if (phase == mpm::ParticlePhase::Liquid) {
    // Smooth pore pressure over particles
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_pore_pressure_smoothing,
                  std::placeholders::_1));
  }
}






////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////                                                                ////////
////////                  THM-MPM Explicit Solver                       ////////
////////                                                                ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//! Thermo-hydro-mechncial MPM Explicit solver
template <unsigned Tdim>
bool mpm::ThermoMPMExplicitTwoPhase<Tdim>::solve() {
  bool status = true;

  console_->info("MPM analysis type {}", io_->analysis_type());

  // Initialise MPI rank and size
  int mpi_rank = 0;
  int mpi_size = 1;

#ifdef USE_MPI
  // Get MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  // Two phases (soil skeleton and pore liquid)
  const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
  const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
  const unsigned mixture = mpm::ParticlePhase::Mixture;

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Pressure smoothing
  if (analysis_.find("pressure_smoothing") != analysis_.end())
    pressure_smoothing_ = analysis_["pressure_smoothing"].template get<bool>();

  // Interface
  if (analysis_.find("interface") != analysis_.end())
    interface_ = analysis_.at("interface").template get<bool>();

  // Free surface
  if (analysis_.find("free_surface") != analysis_.end()) {
    free_surface_particle_ = analysis_["free_surface"]["free_surface_particle"]
                                 .template get<std::string>();
  }

  // Variable timestep
  if (analysis_.find("variable_timestep") != analysis_.end())
    variable_timestep_ = analysis_["variable_timestep"].template get<bool>();

  // Log output steps
  if (post_process_.find("log_output_steps") != post_process_.end())
    log_output_steps_ =post_process_["log_output_steps"].template get<bool>();

  // Initialise materials
  bool mat_status = this->initialise_materials();
  if (!mat_status) {
    status = false;
    throw std::runtime_error("Initialisation of materials failed");
  }

  // Initialise mesh
  bool mesh_status = this->initialise_mesh();
  if (!mesh_status) {
    status = false;
    throw std::runtime_error("Initialisation of mesh failed");
  }

  // Initialise particles
  bool particle_status = this->initialise_particles();
  if (!particle_status) {
    status = false;
    throw std::runtime_error("Initialisation of particles failed");
  }

  // Initialise loading conditions
  bool loading_status = this->initialise_loads();
  if (!loading_status) {
    status = false;
    throw std::runtime_error("Initialisation of loads failed");
  }

  // Initialise vtk output for two phase
  this->initialise_vtk_twophase();

  // Compute mass for each phase
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1));

  // Assign initial particles at free surface
  if (free_surface_particle_ == "assign") {
    bool assign_status = mesh_->assign_free_surface_particles(io_);
    if (!assign_status) {
      status = false;
      throw std::runtime_error("Initialisation free surface particles failed");
    }
  }
  // Check point resume
  if (resume) {
    this->checkpoint_resume();
    this->current_time_ = analysis_["resume"]["current_time"].template get<double>();
  }

  solver_begin = std::chrono::steady_clock::now();

  // Creat empty reaction force file
  this->write_reaction_force(true, this->step_, this->nsteps_); 

  // Main loop
  for (step_ = 0; step_ <= nsteps_; ++step_) {

    if (variable_timestep_) {
      if (dt_matrix_size == 2) {
        if (step_ <= nsteps_matrix_[0]) {
          dt_ += (dt_matrix_[1] - dt_matrix_[0]) / nsteps_matrix_[0];
        } else {
          dt_ = dt_matrix_[1];
        };
      }
      else if (dt_matrix_size == 3){
        if (step_ <= nsteps_matrix_[0]) {
          dt_ += (dt_matrix_[1] - dt_matrix_[0]) / nsteps_matrix_[0];
        } else if (step_ > nsteps_matrix_[0] && step_ <= nsteps_matrix_[1]){
          dt_ += (dt_matrix_[2] - dt_matrix_[1]) / (nsteps_matrix_[1] - nsteps_matrix_[0]);
        } else {
          dt_ = dt_matrix_[2];
        };
      }
    }

    current_time_ += dt_;

    // Record current time
    mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::record_time, std::placeholders::_1, current_time_));

    if (mpi_rank == 0) console_->info("uuid : [{}], Step: {} of {}, timestep = {}, time = {}.\n", 
                                       uuid_, step_, nsteps_, dt_, current_time_);

    // Create a TBB task group
    tbb::task_group task_group;

    // Spawn a task for initialising nodes and cells
    task_group.run([&] {
      // Initialise nodes
      mesh_->iterate_over_nodes(
          std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

      mesh_->iterate_over_cells(
          std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));
    });

    // Spawn a task for particles
    task_group.run([&] {
      // Iterate over each particle to compute shapefn
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));
    });

    task_group.wait();

    // Assign material ids to node
    if (interface_)
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::append_material_id_to_nodes,
                    std::placeholders::_1));
   
    // Assign mass and momentum to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1));
    
    // Assign heat capacity and heat to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_heat_to_nodes,
                  std::placeholders::_1));

    // Compute nodal velocity at the begining of time step
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity, 
                  std::placeholders::_1, dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

        // Apply particle and nodal velocity constraints
    mesh_->apply_velocity_constraints(current_time_); 

    // Compute nodal temperature
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_temperature,
                  std::placeholders::_1, soil_skeleton),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));  

    // Apply nodal temperature constraints
    mesh_->apply_nodal_temperature_constraints(soil_skeleton, current_time_);

    // Update stress first
    if (this->stress_update_ == mpm::StressUpdate::USF){
      this->compute_stress_strain();
    }


    // Spawn a task for external force
    task_group.run([&] {
      // Iterate over particles to compute nodal body force
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_external_force,
                    std::placeholders::_1, this->gravity_));

      // Apply particle traction and map to nodes
      mesh_->apply_traction_on_particles(current_time_);
    });

    // Spawn a task for internal force
    task_group.run([&] {
      // Iterate over particles to compute nodal mixture and fluid internal
      // force
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_internal_force, std::placeholders::_1));

      // Iterate over particles to compute nodal drag force coefficient
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_drag_force_coefficient,
                    std::placeholders::_1));
    });
    task_group.wait();

    // Spawn a task for heat conduction
    task_group.run([&] {
      // Iterate over each particle to compute nodal heat conduction
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_heat_conduction, std::placeholders::_1));
    });
    task_group.wait();

    // // Spawn a task for heat convection
    // task_group.run([&] {
    //   // Iterate over each particle to compute nodal heat convection
    //   mesh_->iterate_over_particles(std::bind(
    //       &mpm::ParticleBase<Tdim>::map_heat_convection, std::placeholders::_1));
    // });
    // task_group.wait();

    // // Spawn a task for plastic work
    // task_group.run([&] {
    //   // Iterate over each particle to compute nodal heat convection
    //   mesh_->iterate_over_particles(std::bind(
    //       &mpm::ParticleBase<Tdim>::map_plastic_work, std::placeholders::_1, dt_));
    // });
    // task_group.wait();

    // Apply particle heat source and map to nodes
    mesh_->apply_heat_source_on_particles(current_time_, dt_);
    
    // Apply heat source on nodes
    mesh_->apply_heat_source_on_nodes(soil_skeleton, current_time_); 

// #ifdef USE_MPI
//     // Run if there is more than a single MPI task
//     if (mpi_size > 1) {
//       // MPI all reduce external force of mixture
//       mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
//           std::bind(&mpm::NodeBase<Tdim>::external_force, std::placeholders::_1,
//                     mixture),
//           std::bind(&mpm::NodeBase<Tdim>::update_external_force,
//                     std::placeholders::_1, false, mixture,
//                     std::placeholders::_2));
//       // MPI all reduce external force of pore fluid
//       mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
//           std::bind(&mpm::NodeBase<Tdim>::external_force, std::placeholders::_1,
//                     pore_liquid),
//           std::bind(&mpm::NodeBase<Tdim>::update_external_force,
//                     std::placeholders::_1, false, pore_liquid,
//                     std::placeholders::_2));

//       // MPI all reduce internal force of mixture
//       mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
//           std::bind(&mpm::NodeBase<Tdim>::internal_force, std::placeholders::_1,
//                     mixture),
//           std::bind(&mpm::NodeBase<Tdim>::update_internal_force,
//                     std::placeholders::_1, false, mixture,
//                     std::placeholders::_2));
//       // MPI all reduce internal force of pore liquid
//       mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
//           std::bind(&mpm::NodeBase<Tdim>::internal_force, std::placeholders::_1,
//                     pore_liquid),
//           std::bind(&mpm::NodeBase<Tdim>::update_internal_force,
//                     std::placeholders::_1, false, pore_liquid,
//                     std::placeholders::_2));

//       // MPI all reduce drag force
//       mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
//           std::bind(&mpm::NodeBase<Tdim>::drag_force_coefficient,
//                     std::placeholders::_1),
//           std::bind(&mpm::NodeBase<Tdim>::update_drag_force_coefficient,
//                     std::placeholders::_1, false, std::placeholders::_2));
      
//       // MPI all reduce heat conduction
//       mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
//           std::bind(&mpm::NodeBase<Tdim>::heat_conduction, std::placeholders::_1),
//           std::bind(&mpm::NodeBase<Tdim>::update_heat_conduction,
//                     std::placeholders::_1, false, 
//                     std::placeholders::_2));

//       // MPI all reduce heat convection
//       mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
//           std::bind(&mpm::NodeBase<Tdim>::heat_convection, std::placeholders::_1),
//           std::bind(&mpm::NodeBase<Tdim>::update_heat_convection,
//                     std::placeholders::_1, false, 
//                     std::placeholders::_2));
//     }
// #endif

    // Pressure smoothing
    if (pressure_smoothing_) {
      this->pressure_smoothing(pore_liquid);
    }

    // Compute free surface cells, nodes, and particles
    mesh_->compute_free_surface(free_surface_particle_, volume_tolerance_);

    // Compute nodal acceleration and update nodal velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<
                      Tdim>::compute_acc_vel_twophase_explicit,
                  std::placeholders::_1, soil_skeleton, pore_liquid, mixture,
                  this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Apply nodal velocity constraints
    mesh_->apply_velocity_constraints(current_time_); 

    // Compute nodal temperature acceleration and update nodal temperature
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_temperature,
                  std::placeholders::_1, soil_skeleton, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Apply nodal temperature constraints
    mesh_->apply_nodal_temperature_constraints(soil_skeleton, current_time_);

    // Iterate over each particle to compute updated position
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_updated_velocity,
        std::placeholders::_1, this->dt_, this->pic_, damping_factor_));

    // Apply particle velocity constraints
    mesh_->apply_velocity_constraints(current_time_);

    // Iterate over each particle to compute updated temperature
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_temperature,
        std::placeholders::_1, this->dt_, this->pic_t_));
        
    // Apply particle temperature constraints
    mesh_->apply_particle_temperature_constraints(current_time_);   

    // Update Stress Last
    if (this->stress_update_ == mpm::StressUpdate::USL) {
      this->compute_stress_strain();
        // Pressure smoothing
      if (pressure_smoothing_) { this->pressure_smoothing(pore_liquid);}
      // Compute free surface cells, nodes, and particles
      mesh_->compute_free_surface(free_surface_particle_, volume_tolerance_);  
    }

    // Locate particles
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

    // Fixed timestep, and data output linearly (every const steps = output_steps)
    if ((!variable_timestep_) & (!log_output_steps_)){
      if (step_ % output_steps_ == 0) {
        // HDF5 outputs
        if (write_hdf5_) this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
        // VTK outputs
        this->write_vtk(this->step_, this->nsteps_);
#endif
        No_output++;
        std::cout << "Number output = " << No_output << "\n";
      }
    }

    // Fixed timestep, and data output logarithmically
    // e.g., output_steps = 100 from 1 sec to 10 sec
    else if ((!variable_timestep_) & log_output_steps_){
      if (int(log10(current_time_) * output_steps_) == 
               int(output_steps_ * log10(dt_ * output_steps_) + No_output + 1)) {
        // HDF5 outputs
        if (write_hdf5_) this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
        // VTK outputs
        this->write_vtk(this->step_, this->nsteps_);
#endif
        No_output++;
        std::cout << "Number output = " << No_output << "\n";
      }
    }

    // Varied timestep, and data output logarithmically
    // e.g., output_steps = 100 from 1 sec to 10 sec
    else if (log_output_steps_ & variable_timestep_) {
      if (int(log10(current_time_) * output_steps_) == 
               int(output_steps_ * log10(dt_matrix_[0] * output_steps_) + No_output + 1)) {
        // HDF5 outputs
        if (write_hdf5_) this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
        // VTK outputs
        this->write_vtk(this->step_, this->nsteps_);
#endif
        No_output++;
        std::cout << "Number output = " << No_output << "\n";
      }
    }

  }
  auto solver_end = std::chrono::steady_clock::now();
  console_->info(
      "Rank {}, Explicit {} solver duration: {} ms", mpi_rank,
      (this->stress_update_ == mpm::StressUpdate::USL ? "USL" : "USF"),
      std::chrono::duration_cast<std::chrono::milliseconds>(solver_end -
                                                            solver_begin)
          .count());

  return status;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////                                                                ////////
////////                       Multiscale solver                        ////////
////////                                                                ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//! Thermo-mechancial MPM explicit solver for MPM-DEM
template <unsigned Tdim>
bool mpm::ThermoMPMExplicitTwoPhase<Tdim>::pre_process() {
  bool status = true;
  console_->info("MPM analysis type {}", io_->analysis_type());

  // Two phases (soil skeleton and pore liquid)
  const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
  const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
  const unsigned mixture = mpm::ParticlePhase::Mixture;

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Pressure smoothing
  if (analysis_.find("pressure_smoothing") != analysis_.end())
    pressure_smoothing_ =
        analysis_.at("pressure_smoothing").template get<bool>();

  // Interface
  if (analysis_.find("interface") != analysis_.end())
    interface_ = analysis_.at("interface").template get<bool>();

  // Free surface
  if (analysis_.find("free_surface") != analysis_.end()) {
    free_surface_particle_ = analysis_["free_surface"]["free_surface_particle"]
                                 .template get<std::string>();
  }

    // Variable timestep
  if (analysis_.find("variable_timestep") != analysis_.end())
    variable_timestep_ = analysis_["variable_timestep"].template get<bool>();

  // Log output steps
  if (post_process_.find("log_output_steps") != post_process_.end())
    log_output_steps_ =post_process_["log_output_steps"].template get<bool>(); 

  // Initialise material
  bool mat_status = this->initialise_materials();
  if (!mat_status) {
    status = false;
    throw std::runtime_error("Initialisation of materials failed");
  }

  // Initialise mesh
  bool mesh_status = this->initialise_mesh();
  if (!mesh_status) {
    status = false;
    throw std::runtime_error("Initialisation of mesh failed");
  }

  // Initialise particles
  bool particle_status = this->initialise_particles();
  if (!particle_status) {
    status = false;
    throw std::runtime_error("Initialisation of particles failed");
  }

  // Initialise loading conditions
  bool loading_status = this->initialise_loads();
  if (!loading_status) {
    status = false;
    throw std::runtime_error("Initialisation of loading failed");
  }

  // Initialise vtk output for two phase
  this->initialise_vtk_twophase();

  // Compute mass
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1));

  // Check point resume
  if (resume) this->checkpoint_resume();

  // Get total reaction force
  this->write_reaction_force(true, this->step_, this->nsteps_);
  solver_begin = std::chrono::steady_clock::now();

  // Creat empty reaction force file
  this->write_reaction_force(true, this->step_, this->nsteps_); 

  return status;
}

// Main loop
template <unsigned Tdim>
bool mpm::ThermoMPMExplicitTwoPhase<Tdim>::get_deformation_task() {
  bool status = true;
  if (step_ < nsteps_) {

    if (variable_timestep_) {
      if (dt_matrix_size == 2) {
        if (step_ <= nsteps_matrix_[0]) {
          dt_ += (dt_matrix_[1] - dt_matrix_[0]) / nsteps_matrix_[0];
        } else {
          dt_ = dt_matrix_[1];
        };
      }
      else if (dt_matrix_size == 3){
        if (step_ <= nsteps_matrix_[0]) {
          dt_ += (dt_matrix_[1] - dt_matrix_[0]) / nsteps_matrix_[0];
        } else if (step_ > nsteps_matrix_[0] && step_ <= nsteps_matrix_[1]){
          dt_ += (dt_matrix_[2] - dt_matrix_[1]) / (nsteps_matrix_[1] - nsteps_matrix_[0]);
        } else {
          dt_ = dt_matrix_[2];
        };
      }
    }

    current_time_ += dt_;

    // Record current time
    mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::record_time, std::placeholders::_1, current_time_));

  // Two phases and its mixture (soil skeleton and pore liquid)
  // NOTE: Mixture nodal variables are stored at the same memory index as the
  // solid phase
  const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
  const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
  const unsigned mixture = mpm::ParticlePhase::Mixture;

    console_->info("Step: {} of {}, timestep = {}, time = {}.\n", 
                                       step_, nsteps_, dt_, current_time_);  

    // Create a TBB task group
    tbb::task_group task_group;

    // Spawn a task for initialising nodes and cells
    task_group.run([&] {
      // Initialise nodes
      mesh_->iterate_over_nodes(
          std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

      mesh_->iterate_over_cells(
          std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));
    });

    // Spawn a task for particles
    task_group.run([&] {
      // Iterate over each particle to compute shapefn
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));
    });

    task_group.wait();

    // Assign material ids to node
    if (interface_)
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::append_material_id_to_nodes,
                    std::placeholders::_1));

    // Assign mass and momentum to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1));

    // Assign heat capacity and heat to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_heat_to_nodes,
                  std::placeholders::_1));

    // Compute nodal velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity, 
                  std::placeholders::_1, dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Apply nodal velocity constraints
    mesh_->apply_velocity_constraints(current_time_); 
    
    // Compute nodal temperature
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_temperature,
                  std::placeholders::_1, soil_skeleton),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));  

    // Apply nodal temperature constraints
    mesh_->apply_nodal_temperature_constraints(soil_skeleton, current_time_);  

    // Update stress first
    if (this->stress_update_ == mpm::StressUpdate::USF) {
      // Apply particle velocity constraints
      mesh_->apply_moving_rigid_boundary(current_time_, dt_);

      // Iterate over each particle to calculate strain of soil_skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_displacement_gradient,
                    std::placeholders::_1, dt_, true));

      // Iterate over each particle to calculate mechancial strain
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::update_particle_strain, std::placeholders::_1, dt_));

      // Iterate over each particle to calculate thermal strain
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::update_particle_thermal_strain, std::placeholders::_1));

      // Iterate over each particle to update particle volume
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::update_particle_volume, std::placeholders::_1));

      // Iterate over each particle to update material density of particle
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::update_particle_density, std::placeholders::_1, dt_));

      // Iterate over each particle to calculate particle porosity
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::update_particle_porosity, std::placeholders::_1, dt_));

      // Iterate over each particle to compute pore pressure
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_thermal_pore_pressure,
                    std::placeholders::_1, dt_)); 

      // Pressure smoothing
      if (pressure_smoothing_) {
      this->pressure_smoothing(pore_liquid);
      }
        
    }
  }
  return status;
}

template <unsigned Tdim>
void mpm::ThermoMPMExplicitTwoPhase<Tdim>::get_info(unsigned& dim, bool& resume,
                                      unsigned& checkpoint_step) {
  dim = Tdim;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  if (resume)
    checkpoint_step = analysis_["resume"]["step"].template get<mpm::Index>();
  else
    checkpoint_step = 0;
}

template <unsigned Tdim>
void mpm::ThermoMPMExplicitTwoPhase<Tdim>::get_status(double& dt, unsigned& step,
                                        unsigned& nsteps,
                                        unsigned& output_steps) {
  dt = dt_;
  step = step_;
  nsteps = nsteps_;
  output_steps = output_steps_;
}

template <unsigned Tdim>
bool mpm::ThermoMPMExplicitTwoPhase<Tdim>::send_deformation_task(
    std::vector<unsigned>& id,
    std::vector<Eigen::MatrixXd>& displacement_gradients) {
  bool status = true;
  mesh_->get_displacement_gradient(id, displacement_gradients);
  return status;
}

template <unsigned Tdim>
bool mpm::ThermoMPMExplicitTwoPhase<Tdim>::set_stress_task(const Eigen::MatrixXd& stresses,
                                             bool increment) {
  bool status = true;
  if (step_ < nsteps_) {
    // Iterate over each particle to compute stress of soil skeleton
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::set_stress, std::placeholders::_1,
                  stresses, increment));
  }
  return status;
}

template <unsigned Tdim>
bool mpm::ThermoMPMExplicitTwoPhase<Tdim>::set_porosity_task(
    const Eigen::MatrixXd& porosities) {
  bool status = true;
  if (step_ < nsteps_) {
    // Iterate over each particle to set porosity
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::set_porosity, std::placeholders::_1,
                  porosities));
  }
  return status;
}

template <unsigned Tdim>
bool mpm::ThermoMPMExplicitTwoPhase<Tdim>::set_fabric_task(std::string fabric_type,
                                             const Eigen::MatrixXd& fabrics) {
  bool status = true;
  if (step_ < nsteps_) {
    // Iterate over each particle to set fabrics
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::set_fabric, std::placeholders::_1,
                  fabric_type, fabrics));
  }
  return status;
}

template <unsigned Tdim>
bool mpm::ThermoMPMExplicitTwoPhase<Tdim>::set_rotation_task(
    const Eigen::MatrixXd& rotations) {
  bool status = true;
  if (step_ < nsteps_) {
    // Iterate over each particle to set rotation
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::set_rotation, std::placeholders::_1,
                  rotations));
  }
  return status;
}

template <unsigned Tdim>
bool mpm::ThermoMPMExplicitTwoPhase<Tdim>::update_state_task() {
  bool status = true;
  if (step_ < nsteps_) {

  // Two phases and its mixture (soil skeleton and pore liquid)
  // NOTE: Mixture nodal variables are stored at the same memory index as the
  // solid phase
  const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
  const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
  const unsigned mixture = mpm::ParticlePhase::Mixture;

    // Create a TBB task group
    tbb::task_group task_group;

    // Spawn a task for external force
    task_group.run([&] {
      // Iterate over each particle to compute nodal body force
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_external_force,
                    std::placeholders::_1, this->gravity_));

      // Apply particle traction and map to nodes
      mesh_->apply_traction_on_particles(current_time_);

    });

    // Spawn a task for internal force
    task_group.run([&] {
      // Iterate over particles to compute nodal mixture and fluid internal
      // force
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_internal_force, std::placeholders::_1));

      // Iterate over particles to compute nodal drag force coefficient
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_drag_force_coefficient,
                    std::placeholders::_1));
    });
    task_group.wait();

    // Spawn a task for heat conduction
    task_group.run([&] {
      // Iterate over each particle to compute nodal heat conduction
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_heat_conduction, std::placeholders::_1));
    });
    task_group.wait();

    // Spawn a task for heat convection
    task_group.run([&] {
      // Iterate over each particle to compute nodal heat convection
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_heat_convection, std::placeholders::_1));
    });
    task_group.wait();    

    // Spawn a task for plastic work
    task_group.run([&] {
      // Iterate over each particle to compute nodal heat convection
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_plastic_work, std::placeholders::_1, dt_));
    });
    task_group.wait();

    // Apply particle heat source and map to nodes
    mesh_->apply_heat_source_on_particles(current_time_, dt_);
    
    // Apply heat source on nodes
    mesh_->apply_heat_source_on_nodes(soil_skeleton, current_time_); 

    // // Pressure smoothing
    // if (pressure_smoothing_) {
    //   this->pressure_smoothing(pore_liquid);
    // }

    // Compute free surface cells, nodes, and particles
    mesh_->compute_free_surface(free_surface_particle_, volume_tolerance_);

    // Compute nodal acceleration and update nodal velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<
                      Tdim>::compute_acc_vel_twophase_explicit,
                  std::placeholders::_1, soil_skeleton, pore_liquid, mixture,
                  this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Apply nodal velocity constraints
    mesh_->apply_velocity_constraints(current_time_); 

    // Compute nodal temperature acceleration and update nodal temperature
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_temperature,
                  std::placeholders::_1, soil_skeleton, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Apply nodal temperature constraints
    mesh_->apply_nodal_temperature_constraints(soil_skeleton, current_time_);

    // Iterate over each particle to compute updated position
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_updated_velocity,
        std::placeholders::_1, this->dt_, this->pic_, damping_factor_));

    // Apply particle velocity constraints
    mesh_->apply_velocity_constraints(current_time_);

    // Iterate over each particle to compute updated temperature
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_temperature,
        std::placeholders::_1, this->dt_, this->pic_t_));
        
    // Apply particle temperature constraints
    mesh_->apply_particle_temperature_constraints(current_time_); 

    // Locate particles
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty()) {
      status = false;
      throw std::runtime_error("Particle outside the mesh domain");
    }

    // Fixed timestep, and data output linearly (every const steps = output_steps)
    if ((!variable_timestep_) & (!log_output_steps_)){
      if (step_ % output_steps_ == 0) {
        // HDF5 outputs
        if (write_hdf5_) this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
        // VTK outputs
        this->write_vtk(this->step_, this->nsteps_);
#endif
        No_output++;
        std::cout << "Number output = " << No_output << "\n";
      }
    }

    // Fixed timestep, and data output logarithmically
    // e.g., output_steps = 100 from 1 sec to 10 sec
    else if ((!variable_timestep_) & log_output_steps_){
      if (int(log10(current_time_) * output_steps_) == 
               int(output_steps_ * log10(dt_ * output_steps_) + No_output + 1)) {
        // HDF5 outputs
        if (write_hdf5_) this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
        // VTK outputs
        this->write_vtk(this->step_, this->nsteps_);
#endif
        No_output++;
        std::cout << "Number output = " << No_output << "\n";
      }
    }

    // Varied timestep, and data output logarithmically
    // e.g., output_steps = 100 from 1 sec to 10 sec
    else if (log_output_steps_ & variable_timestep_) {
      if (int(log10(current_time_) * output_steps_) == 
               int(output_steps_ * log10(dt_matrix_[0] * output_steps_) + No_output + 1)) {
        // HDF5 outputs
        if (write_hdf5_) this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
        // VTK outputs
        this->write_vtk(this->step_, this->nsteps_);
#endif
        No_output++;
        std::cout << "Number output = " << No_output << "\n";
      }
    } 

    step_++;
    if (step_ == nsteps_) {
      auto solver_end = std::chrono::steady_clock::now();
      console_->info(
          "Rank {}, Explicit {} solver duration: {} ms", 0,
          (this->stress_update_ == mpm::StressUpdate::USL ? "USL" : "USF"),
          std::chrono::duration_cast<std::chrono::milliseconds>(solver_end -
                                                                solver_begin)
              .count());
    }
  }
  return status;
}