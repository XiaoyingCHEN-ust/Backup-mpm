//! Constructor
template <unsigned Tdim>
mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::ThermoMPMSemiImplicitTwoPhase(
    const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("ThermoMPMSemiImplicitTwoPhase");
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////                                                                ////////
////////                THM-MPM Semi-implicit Solver                    ////////
////////                                                                ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//! MPM semi-implicit two phase solver
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::solve() {
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

  // Two phases and its mixture (soil skeleton and pore liquid)
  // NOTE: Mixture nodal variables are stored at the same memory index as the
  // solid phase
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

  // Projection method paramter (beta)
  if (analysis_.find("semi_implicit") != analysis_.end()) {
    beta_ = analysis_["semi_implicit"]["beta"].template get<double>();
    free_surface_particle_ = analysis_["semi_implicit"]["free_surface_particle"]
                                 .template get<std::string>();
    implicit_drag_force_ = analysis_["semi_implicit"]["implicit_drag_force"]
                                .template get<bool>();                                      
    alpha_ = analysis_["semi_implicit"]["alpha"].template get<double>();                               
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
    throw std::runtime_error("Initialisation of loads failed");
  }

  // Initialise matrix
  bool matrix_status = this->initialise_matrix();
  if (!matrix_status) {
    status = false;
    throw std::runtime_error("Initialisation of matrix failed");
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

  this->compute_critical_timestep_size(dt_);

  auto solver_begin = std::chrono::steady_clock::now();

  // Main loop
  for (step_ = 0; step_ <= nsteps_; ++step_) {
auto onestep_begin = std::chrono::steady_clock::now();
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
// auto branch_begin = std::chrono::steady_clock::now();

        // Initialise nodes
        mesh_->iterate_over_nodes(
            std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

        mesh_->iterate_over_cells(
            std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));

    // if (step_ == 0) {
        // Iterate over each particle to compute shapefn
        mesh_->iterate_over_particles(std::bind(
            &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));
    // }

/////////////////////////////////////////////////////////////////////////////////
////////            Compute nodal velocities and temperatures            ////////
/////////////////////////////////////////////////////////////////////////////////

    // Apply particle and nodal velocity constraints
    mesh_->apply_velocity_constraints(current_time_);

    // Assign mass and momentum to nodes
    // Solid skeleton
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1));

    // Compute nodal velocity at the begining of time step
    // Nodal velocity constraint are also applied
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity, 
                  std::placeholders::_1, dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // // Apply particle velocity constraints
    // mesh_->apply_moving_rigid_boundary(current_time_, dt_);             

    // ! TODO:REMOVE APPLY WATER TABLE
    // Apply nodal water table
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::apply_water_table,
                  std::placeholders::_1, current_time_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

////////////////////////////////////////////////////////////////////////////////
////////                      Update stress first                       ////////
////////////////////////////////////////////////////////////////////////////////

    // Update stress first
    if (this->stress_update_ == mpm::StressUpdate::USF) {
     
      // Iterate over each particle to calculate strain of soil_skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_particle_strain,
                    std::placeholders::_1, dt_));

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

      // Iterate over each particle to compute stress of soil skeleton
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::update_particle_stress, std::placeholders::_1));

      // // Iterate over each particle to update permeability
      // mesh_->iterate_over_particles(
      //   std::bind(&mpm::ParticleBase<Tdim>::update_permeability,
      //             std::placeholders::_1));

    }

////////////////////////////////////////////////////////////////////////////////
////////               Update nodal and particle temperatures           ////////
////////////////////////////////////////////////////////////////////////////////

    // Assign heat capacity and heat to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_heat_to_nodes, std::placeholders::_1));   
   
    // Compute nodal temperature
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_temperature,
                  std::placeholders::_1, soil_skeleton),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1)); 

    // Apply nodal temperature constraints
    mesh_->apply_nodal_temperature_constraints(soil_skeleton, current_time_);
       
    // Iterate over each particle to compute nodal heat conduction
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_heat_conduction, std::placeholders::_1));

    // // Iterate over each particle to compute nodal heat convection
    // mesh_->iterate_over_particles(std::bind(
    //     &mpm::ParticleBase<Tdim>::map_heat_convection, std::placeholders::_1));

    // // Iterate over each particle to compute nodal plastic work
    // mesh_->iterate_over_particles(std::bind(
    //     &mpm::ParticleBase<Tdim>::map_plastic_work, std::placeholders::_1, dt_));

    // Apply particle heat source and map to nodes
    mesh_->apply_heat_source_on_particles(current_time_, dt_);
    
    // // Apply heat source on nodes
    // mesh_->apply_heat_source_on_nodes(soil_skeleton, current_time_);

    // Compute nodal temperature acceleration and update nodal temperature
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_temperature,
                  std::placeholders::_1, soil_skeleton, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Apply nodal temperature constraints
    mesh_->apply_nodal_temperature_constraints(soil_skeleton, current_time_);

/////////////////////////////////////////////////////////////////////////////////
////////              Update nodal accelerations and pressure            ////////
/////////////////////////////////////////////////////////////////////////////////

    // Iterate over particles to compute nodal body force of soil skeleton
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_external_force,
                  std::placeholders::_1, this->gravity_));

    // Apply particle traction and map to nodes
    mesh_->apply_traction_on_particles(current_time_);
  
    // Iterate over particles to compute nodal mixture internal force
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_internal_force_semi,
                  std::placeholders::_1, beta_));

    // Iterate over particles to compute nodal drag force coefficient
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_drag_force_coefficient,
                  std::placeholders::_1));  

    // Reinitialise system matrix
    bool matrix_reinitialization_status = this->reinitialise_matrix();
    if (!matrix_reinitialization_status) {
      status = false;
      throw std::runtime_error("Reinitialisation of matrix failed");
    }

    // Compute free surface cells, nodes, and particles
    mesh_->compute_free_surface(free_surface_particle_, volume_tolerance_);

    // // map_mass_pressure_to_nodes
    // mesh_->iterate_over_particles(
    //     std::bind(&mpm::ParticleBase<Tdim>::map_pore_pressure_to_nodes,
    //               std::placeholders::_1, current_time_));


    // Iterate over active nodes to compute acceleratation and velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<
                      Tdim>::compute_inter_acc_vel_twophase_semi,
                  std::placeholders::_1, soil_skeleton, pore_liquid, mixture,
                  this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));


    // // Compute intermediate velocity
    // this->compute_intermediate_velocity();

    // // Update intermediate velocity of solid phase
    // mesh_->iterate_over_nodes_predicate(
    //     std::bind(&mpm::NodeBase<Tdim>::update_intermediate_velocity,
    //               std::placeholders::_1, soil_skeleton,
    //               matrix_assembler_->acceleration_inter().topRows(
    //                   matrix_assembler_->active_dof()),
    //               dt_),
    //     std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // // Update intermediate velocity of water phase
    // mesh_->iterate_over_nodes_predicate(
    //     std::bind(&mpm::NodeBase<Tdim>::update_intermediate_velocity,
    //               std::placeholders::_1, pore_liquid,
    //               matrix_assembler_->acceleration_inter().bottomRows(
    //                   matrix_assembler_->active_dof()),
    //               dt_),
    //     std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));


    // // Apply particle velocity constraints
    mesh_->apply_moving_rigid_boundary(current_time_, dt_);

    // // Assign contact velocity
    // mesh_->iterate_over_nodes_predicate(
    //     std::bind(&mpm::NodeBase<Tdim>::compute_rigid_velocity,
    //               std::placeholders::_1, dt_),
    //     std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Assign contact velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::assign_intermediate_velocity_from_rigid,
                  std::placeholders::_1, dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));   

    // Compute poisson equation
    this->compute_poisson_equation();

    // Assign pore pressure to nodes
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::update_pore_pressure_increment,
                  std::placeholders::_1,
                  matrix_assembler_->pore_pressure_increment(),
                  current_time_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Use nodal pore pressure to update particle pore pressure
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_updated_pore_pressure,
                  std::placeholders::_1, beta_));

    // Apply particle pore pressusre constraints
    mesh_->apply_particle_pore_pressure_constraints(current_time_); 

    // Pressure smoothing
    if (pressure_smoothing_) {
      this->pressure_smoothing(pore_liquid);     
    } 

    // Compute corrected force
    this->compute_corrected_force();

/////////////////////////////////////////////////////////////////////////////////
////////              Update particle velocities, pressures              ////////
/////////////////////////////////////////////////////////////////////////////////
    // Iterate over active nodes to compute acceleratation and velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(
          &mpm::NodeBase<Tdim>::compute_acc_vel_twophase_semi,
            std::placeholders::_1, soil_skeleton, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    mesh_->iterate_over_nodes_predicate(
        std::bind(
            &mpm::NodeBase<Tdim>::compute_acc_vel_twophase_semi,
            std::placeholders::_1, pore_liquid, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Assign contact velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::assign_corrected_velocity_from_rigid,
                  std::placeholders::_1, dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1)); 

    // Apply particle and nodal velocity constraints
    mesh_->apply_velocity_constraints(current_time_);

    // Update particle position and kinematics
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_updated_velocity,
        std::placeholders::_1, this->dt_, this->pic_, damping_factor_));
    
    // Iterate over each particle to compute updated temperature
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_temperature,
        std::placeholders::_1, this->dt_, this->pic_t_));
        
    // Apply particle temperature constraints
    mesh_->apply_particle_temperature_constraints(current_time_); 

    // Assign nodal pressure increment constraints
    mesh_->iterate_over_nodes_predicate(
        std::bind(
            &mpm::NodeBase<Tdim>::assign_nodal_pressure_increment_constraints,
            std::placeholders::_1, this->step_, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

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
      "Rank {}, SemiImplicit_Twophase {} solver duration: {} ms", mpi_rank,
      (this->stress_update_ == mpm::StressUpdate::USL ? "USL" : "USF"),
      std::chrono::duration_cast<std::chrono::milliseconds>(solver_end -
                                                            solver_begin)
          .count());
  return status;
}

// Compute time step size
template <unsigned Tdim>
void mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::compute_critical_timestep_size(double dt) {
  const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
  const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
  const unsigned mixture = mpm::ParticlePhase::Mixture;

  // cell minimum size
  auto mesh_props = io_->json_object("mesh");
  // Get Mesh reader from JSON object
  double cellsize_min = mesh_props.at("cellsize_min").template get<double>();
  // Solid Material parameters 
  auto materials =  materials_.at(soil_skeleton);
  double porosity = materials->template property<double>(std::string("porosity"));
  double youngs_modulus = materials->template property<double>(std::string("youngs_modulus"));
  double poisson_ratio = materials->template property<double>(std::string("poisson_ratio"));
  double density = materials->template property<double>(std::string("density"));
  double specific_heat = materials->template property<double>(std::string("specific_heat"));
  double thermal_conductivity = materials->template property<double>(std::string("thermal_conductivity"));
  // Compute timestep fpor one phase MPM                              
  double critical_dt = cellsize_min / std::pow(youngs_modulus/density/(1 - porosity), 0.5);
  console_->info("Critical time step size is {} s", critical_dt);
  // Liquid Material parameters 
  auto liquid_materials =  materials_.at(pore_liquid);
  double liquid_density = liquid_materials->template property<double>(std::string("density"));
  double liquid_specific_heat = liquid_materials->template property<double>(std::string("liquid_specific_heat"));
  double liquid_thermal_conductivity = liquid_materials->template property<double>(std::string("liquid_thermal_conductivity"));

  // Compute timestep for momentum eqaution 
  double density_mixture1 = (1 - porosity) * density;
  double density_mixture2 = (1 - porosity) * density + porosity * liquid_density;
  double critical_dt11 = cellsize_min / std::pow(youngs_modulus/density_mixture1, 0.5);
  double critical_dt12 = cellsize_min / std::pow(youngs_modulus/density_mixture2, 0.5);
  console_->info("Critical time step size for elastic wave propagation (solid base) is {} s", critical_dt11);
  console_->info("Critical time step size for elastic wave propagation (liquid base) is {} s", critical_dt12);

  // Compute timestep for heat transfer eqaution - pure liquid
  double k_mixture1 = (1 - porosity) * thermal_conductivity + porosity * liquid_thermal_conductivity;
  double c_mixture1 = (1 - porosity) * density * specific_heat + porosity * liquid_density * liquid_specific_heat;
  double critical_dt21 = cellsize_min * cellsize_min * c_mixture1 / k_mixture1;                                  
  console_->info("Critical time step size for thermal conduction equation (liquid base) is {} s", critical_dt21);

  if (dt >= std::min(critical_dt11, critical_dt12) || dt >= critical_dt21)
      throw std::runtime_error("Time step size is too large");                             
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////                                                                ////////
////////                       Multiscale solver                        ////////
////////                                                                ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Pre process for MPM-DEM
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::pre_process() {
  bool status = true;
  console_->info("MPM analysis type {}", io_->analysis_type());

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Pressure smoothing
  if (analysis_.find("pressure_smoothing") != analysis_.end())
    pressure_smoothing_ = analysis_["pressure_smoothing"].template get<bool>();

  // Variable timestep
  if (analysis_.find("variable_timestep") != analysis_.end())
    variable_timestep_ = analysis_["variable_timestep"].template get<bool>();

  // Log output steps
  if (post_process_.find("log_output_steps") != post_process_.end())
    log_output_steps_ =post_process_["log_output_steps"].template get<bool>();

  // Projection method paramter (beta)
  if (analysis_.find("semi_implicit") != analysis_.end()) {
    beta_ = analysis_["semi_implicit"]["beta"].template get<double>();
    free_surface_particle_ = analysis_["semi_implicit"]["free_surface_particle"]
                                 .template get<std::string>();
  }

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
    throw std::runtime_error("Initialisation of loads failed");
  }

  // Initialise matrix
  bool matrix_status = this->initialise_matrix();
  if (!matrix_status) {
    status = false;
    throw std::runtime_error("Initialisation of matrix failed");
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
  if (resume) this->checkpoint_resume();

  auto solver_begin = std::chrono::steady_clock::now();

  return status;
}

// Main loop
// (1) Compute deformation gradient for MPM-DEM
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::get_deformation_task() {
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

    // Two phases and its mixture (soil skeleton and pore liquid)
    // NOTE: Mixture nodal variables are stored at the same memory index as the
    // solid phase
    const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
    const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
    const unsigned mixture = mpm::ParticlePhase::Mixture;

    current_time_ += dt_;

    // Record current time
    mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::record_time, std::placeholders::_1, current_time_));

    console_->info("uuid : [{}], Step: {} of {}, timestep = {}, time = {}.\n", 
                                       uuid_, step_, nsteps_, dt_, current_time_);

    // Initialise nodes
    mesh_->iterate_over_nodes(
        std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));
    mesh_->iterate_over_cells(
        std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));
    // Iterate over each particle to compute shapefn
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));

////////////////////////////////////////////////////////////////////////////////
////////          Compute nodal velocities and temperatures             ////////
////////////////////////////////////////////////////////////////////////////////

    // Apply particle and nodal velocity constraints
    mesh_->apply_velocity_constraints(current_time_);

    // Assign mass and momentum to nodes
    // Solid skeleton
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1));

    // Compute nodal velocity at the begining of time step
    // Nodal velocity constraint are also applied
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity, 
                  std::placeholders::_1, dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Apply particle velocity constraints
    mesh_->apply_moving_rigid_boundary(current_time_, dt_);             

    //! TODO:REMOVE APPLY WATER TABLE
    // // Apply nodal water table
    // mesh_->iterate_over_nodes_predicate(
    //     std::bind(&mpm::NodeBase<Tdim>::apply_water_table,
    //               std::placeholders::_1, current_time_),
    //     std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1)); 

////////////////////////////////////////////////////////////////////////////////
////////                      Update stress first                       ////////
////////////////////////////////////////////////////////////////////////////////

    // Update stress first
    if (this->stress_update_ == mpm::StressUpdate::USF) {
      // Iterate over each particle to calculate strain of soil_skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_displacement_gradient,
                    std::placeholders::_1, dt_, true));

      // Iterate over each particle to calculate strain of soil_skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_particle_strain,
                    std::placeholders::_1, dt_));

        // Iterate over each particle to calculate thermal strain
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::update_particle_thermal_strain, std::placeholders::_1));

      // Iterate over each particle to update particle volume
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::update_particle_volume, std::placeholders::_1));

      // // Iterate over each particle to update material density of particle
      // mesh_->iterate_over_particles(std::bind(
      //     &mpm::ParticleBase<Tdim>::update_particle_density, std::placeholders::_1, dt_));

      // // Iterate over each particle to calculate particle porosity
      // mesh_->iterate_over_particles(std::bind(
      //     &mpm::ParticleBase<Tdim>::update_particle_porosity, std::placeholders::_1, dt_));

      // // Iterate over each particle to update permeability
      // mesh_->iterate_over_particles(
      //     std::bind(&mpm::ParticleBase<Tdim>::update_permeability,
      //               std::placeholders::_1));
    }
      // Assign heat capacity and heat to nodes
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_heat_to_nodes, std::placeholders::_1));   
    
      // Compute nodal temperature
      mesh_->iterate_over_nodes_predicate(
          std::bind(&mpm::NodeBase<Tdim>::compute_temperature,
                    std::placeholders::_1, soil_skeleton),
          std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1)); 

      // Apply nodal temperature constraints
      mesh_->apply_nodal_temperature_constraints(soil_skeleton, current_time_);

////////////////////////////////////////////////////////////////////////////////
////////               Update nodal and particle temperatures           ////////
////////////////////////////////////////////////////////////////////////////////
 
      // Iterate over each particle to compute nodal heat conduction
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_heat_conduction, std::placeholders::_1));

      // // // Iterate over each particle to compute nodal heat convection
      // mesh_->iterate_over_particles(std::bind(
      //     &mpm::ParticleBase<Tdim>::map_heat_convection, std::placeholders::_1));

      // // Iterate over each particle to compute nodal plastic work
      // mesh_->iterate_over_particles(std::bind(
      //     &mpm::ParticleBase<Tdim>::map_plastic_work, std::placeholders::_1, dt_));

      // // Apply particle heat source and map to nodes
      // mesh_->apply_heat_source_on_particles(current_time_, dt_);
      
      // // Apply heat source on nodes
      // mesh_->apply_heat_source_on_nodes(soil_skeleton, current_time_);

      // Compute nodal temperature acceleration and update nodal temperature
      mesh_->iterate_over_nodes_predicate(
          std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_temperature,
                    std::placeholders::_1, soil_skeleton, this->dt_),
          std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

      // Apply nodal temperature constraints
      mesh_->apply_nodal_temperature_constraints(soil_skeleton, current_time_);    
  }
  return status;
}

// (2) Get analysis information
template <unsigned Tdim>
void mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::get_info(unsigned& dim, bool& resume,
                                                  unsigned& checkpoint_step) {
  dim = Tdim;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  if (resume)
    checkpoint_step = analysis_["resume"]["step"].template get<mpm::Index>();
  else
    checkpoint_step = 0;
}

// (3) Get steps, timesteps, and output steps
template <unsigned Tdim>
void mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::get_status(double& dt, unsigned& step,
                                                    unsigned& nsteps,
                                                    unsigned& output_steps) {
  dt = dt_;
  step = step_;
  nsteps = nsteps_;
  output_steps = output_steps_;
}

// (4) Get deformation for DEM
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::send_deformation_task(
    std::vector<unsigned>& id,
    std::vector<Eigen::MatrixXd>& displacement_gradients) {
  bool status = true;
  mesh_->get_displacement_gradient(id, displacement_gradients);
  return status;
}

template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::send_temperature_task(
    std::vector<unsigned>& id,
    std::vector<double>& particle_temperature) {
  bool status = true;
  mesh_->get_particle_temperature(id, particle_temperature);
  return status;
}

// Set particle stress
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::set_stress_task(
    const Eigen::MatrixXd& stresses, bool increment) {
  bool status = true;
  if (step_ < nsteps_) {
    // Iterate over each particle to compute stress of soil skeleton
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::set_stress, std::placeholders::_1,
                  stresses, increment));
  }
  return status;
}

// (5) Set particle porosity
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::set_porosity_task(
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

// (6) Set particle fabric
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::set_fabric_task(
    std::string fabric_type, const Eigen::MatrixXd& fabrics) {
  bool status = true;
  if (step_ < nsteps_) {
    // Iterate over each particle to set fabrics
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::set_fabric, std::placeholders::_1,
                  fabric_type, fabrics));
  }
  return status;
}

// (7) Set particle rotation
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::set_rotation_task(
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

// (8) Update particle state, e.g., position, velocity
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::update_state_task() {
  // Two phases and its mixture (soil skeleton and pore liquid)
  // NOTE: Mixture nodal variables are stored at the same memory index as the
  // solid phase

  bool status = true;
  if (step_ < nsteps_) {
    const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
    const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
    const unsigned mixture = mpm::ParticlePhase::Mixture;

/////////////////////////////////////////////////////////////////////////////////
////////              Update nodal accelerations and pressure            ////////
/////////////////////////////////////////////////////////////////////////////////

      // Iterate over particles to compute nodal body force of soil skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_external_force,
                    std::placeholders::_1, this->gravity_));

      // Apply particle traction and map to nodes
      mesh_->apply_traction_on_particles(current_time_);

      // Iterate over particles to compute nodal mixture internal force
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_internal_force_semi,
                    std::placeholders::_1, beta_));

    // Reinitialise system matrix
    bool matrix_reinitialization_status = this->reinitialise_matrix();
    if (!matrix_reinitialization_status) {
      status = false;
      throw std::runtime_error("Reinitialisation of matrix failed");
    }

    // Compute free surface cells, nodes, and particles
    mesh_->compute_free_surface(free_surface_particle_, volume_tolerance_);

    // Compute intermediate velocity
    this->compute_intermediate_velocity();

    // Update intermediate velocity of solid phase
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::update_intermediate_velocity,
                  std::placeholders::_1, soil_skeleton,
                  matrix_assembler_->acceleration_inter().topRows(
                      matrix_assembler_->active_dof()),
                  dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Update intermediate velocity of water phase
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::update_intermediate_velocity,
                  std::placeholders::_1, pore_liquid,
                  matrix_assembler_->acceleration_inter().bottomRows(
                      matrix_assembler_->active_dof()),
                  dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Compute poisson equation
    this->compute_poisson_equation();

    // Assign pore pressure to nodes
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::update_pore_pressure_increment,
                  std::placeholders::_1,
                  matrix_assembler_->pore_pressure_increment(),
                  current_time_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

/////////////////////////////////////////////////////////////////////////////////
////////              Update particle velocities, pressures              ////////
/////////////////////////////////////////////////////////////////////////////////

    // Use nodal pore pressure to update particle pore pressure
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_updated_pore_pressure,
                  std::placeholders::_1, beta_));

    // Pressure smoothing
    if (pressure_smoothing_) {
      this->pressure_smoothing(pore_liquid);
    } 

    // Compute corrected force
    this->compute_corrected_force();

    // Iterate over active nodes to compute acceleratation and velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(
            &mpm::NodeBase<
                Tdim>::compute_acc_vel_twophase_semi,
            std::placeholders::_1, soil_skeleton, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    mesh_->iterate_over_nodes_predicate(
        std::bind(
            &mpm::NodeBase<
                Tdim>::compute_acc_vel_twophase_semi,
            std::placeholders::_1, pore_liquid, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Apply particle and nodal velocity constraints
    mesh_->apply_velocity_constraints(current_time_);

    // Update particle position and kinematics
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_updated_velocity,
        std::placeholders::_1, this->dt_, this->pic_, damping_factor_));
    
    // Iterate over each particle to compute updated temperature
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_temperature,
        std::placeholders::_1, this->dt_, this->pic_t_));
        
    // Apply particle temperature constraints
    mesh_->apply_particle_temperature_constraints(current_time_); 

    // Assign nodal pressure increment constraints
    mesh_->iterate_over_nodes_predicate(
        std::bind(
            &mpm::NodeBase<Tdim>::assign_nodal_pressure_increment_constraints,
            std::placeholders::_1, this->step_, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

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
          "Rank {}, ThermoMPMSemiImplicitTwoPhase {} solver duration: {} ms", 0,
          (this->stress_update_ == mpm::StressUpdate::USL ? "USL" : "USF"),
          std::chrono::duration_cast<std::chrono::milliseconds>(solver_end -
                                                                solver_begin)
              .count());
    }
  }
  return status;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////                                                                ////////
////////                 Semi-implicit functions                        ////////
////////                                                                ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Initialise matrix
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::initialise_matrix() {
  bool status = true;
  try {
    // Max iteration steps
    unsigned max_iter =
        analysis_["matrix"]["max_iter"].template get<unsigned>();
    // Tolerance
    double tolerance = analysis_["matrix"]["tolerance"].template get<double>();
    // Get matrix assembler type
    std::string assembler_type =
        analysis_["matrix"]["assembler_type"].template get<std::string>();
    // Get matrix solver type
    std::string solver_type =
        analysis_["matrix"]["solver_type"].template get<std::string>();
    // Get volume tolerance for free surface
    volume_tolerance_ =
        analysis_["matrix"]["volume_tolerance"].template get<double>();
    // Get thread used for cg parallel computing
    num_threads = 
        analysis_["matrix"]["num_threads"].template get<int>();
    // Get entries number to speed assembly matrix
    K_entries_number_ =
        analysis_["matrix"]["entries_number"]["K_matrix"].template get<int>();
    L_entries_number_ =
        analysis_["matrix"]["entries_number"]["L_matrix"].template get<int>();
    F_entries_number_ =
        analysis_["matrix"]["entries_number"]["F_matrix"].template get<int>();
    T_entries_number_ =
        analysis_["matrix"]["entries_number"]["T_matrix"].template get<int>();
    K_cor_entries_number_ =
        analysis_["matrix"]["entries_number"]["K_cor_matrix"].template get<int>();      
    // Create matrix assembler
    matrix_assembler_ =
        Factory<mpm::AssemblerBase<Tdim>>::instance()->create(assembler_type);
    // Create matrix solver
    matrix_solver_ =
        Factory<mpm::SolverBase<Tdim>, unsigned, double>::instance()->create(
            solver_type, std::move(max_iter), std::move(tolerance));
    // Transfer the mesh pointer
    matrix_assembler_->assign_mesh_pointer(mesh_);

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Reinitialise and resize matrices at the beginning of every time step
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::reinitialise_matrix() {
  bool status = true;
  try {
    const auto active_dof = mesh_->assign_active_node_id();

    // Assign global node indice
    matrix_assembler_->assign_global_node_indices(active_dof);

    // Assign pressure cpnstraints
    matrix_assembler_->assign_pressure_constraints(this->beta_,
                                                   current_time_);

    // Assign velocity constraints
    matrix_assembler_->assign_velocity_constraints();

    // Initialise element matrix
    mesh_->iterate_over_cells(std::bind(
        &mpm::Cell<Tdim>::initialise_element_matrix, std::placeholders::_1));

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// TODO: This is a copy of pressure_smoothing in explicit two-phase
//! MPM Explicit two-phase pressure smoothing
template <unsigned Tdim>
void mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::pressure_smoothing(unsigned phase) {

  // Map pressures to nodes
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

  // Map Pressure back to particles
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

//! Compute intermediate velocity
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::compute_intermediate_velocity(
    std::string solver_type) {
  bool status = true;
  try {
    // Map K_inter to cell
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_K_inter_to_cell, std::placeholders::_1));

    // Assemble stiffness matrix A
    for (unsigned i = 0; i < Tdim; ++i) {
      matrix_assembler_->assemble_stiffness_matrix(i, dt_, implicit_drag_force_, K_entries_number_);
    } 
    // Assemble force vector b
    matrix_assembler_->assemble_force_vector(dt_);
    // Apply velocity constraints to A and b
    matrix_assembler_->apply_velocity_constraints();
  
    // Compute matrix equation of each direction
    for (unsigned i = 0; i < Tdim; ++i) { 
      // Solve equation 1 to compute intermediate acceleration
      matrix_assembler_->assign_intermediate_acceleration(
          i, matrix_solver_->solve(matrix_assembler_->stiffness_matrix(i),
                                   matrix_assembler_->force_inter().col(i),
                                   solver_type, num_threads));
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute poisson equation
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::compute_poisson_equation(
    std::string solver_type) {
  bool status = true;
  try {
    // Map Laplacian elements
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_L_to_cell, std::placeholders::_1, dt_, alpha_));         
    // Assemble laplacian matrix
    matrix_assembler_->assemble_laplacian_matrix(dt_, L_entries_number_);

    // Map Fs & Fm matrix
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_F_to_cell, std::placeholders::_1));
  
    // Map Ts & Tw matrix
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_T_to_cell, std::placeholders::_1));       
  
    // Assemble force vector
    matrix_assembler_->assemble_poisson_right(mesh_, dt_, F_entries_number_);
    matrix_assembler_->assemble_poisson_right_thermal(mesh_, dt_, T_entries_number_);
  
    // Assign free surface
    matrix_assembler_->assign_free_surface(mesh_->free_surface_nodes());
  
    // Apply constraints
    matrix_assembler_->apply_pressure_constraints();

    // Solve matrix equation (compute pore pressure)
    matrix_assembler_->assign_pore_pressure_increment(matrix_solver_->solve(
        matrix_assembler_->laplacian_matrix(),
        matrix_assembler_->force_laplacian_matrix(), solver_type, num_threads));

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Compute corrected force
template <unsigned Tdim>
bool mpm::ThermoMPMSemiImplicitTwoPhase<Tdim>::compute_corrected_force() {
  bool status = true;
  try {
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_K_cor_to_cell, std::placeholders::_1, 
                                                              dt_, alpha_));
    // Assemble corrected force matrix
    matrix_assembler_->assemble_K_cor_matrix(mesh_, dt_, K_cor_entries_number_);
      // Assign corrected force
    mesh_->compute_nodal_corrected_force(
        matrix_assembler_->K_cor_matrix(),
        matrix_assembler_->pore_pressure_increment(), dt_);

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}
