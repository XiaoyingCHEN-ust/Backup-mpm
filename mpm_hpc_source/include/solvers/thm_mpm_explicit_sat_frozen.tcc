//! Constructor
template <unsigned Tdim>
mpm::THMMPMExplicitSatFrozen<Tdim>::THMMPMExplicitSatFrozen(
    const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("THMMPMExplicitSatFrozen");
}

//! THM-MPM semi-implicit two phase with phase change solver
template <unsigned Tdim>
bool mpm::THMMPMExplicitSatFrozen<Tdim>::solve() {
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
  const unsigned pore_ice = mpm::ParticlePhase::Solid;
  const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
  const unsigned mixture = mpm::ParticlePhase::Mixture;

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Pressure smoothing
  if (analysis_.find("pressure_smoothing") != analysis_.end())
    pressure_smoothing_ = analysis_["pressure_smoothing"].template get<bool>();

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
    throw std::runtime_error("Initialisation of loads failed");
  }

  // Initialise vtk output for two phase
  this->initialise_vtk_twophase();

  // Assign porosity
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_liquid_saturation_degree, std::placeholders::_1));      
      
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
    std::cout << "current_time" << this->current_time_ << "\n";
  }
      
  auto solver_begin = std::chrono::steady_clock::now();

  this->compute_critical_timestep_size(dt_);

  // Main loop
  for (; step_ <= nsteps_; ++step_) {

    if (variable_timestep_) {
      if (dt_matrix_size == 2) {
        if (step_ <= nsteps_matrix_[0]) 
          dt_ += (dt_matrix_[1] - dt_matrix_[0]) / nsteps_matrix_[0];
        else 
          dt_ = dt_matrix_[1];
      }
      else if (dt_matrix_size == 3){
        if (step_ <= nsteps_matrix_[0]) 
          dt_ += (dt_matrix_[1] - dt_matrix_[0]) / nsteps_matrix_[0];
        else if (step_ > nsteps_matrix_[0] && step_ <= nsteps_matrix_[1])
          dt_ += (dt_matrix_[2] - dt_matrix_[1]) / (nsteps_matrix_[1] - nsteps_matrix_[0]);
        else 
          dt_ = dt_matrix_[2];
      }
    }

    current_time_ += dt_;
    // Record current time
    mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::record_time, std::placeholders::_1, current_time_));

    if (mpi_rank == 0) console_->info("uuid : [{}], Step: {} of {}, timestep = {}, time = {}.\n", 
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
////////            Compute nodal velocities and temperatures           ////////
////////////////////////////////////////////////////////////////////////////////

    // Assign mass and momentum to nodes
    // Solid skeleton
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1));

    // Apply particle velocity constraints
    mesh_->apply_moving_rigid_boundary(current_time_, dt_);

    // Compute nodal velocity at the begining of time step
    // Nodal velocity constraint are also applied
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity, 
                  std::placeholders::_1, dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
        
    // Apply particle and nodal velocity constraints
    mesh_->apply_velocity_constraints(current_time_);

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

      // Iterate over each particle to update particle volume
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::update_particle_volume, std::placeholders::_1)); 
    
      // // // Iterate over each particle to update material density of particle
      // mesh_->iterate_over_particles(std::bind(
      //     &mpm::ParticleBase<Tdim>::update_particle_density, std::placeholders::_1, dt_));

      // // Iterate over each particle to calculate particle porosity
      // mesh_->iterate_over_particles(std::bind(
      //     &mpm::ParticleBase<Tdim>::update_particle_porosity, std::placeholders::_1, dt_));
 
      // Iterate over each particle to update liquid water saturation
      mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::update_liquid_saturation,
                  std::placeholders::_1, dt_));  

      // Iterate over each particle to map particle scalers
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_scalers_to_nodes, std::placeholders::_1));

      // Compute nodal temperature
      mesh_->iterate_over_nodes_predicate(
          std::bind(&mpm::NodeBase<Tdim>::compute_scalers,
                    std::placeholders::_1),
          std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
          
      // Iterate over each particle to calculate particle PIC porosity/saturation
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_scalers_smoothing, std::placeholders::_1));

      // // Iterate over each particle to update permeability
      // mesh_->iterate_over_particles(
      //   std::bind(&mpm::ParticleBase<Tdim>::update_viscosity,
      //             std::placeholders::_1));

      // // Iterate over each particle to update permeability
      // mesh_->iterate_over_particles(
      //   std::bind(&mpm::ParticleBase<Tdim>::update_permeability,
      //             std::placeholders::_1));

      // Iterate over each particle to compute stress of soil skeleton
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::update_particle_stress, std::placeholders::_1));

      // Iterate over each particle to compute pore pressure
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_pore_pressure,
                    std::placeholders::_1, dt_));

      mesh_->apply_particle_pore_pressure_constraints(current_time_); 

      // Pressure smoothing
      if (pressure_smoothing_) { this->pressure_smoothing(pore_liquid);}

      mesh_->apply_particle_pore_pressure_constraints(current_time_); 
    }

     // Compute free surface cells, nodes, and particles
    mesh_->compute_free_surface(free_surface_particle_, volume_tolerance_);
       
////////////////////////////////////////////////////////////////////////////////
////////               Update nodal temperatures           ////////
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

    // Iterate over each particle to compute nodal latent heat
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_latent_heat, std::placeholders::_1, this->dt_));

    // Iterate over each particle to compute nodal heat conduction
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_heat_conduction, std::placeholders::_1));

    // Iterate over each particle to compute nodal heat convection
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_heat_convection, std::placeholders::_1));

    // Iterate over each particle to compute nodal plastic work
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_plastic_work, std::placeholders::_1, dt_));

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

    Eigen::Matrix<double, Tdim, 1> gravity;
    if (current_time_<= 1.) 
      gravity = this->gravity_ * current_time_ / 1.;
    else gravity = this->gravity_;

    // std:cout << "this->gravity_" << gravity << "\n";
    // Iterate over particles to compute nodal body force of soil skeleton
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_external_force,
                  std::placeholders::_1, this->gravity_));
    // Apply particle traction and map to nodes
    mesh_->apply_traction_on_particles(current_time_);
  
    // Iterate over particles to compute nodal mixture internal force
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_internal_force_semi,
                  std::placeholders::_1, 1.0));

    // Iterate over particles to compute nodal drag force coefficient
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_drag_force_coefficient,
                  std::placeholders::_1));

    // Iterate over active nodes to compute acceleratation and velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<
                      Tdim>::compute_acc_vel_twophase_explicit,
                  std::placeholders::_1, soil_skeleton, pore_liquid, mixture,
                  this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1)); 

    // Apply particle and nodal velocity constraints
    mesh_->apply_velocity_constraints(current_time_);

// /////////////////////////////////////////////////////////////////////////////////
// ////////              Update particle velocities, pressures              ////////
// /////////////////////////////////////////////////////////////////////////////////
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

    // Write reaction force file
    this->write_reaction_force(false, this->step_, this->nsteps_);

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

// TODO: This is a copy of pressure_smoothing in explicit two-phase
//! MPM Explicit two-phase pressure smoothing
template <unsigned Tdim>
void mpm::THMMPMExplicitSatFrozen<Tdim>::pressure_smoothing(unsigned phase) {

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

// Compute time step size
template <unsigned Tdim>
void mpm::THMMPMExplicitSatFrozen<Tdim>::compute_critical_timestep_size(double dt) {
  const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
  const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
  const unsigned mixture = mpm::ParticlePhase::Mixture;

  // cell minimum size
  auto mesh_props = io_->json_object("mesh");
  // Get Mesh reader from JSON object
  double cellsize_min = mesh_props.at("cellsize_min").template get<double>();
  double cell_size = 0.1;
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
  double liquid_specific_heat = liquid_materials->template property<double>(std::string("specific_heat"));
  double liquid_thermal_conductivity = liquid_materials->template property<double>(std::string("thermal_conductivity"));
  double ice_density = liquid_materials->template property<double>(std::string("ice_density"));
  double ice_specific_heat = liquid_materials->template property<double>(std::string("ice_specific_heat"));
  double ice_thermal_conductivity = liquid_materials->template property<double>(std::string("ice_thermal_conductivity"));

  // Compute timestep for momentum eqaution 
  double density_mixture1 = (1 - porosity) * density;
  double density_mixture2 = (1 - porosity) * density + porosity * liquid_density;
  double density_mixture3 = (1 - porosity) * density + porosity * ice_density;
  double critical_dt11 = cellsize_min / std::pow(youngs_modulus/density_mixture1, 0.5);
  double critical_dt12 = cellsize_min / std::pow(youngs_modulus/density_mixture2, 0.5);
  double critical_dt13 = cellsize_min / std::pow(youngs_modulus/density_mixture3, 0.5);
  console_->info("Critical time step size for elastic wave propagation (solid base) is {} s", critical_dt11);
  console_->info("Critical time step size for elastic wave propagation (liquid base) is {} s", critical_dt12);
  console_->info("Critical time step size for elastic wave propagation (ice base) is {} s", critical_dt13);  

  // Compute timestep for heat transfer eqaution - pure liquid
  double k_mixture1 = (1 - porosity) * thermal_conductivity + porosity * liquid_thermal_conductivity;
  double c_mixture1 = (1 - porosity) * density * specific_heat + porosity * liquid_density * liquid_specific_heat;
  double critical_dt21 = cellsize_min * cellsize_min * c_mixture1 / k_mixture1;                                  
  console_->info("Critical time step size for thermal conduction equation (liquid base) is {} s", critical_dt21);
  
  // Compute timestep for heat transfer eqaution - pure ice
  double k_mixture2 = (1 - porosity) * thermal_conductivity + porosity * ice_thermal_conductivity;
  double c_mixture2 = (1 - porosity) * density * specific_heat + porosity * ice_density * ice_specific_heat;
  double critical_dt22 = cellsize_min * cellsize_min * c_mixture2 / k_mixture2;                                
  console_->info("Critical time step size for thermal conduction equation (ice base) is {} s", critical_dt22);

  if (dt >= std::min(critical_dt11, critical_dt12) || dt >= critical_dt13 || dt >= std::min(critical_dt21, critical_dt22))
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
bool mpm::THMMPMExplicitSatFrozen<Tdim>::pre_process() {
  bool status = true;
  console_->info("MPM analysis type {}", io_->analysis_type());

    // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Pressure smoothing
  if (analysis_.find("pressure_smoothing") != analysis_.end())
    pressure_smoothing_ = analysis_["pressure_smoothing"].template get<bool>();

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
    throw std::runtime_error("Initialisation of loads failed");
  }

  // Initialise vtk output for two phase
  this->initialise_vtk_twophase();

  // Assign porosity
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_liquid_saturation_degree, std::placeholders::_1));      
      
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
    std::cout << "current_time" << this->current_time_ << "\n";
  }     
  this->compute_critical_timestep_size(dt_);
  auto solver_begin = std::chrono::steady_clock::now();

  return status;
}

// Main loop
// (1) Compute deformation gradient for MPM-DEM
template <unsigned Tdim>
bool mpm::THMMPMExplicitSatFrozen<Tdim>::get_deformation_task() {
  bool status = true;

  // Main loop
  if (step_ < nsteps_) {

    if (variable_timestep_) {
      if (dt_matrix_size == 2) {
        if (step_ <= nsteps_matrix_[0]) 
          dt_ += (dt_matrix_[1] - dt_matrix_[0]) / nsteps_matrix_[0];
        else 
          dt_ = dt_matrix_[1];
      }
      else if (dt_matrix_size == 3){
        if (step_ <= nsteps_matrix_[0]) 
          dt_ += (dt_matrix_[1] - dt_matrix_[0]) / nsteps_matrix_[0];
        else if (step_ > nsteps_matrix_[0] && step_ <= nsteps_matrix_[1])
          dt_ += (dt_matrix_[2] - dt_matrix_[1]) / (nsteps_matrix_[1] - nsteps_matrix_[0]);
        else 
          dt_ = dt_matrix_[2];
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
////////            Compute nodal velocities and temperatures           ////////
////////////////////////////////////////////////////////////////////////////////

    // Apply particle and nodal velocity constraints
    mesh_->apply_velocity_constraints(current_time_);

    // Assign mass and momentum to nodes
    // Solid skeleton
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1));               

    // Apply particle velocity constraints
    mesh_->apply_moving_rigid_boundary(current_time_, dt_);

    // Compute nodal velocity at the begining of time step
    // Nodal velocity constraint are also applied
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity, 
                  std::placeholders::_1, dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));  

    // Apply particle velocity constraints
    mesh_->apply_moving_rigid_boundary(current_time_, dt_);    

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

      // Iterate over each particle to update particle volume
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::update_particle_volume, std::placeholders::_1)); 
    
      // // // Iterate over each particle to update material density of particle
      // mesh_->iterate_over_particles(std::bind(
      //     &mpm::ParticleBase<Tdim>::update_particle_density, std::placeholders::_1, dt_));

      // Iterate over each particle to calculate particle porosity
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::update_particle_porosity, std::placeholders::_1, dt_));
 
      // Iterate over each particle to update liquid water saturation
      mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::update_liquid_saturation,
                  std::placeholders::_1, dt_));  

      // Iterate over each particle to map particle scalers
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_scalers_to_nodes, std::placeholders::_1));

      // Compute nodal temperature
      mesh_->iterate_over_nodes_predicate(
          std::bind(&mpm::NodeBase<Tdim>::compute_scalers,
                    std::placeholders::_1),
          std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
          
      // Iterate over each particle to calculate particle PIC porosity/saturation
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_scalers_smoothing, std::placeholders::_1));

      // // Iterate over each particle to update permeability
      // mesh_->iterate_over_particles(
      //   std::bind(&mpm::ParticleBase<Tdim>::update_viscosity,
      //             std::placeholders::_1));

      // Iterate over each particle to update permeability
      mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::update_permeability,
                  std::placeholders::_1));

      // Iterate over each particle to compute pore pressure
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_pore_pressure,
                    std::placeholders::_1, dt_));

      mesh_->apply_particle_pore_pressure_constraints(current_time_); 

      // Pressure smoothing
      if (pressure_smoothing_) { this->pressure_smoothing(pore_liquid);}

      mesh_->apply_particle_pore_pressure_constraints(current_time_); 
    }

     // Compute free surface cells, nodes, and particles
    mesh_->compute_free_surface(free_surface_particle_, volume_tolerance_); 

////////////////////////////////////////////////////////////////////////////////
////////               Update nodal temperatures           ////////
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

    // Iterate over each particle to compute nodal latent heat
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_latent_heat, std::placeholders::_1, this->dt_));        

    // Iterate over each particle to compute nodal heat conduction
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_heat_conduction, std::placeholders::_1));       

    // Iterate over each particle to compute nodal heat convection
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_heat_convection, std::placeholders::_1));

    // // Iterate over each particle to compute nodal heat convection
    // mesh_->iterate_over_particles(std::bind(
    //     &mpm::ParticleBase<Tdim>::map_covective_heat_flux, std::placeholders::_1, current_time_));    

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
void mpm::THMMPMExplicitSatFrozen<Tdim>::get_info(unsigned& dim, bool& resume,
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
void mpm::THMMPMExplicitSatFrozen<Tdim>::get_status(double& dt, unsigned& step,
                                                    unsigned& nsteps,
                                                    unsigned& output_steps) {
  dt = dt_;
  step = step_;
  nsteps = nsteps_;
  output_steps = output_steps_;
}

// (4) Get deformation for DEM
template <unsigned Tdim>
bool mpm::THMMPMExplicitSatFrozen<Tdim>::send_deformation_task(
    std::vector<unsigned>& id,
    std::vector<Eigen::MatrixXd>& displacement_gradients) {
  bool status = true;
  mesh_->get_displacement_gradient(id, displacement_gradients);
  return status;
}

template <unsigned Tdim>
bool mpm::THMMPMExplicitSatFrozen<Tdim>::send_temperature_task(
    std::vector<unsigned>& id,
    std::vector<double>& particle_temperature) {
  bool status = true;
  mesh_->get_particle_temperature(id, particle_temperature);
  return status;
}

// Set particle stress
template <unsigned Tdim>
bool mpm::THMMPMExplicitSatFrozen<Tdim>::set_stress_task(
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
bool mpm::THMMPMExplicitSatFrozen<Tdim>::set_porosity_task(
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
bool mpm::THMMPMExplicitSatFrozen<Tdim>::set_fabric_task(
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
bool mpm::THMMPMExplicitSatFrozen<Tdim>::set_rotation_task(
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
bool mpm::THMMPMExplicitSatFrozen<Tdim>::update_state_task() {
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

// if (current_time_<= 0.1) this->gravity_(1) = -981 * current_time_;
// else this->gravity_(1) = -98.1;

// std:cout << "this->gravity_" << this->gravity_ << "\n";

    // Iterate over particles to compute nodal body force of soil skeleton
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_external_force,
                  std::placeholders::_1, this->gravity_));
    // Apply particle traction and map to nodes
    mesh_->apply_traction_on_particles(current_time_);
  
    // Iterate over particles to compute nodal mixture internal force
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_internal_force,
                  std::placeholders::_1));

    // Iterate over particles to compute nodal drag force coefficient
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_drag_force_coefficient,
                  std::placeholders::_1));

    // Iterate over active nodes to compute acceleratation and velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<
                      Tdim>::compute_acc_vel_twophase_explicit,
                  std::placeholders::_1, soil_skeleton, pore_liquid, mixture,
                  this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1)); 

    // Apply particle and nodal velocity constraints
    mesh_->apply_velocity_constraints(current_time_);   

// /////////////////////////////////////////////////////////////////////////////////
// ////////              Update particle velocities, pressures              ////////
// /////////////////////////////////////////////////////////////////////////////////
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

    // Write reaction force file
    this->write_reaction_force(false, this->step_, this->nsteps_);

    // Locate particles
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty()) {
      status = false;
      throw std::runtime_error("Particle outside the mesh domain");
    }

    // Fixed timestep, and data output linearly (every const steps = output_steps)
    if ((!variable_timestep_) & (!log_output_steps_)){
      if (this->step_ % this->output_steps_ == 0) {
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
          "Rank {}, THMMPMExplicitSatFrozen {} solver duration: {} ms", 0,
          (this->stress_update_ == mpm::StressUpdate::USL ? "USL" : "USF"),
          std::chrono::duration_cast<std::chrono::milliseconds>(solver_end -
                                                                solver_begin)
              .count());
    }
  }
  return status;
}
