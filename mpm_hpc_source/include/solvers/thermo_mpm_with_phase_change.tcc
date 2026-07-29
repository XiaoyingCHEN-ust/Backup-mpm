//! Constructor
template <unsigned Tdim>
mpm::ThermoMPMPhaseChange<Tdim>::ThermoMPMPhaseChange(
    const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("ThermoMPMPhaseChange");
}

//! THM-MPM semi-implicit two phase with phase change solver
template <unsigned Tdim>
bool mpm::ThermoMPMPhaseChange<Tdim>::solve() {
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

  // Projection method paramter (beta)
  if (analysis_.find("free_surface") != analysis_.end()) {
    free_surface_particle_ = analysis_["free_surface"]["free_surface_particle"]
                                 .template get<std::string>();                                      
    volume_tolerance_ =
          analysis_["free_surface"]["volume_tolerance"].template get<double>(); 
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

    if (mpi_rank == 0) console_->info("Step: {} of {}, timestep = {}, time = {}.\n", 
                                       step_, nsteps_, dt_, current_time_);

    // Initialise nodes
    mesh_->iterate_over_nodes(
        std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

    mesh_->iterate_over_cells(
        std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));

    // Iterate over each particle to compute shapefn
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));    

////////////////////////////////////////////////////////////////////////////////
////////                   Update liquid water saturation               ////////
////////////////////////////////////////////////////////////////////////////////

    // Update stress first
    if (this->stress_update_ == mpm::StressUpdate::USF) {            

      // Iterate over each particle to update particle volume
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::update_particle_volume, std::placeholders::_1));        
 
      // Iterate over each particle to update liquid water saturation
      mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::update_liquid_saturation,
                  std::placeholders::_1, dt_));          
    }

     // Compute free surface cells, nodes, and particles
    mesh_->compute_free_surface(free_surface_particle_, volume_tolerance_);

////////////////////////////////////////////////////////////////////////////////
////////                     Update nodal temperatures                  ////////
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

    // Iterate over each particle to compute nodal heat convection
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_covective_heat_flux, std::placeholders::_1, current_time_));    

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


/////////////////////////////////////////////////////////////////////////////////
////////                     Update particle temperature                 ////////
/////////////////////////////////////////////////////////////////////////////////
    
    // Iterate over each particle to compute updated temperature
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_temperature,
        std::placeholders::_1, this->dt_, this->pic_t_));
        
    // Apply particle temperature constraints
    mesh_->apply_particle_temperature_constraints(current_time_); 

    // Locate particles
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

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
void mpm::ThermoMPMPhaseChange<Tdim>::compute_critical_timestep_size(double dt) {
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
  double density = materials->template property<double>(std::string("density"));
  double specific_heat = materials->template property<double>(std::string("specific_heat"));
  double thermal_conductivity = materials->template property<double>(std::string("thermal_conductivity"));
  // Liquid Material parameters 
  auto liquid_materials =  materials_.at(pore_liquid);
  double liquid_density = liquid_materials->template property<double>(std::string("density"));
  double liquid_specific_heat = liquid_materials->template property<double>(std::string("specific_heat"));
  double liquid_thermal_conductivity = liquid_materials->template property<double>(std::string("thermal_conductivity"));
  double ice_density = liquid_materials->template property<double>(std::string("ice_density"));
  double ice_specific_heat = liquid_materials->template property<double>(std::string("ice_specific_heat"));
  double ice_thermal_conductivity = liquid_materials->template property<double>(std::string("ice_thermal_conductivity"));
  // Compute timestep for three phase MPM - pure liquid
  double k_mixture1 = (1 - porosity) * thermal_conductivity + porosity * liquid_thermal_conductivity;
  double c_mixture1 = (1 - porosity) * density * specific_heat + porosity * liquid_density * liquid_specific_heat;
  double critical_dt1 = cellsize_min * cellsize_min * c_mixture1 / k_mixture1;                                  
  console_->info("Critical time step size for thermal conduction equation (liquid base) is {} s", critical_dt1);
  
  // Compute timestep for three phase MPM - pure ice
  double k_mixture2 = (1 - porosity) * thermal_conductivity + porosity * ice_thermal_conductivity;
  double c_mixture2 = (1 - porosity) * density * specific_heat + porosity * ice_density * ice_specific_heat;
  double critical_dt2 = cellsize_min * cellsize_min * c_mixture2 / k_mixture2;                                
  console_->info("Critical time step size for thermal conduction equation (ice base) is {} s", critical_dt2); 
  
  if (dt > std::min(critical_dt1, critical_dt2))
      throw std::runtime_error("Time step size is too large");                             
}
