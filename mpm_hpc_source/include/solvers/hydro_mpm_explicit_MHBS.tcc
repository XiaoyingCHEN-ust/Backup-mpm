//! Constructor
template <unsigned Tdim>
mpm::HydroMPMExplicitMHBS<Tdim>::HydroMPMExplicitMHBS(
    const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("HydroMPMExplicitMHBS");
}

//! Thermo-hydro-mechncial MPM Explicit solver
template <unsigned Tdim>
bool mpm::HydroMPMExplicitMHBS<Tdim>::solve() {
  bool status = true;

  console_->info("MPM analysis type {}", io_->analysis_type());

  // Initialise MPI rank and size
  int mpi_rank = 0;
  int mpi_size = 1;

  // Two phases (soil skeleton and pore liquid)
  const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
  const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
  const unsigned mixture = mpm::ParticlePhase::Mixture;
  const unsigned pore_gas = mpm::ParticlePhase::Gas;  

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

  this->compute_critical_timestep_size(dt_);

  // Main loop
  for (step_ = 0; step_ <= nsteps_; ++step_) {

    current_time_ += dt_;

    // Record current time
    mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::record_time, std::placeholders::_1, current_time_));

    console_->info("uuid : [{}], Step: {} of {}, timestep = {}, time = {}.\n", 
                                       uuid_, step_, nsteps_, dt_, current_time_);

    // Fixed timestep, and data output linearly (every const steps = output_steps)

    if (step_ == 0) {
      // HDF5 outputs
      if (write_hdf5_) this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
      // VTK outputs
      this->write_vtk(this->step_, this->nsteps_);
#endif

      std::cout << "Number output = " << No_output << "\n";
      No_output++;
    }

    // Initialise nodes
    mesh_->iterate_over_nodes(
        std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

    mesh_->iterate_over_cells(
        std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));

    // Iterate over each particle to compute shapefn
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));  

  mesh_->iterate_over_cells(std::bind(
            &mpm::Cell<Tdim>::map_cell_volume_to_nodes, std::placeholders::_1, 0));

    // Iterate over each particle to update material density of particle
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_permeability, std::placeholders::_1));

    // Iterate over each particle to update material density of particle
    // mesh_->iterate_over_particles(std::bind(
    //     &mpm::ParticleBase<Tdim>::update_source_term, std::placeholders::_1));

    // Assign mass and momentum to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1));

    // Assign heat capacity and heat to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_heat_to_nodes,
                  std::placeholders::_1));

    // Assign heat capacity and heat to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_pressure_to_nodes,
                  std::placeholders::_1));

    // Assign heat capacity and heat to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_density_to_nodes,
                  std::placeholders::_1));

    // Compute nodal temperature
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_temperature,
                  std::placeholders::_1, soil_skeleton),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));  

    // Apply nodal temperature constraints
    mesh_->apply_nodal_temperature_constraints(soil_skeleton, current_time_);

    // Apply nodal temperature constraints
    mesh_->apply_nodal_convective_heat_constraints(soil_skeleton, current_time_);    
    
    // Compute nodal pressure
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_pressure,
                  std::placeholders::_1, pore_liquid),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Compute nodal pressure
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_pressure,
                  std::placeholders::_1, pore_gas),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Apply nodal pressure constraints
    mesh_->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Tdim>::apply_pressure_constraints,
                std::placeholders::_1, pore_liquid, current_time_),
      std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Iterate over each particle to compute nodal heat conduction
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_heat_conduction, std::placeholders::_1));

    // Iterate over each particle to compute nodal heat convection
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_heat_source, std::placeholders::_1));

    // Iterate over each particle to compute nodal heat conduction
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_hydraulic_conduction, std::placeholders::_1));

    // Iterate over each particle to compute nodal heat convection
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_heat_convection, std::placeholders::_1));

    // Iterate over each particle to compute nodal heat convection
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_mass_convection, std::placeholders::_1));                

    // Compute nodal temperature acceleration and update nodal temperature
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_temperature,
                  std::placeholders::_1, soil_skeleton, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Apply nodal temperature constraints
    mesh_->apply_nodal_temperature_constraints(soil_skeleton, current_time_);

    // Compute nodal pressure acceleration and update nodal pressure
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_pressure_threephase,
                  std::placeholders::_1, pore_liquid, pore_gas, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
    
    // Apply nodal pressure constraints
    mesh_->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Tdim>::apply_pressure_constraints,
                std::placeholders::_1, pore_liquid, current_time_),
      std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Iterate over each particle to compute updated temperature
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_temperature,
        std::placeholders::_1, this->dt_, this->pic_t_));

    // Iterate over each particle to compute updated temperature
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_pore_pressure,
        std::placeholders::_1, this->dt_, pic_t_)); 

    // Iterate over each particle to calculate particle porosity
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_liquid_saturation, std::placeholders::_1, dt_));

    // Iterate over each particle to calculate particle porosity
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_hydrate_saturation, std::placeholders::_1, dt_));        

    // Iterate over each particle to update material density of particle
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_density, std::placeholders::_1, dt_));

    // Iterate over each particle to update material density of particle
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_volume, std::placeholders::_1)); 

    // Locate particles
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

    if (step_ % output_steps_ == 0) {
      // HDF5 outputs
      if (write_hdf5_) this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
      // VTK outputs
      this->write_vtk(this->step_ + 1, this->nsteps_);
#endif
      std::cout << "Number output = " << No_output << "\n";
      No_output++;
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

// Compute time step size
template <unsigned Tdim>
void mpm::HydroMPMExplicitMHBS<Tdim>::compute_critical_timestep_size(double dt) {
  const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
  const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
  const unsigned pore_gas = mpm::ParticlePhase::Gas;
  const unsigned mixture = mpm::ParticlePhase::Mixture;

  const unsigned phase = 0;
  // cell minimum size
  auto mesh_props = io_->json_object("mesh");
  // Get Mesh reader from JSON object
  double cellsize_min = mesh_props.at("cellsize_min").template get<double>();
  // Material parameters 
  auto materials =  materials_.at(0);
  auto liquid_materials =  materials_.at(1);  
  double liquid_modulus = liquid_materials->template property<double>(std::string("bulk_modulus"));
  double density = liquid_materials->template property<double>(std::string("density"));
  double porosity = materials->template property<double>(std::string("porosity"));
  double permeability = materials->template property<double>(std::string("intrinsic_permeability"));
  double viscosity = liquid_materials->template property<double>(std::string("liquid_viscosity")); 

  // Compute timestep for heat transfer eqaution - pure liquid
  double k_mixture = permeability / viscosity * liquid_modulus * density;
  double c_mixture = (1 - porosity) * density;
  double critical_dt = cellsize_min * cellsize_min * c_mixture / k_mixture;
  console_->info("Critical time step size for liquid conduction is {} s", critical_dt);
}


