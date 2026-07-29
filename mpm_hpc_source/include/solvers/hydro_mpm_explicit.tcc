//! Constructor
template <unsigned Tdim>
mpm::HydroMPMExplicit<Tdim>::HydroMPMExplicit(const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("HydroMPMExplicit");
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////                                                                ////////
////////                 TM-MPM Explicit functions                      ////////
////////                                                                ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//! Domain decomposition
template <unsigned Tdim>
void mpm::HydroMPMExplicit<Tdim>::mpi_domain_decompose() {
#ifdef USE_MPI
  // Initialise MPI rank and size
  int mpi_rank = 0;
  int mpi_size = 1;

  // Get MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  if (mpi_size > 1 && mesh_->ncells() > 1) {

    // Initialise MPI
    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    auto mpi_domain_begin = std::chrono::steady_clock::now();
    console_->info("Rank {}, Domain decomposition started\n", mpi_rank);

    // Check if mesh has cells to partition
    if (mesh_->ncells() == 0)
      throw std::runtime_error("Container of cells is empty");

#ifdef USE_GRAPH_PARTITIONING
    // Create graph
    graph_ = std::make_shared<Graph<Tdim>>(mesh_->cells(), mpi_size, mpi_rank);

    // Graph partitioning mode
    int mode = 4;  // FAST
    // Create graph partition
    bool graph_partition = graph_->create_partitions(&comm, mode);
    // Collect the partitions
    graph_->collect_partitions(mpi_size, mpi_rank, &comm);

    // Delete all the particles which is not in local task parititon
    mesh_->remove_all_nonrank_particles();
    // Identify shared nodes across MPI domains
    mesh_->find_domain_shared_nodes();
    // Identify ghost boundary cells
    mesh_->find_ghost_boundary_cells();
#endif
    auto mpi_domain_end = std::chrono::steady_clock::now();
    console_->info("Rank {}, Domain decomposition: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       mpi_domain_end - mpi_domain_begin)
                       .count());
  }
#endif  // MPI
}

//! MPM Explicit pressure smoothing
template <unsigned Tdim>
void mpm::HydroMPMExplicit<Tdim>::pressure_smoothing(unsigned phase) {
  // Assign pressure to nodes
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::map_pressure_to_nodes,
                std::placeholders::_1, current_time_));

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

  // Smooth pressure over particles
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_smoothing,
                std::placeholders::_1));
}

//! MPM Explicit compute stress strain
template <unsigned Tdim>
void mpm::HydroMPMExplicit<Tdim>::compute_stress_strain(unsigned phase) {
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

  // Pressure smoothing
  if (pressure_smoothing_) this->pressure_smoothing(phase);
  
  // Iterate over each particle to compute stress
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::update_particle_stress, std::placeholders::_1));
}

// Compute time step size
template <unsigned Tdim>
void mpm::HydroMPMExplicit<Tdim>::compute_critical_timestep_size(double dt) {
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
  console_->info("Critical time step size for hydraulic conduction is {} s", critical_dt);

  // if (dt >= critical_dt)
  //     throw std::runtime_error("Time step size is too large");  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////                                                                ////////
////////                  TM-MPM Explicit Solver                        ////////
////////                                                                ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//! Thermo-mechanical MPM Explicit solver
template <unsigned Tdim>
bool mpm::HydroMPMExplicit<Tdim>::solve() {
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

  // Phase
  const unsigned phase = mpm::ParticlePhase::Liquid;

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Pressure smoothing
  if (analysis_.find("pressure_smoothing") != analysis_.end())
    pressure_smoothing_ =
        analysis_.at("pressure_smoothing").template get<bool>();

  // Variable timestep
  if (analysis_.find("variable_timestep") != analysis_.end())
    variable_timestep_ = analysis_["variable_timestep"].template get<bool>();

  // Log output steps
  if (post_process_.find("log_output_steps") != post_process_.end())
    log_output_steps_ =post_process_["log_output_steps"].template get<bool>();  

  // Interface
  if (analysis_.find("interface") != analysis_.end())
    interface_ = analysis_.at("interface").template get<bool>();

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


  auto solver_begin = std::chrono::steady_clock::now();

  this->compute_critical_timestep_size(dt_);

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

    if (mpi_rank == 0) console_->info("Step: {} of {}, timestep = {}, time = {}.\n", 
                                       step_, nsteps_, dt_, current_time_);

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

    // Spawn a task for initialising nodes and cells
      // Initialise nodes
      mesh_->iterate_over_nodes(
          std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

      mesh_->iterate_over_cells(
          std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));

    // Spawn a task for particles
      // Iterate over each particle to compute shapefn
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));

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
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_pressure_to_nodes,
                  std::placeholders::_1));
    
    // Compute nodal pressure
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_pressure,
                  std::placeholders::_1, phase),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));  

    // Apply nodal pressure constraints
    mesh_->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Tdim>::apply_pressure_constraints,
                std::placeholders::_1, phase, current_time_),
      std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Iterate over each particle to compute nodal heat conduction
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_hydraulic_conduction, std::placeholders::_1));

    // Compute nodal pressure acceleration and update nodal pressure
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_pressure,
                  std::placeholders::_1, phase, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
    
    // Apply nodal pressure constraints
    mesh_->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Tdim>::apply_pressure_constraints,
                std::placeholders::_1, phase, current_time_),
      std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Iterate over each particle to compute updated temperature
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_pore_pressure,
        std::placeholders::_1, this->dt_, pic_t_));

    // Locate particles
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
    mesh_->transfer_nonrank_particles();
#endif
#endif

    if (step_ % output_steps_ == 0) {
      // HDF5 outputs
      if (write_hdf5_) this->write_hdf5(this->step_ + 1, this->nsteps_);
#ifdef USE_VTK
      // VTK outputs
      this->write_vtk(this->step_, this->nsteps_);
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