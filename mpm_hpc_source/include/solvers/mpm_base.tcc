//! Constructor
template <unsigned Tdim>
mpm::MPMBase<Tdim>::MPMBase(const std::shared_ptr<IO>& io) : mpm::MPM(io) {
  //! Logger
  console_ = spdlog::get("MPMBase");

  // Create a mesh with global id 0
  const mpm::Index id = 0;

  // Set analysis step to start at 0
  step_ = 0;

  // Set mesh as isoparametric
  bool isoparametric = is_isoparametric();

  mesh_ = std::make_unique<mpm::Mesh<Tdim>>(id, isoparametric);

  // Empty all materials
  materials_.clear();

  try {
    analysis_ = io_->analysis();
    // Time-step size
    dt_ = analysis_["dt"].template get<double>();
    // Number of time steps
    nsteps_ = analysis_["nsteps"].template get<mpm::Index>();

    // Variable time steps
    try {
        if (analysis_.find("dt_matrix") != analysis_.end() &&
            analysis_.find("nsteps_matrix") != analysis_.end()) {
            dt_matrix_size = analysis_.at("dt_matrix").size();
            dt_matrix_.resize(dt_matrix_size);
            nsteps_matrix_.resize(dt_matrix_size);
            dt_matrix_.setZero();
            nsteps_matrix_.setZero();
            for (unsigned i = 0; i < analysis_.at("dt_matrix").size(); ++i) {
                dt_matrix_[i] = analysis_.at("dt_matrix").at(i);
                nsteps_matrix_[i] = analysis_.at("nsteps_matrix").at(i);
                nsteps_ = analysis_.at("nsteps_matrix").back();
      }
  } 
    } catch (std::exception& exception) {
      console_->warn(
          "{} #{}: {}. Stress update method is not specified, using USF as "
          "default\n",
          __FILE__, __LINE__, exception.what());
    }

    // Stress update method (USF/USL/MUSL)
    try {
      if (analysis_.find("stress_update") != analysis_.end())
        stress_update_ = mpm::stress_update.at(
            analysis_["stress_update"].template get<std::string>());
    } catch (std::exception& exception) {
      console_->warn(
          "{} #{}: {}. Stress update method is not specified, using USF as "
          "default\n",
          __FILE__, __LINE__, exception.what());
    }

    // Velocity update
    try {
      pic_ = analysis_["PIC"].template get<double>();
    } catch (std::exception& exception) {
      console_->warn(
          "{} #{}: PIC  parameter is not specified, using default "
          "as FLIP PIC=0",
          __FILE__, __LINE__, exception.what());
      pic_ = 0.;
    }

    // Temperature update
    try {
      pic_t_ = analysis_["PIC_T"].template get<double>();
    } catch (std::exception& exception) {
      console_->warn(
          "{} #{}: PIC_T  parameter is not specified, using default "
          "as FLIP PIC_T=0",
          __FILE__, __LINE__, exception.what());
      pic_t_ = 0.;
    }

    // Solid affine particle-grid velocity transfer
    try {
      if (analysis_.find("affine_mpm") != analysis_.end())
        affine_mpm_ = analysis_["affine_mpm"].template get<bool>();
      else if (analysis_.find("APIC") != analysis_.end())
        affine_mpm_ = analysis_["APIC"].template get<bool>();
    } catch (std::exception& exception) {
      console_->warn(
          "{} #{}: {}. APIC parameter is invalid, using affine_mpm=false",
          __FILE__, __LINE__, exception.what());
      affine_mpm_ = false;
    }

    // Damping
    try {
      if (analysis_.find("damping") != analysis_.end()) {
        if (!initialise_damping(analysis_.at("damping")))
          throw std::runtime_error("Damping parameters are not defined");
      }
    } catch (std::exception& exception) {
      console_->warn("{} #{}: Damping is not specified, using none as default",
                     __FILE__, __LINE__, exception.what());
    }

    // Math functions
    try {
      // Get materials properties
      auto math_functions = io_->json_object("math_functions");
      if (!math_functions.empty())
        this->initialise_math_functions(math_functions);
    } catch (std::exception& exception) {
      console_->warn("{} #{}: No math functions are defined", __FILE__,
                     __LINE__, exception.what());
    }

    post_process_ = io_->post_processing();
    // Output steps
    output_steps_ = post_process_["output_steps"].template get<mpm::Index>();
    // Write hdf5 or vtk file or not
    write_hdf5_ = post_process_["write_hdf5"].template get<bool>();
    write_vtk_ = post_process_["write_vtk"].template get<bool>();

  } catch (std::domain_error& domain_error) {
    console_->error("{} {} Get analysis object: {}", __FILE__, __LINE__,
                    domain_error.what());
    abort();
  }
  // Initialise vtk output
  this->initialise_vtk();
}

// Initialise mesh
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_mesh() {

  bool status = true;

  try {
    // Initialise MPI rank and size
    int mpi_rank = 0;
    int mpi_size = 1;

#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    // Get number of MPI ranks
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

    // Get mesh properties
    auto mesh_props = io_->json_object("mesh");
    // Get Mesh reader from JSON object
    const std::string io_type =
        mesh_props["io_type"].template get<std::string>();

    bool check_duplicates = true;
    try {
      check_duplicates = mesh_props["check_duplicates"].template get<bool>();
    } catch (std::exception& exception) {
      console_->warn(
          "{} #{}: Check duplicates, not specified setting default as true",
          __FILE__, __LINE__, exception.what());
      check_duplicates = true;
    }

    // Create a mesh reader
    auto mesh_io = Factory<mpm::IOMesh<Tdim>>::instance()->create(io_type);

    auto nodes_begin = std::chrono::steady_clock::now();
    // Global Index
    mpm::Index gid = 0;
    // Node type
    const auto node_type = mesh_props["node_type"].template get<std::string>();

    // Mesh file
    std::string mesh_file =
        io_->file_name(mesh_props["mesh"].template get<std::string>());

    // Create nodes from file
    bool node_status =
        mesh_->create_nodes(gid,                                  // global id
                            node_type,                            // node type
                            mesh_io->read_mesh_nodes(mesh_file),  // coordinates
                            check_duplicates);                    // check dups

    if (!node_status) {
      status = false;
      throw std::runtime_error("Addition of nodes to mesh failed");
    }

    auto nodes_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Read nodes: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       nodes_end - nodes_begin)
                       .count());

    // Read and assign node sets
    this->node_entity_sets(mesh_props, check_duplicates);

    // Read nodal euler angles and assign rotation matrices
    this->node_euler_angles(mesh_props, mesh_io);

    // Read and assign velocity constraints
    this->nodal_velocity_constraints(mesh_props, mesh_io);

    // Read and assign friction constraints
    this->nodal_frictional_constraints(mesh_props, mesh_io);

    // Read and assign pore pressure constraints
    this->nodal_pore_pressure_constraints(mesh_props, mesh_io);

    // Read and assign temperature constraints
    this->nodal_temperature_constraints(mesh_props, mesh_io);

    // Read and assign temperature constraints
    this->nodal_convective_heat_constraints(mesh_props, mesh_io);    

    // Read and assign temperature constraints
    this->nodal_heat_source(mesh_props, mesh_io);

    //! TODO FIX water table
    // // Initialise water table
    // try {
    //   // Get water table properties
    //   auto water_table = io_->json_object("water_table");
    //   if (!water_table.empty())
    //   this->initialise_nodal_water_table(water_table);
    // } catch (std::exception& exception) {
    //   console_->warn("#{}: No water tables are defined ", __LINE__);
    // }

    // Initialise cell
    auto cells_begin = std::chrono::steady_clock::now();
    // Shape function name
    const auto cell_type = mesh_props["cell_type"].template get<std::string>();
    // Shape function
    std::shared_ptr<mpm::Element<Tdim>> element =
        Factory<mpm::Element<Tdim>>::instance()->create(cell_type);

    // Create cells from file
    bool cell_status =
        mesh_->create_cells(gid,      // global id
                            element,  // element type
                            mesh_io->read_mesh_cells(mesh_file),  // Node ids
                            check_duplicates);                    // Check dups

    if (!cell_status) {
      status = false;
      throw std::runtime_error("Addition of cells to mesh failed");
    }

    // Compute cell neighbours
    mesh_->compute_cell_neighbours();

    // Read and assign cell sets
    this->cell_entity_sets(mesh_props, check_duplicates);

    auto cells_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Read cells: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       cells_end - cells_begin)
                       .count());
  } catch (std::exception& exception) {
    console_->error("#{}: Reading mesh: {}", __LINE__, exception.what());
  }

  // Terminate if mesh creation failed
  if (!status) throw std::runtime_error("Initialisation of mesh failed");
  return status;
}

// Initialise particles
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_particles() {

  bool status = true;

  try {
    // Initialise MPI rank and size
    int mpi_rank = 0;
    int mpi_size = 1;

#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    // Get number of MPI ranks
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

    // Get mesh properties
    auto mesh_props = io_->json_object("mesh");
    // Get Mesh reader from JSON object
    const std::string io_type =
        mesh_props["io_type"].template get<std::string>();

    bool check_duplicates = true;
    try {
      check_duplicates = mesh_props["check_duplicates"].template get<bool>();
    } catch (std::exception& exception) {
      console_->warn(
          "{} #{}: Check duplicates, not specified setting default as true",
          __FILE__, __LINE__, exception.what());
      check_duplicates = true;
    }

    auto particles_gen_begin = std::chrono::steady_clock::now();

    // Get particles properties
    auto json_particles = io_->json_object("particles");

    for (const auto& json_particle : json_particles) {
      // Generate particles
      bool gen_status =
          mesh_->generate_particles(io_, json_particle["generator"]);
      if (!gen_status) status = false;
    }

    auto particles_gen_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Generate particles: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       particles_gen_end - particles_gen_begin)
                       .count());

    auto particles_locate_begin = std::chrono::steady_clock::now();

    // Create a mesh reader
    auto particle_io = Factory<mpm::IOMesh<Tdim>>::instance()->create(io_type);

    // Read and assign particles cells
    this->particles_cells(mesh_props, particle_io);

    // Locate particles in cell
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

    // Write particles and cells to file
    particle_io->write_particles_cells(
        io_->output_file("particles-cells", ".txt", uuid_, 0, 0).string(),
        mesh_->particles_cells());

    auto particles_locate_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Locate particles: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       particles_locate_end - particles_locate_begin)
                       .count());

    auto particles_volume_begin = std::chrono::steady_clock::now();
    

    // Read and assign particles volumes
    this->particles_volumes(mesh_props, particle_io);



    // Check if axisymmetric
    if (analysis_.find("isaxisymmetric") != analysis_.end())
    is_axisymmetric_ = analysis_["isaxisymmetric"].template get<bool>();
     
    // Compute volume
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_volume, std::placeholders::_1, 
        is_axisymmetric_));

    // Read and assign particles stresses
    this->particles_stresses(mesh_props, particle_io);

    // Read and assign particles initial pore pressure
    this->particles_pore_pressures(mesh_props, particle_io);

    // Read and assign particles initial temperature
    this->particles_temperatures(mesh_props, particle_io);  

    auto particles_volume_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Read volume, velocity and stresses: {} ms",
                   mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       particles_volume_end - particles_volume_begin)
                       .count());

    // Particle entity sets
    auto particles_sets_begin = std::chrono::steady_clock::now();
    this->particle_entity_sets(mesh_props, check_duplicates);
    this->particle_materials(mesh_props);
    auto particles_sets_end = std::chrono::steady_clock::now();

    // Read and assign particles velocity constraints
    this->particle_velocity_constraints(mesh_props, particle_io);

    // Read and assign particles temperature constraints
    this->particle_temperature_constraints(mesh_props, particle_io);

    // Read and assign particles pore pressure constraints
    this->particle_pore_pressure_constraints(mesh_props, particle_io);        

    console_->info("Rank {} Create particle sets: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       particles_volume_end - particles_volume_begin)
                       .count());

  } catch (std::exception& exception) {
    console_->error("#{}: MPM Base generating particles: {}", __LINE__,
                    exception.what());
    status = false;
  }
  if (!status) throw std::runtime_error("Initialisation of particles failed");
  return status;
}

// Initialise materials
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_materials() {
  bool status = true;
  try {
    // Get materials properties
    auto materials = io_->json_object("materials");

    for (const auto material_props : materials) {
      // Get material type
      const std::string material_type =
          material_props["type"].template get<std::string>();

      // Get material id
      auto material_id = material_props["id"].template get<unsigned>();

      // Create a new material from JSON object
      auto mat =
          Factory<mpm::Material<Tdim>, unsigned, const Json&>::instance()
              ->create(material_type, std::move(material_id), material_props);

      // Add material to list
      auto result = materials_.insert(std::make_pair(mat->id(), mat));

      // If insert material failed
      if (!result.second) {
        status = false;
        throw std::runtime_error(
            "New material cannot be added, insertion failed");
      }
    }
    // Copy materials to mesh
    mesh_->initialise_material_models(this->materials_);
  } catch (std::exception& exception) {
    console_->error("#{}: Reading materials: {}", __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Initialise vtk output for single phase
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_vtk() {

  bool status = true;

  // Default VTK attributes for single phase
  std::vector<std::string> vtk = {"ids",
                                  "current_time",
                                  "porosities",
                                  "PIC_porosities",
                                  "volumes",
                                  "PIC_volumes",
                                  "masses",
                                  "solid_fractions",        
                                  "densities",         
                                  "temperatures", 
                                  "temperature_gradients",
                                  "mass_gradients",
                                  "outward_normals",
                                  "PIC_temperatures",  
                                  "temperature_accelerations",
                                  "volumetric_strains", 
                                  "PIC_volumetric_strains",
                                  "dvolumetric_strains",
                                  "thermal_volumetric_strains",
                                  "dthermal_volumetric_strains",
                                  "deviatoric_strains",
                                  "PIC_deviatoric_strains",
                                  "deviatoric_stresses",
                                  "mean_stresses",
                                  "displacements",         
                                  "rotations",             
                                  "velocities",
                                  "accelerations",            
                                  "stresses",              
                                  "strains",               
                                  "thermal_strains",               
                                  "fabric_CNs",            
                                  "fabric_POs",
                                  "velocity_gradients",            
                                  "displacement_gradients",
                                  "deformation_gradients",
                                  "free_surfaces",
                                  "grad_shapefns",
                                  "K_matrix",
                                  "AlphaXX", "AlphaYY", "AlphaZZ",
                                  "AlphaZY", "AlphaZX", "AlphaXY",
                                  "AlphaInitialXX", "AlphaInitialYY",
                                  "AlphaInitialZZ", "AlphaInitialZY",
                                  "AlphaInitialZX", "AlphaInitialXY",
                                  "ZXX", "ZYY", "ZZZ", "ZZY", "ZZX",
                                  "ZXY", "eps_p_q", "void_ratio"};

  try {
    if (post_process_.at("vtk").is_array() &&
        post_process_.at("vtk").size() > 0) {
      for (unsigned i = 0; i < post_process_.at("vtk").size(); ++i) {
        std::string attribute =
            post_process_["vtk"][i].template get<std::string>();
        if (std::find(vtk.begin(), vtk.end(), attribute) != vtk.end())
          vtk_attributes_.emplace_back(attribute);
        else
          throw std::runtime_error("Specificed VTK argument is incorrect :" +
                                   attribute);
      }
    }
  } catch (std::exception& exception) {
    status = false;
    console_->error("{} {}: {}", __FILE__, __LINE__, exception.what());
    abort();
  }
  return status;
}

// Initialise vtk output for two phase
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_vtk_twophase() {

  bool status = true;

  // Default VTK attributes for water phase
  std::vector<std::string> liquid_vtk = {"pore_pressures",
                                        "PIC_pore_pressure_excess",
                                        "PIC_ru",
                                        "initial_vertical_effective_stresses",
                                        "dynamic_vertical_effective_stresses",
                                        "liquefaction_potentials",
                                        "momentary_liquefied",
                                        "gamma_sub",
                                        "PIC_pore_pressures",
                                        "permeabilities",
                                        "suction_pressures",
                                        "liquid_volumetric_strains",
                                        "liquid_dvolumetric_strains",
                                        "liquid_velocities", 
                                        "liquid_accelerations",
                                        "liquid_strains",
                                        "liquid_fluxes",
                                         "liquid_pressures",
                                         "PIC_liquid_pressures",      
                                         "liquid_pressure_accelerations",
                                         "liquid_saturations",    
                                        "liquid_fractions",      
                                        "liquid_chis",           
                                        "liquid_densities",      
                                        "liquid_sources",        
                                        "liquid_permeabilities",
                                        "liquid_volumes",        
                                        "liquid_masses", 
                                        "liquid_mass_densities",        
                                        "liquid_vol_strains",    
                                        "liquid_viscosities",
                                        "liquid_pressure_gradients",
                                        "liquid_critical_times",
                                        "ice_fractions",
                                        "ice_masses",
                                        "ice_mass_densities",
                                        "ice_densities",
                                        "ice_saturations",
                                        "ice_saturation_rates",
                                        "hydrate_saturations",
                                        "PIC_hydrate_saturations",
                                        "hydrate_fractions",     
                                        "hydrate_densities",     
                                        "hydrate_sources",       
                                        "hydrate_volumes",       
                                        "hydrate_masses",        
                                        "hydrate_mass_densities",
                                        "gas_velocities",
                                        "gas_accelerations",
                                        "gas_strains",
                                        "gas_fluxes",
                                        "gas_pressure_gradients",
                                        "gas_pressures",
                                        "PIC_gas_pressures",         
                                        "gas_saturations",       
                                        "gas_densities", 
                                        "gas_fractions",        
                                        "gas_sources",           
                                        "gas_permeabilities",     
                                        "gas_volumes",       
                                        "gas_masses",            
                                        "gas_mass_densities",    
                                        "gas_vol_strains",
                                        "gas_critical_times",
                                        "gas_viscosities",
                                        "viscosities",
                                        "permeabilities"};

  std::vector<std::string> liquid_vtk_allowed = liquid_vtk;
  // Force-only pressures are opt-in three-phase diagnostics.  Keep them out
  // of the default list because other particle types do not register them.
  liquid_vtk_allowed.emplace_back("force_liquid_pressures");
  liquid_vtk_allowed.emplace_back("force_gas_pressures");
  liquid_vtk_allowed.emplace_back("liquid_seepage_velocities");
  liquid_vtk_allowed.emplace_back("seepage_velocities");
  liquid_vtk_allowed.emplace_back("liquid_seepage_forces");
  liquid_vtk_allowed.emplace_back("seepage_forces");

  try {
    if (post_process_.at("liquid_vtk").is_array() &&
        post_process_.at("liquid_vtk").size() > 0) {
      for (unsigned i = 0; i < post_process_.at("liquid_vtk").size(); ++i) {
        std::string liquid_attribute =
            post_process_["liquid_vtk"][i].template get<std::string>();
        if (std::find(liquid_vtk_allowed.begin(), liquid_vtk_allowed.end(),
                      liquid_attribute) != liquid_vtk_allowed.end())
          vtk_attributes_.emplace_back(liquid_attribute);
        else
          throw std::runtime_error(
              "Specificed VTK argument for water phase is incorrect");
      }
    } else {
      throw std::runtime_error(
          "Specificed VTK arguments for water phase are incorrect, using "
          "defaults");
    }
  } catch (std::exception& exception) {
    status = false;
    for (auto& lvtk : liquid_vtk) vtk_attributes_.emplace_back(lvtk);
    console_->error("{} {}: {}", __FILE__, __LINE__, exception.what());
  }

  return status;
}

//! Checkpoint resume
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::checkpoint_resume() {
  bool checkpoint = true;
  try {
    // TODO: Set phase
    const unsigned phase = 0;

    if (!analysis_["resume"]["resume"].template get<bool>())
      throw std::runtime_error("Resume analysis option is disabled!");

    // Get unique analysis id
    std::string resume_uuid = analysis_["resume"]["uuid"].template get<std::string>();
    // Get step
    int resume_step = analysis_["resume"]["step"].template get<mpm::Index>();
    // Get nsteps
    int resume_nsteps = analysis_["resume"]["nsteps"].template get<mpm::Index>();
    // Start from this step
    bool start_from_this_step = analysis_["resume"]["start_from_this_step"].template get<bool>();

    if (start_from_this_step) this->step_ = analysis_["resume"]["this_step"].template get<mpm::Index>();
  
    // Input particle h5 file for resume
    std::string attribute = "particles";
    std::string extension = ".h5";

    auto particles_file =
        io_->output_file(attribute, extension, resume_uuid, resume_step, resume_nsteps)
            .string();
    // Load particle information from file
    mesh_->read_particles_hdf5(phase, particles_file);

    // Clear all particle ids
    mesh_->iterate_over_cells(
        std::bind(&mpm::Cell<Tdim>::clear_particle_ids, std::placeholders::_1));

    // Locate particles
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

    // Increament step
    ++this->step_;

    console_->info("Checkpoint resume at step {} of {}", this->step_,
                   this->nsteps_);

  } catch (std::exception& exception) {
    console_->info("{} {} Resume failed, restarting analysis: {}", __FILE__,
                   __LINE__, exception.what());
    this->step_ = 0;
    checkpoint = false;
  }
  return checkpoint;
}

//! Write HDF5 files
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::write_hdf5(mpm::Index step, mpm::Index max_steps) {
  // Write input geometry to vtk file
  std::string attribute = "particles";
  std::string extension = ".h5";

  auto particles_file =
      io_->output_file(attribute, extension, uuid_, step, max_steps).string();

  const unsigned phase = 0;
  mesh_->write_particles_hdf5(phase, particles_file);
}

//! Write VTK files
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::write_vtk(mpm::Index step, mpm::Index max_steps) {

  // VTK PolyData writer
  auto vtk_writer = std::make_unique<VtkWriter>(mesh_->particle_coordinates());
  const std::string extension = ".vtp";

  // Write mesh on step 0
  if (step == 0 | step == 100) {
    vtk_writer->write_mesh(
        io_->output_file("mesh", extension, uuid_, step, max_steps).string(),
        mesh_->nodal_coordinates(), mesh_->node_pairs());

    // Write input geometry to vtk file
    auto meshfile =
        io_->output_file("geometry", extension, uuid_, step, max_steps)
            .string();
    vtk_writer->write_geometry(meshfile);
  }

  // TODO: Generalize this with the one phase code
  // Solid phase output
  const unsigned soil_skeleton = mpm::ParticlePhase::Solid;

  //! VTK scalar variable, 1 component
  std::vector<std::string> vtk_scalar_data = {"ids",
                                              "porosities",
                                              "PIC_porosities",
                                              "volumetric_strains",
                                              "dvolumetric_strains",
                                              "liquid_volumetric_strains",
                                              "liquid_dvolumetric_strains",
                                              "thermal_volumetric_strains",
                                              "dthermal_volumetric_strains",
                                              "PIC_volumetric_strains",
                                              "PIC_deviatoric_strains",
                                              "mean_stresses",
                                              "deviatoric_stresses",
                                              "densities",
                                              "volumes",
                                              "PIC_volumes",
                                              "masses",
                                              "pore_pressures",
                                              "PIC_pore_pressure_excess",
                                              "PIC_ru",
                                              "initial_vertical_effective_stresses",
                                              "dynamic_vertical_effective_stresses",
                                              "liquefaction_potentials",
                                              "momentary_liquefied",
                                              "gamma_sub",
                                              "PIC_pore_pressures",
                                              "temperatures",
                                              "PIC_temperatures",
                                              "temperature_accelerations",
                                              "solid_fractions",
                                              "current_time",
                                              "deviatoric_strains",
                                              "suction_pressures",
                                               "liquid_pressures",
                                               "PIC_liquid_pressures",
                                               "force_liquid_pressures",
                                               "liquid_pressure_accelerations",
                                               "liquid_saturations",
                                              "liquid_fractions",
                                              "liquid_chis",
                                              "liquid_densities",
                                              "liquid_sources",
                                              "liquid_permeabilities",
                                              "liquid_volumes",
                                              "liquid_masses",
                                              "liquid_mass_densities",
                                              "liquid_vol_strains",
                                              "liquid_viscosities",
                                              "liquid_critical_times",
                                              "ice_fractions",
                                              "ice_densities",
                                              "ice_masses",
                                              "ice_mass_densities",
                                              "ice_saturations",
                                              "ice_saturation_rates",
                                              "hydrate_saturations",
                                              "PIC_hydrate_saturations",
                                              "hydrate_fractions",
                                              "hydrate_densities",
                                              "hydrate_sources",
                                              "hydrate_volumes",
                                              "hydrate_masses",
                                              "hydrate_mass_densities",
                                              "gas_pressures",
                                              "PIC_gas_pressures",
                                              "force_gas_pressures",
                                              "gas_saturations",
                                              "gas_densities",
                                              "gas_fractions",
                                              "gas_sources",
                                              "gas_permeabilities",
                                              "gas_volumes",
                                              "gas_masses",
                                              "gas_mass_densities",
                                              "gas_vol_strains",
                                              "gas_viscosities",
                                              "viscosities",
                                              "gas_critical_times",
                                              "free_surfaces",
                                              "permeabilities"};

  //! VTK vector variable, 3 component
  std::vector<std::string> vtk_vector_data = {"displacements",
                                              "temperature_gradients",
                                              "mass_gradients",
                                              "outward_normals", 
                                              "velocities",
                                              "accelerations",
                                              "liquid_velocities",
                                              "liquid_accelerations",
                                              "liquid_fluxes",
                                              "liquid_pressure_gradients",
                                              "liquid_seepage_velocities",
                                              "seepage_velocities",
                                              "liquid_seepage_forces",
                                              "seepage_forces",
                                              "rotations",
                                              "gas_velocities",
                                              "gas_accelerations",
                                              "gas_fluxes",
                                              "gas_pressure_gradients",
                                              };

  //! VTK tensor variable, 6 component
  std::vector<std::string> vtk_tensor_vector_data = {
                                              "stresses",
                                              "strains",
                                              "liquid_strains",
                                              "gas_strains",
                                              "thermal_strains",
                                              "K_matrix"};

  //! VTK tensor variable, 9 component
  std::vector<std::string> vtk_tensor_data = {"fabric_CNs",
                                              "fabric_POs",
                                              "velocity_gradients",
                                              "deformation_gradients",
                                              "displacement_gradients",
                                              "grad_shapefns"};

  //! Material state variable, 1 component
  const std::vector<std::string> vtk_state_scalar_data = {
      "AlphaXX",        "AlphaYY",        "AlphaZZ",
      "AlphaZY",        "AlphaZX",        "AlphaXY",
      "AlphaInitialXX", "AlphaInitialYY", "AlphaInitialZZ",
      "AlphaInitialZY", "AlphaInitialZX", "AlphaInitialXY",
      "ZXX",            "ZYY",            "ZZZ",
      "ZZY",            "ZZX",            "ZXY",
      "eps_p_q",        "void_ratio"};

  vtk_writer->create_new_dataset();

  for (auto& attribute : vtk_attributes_) {
    // Write scalar
    if (std ::find(vtk_scalar_data.begin(), vtk_scalar_data.end(), 
                                        attribute) != vtk_scalar_data.end()) {
      vtk_writer->write_scalar_point_data(
          mesh_->particles_scalar_data(attribute), attribute);
    }
    // Write constitutive state variable
    else if (std::find(vtk_state_scalar_data.begin(),
                       vtk_state_scalar_data.end(), attribute) !=
             vtk_state_scalar_data.end()) {
      vtk_writer->write_scalar_point_data(
          mesh_->particles_statevars_data(attribute), attribute);
    }
    // Write vector    
    else if (std ::find(vtk_vector_data.begin(), vtk_vector_data.end(),
                                        attribute) != vtk_vector_data.end()) {
      vtk_writer->write_vector_point_data(
          mesh_->template particles_vector_data<3>(attribute), attribute);
    }
    // Write tensor vector
    else if (std ::find(vtk_tensor_vector_data.begin(), vtk_tensor_vector_data.end(),
                                        attribute) != vtk_tensor_vector_data.end()) {
      vtk_writer->write_tensor_vector_point_data(
          mesh_->template particles_vector_data<6>(attribute), attribute);
    }
    // Write tensor
    else if (std ::find(vtk_tensor_data.begin(), vtk_tensor_data.end(),
                                        attribute) != vtk_tensor_data.end()) {
      vtk_writer->write_tensor_point_data(
          mesh_->template particles_vector_data<9>(attribute), attribute);
    }

  }
  auto file = io_->output_file("particle", extension, uuid_, step, max_steps)
                      .string();
  vtk_writer->write_file(file);
  

  // VTK PolyData writer for nodal properties
  auto nodal_vtk_writer = std::make_unique<VtkWriter>(mesh_->nodal_coordinates());

  //! VTK scalar variable, 1 component
  std::vector<std::string> nodal_vtk_scalar_data = {"temperature",
                                                    "temperature_acc",
                                                    "mixture_mass",
                                                    "liquid_mass", 
                                                    "gas_mass",
                                                    "free_surface",
                                                    "liquid_pressure",
                                                    "liquid_pressure_acc",
                                                    "gas_pressure",
                                                    "gas_pressure_acc",
                                                    "gas_mass_source",
                                                    "liquid_mass_source",
                                                    "gas_hydraulic_conduction",
                                                    "liquid_hydraulic_conduction"};

  //! VTK vector variable, 3 component
  std::vector<std::string> nodal_vtk_vector_data = {"solid_velocity",
                                                    "liquid_velocity",
                                                    "gas_velocity",
                                                    "solid_acceleration",
                                                    "liquid_acceleration",
                                                    "gas_acceleration",
                                                    "mix_ext_force",   
                                                    "liquid_ext_force",
                                                    "gas_ext_force",
                                                    "mix_int_force",   
                                                    "liquid_int_force",
                                                    "gas_int_force",                                                    
                                                    "drag_force_liquid",
                                                    "drag_force_gas"
                                                    }; 

  std::vector<std::string> nodal_vtk_attributes = {"temperature",
                                                    "temperature_acc",
                                                    "mixture_mass",
                                                    "solid_velocity",
                                                    "solid_acceleration",
                                                    "mix_ext_force",   
                                                    "mix_int_force", 
                                                    "free_surface",
                                                    "reaction_force"
                                                    // "liquid_mass",
                                                    // "liquid_velocity",
                                                    // "liquid_acceleration",
                                                    // "liquid_ext_force",
                                                    // "liquid_int_force",
                                                    // "liquid_pressure",
                                                    // "liquid_pressure_acc",
                                                    // "liquid_mass_source",
                                                    // "liquid_hydraulic_conduction",
                                                    // "gas_mass",
                                                    // "gas_velocity",                                                    
                                                    // "gas_acceleration",
                                                    // "gas_ext_force",
                                                    // "gas_int_force",
                                                    // "drag_force_liquid",
                                                    // "drag_force_gas",
                                                    // "gas_pressure",
                                                    // "gas_pressure_acc",
                                                    // "gas_mass_source",
                                                    // "gas_hydraulic_conduction"
                                                    };
  nodal_vtk_writer->create_new_dataset();

for (auto& attribute : nodal_vtk_attributes) {
    // Write scalar
    if (std ::find(nodal_vtk_scalar_data.begin(), nodal_vtk_scalar_data.end(), 
                                        attribute) != nodal_vtk_scalar_data.end()) {
      nodal_vtk_writer->write_scalar_point_data(
          mesh_->nodal_scalar_data(attribute), attribute);
    }
    // Write vector    
    else if (std ::find(nodal_vtk_vector_data.begin(), nodal_vtk_vector_data.end(),
                                        attribute) != nodal_vtk_vector_data.end()) {
      nodal_vtk_writer->write_vector_point_data(
          mesh_->template nodal_vector_data<3>(attribute), attribute);
    }
  }
  auto file_node = io_->output_file("node", extension, uuid_, step, max_steps)
                      .string();
  nodal_vtk_writer->write_file(file_node);

}

#ifdef USE_PARTIO
//! Write Partio files
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::write_partio(mpm::Index step, mpm::Index max_steps) {

  // MPI parallel partio file
  int mpi_rank = 0;
  int mpi_size = 1;
#ifdef USE_MPI
  // Get MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  // Get Partio file extensions
  const std::string extension = ".bgeo";
  const std::string attribute = "partio";
  // Create filename
  auto file =
      io_->output_file(attribute, extension, uuid_, step, max_steps).string();
  // Write partio file
  mpm::partio::write_particles(file, mesh_->particles_hdf5());
}
#endif  // USE_PARTIO

//! Return if a mesh is isoparametric
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::is_isoparametric() {
  bool isoparametric = true;

  try {
    const auto mesh_props = io_->json_object("mesh");
    isoparametric = mesh_props.at("isoparametric").template get<bool>();
  } catch (std::exception& exception) {
    console_->warn(
        "{} {} Isoparametric status of mesh: {}\n Setting mesh as "
        "isoparametric.",
        __FILE__, __LINE__, exception.what());
    isoparametric = true;
  }
  return isoparametric;
}

//! Initialise loads
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_loads() {
  bool status = true;
  try {
    auto loads = io_->json_object("external_loading_conditions");
    // Initialise gravity loading
    if (loads.at("gravity").is_array() &&
        loads.at("gravity").size() == gravity_.size()) {
      for (unsigned i = 0; i < gravity_.size(); ++i) {
        gravity_[i] = loads.at("gravity").at(i);
      }
    } else {
      throw std::runtime_error("Specified gravity dimension is invalid");
    }

    // Create a file reader
    const std::string io_type =
        io_->json_object("mesh")["io_type"].template get<std::string>();
    auto reader = Factory<mpm::IOMesh<Tdim>>::instance()->create(io_type);

    // Read and assign particles surface tractions
    if (loads.find("particle_surface_traction") != loads.end()) {
      for (const auto& ptraction : loads["particle_surface_traction"]) {
        // Get the math function
        std::shared_ptr<FunctionBase> Tfunction = nullptr;
        // If a math function is defined set to function or use scalar
        if (ptraction.find("math_function_id") != ptraction.end())
          Tfunction = math_functions_.at(
              ptraction.at("math_function_id").template get<unsigned>());
        // Set id
        int pset_id = ptraction.at("pset_id").template get<int>();
        // Facet
        unsigned facet = ptraction.at("facet").template get<unsigned>();
        // Direction
        unsigned dir = ptraction.at("dir").template get<unsigned>();
        // Traction
        double traction = ptraction.at("traction").template get<double>();

        // Create particle surface tractions
        bool particles_tractions = mesh_->create_particles_tractions(
            Tfunction, pset_id, facet, dir, traction);
        if (!particles_tractions)
          throw std::runtime_error(
              "Particles tractions are not properly assigned");
      }
    } else
      console_->warn(
          "No particle surface traction is defined for the analysis");

    // Read and assign particles heat sources
    if (loads.find("particle_heat_source") != loads.end()) {
      for (const auto& pheat_source : loads["particle_heat_source"]) {

        // Get the math function
        std::shared_ptr<FunctionBase> Hfunction = nullptr;
        // If a math function is defined set to function or use scalar
        if (pheat_source.find("math_function_id") != pheat_source.end())
          Hfunction = math_functions_.at(
              pheat_source.at("math_function_id").template get<unsigned>());
        
        // Set id
        int pset_id = pheat_source.at("pset_id").template get<int>();
        // Heat source
        double heat_source = pheat_source.at("heat_source").template get<double>();

        // Create particle heat source
        bool particles_heat_sources = mesh_->create_particles_heat_sources(
            Hfunction, pset_id, heat_source);
        if (!particles_heat_sources)
          throw std::runtime_error(
              "Particles heat sources are not properly assigned");
      }
    } else
      console_->warn(
          "No particle heat source is defined for the analysis");

    // Read and assign particles surface tractions
    if (loads.find("particle_surface_contact") != loads.end()) {
      for (const auto& pcontact : loads["particle_surface_contact"]) {
        // Set id
        int pset_id = pcontact.at("pset_id").template get<int>();
        // Components of normal vector
        unsigned dir = pcontact.at("dir").template get<unsigned>();
        double normal = pcontact.at("normal").template get<double>();
        // // Normal vector
        // Eigen::Matrix<double, Tdim, 1> normal.setZero();
        // if (Tdim == 2) {
        //   normal(0) = normal_x;
        //   normal(1) = normal_y;
        // } 
        // if (Tdim == 3) {
        //   normal(0) = normal_x;
        //   normal(1) = normal_y;
        //   normal(2) = normal_z;          
        // }
        // Create particle surface tractions
        bool particles_contacts = mesh_->create_particles_contacts(
            pset_id, dir, normal);

        if (!particles_contacts)
          throw std::runtime_error(
              "Particles contacts are not properly assigned");
      } 
    } else
      console_->warn("No particle surface contacts is defined for the analysis"); 

    // // Read and assign nodal concentrated forces
    // if (loads.find("concentrated_nodal_forces") != loads.end()) {
    //   for (const auto& nforce : loads["concentrated_nodal_forces"]) {
    //     // Forces are specified in a file
    //     if (nforce.find("file") != nforce.end()) {
    //       std::string force_file =
    //           nforce.at("file").template get<std::string>();
    //       bool nodal_forces = mesh_->assign_nodal_concentrated_forces(
    //           reader->read_forces(io_->file_name(force_file)));
    //       if (!nodal_forces)
    //         throw std::runtime_error(
    //             "Nodal force file is invalid, forces are not properly "
    //             "assigned");
    //       set_node_concentrated_force_ = true;
    //     } else {
    //       // Get the math function
    //       std::shared_ptr<FunctionBase> ffunction = nullptr;
    //       if (nforce.find("math_function_id") != nforce.end())
    //         ffunction = math_functions_.at(
    //             nforce.at("math_function_id").template get<unsigned>());
    //       // Set id
    //       int nset_id = nforce.at("nset_id").template get<int>();
    //       // Direction
    //       unsigned dir = nforce.at("dir").template get<unsigned>();
    //       // Force
    //       double force = nforce.at("force").template get<double>();

    //       // Read and assign nodal concentrated forces
    //       bool nodal_force = mesh_->assign_nodal_concentrated_forces(
    //           ffunction, nset_id, dir, force);
    //       if (!nodal_force)
    //         throw std::runtime_error(
    //             "Concentrated nodal forces are not properly assigned");
    //       set_node_concentrated_force_ = true;
    //     }
    //   }  
    // } else
    //   console_->warn("No concentrated nodal force is defined for the analysis");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Initialise math functions
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_math_functions(const Json& math_functions) {
  bool status = true;
  try {
    // Get materials properties
    for (const auto& function_props : math_functions) {

      // Get math function id
      auto function_id = function_props["id"].template get<unsigned>();

      // Get function type
      const std::string function_type =
          function_props["type"].template get<std::string>();

      // Create a new function from JSON object
      auto function =
          Factory<mpm::FunctionBase, unsigned, const Json&>::instance()->create(
              function_type, std::move(function_id), function_props);

      // Add material to list
      auto insert_status =
          math_functions_.insert(std::make_pair(function->id(), function));

      // If insert material failed
      if (!insert_status.second) {
        status = false;
        throw std::runtime_error(
            "Invalid properties for new math function, fn insertion failed");
      }
    }
  } catch (std::exception& exception) {
    console_->error("#{}: Reading math functions: {}", __LINE__,
                    exception.what());
    status = false;
  }
  return status;
}

//! Node entity sets
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::node_entity_sets(const Json& mesh_props,
                                          bool check_duplicates) {
  try {
    if (mesh_props.find("entity_sets") != mesh_props.end()) {
      std::string entity_sets =
          mesh_props["entity_sets"].template get<std::string>();
      if (!io_->file_name(entity_sets).empty()) {
        bool node_sets = mesh_->create_node_sets(
            (io_->entity_sets(io_->file_name(entity_sets), "node_sets")),
            check_duplicates);
        if (!node_sets)
          throw std::runtime_error("Node sets are not properly assigned");
      }
    } else
      throw std::runtime_error("Entity set JSON not found");
  } catch (std::exception& exception) {
    console_->warn("#{}: Entity sets are undefined {} ", __LINE__,
                   exception.what());
  }
}

//! Node Euler angles
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::node_euler_angles(
    const Json& mesh_props, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io) {
  try {
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find("nodal_euler_angles") !=
            mesh_props["boundary_conditions"].end()) {
      std::string euler_angles =
          mesh_props["boundary_conditions"]["nodal_euler_angles"]
              .template get<std::string>();
      if (!io_->file_name(euler_angles).empty()) {
        bool rotation_matrices = mesh_->compute_nodal_rotation_matrices(
            mesh_io->read_euler_angles(io_->file_name(euler_angles)));
        if (!rotation_matrices)
          throw std::runtime_error(
              "Euler angles are not properly assigned/computed");
      }
    } else
      throw std::runtime_error("Euler angles JSON not found");
  } catch (std::exception& exception) {
    console_->warn("#{}: Euler angles are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Nodal velocity constraints
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::nodal_velocity_constraints(
    const Json& mesh_props, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io) {
  try {
    // Read and assign velocity constraints
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find("velocity_constraints") !=
            mesh_props["boundary_conditions"].end()) {
      // Iterate over velocity constraints
      for (const auto& constraints :
           mesh_props["boundary_conditions"]["velocity_constraints"]) {
        // Velocity constraints are specified in a file
        if (constraints.find("file") != constraints.end()) {
          std::string velocity_constraints_file =
              constraints.at("file").template get<std::string>();
          bool velocity_constraints = mesh_->assign_nodal_velocity_constraints(
              mesh_io->read_velocity_constraints(
                  io_->file_name(velocity_constraints_file)));
          if (!velocity_constraints)
            throw std::runtime_error(
                "Velocity constraints are not properly assigned");

        } else {
          // Set id
          int nset_id = constraints.at("nset_id").template get<int>();
          // Direction
          unsigned dir = constraints.at("dir").template get<unsigned>();
          // Velocity
          double velocity = constraints.at("velocity").template get<double>();

          // Get the math function
          std::shared_ptr<FunctionBase> vfunction = nullptr;
          if (constraints.find("math_function_id") != constraints.end()) {
            vfunction = math_functions_.at(
                constraints.at("math_function_id").template get<unsigned>());
          }

          // Add velocity constraint to mesh
          auto velocity_constraint = std::make_shared<mpm::VelocityConstraint>(
              nset_id, dir, velocity, vfunction);

          // mesh_->creat_nodal_velocity_constraint()
          mesh_->create_nodal_velocity_constraint(nset_id, velocity_constraint);
        }
      }
    } else
      throw std::runtime_error("Velocity constraints JSON not found");
  } catch (std::exception& exception) {
    console_->error("#{}: Velocity constraints are undefined {} ", __LINE__,
                    exception.what());
  }
}

// Nodal frictional constraints
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::nodal_frictional_constraints(
    const Json& mesh_props, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io) {
  try {
    // Read and assign friction constraints
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find("friction_constraints") !=
            mesh_props["boundary_conditions"].end()) {
      // Iterate over velocity constraints
      for (const auto& constraints :
           mesh_props["boundary_conditions"]["friction_constraints"]) {
        // Friction constraints are specified in a file
        if (constraints.find("file") != constraints.end()) {
          std::string friction_constraints_file =
              constraints.at("file").template get<std::string>();
          bool friction_constraints = mesh_->assign_nodal_friction_constraints(
              mesh_io->read_friction_constraints(
                  io_->file_name(friction_constraints_file)));
          if (!friction_constraints)
            throw std::runtime_error(
                "Friction constraints are not properly assigned");

        } else {

          // Set id
          int nset_id = constraints.at("nset_id").template get<int>();
          // Direction
          unsigned dir = constraints.at("dir").template get<unsigned>();
          // Sign n
          int sign_n = constraints.at("sign_n").template get<int>();
          // Friction
          double friction = constraints.at("friction").template get<double>();
          // Add friction constraint to mesh
          auto friction_constraint = std::make_shared<mpm::FrictionConstraint>(
              nset_id, dir, sign_n, friction);
          mesh_->assign_nodal_frictional_constraint(nset_id,
                                                    friction_constraint);
        }
      }
    } else
      throw std::runtime_error("Friction constraints JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: Friction conditions are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Nodal pore pressure constraints (Coupled solid-fluid formulation)
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::nodal_pore_pressure_constraints(
    const Json& mesh_props, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io) {
  try {
    // Water phase indice
    unsigned phase = mpm::ParticlePhase::Liquid;
    // Total phase
    const unsigned Tnphases = 2;

    // Read and assign pore pressure constraints
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find("pore_pressure_constraints") !=
            mesh_props["boundary_conditions"].end()) {

      // Iterate over pore pressure constraints
      for (const auto& constraints :
           mesh_props["boundary_conditions"]["pore_pressure_constraints"]) {

        // Check if it is pressure increment constraints
        if (constraints.find("increment_boundary") != constraints.end() &&
            constraints["increment_boundary"]) {
          // Storage location of the increment boundary
          phase += Tnphases;
        }

        // Pore pressure constraints are specified in a file
        if (constraints.find("file") != constraints.end()) {
          std::string pore_pressure_constraints_file =
              constraints.at("file").template get<std::string>();
          bool ppressure_constraints = mesh_->assign_nodal_pressure_constraints(
              phase, mesh_io->read_pressure_constraints(
                         io_->file_name(pore_pressure_constraints_file)));
          if (!ppressure_constraints)
            throw std::runtime_error(
                "Pore pressure constraints are not properly assigned");
        } else {
          // Get the math function
          std::shared_ptr<FunctionBase> pfunction = nullptr;
          if (constraints.find("math_function_id") != constraints.end())
            pfunction = math_functions_.at(
                constraints.at("math_function_id").template get<unsigned>());
          // Set id
          int nset_id = constraints.at("nset_id").template get<int>();
          // Pore Pressure
          double pore_pressure =
              constraints.at("pore_pressure").template get<double>();
          // Add pore pressure constraint to mesh
          mesh_->assign_nodal_pressure_constraint(pfunction, nset_id, phase,
                                                  pore_pressure);
          // if (phase >= Tnphases) {
          //   // Reference step
          //   const Index ref_step =
          //       constraints.at("ref_step").template get<Index>();
          //   // Add reference step to nodes
          //   mesh_->assign_nodal_pressure_reference_step(nset_id, ref_step);
          // }
        }
      }
    } else
      throw std::runtime_error("Pore pressure constraints JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: Pore pressure conditions are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Nodal temperature constraints 
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::nodal_temperature_constraints(
    const Json& mesh_props, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io) {
  try {
    // phase indice
    unsigned phase = mpm::ParticlePhase::Solid;
    // Total phase
    const unsigned Tnphases = 1;
    
    // Read and assign temperature constraints
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find("temperature_constraints") !=
            mesh_props["boundary_conditions"].end()) {

      // Iterate over temperature constraints
      for (const auto& constraints :
           mesh_props["boundary_conditions"]["temperature_constraints"]) {

        // temperature constraints are specified in a file
        if (constraints.find("file") != constraints.end()) {
          std::string temperature_constraints_file =
              constraints.at("file").template get<std::string>();
          bool ptemperature_constraints = mesh_->assign_nodal_temperature_constraints(
              phase, mesh_io->read_temperature_constraints(
                         io_->file_name(temperature_constraints_file)));              
          if (!ptemperature_constraints)
            throw std::runtime_error(
                "Temperature constraints are not properly assigned");
        } else {
          // Get the math function
          std::shared_ptr<FunctionBase> Tfunction = nullptr;
          if (constraints.find("math_function_id") != constraints.end())
            Tfunction = math_functions_.at(
                constraints.at("math_function_id").template get<unsigned>());
          // Set id
          int nset_id = constraints.at("nset_id").template get<int>();
          // Temperature
          double temperature =
              constraints.at("temperature").template get<double>(); 
          // if (nset_id == -2) {
          //   // Add temperature constraint to mesh
          //   mesh_->assign_moving_free_surface_temperature(Tfunction, nset_id, phase,
          //                                             temperature);
          // } else {
            // Add temperature constraint to mesh
            mesh_->assign_nodal_temperature_constraint(Tfunction, nset_id, phase,
                                                      temperature);            
          // }  
        }
      }
    } else
      throw std::runtime_error("temperature constraints JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: temperature conditions are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Nodal temperature constraints 
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::nodal_convective_heat_constraints(
    const Json& mesh_props, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io) {
  try {
    // phase indice
    unsigned phase = mpm::ParticlePhase::Solid;
    // Total phase
    const unsigned Tnphases = 1;
    
    // Read and assign temperature constraints
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find("convective_heat_constraints") !=
            mesh_props["boundary_conditions"].end()) {

      // Iterate over temperature constraints
      for (const auto& constraints :
           mesh_props["boundary_conditions"]["convective_heat_constraints"]) {

          // Get the math function
          std::shared_ptr<FunctionBase> Tfunction = nullptr;
          if (constraints.find("math_function_id") != constraints.end())
            Tfunction = math_functions_.at(
                constraints.at("math_function_id").template get<unsigned>());
          // Set id
          int nset_id = constraints.at("nset_id").template get<int>();
          // Temperature
          double temperature =
              constraints.at("temperature").template get<double>();
          double convective_coeff =
              constraints.at("convective_coeff").template get<double>();

          // Add temperature constraint to mesh
          mesh_->assign_nodal_convective_heat_constraint(Tfunction, nset_id, phase,
                                                    temperature, convective_coeff);
      }
    } else
      throw std::runtime_error("convective heat constraints JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: convective heat conditions are undefined {} ", __LINE__,
                   exception.what());
  }
}


// Nodal heat source
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::nodal_heat_source(
    const Json& mesh_props, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io) {
  try {
    // phase indice
    unsigned phase = mpm::ParticlePhase::Solid;
    // Total phase
    const unsigned Tnphases = 1;

    // Read and assign heat source
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find("heat_source") !=
            mesh_props["boundary_conditions"].end()) {

      // Iterate over heat source
      for (const auto& constraints :
           mesh_props["boundary_conditions"]["heat_source"]) {

        // // temperature constraints are specified in a file
        // if (constraints.find("file") != constraints.end()) {
        //   std::string heat_source_file =
        //       constraints.at("file").template get<std::string>();
        //   bool nheat_source = mesh_->assign_nodal_heat_sources(
        //       phase, mesh_io->read_heat_source(
        //                  io_->file_name(heat_source_file)));
        //   if (!nheat_source)
        //     throw std::runtime_error(
        //         "Heat source is not properly assigned");
        // } else {
          // Get the math function
          std::shared_ptr<FunctionBase> Hfunction = nullptr;
          if (constraints.find("math_function_id") != constraints.end())
            Hfunction = math_functions_.at(
                constraints.at("math_function_id").template get<unsigned>());
          // Set id
          int nset_id = constraints.at("nset_id").template get<int>();
          // Heat source
          double heat_source =
              constraints.at("heat_source").template get<double>();
          // Add heat source constraint to mesh
          mesh_->assign_nodal_heat_source(Hfunction, nset_id, phase,
                                                     heat_source);
        // }
      }
    } else
      throw std::runtime_error("heat source JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: no heat source are defined {} ", __LINE__,
                   exception.what());
  }
}

//! Cell entity sets
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::cell_entity_sets(const Json& mesh_props,
                                          bool check_duplicates) {
  try {
    if (mesh_props.find("entity_sets") != mesh_props.end()) {
      // Read and assign cell sets
      std::string entity_sets =
          mesh_props["entity_sets"].template get<std::string>();
      if (!io_->file_name(entity_sets).empty()) {
        bool cell_sets = mesh_->create_cell_sets(
            (io_->entity_sets(io_->file_name(entity_sets), "cell_sets")),
            check_duplicates);
        if (!cell_sets)
          throw std::runtime_error("Cell sets are not properly assigned");
      }
    } else
      throw std::runtime_error("Cell entity sets JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: Cell entity sets are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Particles cells
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particles_cells(
    const Json& mesh_props,
    const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io) {
  try {
    if (mesh_props.find("particle_cells") != mesh_props.end()) {
      std::string fparticles_cells =
          mesh_props["particle_cells"].template get<std::string>();

      if (!io_->file_name(fparticles_cells).empty()) {
        bool particles_cells =
            mesh_->assign_particles_cells(particle_io->read_particles_cells(
                io_->file_name(fparticles_cells)));
        if (!particles_cells)
          throw std::runtime_error(
              "Particle cells are not properly assigned to particles");
      }
    } else
      throw std::runtime_error("Particle cells JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: Particle cells are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Particles volumes
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particles_volumes(
    const Json& mesh_props,
    const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io) {
  try {
    if (mesh_props.find("particles_volumes") != mesh_props.end()) {
      std::string fparticles_volumes =
          mesh_props["particles_volumes"].template get<std::string>();
      if (!io_->file_name(fparticles_volumes).empty()) {
        bool particles_volumes =
            mesh_->assign_particles_volumes(particle_io->read_particles_volumes(
                io_->file_name(fparticles_volumes)));
        if (!particles_volumes)
          throw std::runtime_error(
              "Particles volumes are not properly assigned");
      }
    } else
      throw std::runtime_error("Particle volumes JSON not found");
  } catch (std::exception& exception) {
    console_->warn("#{}: Particle volumes are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Particle velocity constraints
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particle_velocity_constraints(
    const Json& mesh_props,
    const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io) {
  try {
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find(
            "particles_velocity_constraints") !=
            mesh_props["boundary_conditions"].end()) {

      // Iterate over velocity constraints
      for (const auto& constraints :
           mesh_props["boundary_conditions"]
                     ["particles_velocity_constraints"]) {
        // Set id
        int pset_id = constraints.at("pset_id").template get<int>();
        // Direction
        unsigned dir = constraints.at("dir").template get<unsigned>();
        // Velocity
        double velocity = constraints.at("velocity").template get<double>();

        // Get the math function
        std::shared_ptr<FunctionBase> vfunction = nullptr;
        if (constraints.find("math_function_id") != constraints.end())
          vfunction = math_functions_.at(
              constraints.at("math_function_id").template get<unsigned>());

        // Add velocity constraint to mesh
        auto velocity_constraint = std::make_shared<mpm::VelocityConstraint>(
            pset_id, dir, velocity, vfunction);

        mesh_->create_particle_velocity_constraint(pset_id,
                                                   velocity_constraint);
      }
    } else
      throw std::runtime_error("Particle velocity constraints JSON not found");
  } catch (std::exception& exception) {
    console_->warn("#{}: Particle velocity constraints are undefined {} ",
                   __LINE__, exception.what());
  }
}

// Particle pore pressure constraints
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particle_pore_pressure_constraints(
    const Json& mesh_props,
    const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io) {
  try {
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find(
            "particle_pore_pressure_constraints") !=
            mesh_props["boundary_conditions"].end()) {

      // Iterate over pore_pressure constraints
      for (const auto& constraints :
           mesh_props["boundary_conditions"]
                     ["particle_pore_pressure_constraints"]) {
        // Set id
        int pset_id = constraints.at("pset_id").template get<int>();
        // pore_pressure
        double pore_pressure = constraints.at("pore_pressure").template get<double>();

        // Get the math function
        std::shared_ptr<FunctionBase> Pfunction = nullptr;
        if (constraints.find("math_function_id") != constraints.end())
          Pfunction = math_functions_.at(
              constraints.at("math_function_id").template get<unsigned>());

        // Add pore_pressure constraint to mesh
        auto pore_pressure_constraint = std::make_shared<mpm::PorepressureConstraint>(
            pset_id, pore_pressure, Pfunction);

        mesh_->create_particle_pore_pressure_constraint(pset_id,
                                                   pore_pressure_constraint);
      }
    } else
      throw std::runtime_error("Particle pore pressure constraints JSON not found");
  } catch (std::exception& exception) {
    console_->warn("#{}: Particle pore pressure constraints are undefined {} ",
                   __LINE__, exception.what());
  }
}

// Particle temperature constraints
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particle_temperature_constraints(
    const Json& mesh_props,
    const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io) {
  try {
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find(
            "particles_temperature_constraints") !=
            mesh_props["boundary_conditions"].end()) {

      // Iterate over temperature constraints
      for (const auto& constraints :
           mesh_props["boundary_conditions"]
                     ["particles_temperature_constraints"]) {
        // Set id
        int pset_id = constraints.at("pset_id").template get<int>();
        // temperature
        double temperature = constraints.at("temperature").template get<double>();

        // Get the math function
        std::shared_ptr<FunctionBase> Tfunction = nullptr;
        if (constraints.find("math_function_id") != constraints.end())
          Tfunction = math_functions_.at(
              constraints.at("math_function_id").template get<unsigned>());

        // Add temperature constraint to mesh
        auto temperature_constraint = std::make_shared<mpm::TemperatureConstraint>(
            pset_id, temperature, Tfunction);

        mesh_->create_particle_temperature_constraint(pset_id,
                                                   temperature_constraint);
      }
    } else
      throw std::runtime_error("Particle temperature constraints JSON not found");
  } catch (std::exception& exception) {
    console_->warn("#{}: Particle temperature constraints are undefined {} ",
                   __LINE__, exception.what());
  }
}

// Particles stresses
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particles_stresses(
    const Json& mesh_props,
    const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io) {
  try {
    if (mesh_props.find("particles_stresses") != mesh_props.end()) {
      std::string fparticles_stresses =
          mesh_props["particles_stresses"].template get<std::string>();
      if (!io_->file_name(fparticles_stresses).empty()) {

        // Get stresses of all particles
        const auto all_particles_stresses =
            particle_io->read_particles_stresses(
                io_->file_name(fparticles_stresses));

        // Read and assign particles stresses
        if (!mesh_->assign_particles_stresses(all_particles_stresses))
          throw std::runtime_error(
              "Particles stresses are not properly assigned");
      }
    } else
      throw std::runtime_error("Particle stresses JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: Particle stresses are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Particles pore pressures
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particles_pore_pressures(
    const Json& mesh_props,
    const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io) {
  try {
    if (mesh_props.find("particles_pore_pressures") != mesh_props.end()) {
      // Assign initial pore pressure by file
      if (mesh_props["particles_pore_pressures"].find("file") !=
          mesh_props["particles_pore_pressures"].end()) {
        std::string fparticles_pore_pressures =
            mesh_props["particles_pore_pressures"]["file"]
                .template get<std::string>();
        if (!io_->file_name(fparticles_pore_pressures).empty()) {

          // Get pore pressures of all particles
          const auto all_particles_pore_pressures =
              particle_io->read_particles_pressures(
                  io_->file_name(fparticles_pore_pressures));

          // Read and assign particles pore pressures
          if (!mesh_->assign_particles_pore_pressures(
                  all_particles_pore_pressures))
            throw std::runtime_error(
                "Particles pore pressures are not properly assigned");
        } else
          throw std::runtime_error("Particle pore pressures JSON not found");

      } else {
        // Initialise water tables
        std::map<double, double> refernece_points;
        // Vertical direction
        const unsigned dir_v = mesh_props["particles_pore_pressures"]["dir_v"]
                                    .template get<unsigned>();
        // Horizontal direction
        const unsigned dir_h = mesh_props["particles_pore_pressures"]["dir_h"]
                                    .template get<unsigned>();
        // Iterate over water tables
        for (const auto& water_table :
              mesh_props["particles_pore_pressures"]["water_tables"]) {
          // Position coordinate
          double position = water_table.at("position").template get<double>();
          // Direction
          double h0 = water_table.at("h0").template get<double>();
          // Add refernece points to mesh
          refernece_points.insert(std::make_pair<double, double>(
              static_cast<double>(position), static_cast<double>(h0)));
        }
        // Initialise particles pore pressures by watertable
        mesh_->iterate_over_particles(std::bind(
            &mpm::ParticleBase<Tdim>::initialise_pore_pressure_watertable,
            std::placeholders::_1, dir_v, dir_h, refernece_points));
      }
    }
  } catch (std::exception& exception) {
    console_->warn("#{}: Particle pore pressures are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Particles temperatures
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particles_temperatures(
    const Json& mesh_props,
    const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io) {
  try {
    if (mesh_props.find("particles_temperatures") != mesh_props.end()) {
      // Assign initial temperature by file
        std::string fparticles_temperatures =
            mesh_props["particles_temperatures"].template get<std::string>();
      if (!io_->file_name(fparticles_temperatures).empty()) {

        // Get temperatures of all particles
        const auto all_particles_temperatures =
            particle_io->read_particles_temperatures(
                io_->file_name(fparticles_temperatures));

          // Read and assign particles temperatures
          if (!mesh_->assign_particles_temperatures(all_particles_temperatures))
            throw std::runtime_error(
                "Particles temperatures are not properly assigned");
        } 
      }else
       throw std::runtime_error("Particle temperatures JSON not found");
       
  } catch (std::exception& exception) {
    console_->warn("#{}: Particle temperatures are undefined {} ", __LINE__,
                   exception.what());
  }
}

//! Particle entity sets
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particle_entity_sets(const Json& mesh_props,
                                              bool check_duplicates) {
  // Read and assign particle sets
  try {
    if (mesh_props.find("entity_sets") != mesh_props.end()) {
      std::string entity_sets =
          mesh_props["entity_sets"].template get<std::string>();
      if (!io_->file_name(entity_sets).empty()) {
        bool particle_sets = mesh_->create_particle_sets(
            (io_->entity_sets(io_->file_name(entity_sets), "particle_sets")),
            check_duplicates);
      }
    } else
      throw std::runtime_error("Particle entity set JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: Particle sets are undefined {} ", __LINE__,
                   exception.what());
  }
}

//! Assign particle materials by particle sets
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particle_materials(const Json& mesh_props) {
  const Json* particle_materials = nullptr;

  if (mesh_props.find("boundary_conditions") != mesh_props.end()) {
    const auto& boundary_conditions = mesh_props.at("boundary_conditions");
    if (boundary_conditions.find("particle_materials") !=
        boundary_conditions.end())
      particle_materials = &boundary_conditions.at("particle_materials");
  }

  if (particle_materials == nullptr &&
      mesh_props.find("particle_materials") != mesh_props.end())
    particle_materials = &mesh_props.at("particle_materials");

  if (particle_materials == nullptr) return;

  for (const auto& particle_material : *particle_materials) {
    const int pset_id = particle_material.at("pset_id").template get<int>();
    const unsigned material_id =
        particle_material.at("material_id").template get<unsigned>();

    if (!mesh_->assign_material_to_particle_set(pset_id, material_id))
      throw std::runtime_error("Particle material assignment failed");
  }
}

// Initialise Damping
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_damping(const Json& damping_props) {

  // Read damping JSON object
  bool status = true;
  try {
    // Read damping type
    std::string type = damping_props.at("type").template get<std::string>();
    if (type == "Cundall") damping_type_ = mpm::Damping::Cundall;
    if (type == "Viscous") damping_type_ = mpm::Damping::Viscous;
    // Read damping factor
    damping_factor_ = damping_props.at("damping_factor").template get<double>();

  } catch (std::exception& exception) {
    console_->warn("#{}: Damping parameters are undefined {} ", __LINE__,
                   exception.what());
    status = false;
  }

  return status;
}

// Initialise water table
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_nodal_water_table(
    const Json& water_tables) {
  bool status = true;
  try {
    // Iterate over each reference position
    for (const auto& water_table : water_tables) {
      // Set id
      int nset_id = water_table.at("nset_id").template get<int>();
      // Direction
      unsigned dir = water_table.at("dir").template get<unsigned>();
      // h0
      double h0 = water_table.at("h0").template get<double>();
      // Get the math function
      std::shared_ptr<FunctionBase> wfunction = nullptr;
      if (water_table.find("math_function_id") != water_table.end())
        wfunction = math_functions_.at(
            water_table.at("math_function_id").template get<unsigned>());
      // Add water table to mesh
      mesh_->assign_nodal_water_table(wfunction, nset_id, dir, h0);
    }
  } catch (std::exception& exception) {
    console_->warn("#{}: Water table are undefined {} ", __LINE__,
                   exception.what());
    status = false;
  }

  return status;
}

//! Write reaction force
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::write_reaction_force(bool overwrite, mpm::Index step,
                                              mpm::Index max_steps) {

  auto reaction_file =
      io_->output_file("reaction-force", ".txt", uuid_, 0, 0).string();

  //! total displacement
  Eigen::Matrix<double, Tdim, 1> disp = Eigen::Matrix<double, Tdim, 1>::Zero();

  //! total reactiton force
  Eigen::Matrix<double, Tdim, 1> reaction_force =
      Eigen::Matrix<double, Tdim, 1>::Zero();

  //! store the displacement and reaction in disp and reaction_force para
  mesh_->get_reaction_force(disp, reaction_force);

  std::fstream file;
  if (overwrite) {
    file.open(reaction_file.c_str(), std::fstream::out);
  } else {
    file.open(reaction_file.c_str(), std::fstream::app);
    if (Tdim == 1) {
      file << disp[0] << "\t" << reaction_force[0] << "\n";
    } else if (Tdim == 2) {
      file << disp[0] << "\t" << disp[1] << "\t" << reaction_force[0] << "\t"
           << reaction_force[1] << "\n";
    } else if (Tdim == 3) {
      file << disp[0] << "\t" << disp[1] << "\t" << disp[2] << "\t"
           << reaction_force[0] << "\t" << reaction_force[1] << "\t"
           << reaction_force[2] << "\n";
    }
  }
  file.close();
}
