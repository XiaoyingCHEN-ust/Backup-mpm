//! Constructor
template <unsigned Tdim>
mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::ThermoMPMExplicitThreePhaseLag(
    const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("ThermoMPMExplicitThreePhaseLag");
}

#ifdef USE_VTK
template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::prepare_seepage_velocity_vtk() {
  const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
  const unsigned pore_gas = mpm::ParticlePhase::Gas;

  mesh_->iterate_over_nodes(
      std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));
  mesh_->iterate_over_cells(
      std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::map_mass_pressure_to_nodes,
                std::placeholders::_1));
  mesh_->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Tdim>::compute_pressure,
                std::placeholders::_1, pore_liquid),
      std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
  mesh_->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Tdim>::compute_pressure,
                std::placeholders::_1, pore_gas),
      std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

#ifdef USE_MPI
  int mpi_size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if (mpi_size > 1) {
    mesh_->template nodal_halo_exchange<double, 1>(
        std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1,
                  pore_liquid),
        std::bind(&mpm::NodeBase<Tdim>::assign_pressure,
                  std::placeholders::_1, pore_liquid, std::placeholders::_2));
    mesh_->template nodal_halo_exchange<double, 1>(
        std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1,
                  pore_gas),
        std::bind(&mpm::NodeBase<Tdim>::assign_pressure,
                  std::placeholders::_1, pore_gas, std::placeholders::_2));
  }
#endif

  if (!pressure_smoothing_) {
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_gradient_to_cache,
                  std::placeholders::_1, pore_liquid));
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_gradient_to_cache,
                  std::placeholders::_1, pore_gas));
  }

  if (pressure_smoothing_) {
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_pore_pressure_smoothing_for_pic,
                  std::placeholders::_1));
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_gradient_to_cache,
                  std::placeholders::_1, pore_liquid));
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_gradient_to_cache,
                  std::placeholders::_1, pore_gas));
  }

  // Keep the VTK pressure gradient from the single particle-PIC pressure field.
  // Particle gradients are evaluated by a local least-squares fit, with
  // boundary components extrapolated from interior samples.
}

template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::
    write_particle_vtk_with_seepage_velocity(mpm::Index step,
                                             mpm::Index max_steps) {
  auto vtk_writer = std::make_unique<VtkWriter>(mesh_->particle_coordinates());
  const std::string extension = ".vtp";

  std::vector<std::string> vtk_scalar_data = {
      "porosities",          "PIC_porosities",
      "volumetric_strains",  "dvolumetric_strains",
      "liquid_volumetric_strains",
      "liquid_dvolumetric_strains",
      "thermal_volumetric_strains",
      "dthermal_volumetric_strains",
      "PIC_volumetric_strains",
      "PIC_deviatoric_strains",
      "mean_stresses",       "deviatoric_stresses",
      "densities",           "volumes",
      "PIC_volumes",         "masses",
      "pore_pressures",      "PIC_pore_pressure_excess",
      "PIC_ru",
      "initial_vertical_effective_stresses",
      "dynamic_vertical_effective_stresses",
      "vertical_effective_stress_remaining_ratios",
      "liquefaction_potentials",
      "momentary_liquefied",
      "gamma_sub",
      "PIC_pore_pressures",
      "temperatures",        "PIC_temperatures",
      "temperature_accelerations",
      "solid_fractions",     "current_time",
      "deviatoric_strains",  "suction_pressures",
      "liquid_pressures",    "PIC_liquid_pressures",
      "liquid_pressure_accelerations",
      "liquid_saturations",  "liquid_fractions",
      "liquid_chis",         "liquid_densities",
      "liquid_sources",      "liquid_permeabilities",
      "liquid_volumes",      "liquid_masses",
      "liquid_mass_densities",
      "liquid_vol_strains",  "liquid_viscosities",
      "liquid_critical_times",
      "ice_fractions",       "ice_densities",
      "ice_masses",          "ice_mass_densities",
      "ice_saturations",     "ice_saturation_rates",
      "hydrate_saturations", "PIC_hydrate_saturations",
      "hydrate_fractions",   "hydrate_densities",
      "hydrate_sources",     "hydrate_volumes",
      "hydrate_masses",      "hydrate_mass_densities",
      "gas_pressures",       "PIC_gas_pressures",
      "gas_saturations",     "gas_densities",
      "gas_fractions",       "gas_sources",
      "gas_permeabilities",  "gas_volumes",
      "gas_masses",          "gas_mass_densities",
      "gas_vol_strains",     "gas_viscosities",
      "viscosities",         "gas_critical_times",
      "free_surfaces",       "permeabilities",
      "ids"};

  std::vector<std::string> vtk_vector_data = {
      "displacements",       "temperature_gradients",
      "mass_gradients",      "outward_normals",
      "velocities",          "accelerations",
      "liquid_velocities",   "liquid_accelerations",
      "liquid_fluxes",       "liquid_pressure_gradients",
      "liquid_seepage_velocities", "seepage_velocities",
      "liquid_seepage_forces", "seepage_forces",
      "rotations",           "gas_velocities",
      "gas_accelerations",   "gas_fluxes",
      "gas_pressure_gradients"};

  std::vector<std::string> vtk_tensor_vector_data = {
      "stresses", "strains", "liquid_strains", "gas_strains",
      "thermal_strains", "K_matrix"};

  std::vector<std::string> vtk_tensor_data = {
      "fabric_CNs", "fabric_POs", "velocity_gradients",
      "deformation_gradients", "displacement_gradients", "grad_shapefns"};

  const std::vector<std::string> vtk_state_scalar_data = {
      "AlphaXX",        "AlphaYY",        "AlphaZZ",
      "AlphaZY",        "AlphaZX",        "AlphaXY",
      "AlphaInitialXX", "AlphaInitialYY", "AlphaInitialZZ",
      "AlphaInitialZY", "AlphaInitialZX", "AlphaInitialXY",
      "ZXX",            "ZYY",            "ZZZ",
      "ZZY",            "ZZX",            "ZXY",
      "eps_p_q",        "void_ratio",      "phi",
      "psi",            "cohesion",        "pdstrain"};

  auto particle_vtk_attributes = vtk_attributes_;
  if (std::find(particle_vtk_attributes.begin(), particle_vtk_attributes.end(),
                "liquid_seepage_velocities") == particle_vtk_attributes.end())
    particle_vtk_attributes.emplace_back("liquid_seepage_velocities");
  if (std::find(particle_vtk_attributes.begin(), particle_vtk_attributes.end(),
                "liquid_seepage_forces") == particle_vtk_attributes.end())
    particle_vtk_attributes.emplace_back("liquid_seepage_forces");
  if (std::find(particle_vtk_attributes.begin(), particle_vtk_attributes.end(),
                "gamma_sub") == particle_vtk_attributes.end())
    particle_vtk_attributes.emplace_back("gamma_sub");

  vtk_writer->create_new_dataset();
  for (auto& attribute : particle_vtk_attributes) {
    if (std::find(vtk_scalar_data.begin(), vtk_scalar_data.end(), attribute) !=
        vtk_scalar_data.end()) {
      vtk_writer->write_scalar_point_data(mesh_->particles_scalar_data(attribute),
                                          attribute);
    } else if (std::find(vtk_state_scalar_data.begin(),
                         vtk_state_scalar_data.end(), attribute) !=
               vtk_state_scalar_data.end()) {
      vtk_writer->write_scalar_point_data(
          mesh_->particles_statevars_data(attribute), attribute);
    } else if (std::find(vtk_vector_data.begin(), vtk_vector_data.end(),
                         attribute) != vtk_vector_data.end()) {
      vtk_writer->write_vector_point_data(
          mesh_->template particles_vector_data<3>(attribute), attribute);
    } else if (std::find(vtk_tensor_vector_data.begin(),
                         vtk_tensor_vector_data.end(), attribute) !=
               vtk_tensor_vector_data.end()) {
      vtk_writer->write_tensor_vector_point_data(
          mesh_->template particles_vector_data<6>(attribute), attribute);
    } else if (std::find(vtk_tensor_data.begin(), vtk_tensor_data.end(),
                         attribute) != vtk_tensor_data.end()) {
      vtk_writer->write_tensor_point_data(
          mesh_->template particles_vector_data<9>(attribute), attribute);
    }
  }

  auto file = io_->output_file("particle", extension, uuid_, step, max_steps)
                  .string();
  vtk_writer->write_file(file);
}
#endif

template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::initialise_rigid_pipeline(
    int mpi_rank) {
  if (analysis_.find("rigid_pipeline") == analysis_.end()) return;
  const auto& config = analysis_.at("rigid_pipeline");
  rigid_pipeline_enabled_ = config.at("enable").template get<bool>();
  if (!rigid_pipeline_enabled_) return;
  if (Tdim != 2)
    throw std::runtime_error("The circular rigid pipeline is available only in 2D");

  const auto centre_values =
      config.at("center").template get<std::vector<double>>();
  if (centre_values.size() != 2)
    throw std::runtime_error("rigid_pipeline.center must contain two values");
  const Eigen::Vector2d centre(centre_values[0], centre_values[1]);

  mpm::pipeline::RigidPipelineSection2D section;
  section.outer_radius =
      config.at("outer_radius").template get<double>();
  section.wall_thickness =
      config.at("wall_thickness").template get<double>();
  section.density = config.at("density").template get<double>();
  section.contents_mass_per_length =
      config.at("contents_mass_per_length").template get<double>();
  section.contents_inertia_per_length =
      config.at("contents_inertia_per_length").template get<double>();
  const unsigned surface_markers =
      config.at("surface_markers").template get<unsigned>();

  rigid_pipeline_fixed_ = config.at("fixed").template get<bool>();
  if (config.find("release_time") != config.end())
    rigid_pipeline_release_time_ =
        config.at("release_time").template get<double>();
  rigid_pipeline_allow_rotation_ =
      config.at("allow_rotation").template get<bool>();
  rigid_pipeline_particle_radius_ =
      config.at("particle_radius").template get<double>();
  rigid_pipeline_normal_penalty_ =
      config.at("normal_penalty").template get<double>();
  rigid_pipeline_normal_damping_ =
      config.at("normal_damping").template get<double>();
  rigid_pipeline_tangential_damping_ =
      config.at("tangential_damping").template get<double>();
  rigid_pipeline_friction_ =
      config.at("friction_coefficient").template get<double>();
  rigid_pipeline_fluid_density_ =
      config.at("fluid_density").template get<double>();
  rigid_pipeline_translational_damping_ =
      config.at("translational_damping").template get<double>();
  rigid_pipeline_rotational_damping_ =
      config.at("rotational_damping").template get<double>();
  rigid_pipeline_history_interval_ = output_steps_;
  if (config.find("history_interval") != config.end())
    rigid_pipeline_history_interval_ =
        config.at("history_interval").template get<mpm::Index>();

  if (rigid_pipeline_particle_radius_ < 0. ||
      rigid_pipeline_normal_penalty_ <= 0. ||
      rigid_pipeline_normal_damping_ < 0. ||
      rigid_pipeline_tangential_damping_ < 0. ||
      rigid_pipeline_friction_ < 0. || rigid_pipeline_fluid_density_ < 0. ||
      rigid_pipeline_translational_damping_ < 0. ||
      rigid_pipeline_rotational_damping_ < 0. ||
      rigid_pipeline_release_time_ < 0. ||
      rigid_pipeline_history_interval_ == 0)
    throw std::runtime_error("Invalid rigid_pipeline contact or output parameter");

  rigid_pipeline_ = std::make_unique<mpm::pipeline::RigidPipeline2D>(
      centre, section, surface_markers);

  if (mpi_rank == 0) {
    const auto history_path =
        io_->output_file("pipeline-history", ".csv", uuid_, 0, nsteps_, false)
            .string();
    rigid_pipeline_history_.open(history_path, std::ios::out | std::ios::trunc);
    if (!rigid_pipeline_history_)
      throw std::runtime_error("Failed to create rigid-pipeline history file");
    rigid_pipeline_history_
        << "step,time,fixed,center_x,center_y,displacement_x,displacement_y,"
           "velocity_x,velocity_y,angle,angular_velocity,contact_force_x,"
           "contact_force_y,contact_moment,total_force_x,total_force_y,"
           "total_moment,contacts,max_penetration\n";
    rigid_pipeline_history_ << std::setprecision(16);
    console_->info(
        "Rigid pipeline: fixed={}, release_time={} s, rotation={}, "
        "mass/length={} kg/m, inertia/length={} kg m, friction={}",
        rigid_pipeline_fixed_, rigid_pipeline_release_time_,
        rigid_pipeline_allow_rotation_,
        rigid_pipeline_->mass_per_length(),
        rigid_pipeline_->inertia_per_length(), rigid_pipeline_friction_);
  }
}

template <unsigned Tdim>
mpm::RigidCircleContactResult<Tdim>
mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::map_rigid_pipeline_contact() {
  mpm::RigidCircleContactResult<Tdim> result;
  if (!rigid_pipeline_enabled_) return result;

  Eigen::Matrix<double, Tdim, 1> pipe_centre =
      Eigen::Matrix<double, Tdim, 1>::Zero();
  Eigen::Matrix<double, Tdim, 1> pipe_velocity =
      Eigen::Matrix<double, Tdim, 1>::Zero();
  pipe_centre[0] = rigid_pipeline_->centre().x();
  pipe_centre[1] = rigid_pipeline_->centre().y();
  pipe_velocity[0] = rigid_pipeline_->velocity().x();
  pipe_velocity[1] = rigid_pipeline_->velocity().y();

  std::mutex reduction_mutex;
  mesh_->iterate_over_particles(
      [&](const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle) {
        const auto local = particle->map_rigid_circle_contact_force(
            pipe_centre, pipe_velocity, rigid_pipeline_->angular_velocity(),
            rigid_pipeline_->outer_radius(),
            rigid_pipeline_particle_radius_, rigid_pipeline_normal_penalty_,
            rigid_pipeline_normal_damping_,
            rigid_pipeline_tangential_damping_, rigid_pipeline_friction_);
        if (local.contacts == 0) return;
        std::lock_guard<std::mutex> guard(reduction_mutex);
        result.force += local.force;
        result.moment += local.moment;
        result.contacts += local.contacts;
        result.maximum_penetration =
            std::max(result.maximum_penetration, local.maximum_penetration);
      });

#ifdef USE_MPI
  double sum_values[4] = {result.force[0], result.force[1], result.moment,
                          static_cast<double>(result.contacts)};
  MPI_Allreduce(MPI_IN_PLACE, sum_values, 4, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  result.force[0] = sum_values[0];
  result.force[1] = sum_values[1];
  result.moment = sum_values[2];
  result.contacts = static_cast<unsigned>(sum_values[3]);
  double maximum_penetration = result.maximum_penetration;
  MPI_Allreduce(&maximum_penetration, &result.maximum_penetration, 1,
                MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  return result;
}

template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::write_rigid_pipeline_output(
    const RigidCircleContactResult<Tdim>& contact,
    const Eigen::Vector2d& total_force, double total_moment, int mpi_rank) {
  if (!rigid_pipeline_enabled_ || mpi_rank != 0 ||
      step_ % rigid_pipeline_history_interval_ != 0)
    return;

  const auto displacement = rigid_pipeline_->displacement();
  const bool currently_fixed =
      rigid_pipeline_fixed_ || current_time_ < rigid_pipeline_release_time_;
  rigid_pipeline_history_
      << step_ << ',' << current_time_ << ',' << currently_fixed << ','
      << rigid_pipeline_->centre().x() << ',' << rigid_pipeline_->centre().y()
      << ',' << displacement.x() << ',' << displacement.y() << ','
      << rigid_pipeline_->velocity().x() << ','
      << rigid_pipeline_->velocity().y() << ',' << rigid_pipeline_->angle()
      << ',' << rigid_pipeline_->angular_velocity() << ',' << contact.force[0]
      << ',' << contact.force[1] << ',' << contact.moment << ','
      << total_force.x() << ',' << total_force.y() << ',' << total_moment << ','
      << contact.contacts << ',' << contact.maximum_penetration << '\n';
  rigid_pipeline_history_.flush();

  if (!write_vtk_) return;

  const auto positions = rigid_pipeline_->marker_positions();
  const auto velocities = rigid_pipeline_->marker_velocities();
  const auto file =
      io_->output_file("pipeline", ".vtp", uuid_, step_, nsteps_, false)
          .string();
  std::ofstream stream(file, std::ios::out | std::ios::trunc);
  if (!stream) throw std::runtime_error("Failed to create pipeline VTP output");
  stream << std::setprecision(16)
         << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"PolyData\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n"
         << "  <PolyData>\n"
         << "    <Piece NumberOfPoints=\"" << positions.rows()
         << "\" NumberOfVerts=\"0\" NumberOfLines=\"1\" "
            "NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
         << "      <PointData Vectors=\"velocity\">\n"
         << "        <DataArray type=\"Float64\" Name=\"velocity\" "
            "NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (Eigen::Index i = 0; i < velocities.rows(); ++i)
    stream << "          " << velocities(i, 0) << ' ' << velocities(i, 1)
           << " 0\n";
  stream << "        </DataArray>\n"
         << "      </PointData>\n"
         << "      <Points>\n"
         << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            "format=\"ascii\">\n";
  for (Eigen::Index i = 0; i < positions.rows(); ++i)
    stream << "          " << positions(i, 0) << ' ' << positions(i, 1)
           << " 0\n";
  stream << "        </DataArray>\n"
         << "      </Points>\n"
         << "      <Lines>\n"
         << "        <DataArray type=\"Int32\" Name=\"connectivity\" "
            "format=\"ascii\">\n          ";
  for (Eigen::Index i = 0; i < positions.rows(); ++i) stream << i << ' ';
  stream << "0\n        </DataArray>\n"
         << "        <DataArray type=\"Int32\" Name=\"offsets\" "
            "format=\"ascii\">\n          "
         << positions.rows() + 1 << "\n        </DataArray>\n"
         << "      </Lines>\n"
         << "    </Piece>\n"
         << "  </PolyData>\n"
         << "</VTKFile>\n";
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////                                                                ////////
////////             THM-MPM Explicit ThreePhase Lag Solver             ////////
////////                                                                ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//! Thermo-hydro-mechncial MPM Explicit solver
template <unsigned Tdim>
bool mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::solve() {
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

  // Three phases (soil skeleton and pore liquid)
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
  if (analysis_.find("pressure_smoothing_in_loop") != analysis_.end())
    pressure_smoothing_in_loop_ =
        analysis_["pressure_smoothing_in_loop"].template get<bool>();
  if (analysis_.find("pressure_smoothing_iterations") != analysis_.end())
    pressure_smoothing_iterations_ =
        std::max(1U, analysis_["pressure_smoothing_iterations"].template get<unsigned>());

  // Interface
  if (analysis_.find("interface") != analysis_.end())
    interface_ = analysis_.at("interface").template get<bool>();

  // Free surface
  if (analysis_.find("free_surface") != analysis_.end()) {
    free_surface_particle_ = analysis_["free_surface"]["free_surface_particle"]
                                  .template get<std::string>();
    //Get volume tolerance for free surface
    volume_tolerance_ = analysis_["free_surface"]["volume_tolerance"]
                                  .template get<double>();
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
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_affine_mpm,
                std::placeholders::_1, this->affine_mpm_));

  if (!this->initialise_prescribed_phase_pressures()) {
    status = false;
    throw std::runtime_error("Initialisation of prescribed phase pressures failed");
  }

  // Initialise loading conditions
  bool loading_status = this->initialise_loads();
  if (!loading_status) {
    status = false;
    throw std::runtime_error("Initialisation of loads failed");
  }

  // Initialise vtk output for three phase
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

  // Assign initial particles at surface but not free, like tunnel
  bool assign_nonfree_surface_status = mesh_->assign_nonfree_surface_particles(io_);
  if (!assign_nonfree_surface_status) {
    status = false;
    throw std::runtime_error("Initialisation nonfree surface particles failed");
  }

  // Check point resume
  if (resume) {
    if (!this->checkpoint_resume())
      throw std::runtime_error(
          "Checkpoint resume was requested but could not be completed");
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::assign_affine_mpm,
                  std::placeholders::_1, this->affine_mpm_));
    this->current_time_ = analysis_["resume"]["current_time"].template get<double>();
    std::cout << "current_time" << this->current_time_ << "\n";
  }

  // Resolve pressure-database samples once, after a possible checkpoint has
  // replaced particle state. Subsequent lookup is by particle id only, so a
  // moving particle retains its original pressure history.
  if (write_prescribed_phase_pressures_) {
    this->write_prescribed_pressure_points();
    this->write_prescribed_pressure_binary_header();
  }
  if (prescribed_phase_pressures_) this->bind_prescribed_pressure_samples();

  solver_begin = std::chrono::steady_clock::now();

  this->compute_critical_timestep_size(dt_);

  // The fixed equilibrium and released SANISAND stages use the same pipe
  // implementation; only rigid_pipeline.fixed changes between their inputs.
  this->initialise_rigid_pipeline(mpi_rank);

  if (write_prescribed_phase_pressures_) {
    if (std::abs(current_time_) >
        1.0e-12 * std::max(1.0, std::abs(dt_)))
      throw std::runtime_error(
          "A prescribed pressure database must start at analysis time zero");
    // Store the true initial state. Later frame numbers then correspond to
    // source_step * source_dt, including when frames are written sparsely.
    this->write_prescribed_pressure_frame(0);
  }

////////////////////////////////////////////////////////////////////////////////
////////////////                  MAIN LOOP               //////////////////////
////////////////////////////////////////////////////////////////////////////////

  // Main loop
  // step_ is the number of completed updates. Starting at one keeps checkpoint
  // and pressure-database frame numbers aligned with physical time while still
  // performing exactly nsteps_ updates.
  for (step_ = 1; step_ <= nsteps_; ++step_) {

    mpm::RigidCircleContactResult<Tdim> rigid_pipeline_contact;
    Eigen::Vector2d rigid_pipeline_total_force = Eigen::Vector2d::Zero();
    double rigid_pipeline_total_moment = 0.;

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

    if (mpi_rank == 0 &&
        (step_ <= 10 || step_ % output_steps_ == 0 || step_ == nsteps_))
      console_->info(
          "uuid : [{}], Step: {} of {}, timestep = {}, time = {}.\n",
          uuid_, step_, nsteps_, dt_, current_time_);

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

    // Assign mass and momentum to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1));

    // Assign heat capacity and heat to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_heat_to_nodes,
                  std::placeholders::_1));

    // Apply particle velocity constraints
    mesh_->apply_moving_rigid_boundary(current_time_, dt_);

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

    // Compute free surface cells, nodes, and particles
    // Cell volume was mapped above, before the support-threshold velocity
    // calculation. Reuse that same mapping so the acceleration stage applies
    // an identical density threshold rather than seeing twice the volume.
    mesh_->compute_free_surface(free_surface_particle_, volume_tolerance_,
                                false);

    // // Assign heat capacity and heat to nodes
    // mesh_->iterate_over_particles(
    //     std::bind(&mpm::ParticleBase<Tdim>::map_mass_pressure_to_nodes,
    //               std::placeholders::_1));

    // // Compute nodal pressure
    // mesh_->iterate_over_nodes_predicate(
    //     std::bind(&mpm::NodeBase<Tdim>::compute_pressure,
    //               std::placeholders::_1, pore_liquid),
    //     std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // // Compute nodal pressure
    // mesh_->iterate_over_nodes_predicate(
    //     std::bind(&mpm::NodeBase<Tdim>::compute_pressure,
    //               std::placeholders::_1, pore_gas),
    //     std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Initialise wave conditions
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::build_wave_pf_context, std::placeholders::_1));

    // Apply nodal temperature constraints
    mesh_->apply_nodal_temperature_constraints(soil_skeleton, current_time_);

    // Apply nodal temperature constraints
    mesh_->apply_nodal_convective_heat_constraints(soil_skeleton, current_time_); 

    // Iterate over each particle to calculate mechancial strain
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_strain, std::placeholders::_1, dt_));

    // Iterate over each particle to calculate thermal strain
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_thermal_strain, std::placeholders::_1));

    // Iterate over each particle to compute stress
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_stress, std::placeholders::_1));

    // Iterate over each particle to compute or prescribe pore pressure
    if (prescribed_phase_pressures_) {
      this->assign_prescribed_phase_pressures();
    } else {
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_pore_pressure,
                    std::placeholders::_1, dt_));
    }

    mesh_->apply_particle_pore_pressure_constraints(current_time_);

    // Pressure smoothing with state write-back for the next time step
    if (pressure_smoothing_ && pressure_smoothing_in_loop_) {
      this->pressure_smoothing(pore_liquid);
    }

    if (write_prescribed_phase_pressures_)
      this->write_prescribed_pressure_frame(step_);

    // Iterate over each particle to calculate particle porosity
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_porosity, std::placeholders::_1, dt_));

    // Iterate over each particle to calculate particle saturation
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_liquid_gas_saturation, std::placeholders::_1, dt_));

    // Iterate over each particle to update material density of particle
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_density, std::placeholders::_1, dt_));

    // Iterate over each particle to update material density of particle
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::update_particle_volume, std::placeholders::_1));

    // Apply particle traction and map to nodes
    mesh_->apply_traction_on_particles(current_time_);

    // Iterate over particles to compute nodal body force
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_external_force,
                  std::placeholders::_1, this->gravity_));

    // Iterate over particles to compute nodal mixture and fluid internal
    // force
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_internal_force, std::placeholders::_1));

    // Couple the circular rigid body to the soil-mixture equation. The pipe is
    // kept fixed during equilibrium and integrated after release. Since this
    // case is fully submerged, Archimedes buoyancy is applied explicitly.
    if (rigid_pipeline_enabled_) {
      rigid_pipeline_->clear_forces();
      rigid_pipeline_contact = this->map_rigid_pipeline_contact();
      const Eigen::Vector2d contact_force(rigid_pipeline_contact.force[0],
                                          rigid_pipeline_contact.force[1]);
      rigid_pipeline_->apply_force(contact_force, rigid_pipeline_->centre());
      if (rigid_pipeline_allow_rotation_)
        rigid_pipeline_->apply_moment(rigid_pipeline_contact.moment);

      Eigen::Vector2d pipe_gravity(this->gravity_[0], this->gravity_[1]);
      rigid_pipeline_->apply_body_acceleration(pipe_gravity);
      constexpr double pi = 3.14159265358979323846;
      const Eigen::Vector2d buoyancy =
          -rigid_pipeline_fluid_density_ * pi *
          rigid_pipeline_->outer_radius() * rigid_pipeline_->outer_radius() *
          pipe_gravity;
      rigid_pipeline_->apply_force(buoyancy, rigid_pipeline_->centre());

      rigid_pipeline_total_force = rigid_pipeline_->resultant_force();
      rigid_pipeline_total_moment = rigid_pipeline_contact.moment;
      const bool currently_fixed =
          rigid_pipeline_fixed_ || current_time_ < rigid_pipeline_release_time_;
      if (currently_fixed) {
        rigid_pipeline_->clear_forces();
      } else {
        rigid_pipeline_->advance(dt_, rigid_pipeline_translational_damping_,
                                 rigid_pipeline_rotational_damping_);
      }
    }

    // Iterate over particles to compute nodal drag force coefficient
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_drag_force_coefficient,
                  std::placeholders::_1));

    // // Pressure smoothing
    // if (pressure_smoothing_) { this->pressure_smoothing(pore_liquid);}

    // Compute nodal acceleration and update nodal velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acc_vel_threephase_explicit,
          std::placeholders::_1, soil_skeleton, pore_liquid, pore_gas, mixture,
          this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Apply nodal and particle velocity constraints
    mesh_->apply_velocity_constraints(current_time_);

    // Iterate over each particle to compute updated position
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_updated_velocity,
        std::placeholders::_1, this->dt_, this->pic_, damping_factor_));

    // Apply nodal and particle velocity constraints
    mesh_->apply_velocity_constraints(current_time_);

    // Iterate over each particle to compute nodal heat conduction
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_heat_conduction, std::placeholders::_1));

    // // Iterate over each particle to compute nodal heat convection
    // mesh_->iterate_over_particles(std::bind(
    //     &mpm::ParticleBase<Tdim>::map_heat_convection, std::placeholders::_1));

    // // Iterate over each particle to compute nodal heat convection
    // mesh_->iterate_over_particles(std::bind(
    //     &mpm::ParticleBase<Tdim>::map_plastic_work, std::placeholders::_1, dt_));

    // Compute nodal temperature acceleration and update nodal temperature
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_temperature,
                  std::placeholders::_1, soil_skeleton, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Apply nodal temperature constraints
    mesh_->apply_nodal_temperature_constraints(soil_skeleton, current_time_);

    // // Iterate over each particle to compute updated temperature
    // mesh_->iterate_over_particles(std::bind(
    //     &mpm::ParticleBase<Tdim>::update_particle_temperature,
    //     std::placeholders::_1, this->dt_, this->pic_t_));

    // // Iterate over each particle to output results
    // mesh_->iterate_over_particles(std::bind(
    //     &mpm::ParticleBase<Tdim>::output_results,
    //     std::placeholders::_1, this->step_));

    // // Apply particle temperature constraints
    // mesh_->apply_particle_temperature_constraints(current_time_);

    // Locate particles
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty()) {
      std::ostringstream message;
      message << "Particle outside the mesh domain: count="
              << unlocatable_particles.size() << ", step=" << step_
              << ", time=" << current_time_;
      const std::size_t diagnostic_count =
          std::min<std::size_t>(unlocatable_particles.size(), 8);
      for (std::size_t i = 0; i < diagnostic_count; ++i) {
        const auto& particle = unlocatable_particles[i];
        const auto coordinates = particle->coordinates();
        const auto velocity = particle->vector_data("velocities");
        message << "; particle[id=" << particle->id() << ", coordinates=("
                << coordinates.transpose() << "), velocity=("
                << velocity.transpose() << "), finite="
                << (coordinates.allFinite() && velocity.allFinite() ? "yes"
                                                                     : "no")
                << "]";
      }
      throw std::runtime_error(message.str());
    }

    this->write_rigid_pipeline_output(
        rigid_pipeline_contact, rigid_pipeline_total_force,
        rigid_pipeline_total_moment, mpi_rank);


    // Fixed timestep, and data output linearly (every const steps = output_steps)
    if ((!variable_timestep_) & (!log_output_steps_)){
      if (step_ % output_steps_ == 0) {
        // HDF5 outputs
        if (write_hdf5_) this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
        // VTK outputs
        this->write_vtk(this->step_, this->nsteps_);
        this->prepare_seepage_velocity_vtk();
        this->write_particle_vtk_with_seepage_velocity(this->step_,
                                                       this->nsteps_);
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
        this->prepare_seepage_velocity_vtk();
        this->write_particle_vtk_with_seepage_velocity(this->step_,
                                                       this->nsteps_);
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
        this->prepare_seepage_velocity_vtk();
        this->write_particle_vtk_with_seepage_velocity(this->step_,
                                                       this->nsteps_);
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
////////////////                LOOP ENDING               //////////////////////
////////////////////////////////////////////////////////////////////////////////

//! Initialise externally prescribed phase pressures
template <unsigned Tdim>
bool mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::
    initialise_prescribed_phase_pressures() {
  bool status = true;
  try {
    if (analysis_.find("prescribed_phase_pressures") == analysis_.end()) {
      prescribed_phase_pressures_ = false;
      write_prescribed_phase_pressures_ = false;
      return status;
    }

    const auto& pressure_props = analysis_["prescribed_phase_pressures"];
    if (pressure_props.find("enable") != pressure_props.end())
      prescribed_phase_pressures_ =
          pressure_props["enable"].template get<bool>();

    if (pressure_props.find("write") != pressure_props.end())
      write_prescribed_phase_pressures_ =
          pressure_props["write"].template get<bool>();

    if (!prescribed_phase_pressures_ && !write_prescribed_phase_pressures_)
      return status;

    if (prescribed_phase_pressures_ && write_prescribed_phase_pressures_)
      throw std::runtime_error(
          "Prescribed phase pressures cannot be read and written in the same "
          "analysis");

    if (pressure_props.find("path") == pressure_props.end())
      throw std::runtime_error("prescribed_phase_pressures.path is required");

    prescribed_pressure_path_ =
        pressure_props["path"].template get<std::string>();

    if (pressure_props.find("file_prefix") != pressure_props.end())
      prescribed_pressure_prefix_ =
          pressure_props["file_prefix"].template get<std::string>();

    if (pressure_props.find("mapping") != pressure_props.end())
      prescribed_pressure_mapping_ =
          pressure_props["mapping"].template get<std::string>();

    if (pressure_props.find("source_dt") != pressure_props.end())
      prescribed_pressure_source_dt_ =
          pressure_props["source_dt"].template get<double>();
    else
      prescribed_pressure_source_dt_ = dt_;

    if (pressure_props.find("step_interval") != pressure_props.end())
      prescribed_pressure_step_interval_ =
          pressure_props["step_interval"].template get<mpm::Index>();
    else
      prescribed_pressure_step_interval_ = output_steps_;

    if (prescribed_pressure_source_dt_ <= 0.0 ||
        prescribed_pressure_step_interval_ == 0)
      throw std::runtime_error("Invalid prescribed pressure time settings");

    const double expected_time_interval =
        prescribed_pressure_source_dt_ *
        static_cast<double>(prescribed_pressure_step_interval_);

    if (pressure_props.find("time_interval") != pressure_props.end())
      prescribed_pressure_time_interval_ =
          pressure_props["time_interval"].template get<double>();
    else
      prescribed_pressure_time_interval_ = expected_time_interval;

    const double time_tolerance =
        1.0e-12 * std::max(1.0, std::abs(expected_time_interval));
    if (prescribed_pressure_time_interval_ <= 0.0 ||
        std::abs(prescribed_pressure_time_interval_ - expected_time_interval) >
            time_tolerance)
      throw std::runtime_error(
          "prescribed_phase_pressures.time_interval must equal source_dt * "
          "step_interval");

    if (pressure_props.find("max_step") != pressure_props.end())
      prescribed_pressure_max_step_ =
          pressure_props["max_step"].template get<mpm::Index>();
    else
      prescribed_pressure_max_step_ = nsteps_;

    if (pressure_props.find("coordinate_bin_size") != pressure_props.end())
      prescribed_pressure_coordinate_bin_size_ =
          pressure_props["coordinate_bin_size"].template get<double>();

    const bool map_by_id = prescribed_pressure_mapping_ == "id" ||
                           prescribed_pressure_mapping_ == "particle_id";
    const bool map_by_coordinates =
        prescribed_pressure_mapping_ == "coordinate" ||
        prescribed_pressure_mapping_ == "coordinates" ||
        prescribed_pressure_mapping_ == "nearest";
    if (!map_by_id && !map_by_coordinates)
      throw std::runtime_error(
          "prescribed_phase_pressures.mapping must be id, particle_id, "
          "coordinate, coordinates, or nearest");

    if (write_prescribed_phase_pressures_) {
      mpm::prescribed_pressure::detail::validate_write_max_step(
          prescribed_pressure_max_step_, nsteps_);
      const double dt_tolerance =
          1.0e-12 * std::max(1.0, std::abs(dt_));
      if (variable_timestep_ ||
          std::abs(prescribed_pressure_source_dt_ - dt_) > dt_tolerance)
        throw std::runtime_error(
            "Writing a prescribed pressure database requires a fixed solver "
            "dt equal to source_dt");
    }

    if (prescribed_pressure_max_step_ % prescribed_pressure_step_interval_ !=
        0) {
      const auto aligned_max_step =
          prescribed_pressure_max_step_ -
          prescribed_pressure_max_step_ % prescribed_pressure_step_interval_;
      console_->warn(
          "Prescribed pressure max_step [{}] is not a stored-frame step; "
          "using [{}]",
          prescribed_pressure_max_step_, aligned_max_step);
      prescribed_pressure_max_step_ = aligned_max_step;
    }

    if (write_prescribed_phase_pressures_) {
      std::filesystem::create_directories(prescribed_pressure_path_);
    }

    if (prescribed_phase_pressures_) {
      this->read_prescribed_pressure_points();
      this->validate_prescribed_pressure_binary_header();
    }

    if (prescribed_pressure_max_step_ % prescribed_pressure_step_interval_ !=
        0) {
      const auto aligned_max_step =
          prescribed_pressure_max_step_ -
          prescribed_pressure_max_step_ % prescribed_pressure_step_interval_;
      console_->warn(
          "Prescribed pressure max_step [{}] is not a stored-frame step; "
          "using [{}]",
          prescribed_pressure_max_step_, aligned_max_step);
      prescribed_pressure_max_step_ = aligned_max_step;
    }

    prescribed_pressure_frame_cache_valid_.fill(false);
    prescribed_pressure_frame_cache_use_.fill(0);
    prescribed_pressure_frame_cache_counter_ = 0;

    console_->info(
        "Prescribed phase pressures: read [{}], write [{}], path [{}], mapping [{}], binary values [{}]",
        prescribed_phase_pressures_, write_prescribed_phase_pressures_,
        prescribed_pressure_path_, prescribed_pressure_mapping_,
        this->prescribed_pressure_values_filename());
  } catch (std::exception& exception) {
    console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                    exception.what());
    status = false;
  }
  return status;
}

//! Prescribed pressure points filename
template <unsigned Tdim>
std::string mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::
    prescribed_pressure_points_filename() const {
  const bool has_separator =
      !prescribed_pressure_path_.empty() &&
      (prescribed_pressure_path_.back() == '/' ||
       prescribed_pressure_path_.back() == '\\');

  return prescribed_pressure_path_ + (has_separator ? "" : "/") +
         prescribed_pressure_prefix_ + "_points.txt";
}

//! Prescribed pressure binary values filename
template <unsigned Tdim>
std::string mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::
    prescribed_pressure_values_filename() const {
  const bool has_separator =
      !prescribed_pressure_path_.empty() &&
      (prescribed_pressure_path_.back() == '/' ||
       prescribed_pressure_path_.back() == '\\');

  return prescribed_pressure_path_ + (has_separator ? "" : "/") +
         prescribed_pressure_prefix_ + "_values.bin";
}

//! Byte offset of one binary pressure frame
template <unsigned Tdim>
std::streamoff mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::
    prescribed_pressure_frame_offset(mpm::Index source_step) const {
  const std::streamoff header_size =
      16 + 2 * sizeof(std::uint32_t) + 3 * sizeof(std::uint64_t) +
      sizeof(double);
  const std::streamoff frame_size =
      sizeof(std::uint64_t) +
      static_cast<std::streamoff>(prescribed_pressure_particle_count_) * 2 *
          sizeof(double);
  const auto frame_id = static_cast<std::uint64_t>(
      source_step / prescribed_pressure_step_interval_);
  return header_size + static_cast<std::streamoff>(frame_id) * frame_size;
}

//! Write the fixed particle coordinates used by the pressure database
template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::
    write_prescribed_pressure_points() {
  const auto filename = this->prescribed_pressure_points_filename();
  std::ofstream file(filename.c_str(), std::ios::out);
  if (!file.is_open())
    throw std::runtime_error("Cannot open prescribed pressure points file: " +
                             filename);

  const auto particle_ids = mesh_->particles_scalar_data("ids");
  const auto coordinates = mesh_->particle_coordinates();
  if (particle_ids.size() != coordinates.size())
    throw std::runtime_error("Particle pressure point arrays are inconsistent");

  prescribed_pressure_particle_count_ = particle_ids.size();
  file << "# particle_count " << prescribed_pressure_particle_count_ << '\n';
  file << "# id";
  for (unsigned i = 0; i < Tdim; ++i) file << " x" << i;
  file << '\n';
  file << std::setprecision(16);

  for (std::size_t i = 0; i < particle_ids.size(); ++i) {
    file << static_cast<mpm::Index>(particle_ids[i]);
    for (unsigned j = 0; j < Tdim; ++j) file << ' ' << coordinates[i][j];
    file << '\n';
  }
}

//! Write the binary pressure database header
template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::
    write_prescribed_pressure_binary_header() {
  prescribed_pressure_legacy_step_offset_ = false;
  const auto filename = this->prescribed_pressure_values_filename();
  std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary |
                                        std::ios::trunc);
  if (!file.is_open())
    throw std::runtime_error("Cannot open prescribed pressure binary file: " +
                             filename);

  const char magic[16] = "MPM_PRESSURE_V2";
  const std::uint32_t dim = Tdim;
  const std::uint32_t values_per_particle = 2;
  const std::uint64_t particle_count =
      static_cast<std::uint64_t>(prescribed_pressure_particle_count_);
  const std::uint64_t step_interval =
      static_cast<std::uint64_t>(prescribed_pressure_step_interval_);
  const std::uint64_t max_step =
      static_cast<std::uint64_t>(prescribed_pressure_max_step_);
  const double source_dt = prescribed_pressure_source_dt_;

  file.write(magic, sizeof(magic));
  file.write(reinterpret_cast<const char*>(&dim), sizeof(dim));
  file.write(reinterpret_cast<const char*>(&values_per_particle),
             sizeof(values_per_particle));
  file.write(reinterpret_cast<const char*>(&particle_count),
             sizeof(particle_count));
  file.write(reinterpret_cast<const char*>(&step_interval),
             sizeof(step_interval));
  file.write(reinterpret_cast<const char*>(&max_step), sizeof(max_step));
  file.write(reinterpret_cast<const char*>(&source_dt), sizeof(source_dt));
}

//! Read the fixed particle coordinates used by the pressure database
template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::
    read_prescribed_pressure_points() {
  const auto filename = this->prescribed_pressure_points_filename();
  std::ifstream file(filename.c_str());
  if (!file.is_open())
    throw std::runtime_error("Cannot open prescribed pressure points file: " +
                             filename);

  prescribed_pressure_samples_.clear();
  prescribed_pressure_id_to_sample_.clear();

  std::string line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::replace(line.begin(), line.end(), ',', ' ');

    std::istringstream line_stream(line);
    std::vector<double> values;
    double value = 0.0;
    while (line_stream >> value) values.emplace_back(value);
    if (values.size() < static_cast<std::size_t>(Tdim + 1)) continue;

    PressureSample sample;
    sample.id = static_cast<mpm::Index>(values[0]);
    for (unsigned i = 0; i < Tdim; ++i) sample.coordinates[i] = values[1 + i];
    sample.has_coordinates = true;
    const auto inserted = prescribed_pressure_id_to_sample_.emplace(
        sample.id, prescribed_pressure_samples_.size());
    if (!inserted.second)
      throw std::runtime_error(
          "Duplicate particle id in prescribed pressure points file: " +
          std::to_string(sample.id));
    prescribed_pressure_samples_.emplace_back(sample);
  }

  if (prescribed_pressure_samples_.empty())
    throw std::runtime_error("Prescribed pressure points file has no samples: " +
                             filename);

  prescribed_pressure_particle_count_ = prescribed_pressure_samples_.size();
  this->build_prescribed_pressure_spatial_index();
}

//! Validate the binary pressure database header
template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::
    validate_prescribed_pressure_binary_header() {
  const auto filename = this->prescribed_pressure_values_filename();
  std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);
  if (!file.is_open())
    throw std::runtime_error("Cannot open prescribed pressure binary file: " +
                             filename);

  char magic[16] = {};
  std::uint32_t dim = 0;
  std::uint32_t values_per_particle = 0;
  std::uint64_t particle_count = 0;
  std::uint64_t step_interval = 0;
  std::uint64_t max_step = 0;
  double source_dt = 0.0;

  file.read(magic, sizeof(magic));
  file.read(reinterpret_cast<char*>(&dim), sizeof(dim));
  file.read(reinterpret_cast<char*>(&values_per_particle),
            sizeof(values_per_particle));
  file.read(reinterpret_cast<char*>(&particle_count), sizeof(particle_count));
  file.read(reinterpret_cast<char*>(&step_interval), sizeof(step_interval));
  file.read(reinterpret_cast<char*>(&max_step), sizeof(max_step));
  file.read(reinterpret_cast<char*>(&source_dt), sizeof(source_dt));

  const std::string pressure_format(magic, magic + 15);
  const bool format_v1 = pressure_format == "MPM_PRESSURE_V1";
  const bool format_v2 = pressure_format == "MPM_PRESSURE_V2";
  if ((!format_v1 && !format_v2) || dim != Tdim || values_per_particle != 2)
    throw std::runtime_error("Invalid prescribed pressure binary header: " +
                             filename);
  prescribed_pressure_legacy_step_offset_ = format_v1;
  if (format_v1)
    console_->warn(
        "Reading legacy MPM_PRESSURE_V1 timing: frame zero represents "
        "source_dt rather than time zero");

  if (particle_count != prescribed_pressure_particle_count_)
    throw std::runtime_error(
        "Prescribed pressure points and binary particle counts differ");

  if (step_interval != prescribed_pressure_step_interval_)
    throw std::runtime_error(
        "Prescribed pressure binary step_interval differs from JSON");

  if (std::abs(source_dt - prescribed_pressure_source_dt_) >
      1.0e-12 * std::max(1.0, prescribed_pressure_source_dt_))
    throw std::runtime_error(
        "Prescribed pressure binary source_dt differs from JSON: file=" +
        std::to_string(source_dt) +
        ", JSON=" + std::to_string(prescribed_pressure_source_dt_));

  if (max_step < prescribed_pressure_max_step_)
    prescribed_pressure_max_step_ = static_cast<mpm::Index>(max_step);
}

//! Write one prescribed pressure frame
template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::write_prescribed_pressure_frame(
    mpm::Index source_step) {
  if (source_step % prescribed_pressure_step_interval_ != 0) return;

  const auto filename = this->prescribed_pressure_values_filename();
  std::fstream file(filename.c_str(), std::ios::in | std::ios::out |
                                      std::ios::binary);
  if (!file.is_open())
    throw std::runtime_error("Cannot open prescribed pressure binary file: " +
                             filename);

  const auto particle_ids = mesh_->particles_scalar_data("ids");
  const auto liquid_pressures =
      mesh_->particles_scalar_data("PIC_liquid_pressures");
  const auto gas_pressures = mesh_->particles_scalar_data("PIC_gas_pressures");

  if (particle_ids.size() != prescribed_pressure_particle_count_ ||
      particle_ids.size() != liquid_pressures.size() ||
      particle_ids.size() != gas_pressures.size())
    throw std::runtime_error("Particle pressure output arrays are inconsistent");

  const std::uint64_t stored_step = static_cast<std::uint64_t>(source_step);
  file.seekp(this->prescribed_pressure_frame_offset(source_step));
  file.write(reinterpret_cast<const char*>(&stored_step), sizeof(stored_step));
  file.write(reinterpret_cast<const char*>(liquid_pressures.data()),
             static_cast<std::streamsize>(liquid_pressures.size() *
                                          sizeof(double)));
  file.write(reinterpret_cast<const char*>(gas_pressures.data()),
             static_cast<std::streamsize>(gas_pressures.size() *
                                          sizeof(double)));
}

//! Return a cached prescribed pressure frame
template <unsigned Tdim>
const typename mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::PressureFrame&
mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::prescribed_pressure_frame(
    mpm::Index source_step) {
  if (source_step % prescribed_pressure_step_interval_ != 0)
    throw std::runtime_error(
        "Requested prescribed pressure step is not a stored-frame step");

  for (std::size_t i = 0; i < prescribed_pressure_frame_cache_.size(); ++i) {
    if (prescribed_pressure_frame_cache_valid_[i] &&
        prescribed_pressure_frame_cache_[i].step == source_step) {
      prescribed_pressure_frame_cache_use_[i] =
          ++prescribed_pressure_frame_cache_counter_;
      return prescribed_pressure_frame_cache_[i];
    }
  }

  std::size_t cache_slot = prescribed_pressure_frame_cache_.size();
  for (std::size_t i = 0; i < prescribed_pressure_frame_cache_.size(); ++i) {
    if (!prescribed_pressure_frame_cache_valid_[i]) {
      cache_slot = i;
      break;
    }
  }
  if (cache_slot == prescribed_pressure_frame_cache_.size())
    cache_slot = prescribed_pressure_frame_cache_use_[0] <=
                         prescribed_pressure_frame_cache_use_[1]
                     ? 0
                     : 1;

  prescribed_pressure_frame_cache_[cache_slot] =
      this->read_prescribed_pressure_frame(source_step);
  prescribed_pressure_frame_cache_valid_[cache_slot] = true;
  prescribed_pressure_frame_cache_use_[cache_slot] =
      ++prescribed_pressure_frame_cache_counter_;
  return prescribed_pressure_frame_cache_[cache_slot];
}

//! Read one prescribed pressure frame
template <unsigned Tdim>
typename mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::PressureFrame
mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::read_prescribed_pressure_frame(
    mpm::Index source_step) {
  PressureFrame frame;
  frame.step = source_step;
  frame.liquid_pressures.resize(prescribed_pressure_particle_count_);
  frame.gas_pressures.resize(prescribed_pressure_particle_count_);

  const auto filename = this->prescribed_pressure_values_filename();
  std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);
  if (!file.is_open())
    throw std::runtime_error("Cannot open prescribed pressure binary file: " +
                             filename);

  std::uint64_t stored_step = 0;
  file.seekg(this->prescribed_pressure_frame_offset(source_step));
  file.read(reinterpret_cast<char*>(&stored_step), sizeof(stored_step));
  file.read(reinterpret_cast<char*>(frame.liquid_pressures.data()),
            static_cast<std::streamsize>(frame.liquid_pressures.size() *
                                         sizeof(double)));
  file.read(reinterpret_cast<char*>(frame.gas_pressures.data()),
            static_cast<std::streamsize>(frame.gas_pressures.size() *
                                         sizeof(double)));

  if (!file || stored_step != static_cast<std::uint64_t>(source_step))
    throw std::runtime_error("Cannot read prescribed pressure frame " +
                             std::to_string(source_step) + " from " +
                             filename);

  return frame;
}

//! Build a simple coordinate bin index for nearest-neighbour lookup
template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::
    build_prescribed_pressure_spatial_index() {
  prescribed_pressure_has_spatial_index_ = false;
  prescribed_pressure_coordinate_bins_.clear();

  bool first = true;
  for (const auto& sample : prescribed_pressure_samples_) {
    if (!sample.has_coordinates) continue;
    if (first) {
      prescribed_pressure_min_coordinates_ = sample.coordinates;
      first = false;
    } else {
      for (unsigned i = 0; i < Tdim; ++i)
        prescribed_pressure_min_coordinates_[i] =
            std::min(prescribed_pressure_min_coordinates_[i],
                     sample.coordinates[i]);
    }
  }

  if (first) return;

  const auto bin_size =
      std::max(prescribed_pressure_coordinate_bin_size_, 1.0e-12);
  for (std::size_t i = 0; i < prescribed_pressure_samples_.size(); ++i) {
    if (!prescribed_pressure_samples_[i].has_coordinates) continue;
    const auto key = this->pressure_bin_key(
        prescribed_pressure_samples_[i].coordinates,
        prescribed_pressure_min_coordinates_, bin_size);
    prescribed_pressure_coordinate_bins_[key].emplace_back(i);
  }
  prescribed_pressure_has_spatial_index_ = true;
}

//! Hash coordinate bin to a key
template <unsigned Tdim>
long long mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::pressure_bin_key(
    const Eigen::Matrix<double, Tdim, 1>& coordinates,
    const Eigen::Matrix<double, Tdim, 1>& min_coordinates,
    double bin_size) const {
  unsigned long long key = 1469598103934665603ULL;
  for (unsigned i = 0; i < Tdim; ++i) {
    const long long bin =
        static_cast<long long>(std::floor((coordinates[i] - min_coordinates[i]) /
                                          std::max(bin_size, 1.0e-12)));
    key ^= static_cast<unsigned long long>(bin + 1099511628211ULL);
    key *= 1099511628211ULL;
  }
  return static_cast<long long>(key);
}

//! Bind every dynamic particle to one database sample at the start of the stage
template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::
    bind_prescribed_pressure_samples() {
  const auto particle_ids = mesh_->particles_scalar_data("ids");
  const auto particle_coordinates = mesh_->particle_coordinates();
  if (particle_ids.size() != particle_coordinates.size())
    throw std::runtime_error(
        "Particle id and coordinate arrays differ while binding prescribed "
        "pressures");

  prescribed_pressure_particle_to_sample_.clear();
  prescribed_pressure_particle_to_sample_.reserve(particle_ids.size());

  std::size_t missing_particles = 0;
  mpm::Index first_missing_id = std::numeric_limits<mpm::Index>::max();
  for (std::size_t i = 0; i < particle_ids.size(); ++i) {
    const auto particle_id = static_cast<mpm::Index>(particle_ids[i]);
    const Eigen::Matrix<double, Tdim, 1> initial_coordinates =
        particle_coordinates[i].template head<Tdim>();
    std::size_t sample_index = std::numeric_limits<std::size_t>::max();
    if (!this->initial_pressure_sample_index(
            particle_id, initial_coordinates, &sample_index)) {
      if (missing_particles == 0) first_missing_id = particle_id;
      ++missing_particles;
      continue;
    }

    const auto inserted = prescribed_pressure_particle_to_sample_.emplace(
        particle_id, sample_index);
    if (!inserted.second)
      throw std::runtime_error(
          "Duplicate dynamic particle id while binding prescribed pressures: " +
          std::to_string(particle_id));
  }

  if (missing_particles != 0)
    throw std::runtime_error(
        "Unable to bind " + std::to_string(missing_particles) +
        " particles to the prescribed pressure database; first missing id=" +
        std::to_string(first_missing_id));

  console_->info(
      "Bound [{}] particles to fixed prescribed-pressure samples using [{}]",
      prescribed_pressure_particle_to_sample_.size(),
      prescribed_pressure_mapping_);
}

//! Resolve a database sample once, using the dynamic-stage initial state
template <unsigned Tdim>
bool mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::
    initial_pressure_sample_index(
        mpm::Index particle_id,
        const Eigen::Matrix<double, Tdim, 1>& particle_coordinates,
        std::size_t* sample_index) const {
  if (prescribed_pressure_mapping_ == "id" ||
      prescribed_pressure_mapping_ == "particle_id") {
    const auto iter = prescribed_pressure_id_to_sample_.find(particle_id);
    if (iter != prescribed_pressure_id_to_sample_.end()) {
      *sample_index = iter->second;
      return true;
    }
    return false;
  }

  std::size_t best_index = std::numeric_limits<std::size_t>::max();
  double best_distance = std::numeric_limits<double>::max();

  auto check_sample = [&](std::size_t sample_index) {
    const auto& sample = prescribed_pressure_samples_[sample_index];
    if (!sample.has_coordinates) return;
    const double distance =
        (sample.coordinates - particle_coordinates).squaredNorm();
    if (distance < best_distance) {
      best_distance = distance;
      best_index = sample_index;
    }
  };

  if (prescribed_pressure_has_spatial_index_) {
    Eigen::Matrix<double, Tdim, 1> bin_coordinates =
        prescribed_pressure_min_coordinates_;
    for (int radius = 0; radius <= 3 && best_index == std::numeric_limits<std::size_t>::max();
         ++radius) {
      for (int ix = -radius; ix <= radius; ++ix) {
        for (int iy = -radius; iy <= radius; ++iy) {
          for (int iz = -radius; iz <= radius; ++iz) {
            if (Tdim < 3 && iz != 0) continue;
            bin_coordinates = particle_coordinates;
            bin_coordinates[0] += ix * prescribed_pressure_coordinate_bin_size_;
            if (Tdim > 1)
              bin_coordinates[1] += iy * prescribed_pressure_coordinate_bin_size_;
            if (Tdim > 2)
              bin_coordinates[2] += iz * prescribed_pressure_coordinate_bin_size_;
            const auto key = this->pressure_bin_key(
                bin_coordinates, prescribed_pressure_min_coordinates_,
                std::max(prescribed_pressure_coordinate_bin_size_, 1.0e-12));
            const auto bin_iter =
                prescribed_pressure_coordinate_bins_.find(key);
            if (bin_iter == prescribed_pressure_coordinate_bins_.end())
              continue;
            for (const auto sample_index : bin_iter->second)
              check_sample(sample_index);
          }
        }
      }
    }
  }

  if (best_index == std::numeric_limits<std::size_t>::max()) {
    for (std::size_t i = 0; i < prescribed_pressure_samples_.size(); ++i)
      check_sample(i);
  }

  if (best_index == std::numeric_limits<std::size_t>::max()) return false;

  *sample_index = best_index;
  return true;
}

//! Lookup the immutable database binding by dynamic particle id
template <unsigned Tdim>
bool mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::pressure_sample_index(
    const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle,
    std::size_t* sample_index) const {
  const auto iter =
      prescribed_pressure_particle_to_sample_.find(particle->id());
  if (iter == prescribed_pressure_particle_to_sample_.end()) return false;
  *sample_index = iter->second;
  return true;
}

//! Lookup pressure in a frame through the immutable particle binding
template <unsigned Tdim>
bool mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::pressure_from_frame(
    const PressureFrame& frame,
    const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle,
    double* liquid_pressure, double* gas_pressure) const {
  std::size_t sample_index = std::numeric_limits<std::size_t>::max();
  if (!this->pressure_sample_index(particle, &sample_index)) return false;
  if (sample_index >= frame.liquid_pressures.size() ||
      sample_index >= frame.gas_pressures.size())
    return false;

  *liquid_pressure = frame.liquid_pressures[sample_index];
  *gas_pressure = frame.gas_pressures[sample_index];
  return true;
}

//! Assign externally prescribed phase pressures for the current time
template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::
    assign_prescribed_phase_pressures() {
  // Database frame numbers belong to the source simulation. Convert the
  // actual analysis time to that source grid; solver step numbers are not
  // interchangeable when the two simulations use different dt values.
  double source_step_position = current_time_ / prescribed_pressure_source_dt_;
  if (prescribed_pressure_legacy_step_offset_)
    source_step_position -= 1.0;
  source_step_position = std::max(0.0, source_step_position);
  const double nearest_source_step = std::round(source_step_position);
  if (std::abs(source_step_position - nearest_source_step) <=
      1.0e-10 * std::max(1.0, std::abs(source_step_position)))
    source_step_position = nearest_source_step;
  const double maximum_source_step =
      static_cast<double>(prescribed_pressure_max_step_);
  if (source_step_position >
      maximum_source_step +
          1.0e-10 * std::max(1.0, std::abs(maximum_source_step)))
    throw std::runtime_error(
        "Prescribed pressure history is shorter than the requested analysis "
        "time (source step " +
        std::to_string(source_step_position) + " > max_step " +
        std::to_string(prescribed_pressure_max_step_) + ")");
  source_step_position =
      std::min(source_step_position, maximum_source_step);

  const double source_frame_position =
      source_step_position /
      static_cast<double>(prescribed_pressure_step_interval_);
  const auto frame_id =
      static_cast<mpm::Index>(std::floor(source_frame_position));
  const mpm::Index lower_step =
      frame_id * prescribed_pressure_step_interval_;

  mpm::Index upper_step =
      std::min(lower_step + prescribed_pressure_step_interval_,
               prescribed_pressure_max_step_);

  const double alpha =
      (upper_step > lower_step)
          ? std::min(
                1.0,
                std::max(0.0, (source_step_position - lower_step) /
                                  static_cast<double>(upper_step - lower_step)))
          : 0.0;

  const auto& lower_frame = this->prescribed_pressure_frame(lower_step);
  if (alpha <= 0.0) {
    mesh_->iterate_over_particles([&](const auto& particle) {
      double liquid_pressure = 0.0;
      double gas_pressure = 0.0;
      if (!this->pressure_from_frame(lower_frame, particle, &liquid_pressure,
                                     &gas_pressure))
        return;
      particle->assign_prescribed_phase_pressures(liquid_pressure,
                                                  gas_pressure);
    });
    return;
  }

  const auto& upper_frame = this->prescribed_pressure_frame(upper_step);
  mesh_->iterate_over_particles([&](const auto& particle) {
    double lower_liquid = 0.0;
    double lower_gas = 0.0;
    double upper_liquid = 0.0;
    double upper_gas = 0.0;

    const bool has_lower =
        this->pressure_from_frame(lower_frame, particle, &lower_liquid,
                                  &lower_gas);
    const bool has_upper =
        this->pressure_from_frame(upper_frame, particle, &upper_liquid,
                                  &upper_gas);

    if (!has_lower || !has_upper) return;

    const double liquid_pressure =
        (1.0 - alpha) * lower_liquid + alpha * upper_liquid;
    const double gas_pressure = (1.0 - alpha) * lower_gas + alpha * upper_gas;
    particle->assign_prescribed_phase_pressures(liquid_pressure, gas_pressure);
  });
}

//! MPM Explicit two-phase pressure smoothing
template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::pressure_smoothing(unsigned phase) {

  const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
  const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
  const unsigned mixture = mpm::ParticlePhase::Mixture;
  const unsigned pore_gas = mpm::ParticlePhase::Gas; 

if (phase == mpm::ParticlePhase::Solid) {
    // Assign pressure to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_pressure_to_nodes,
                  std::placeholders::_1, current_time_));
  } else if (phase == mpm::ParticlePhase::Liquid) {
    // // Assign heat capacity and heat to nodes
    // mesh_->iterate_over_particles(
    //     std::bind(&mpm::ParticleBase<Tdim>::map_pore_pressure_to_nodes,
    //               std::placeholders::_1, current_time_));

    // Assign heat capacity and heat to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_pressure_to_nodes,
                  std::placeholders::_1));

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

    // // Apply nodal pressure constraints
    // mesh_->iterate_over_nodes_predicate(
    //   std::bind(&mpm::NodeBase<Tdim>::apply_pressure_constraints,
    //             std::placeholders::_1, pore_liquid, current_time_),
    //   std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
  
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
    if (phase == mpm::ParticlePhase::Liquid) {
      mesh_->template nodal_halo_exchange<double, 1>(
          std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1,
                    pore_gas),
          std::bind(&mpm::NodeBase<Tdim>::assign_pressure,
                    std::placeholders::_1, pore_gas, std::placeholders::_2));
    }
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

    // mesh_->apply_particle_pore_pressure_constraints(current_time_); 

  }
}

// Compute time step size
template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseLag<Tdim>::compute_critical_timestep_size(double dt) {
  const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
  const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
  const unsigned pore_gas = mpm::ParticlePhase::Gas;
  const unsigned mixture = mpm::ParticlePhase::Mixture;

  // cell minimum size
  auto mesh_props = io_->json_object("mesh");
  // Get Mesh reader from JSON object
  double cellsize_min = mesh_props.at("cellsize_min").template get<double>();
  // Solid Material parameters 
  auto materials =  materials_.at(soil_skeleton);
  double porosity = materials->template property<double>(std::string("porosity"));
  double critical_timestep_modulus =
      materials->template property_or<double>(
          std::string("critical_timestep_modulus"), 0.);
  if (critical_timestep_modulus <= 0.)
    critical_timestep_modulus =
        materials->template property<double>(std::string("youngs_modulus"));
  if (!std::isfinite(critical_timestep_modulus) ||
      critical_timestep_modulus <= 0.)
    throw std::runtime_error("Critical timestep modulus must be finite and positive");
  double density = materials->template property<double>(std::string("density"));
  const double minimum_support_fraction =
      analysis_.value("minimum_nodal_support_fraction", 0.0);
  if (!std::isfinite(minimum_support_fraction) ||
      minimum_support_fraction < 0. || minimum_support_fraction > 0.05)
    throw std::runtime_error(
        "minimum_nodal_support_fraction must lie in [0, 0.05]");
  const double minimum_nodal_density =
      minimum_support_fraction * (1. - porosity) * density;
  mesh_->iterate_over_nodes(std::bind(
      &mpm::NodeBase<Tdim>::assign_minimum_nodal_density,
      std::placeholders::_1, minimum_nodal_density));
  if (minimum_nodal_density > 0.)
    console_->info(
        "Minimum nodal support density is {} kg/m^3 (fraction={})",
        minimum_nodal_density, minimum_support_fraction);
  double specific_heat = materials->template property<double>(std::string("specific_heat"));
  double thermal_conductivity = materials->template property<double>(std::string("thermal_conductivity"));
  // Compute timestep fpor one phase MPM                              
  double critical_dt = cellsize_min / std::pow(critical_timestep_modulus/density/(1 - porosity), 0.5);
  console_->info("Critical time step size is {} s", critical_dt);
  // Liquid Material parameters 
  auto liquid_materials =  materials_.at(pore_liquid);
  double liquid_density = liquid_materials->template property<double>(std::string("density"));
  double liquid_specific_heat = liquid_materials->template property<double>(std::string("liquid_specific_heat"));
  double liquid_thermal_conductivity = liquid_materials->template property<double>(std::string("liquid_thermal_conductivity"));

  // Compute timestep for momentum eqaution 
  double density_mixture1 = (1 - porosity) * density;
  double density_mixture2 = (1 - porosity) * density + porosity * liquid_density;
  double critical_dt11 = cellsize_min / std::pow(critical_timestep_modulus/density_mixture1, 0.5);
  double critical_dt12 = cellsize_min / std::pow(critical_timestep_modulus/density_mixture2, 0.5);
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
