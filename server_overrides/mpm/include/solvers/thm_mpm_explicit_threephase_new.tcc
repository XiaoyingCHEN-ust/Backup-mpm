//! Constructor
template <unsigned Tdim>
mpm::ThermoMPMExplicitThreePhaseNew<Tdim>::ThermoMPMExplicitThreePhaseNew(
    const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("ThermoMPMExplicitThreePhaseNew");
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////                                                                ////////
////////             THM-MPM Explicit ThreePhase New Solver             ////////
////////                                                                ////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//! Thermo-hydro-mechncial MPM Explicit solver
template <unsigned Tdim>
bool mpm::ThermoMPMExplicitThreePhaseNew<Tdim>::solve() {
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

  // Pressure integration.  The default preserves the historical explicit
  // update; "semi_implicit" activates a global backward-Euler liquid/gas
  // pressure solve while retaining explicit solid mechanics.
  if (analysis_.find("pressure_integration") != analysis_.end()) {
    const auto pressure_integration =
        analysis_["pressure_integration"].template get<std::string>();
    if (pressure_integration == "semi_implicit") {
      semi_implicit_pressure_ = true;
    } else if (pressure_integration != "explicit") {
      throw std::runtime_error(
          "pressure_integration must be 'explicit' or 'semi_implicit'");
    }
  }
  if (analysis_.find("semi_implicit_pressure") != analysis_.end()) {
    const auto& pressure_options = analysis_["semi_implicit_pressure"];
    if (pressure_options.find("boundary_penalty") != pressure_options.end())
      pressure_boundary_penalty_ =
          pressure_options["boundary_penalty"].template get<double>();
    if (pressure_options.find("log_solver") != pressure_options.end())
      log_pressure_solver_ =
          pressure_options["log_solver"].template get<bool>();
  }
  if (!(pressure_boundary_penalty_ > 0.0))
    throw std::runtime_error(
        "semi_implicit_pressure.boundary_penalty must be positive");
  console_->info("Three-phase pressure integration: {}",
                 semi_implicit_pressure_ ? "semi_implicit" : "explicit");

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
  mpm::Index start_step = 0;
  if (resume) {
    if (!this->checkpoint_resume())
      throw std::runtime_error("Checkpoint resume failed");
    this->current_time_ = analysis_["resume"]["current_time"].template get<double>();
    start_step = this->step_;
    std::cout << "current_time" << this->current_time_ << "\n";
  }

  solver_begin = std::chrono::steady_clock::now();

  this->compute_critical_timestep_size(dt_);

////////////////////////////////////////////////////////////////////////////////
////////////////                  MAIN LOOP               //////////////////////
////////////////////////////////////////////////////////////////////////////////

  // Main loop
  for (step_ = start_step; step_ <= nsteps_; ++step_) {

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
    mesh_->compute_free_surface(free_surface_particle_, volume_tolerance_);

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

    // Update liquid and gas pressures.  The semi-implicit path solves the
    // coupled spatial pressure problem; the explicit path is unchanged.
    if (semi_implicit_pressure_) {
      if (!this->solve_semi_implicit_pressure(dt_))
        throw std::runtime_error("Semi-implicit pressure solve failed");
    } else {
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_pore_pressure,
                    std::placeholders::_1, dt_));
    }

    mesh_->apply_particle_pore_pressure_constraints(current_time_);

    // Pressure smoothing
    if (pressure_smoothing_ & (step_ % 100 == 0)) { this->pressure_smoothing(pore_liquid);}

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
////////////////                LOOP ENDING               //////////////////////
////////////////////////////////////////////////////////////////////////////////

//! Coupled liquid/gas backward-Euler pressure solve
template <unsigned Tdim>
bool mpm::ThermoMPMExplicitThreePhaseNew<Tdim>::solve_semi_implicit_pressure(
    double dt) {
  using PressureState = mpm::ThreePhasePressureState<Tdim>;
  using Triplet = Eigen::Triplet<double>;
  if (!(dt > 0.0)) return false;

  try {
    const unsigned mesh_active_dof = mesh_->assign_active_node_id();
    if (mesh_active_dof == 0)
      throw std::runtime_error("No active pressure nodes");

    std::vector<PressureState> states;
    states.reserve(mesh_->nparticles());
    std::mutex states_mutex;
    mesh_->iterate_over_particles(
        [&](const std::shared_ptr<mpm::ParticleBase<Tdim>>& base_particle) {
          const auto particle = std::dynamic_pointer_cast<
              mpm::ThreePhasePressureParticle<Tdim>>(base_particle);
          if (!particle) return;
          auto state = particle->semi_implicit_pressure_state();
          if (state.active_node_ids.empty()) return;
          std::lock_guard<std::mutex> lock(states_mutex);
          states.emplace_back(std::move(state));
        });
    if (states.empty())
      throw std::runtime_error("No three-phase particles in pressure solve");

    const bool fixed_gas_pressure = states.front().fixed_gas_pressure;
    for (const auto& state : states) {
      if (state.fixed_gas_pressure != fixed_gas_pressure)
        throw std::runtime_error(
            "Mixed fixed and variable gas-pressure materials are unsupported");
      if (state.active_node_ids.size() !=
              static_cast<std::size_t>(state.shapefn.size()) ||
          state.dn_dx.rows() != state.shapefn.size())
        throw std::runtime_error("Invalid particle pressure interpolation data");
      for (const auto node_id : state.active_node_ids)
        if (node_id >= mesh_active_dof)
          throw std::runtime_error("Inactive node in pressure interpolation");
    }

    // A GIMP cell can mark support nodes as mechanically active even when a
    // node has zero pressure interpolation weight.  Build a compact pressure
    // DOF map from nodes that contribute through either N or grad(N), instead
    // of treating every mechanically active node as a pressure unknown.
    const Index invalid_dof = std::numeric_limits<Index>::max();
    std::vector<bool> pressure_support(mesh_active_dof, false);
    for (const auto& state : states) {
      for (unsigned i = 0; i < state.active_node_ids.size(); ++i) {
        if (std::abs(state.shapefn(i)) > 1.0e-14 ||
            state.dn_dx.row(i).squaredNorm() > 1.0e-28)
          pressure_support[state.active_node_ids[i]] = true;
      }
    }
    std::vector<Index> pressure_dof(mesh_active_dof, invalid_dof);
    unsigned pressure_active_dof = 0;
    for (unsigned node_id = 0; node_id < mesh_active_dof; ++node_id) {
      if (pressure_support[node_id])
        pressure_dof[node_id] = pressure_active_dof++;
    }
    if (pressure_active_dof == 0)
      throw std::runtime_error("No supported pressure degrees of freedom");

    // Project old particle pressures to the compact pressure nodes for the
    // storage term.  Absolute interpolation weights are robust to generalised
    // interpolation functions.  A gradient-only pressure node has no storage
    // contribution, so initialise it with the volume-weighted particle mean.
    Eigen::VectorXd nodal_weight =
        Eigen::VectorXd::Zero(pressure_active_dof);
    Eigen::VectorXd old_liquid_pressure =
        Eigen::VectorXd::Zero(pressure_active_dof);
    Eigen::VectorXd old_gas_pressure =
        Eigen::VectorXd::Zero(pressure_active_dof);
    double particle_volume = 0.0;
    double mean_liquid_pressure = 0.0;
    double mean_gas_pressure = 0.0;
    for (const auto& state : states) {
      particle_volume += state.volume;
      mean_liquid_pressure += state.volume * state.liquid_pressure;
      mean_gas_pressure += state.volume * state.gas_pressure;
      for (unsigned i = 0; i < state.active_node_ids.size(); ++i) {
        const auto dof = pressure_dof[state.active_node_ids[i]];
        if (dof == invalid_dof) continue;
        const double weight = state.volume * std::abs(state.shapefn(i));
        nodal_weight(dof) += weight;
        old_liquid_pressure(dof) += weight * state.liquid_pressure;
        old_gas_pressure(dof) += weight * state.gas_pressure;
      }
    }
    if (!(particle_volume > 0.0))
      throw std::runtime_error("Non-positive particle volume in pressure solve");
    mean_liquid_pressure /= particle_volume;
    mean_gas_pressure /= particle_volume;
    for (unsigned i = 0; i < pressure_active_dof; ++i) {
      if (nodal_weight(i) > 0.0) {
        old_liquid_pressure(i) /= nodal_weight(i);
        old_gas_pressure(i) /= nodal_weight(i);
      } else {
        old_liquid_pressure(i) = mean_liquid_pressure;
        old_gas_pressure(i) = mean_gas_pressure;
      }
    }

    const unsigned system_size =
        fixed_gas_pressure ? pressure_active_dof : 2 * pressure_active_dof;
    Eigen::VectorXd right = Eigen::VectorXd::Zero(system_size);
    std::vector<Triplet> coefficients;
    const unsigned nodes_per_particle =
        states.front().active_node_ids.size();
    coefficients.reserve(states.size() * nodes_per_particle *
                         nodes_per_particle * (fixed_gas_pressure ? 1 : 4));

    for (const auto& state : states) {
      const double c_ww = state.porosity *
          (state.saturation_pressure_tangent +
           state.liquid_saturation * state.liquid_compressibility);
      const double c_wg =
          -state.porosity * state.saturation_pressure_tangent;
      const double c_gw = c_wg;
      const double c_gg = state.porosity *
          (state.gas_compressibility +
           state.saturation_pressure_tangent);
      const double liquid_conductivity =
          std::max(state.liquid_conductivity, 0.0);
      const double gas_conductivity =
          std::max(state.gas_conductivity, 0.0);

      double gradient_scale = 0.0;
      for (int i = 0; i < state.dn_dx.rows(); ++i)
        gradient_scale =
            std::max(gradient_scale, state.dn_dx.row(i).squaredNorm());
      const double liquid_penalty = pressure_boundary_penalty_ *
          std::max({c_ww / dt, liquid_conductivity * gradient_scale,
                    1.0e-20});
      const double gas_penalty = pressure_boundary_penalty_ *
          std::max({c_gg / dt, gas_conductivity * gradient_scale,
                    1.0e-20});

      for (unsigned i = 0; i < state.active_node_ids.size(); ++i) {
        const auto row_w = pressure_dof[state.active_node_ids[i]];
        if (row_w == invalid_dof) continue;
        const auto row_g = row_w + pressure_active_dof;
        const double shape_i = state.shapefn(i);
        const auto gradient_i = state.dn_dx.row(i).transpose();

        right(row_w) +=
            state.volume * shape_i * state.liquid_mass_source;
        right(row_w) += state.volume * liquid_conductivity *
                        gradient_i.dot(state.liquid_density * gravity_);
        if (!fixed_gas_pressure) {
          right(row_g) +=
              state.volume * shape_i * state.gas_mass_source;
          right(row_g) += state.volume * gas_conductivity *
                          gradient_i.dot(state.gas_density * gravity_);
        }

        // Row-sum lump the pressure-storage matrix.  A consistent mass
        // matrix is non-monotone at the sharp dry/wet discontinuity and
        // produced saturated bands disconnected from the imposed surface.
        const double lumped_mass =
            state.volume * std::max(shape_i, 0.0);
        coefficients.emplace_back(
            row_w, row_w, c_ww * lumped_mass / dt);
        right(row_w) +=
            c_ww * lumped_mass / dt * old_liquid_pressure(row_w);
        if (!fixed_gas_pressure) {
          coefficients.emplace_back(
              row_w, row_g, c_wg * lumped_mass / dt);
          coefficients.emplace_back(
              row_g, row_w, c_gw * lumped_mass / dt);
          coefficients.emplace_back(
              row_g, row_g, c_gg * lumped_mass / dt);
          right(row_w) +=
              c_wg * lumped_mass / dt * old_gas_pressure(row_w);
          right(row_g) +=
              c_gw * lumped_mass / dt * old_liquid_pressure(row_w) +
              c_gg * lumped_mass / dt * old_gas_pressure(row_w);
        }

        for (unsigned j = 0; j < state.active_node_ids.size(); ++j) {
          const auto col_w = pressure_dof[state.active_node_ids[j]];
          if (col_w == invalid_dof) continue;
          const auto col_g = col_w + pressure_active_dof;
          const double diffusion = state.volume *
              gradient_i.dot(state.dn_dx.row(j).transpose());

          coefficients.emplace_back(
              row_w, col_w, liquid_conductivity * diffusion);

          if (!fixed_gas_pressure) {
            coefficients.emplace_back(
                row_g, col_g, gas_conductivity * diffusion);
          }
        }

        if (state.pressure_boundary) {
          coefficients.emplace_back(
              row_w, row_w, liquid_penalty * lumped_mass);
          right(row_w) += liquid_penalty * lumped_mass *
                          state.boundary_liquid_pressure;
          if (!fixed_gas_pressure) {
            coefficients.emplace_back(
                row_g, row_g, gas_penalty * lumped_mass);
            right(row_g) += gas_penalty * lumped_mass *
                            state.boundary_gas_pressure;
          }
        }
      }
    }

    Eigen::SparseMatrix<double> pressure_matrix(system_size, system_size);
    pressure_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
    pressure_matrix.makeCompressed();

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> pressure_solver;
    pressure_solver.compute(pressure_matrix);
    if (pressure_solver.info() != Eigen::Success)
      throw std::runtime_error("Pressure matrix factorisation failed");
    const Eigen::VectorXd pressure = pressure_solver.solve(right);
    if (pressure_solver.info() != Eigen::Success || !pressure.allFinite())
      throw std::runtime_error("Pressure linear solve failed");

    const double relative_residual =
        (pressure_matrix * pressure - right).norm() /
        std::max(right.norm(), 1.0);
    if (!std::isfinite(relative_residual) || relative_residual > 1.0e-7)
      throw std::runtime_error("Pressure solve residual exceeded tolerance");

    // Transfer only the pressure increment back to particles.  Replacing the
    // particle pressure with the projected absolute nodal pressure would add
    // one artificial smoothing operation per time step and move a sharp
    // wetting front faster when a smaller time step is used.
    const Eigen::VectorXd compact_liquid_pressure_increment =
        pressure.head(pressure_active_dof) - old_liquid_pressure;
    Eigen::VectorXd compact_gas_pressure_increment =
        Eigen::VectorXd::Zero(pressure_active_dof);
    if (!fixed_gas_pressure)
      compact_gas_pressure_increment =
          pressure.tail(pressure_active_dof) - old_gas_pressure;
    Eigen::VectorXd liquid_pressure_increment =
        Eigen::VectorXd::Zero(mesh_active_dof);
    Eigen::VectorXd gas_pressure_increment =
        Eigen::VectorXd::Zero(mesh_active_dof);
    Eigen::VectorXd nodal_liquid_pressure =
        Eigen::VectorXd::Constant(mesh_active_dof, mean_liquid_pressure);
    Eigen::VectorXd nodal_gas_pressure =
        Eigen::VectorXd::Constant(mesh_active_dof, mean_gas_pressure);
    for (unsigned node_id = 0; node_id < mesh_active_dof; ++node_id) {
      const auto dof = pressure_dof[node_id];
      if (dof == invalid_dof) continue;
      liquid_pressure_increment(node_id) =
          compact_liquid_pressure_increment(dof);
      gas_pressure_increment(node_id) =
          compact_gas_pressure_increment(dof);
      nodal_liquid_pressure(node_id) = pressure(dof);
      nodal_gas_pressure(node_id) =
          fixed_gas_pressure ? old_gas_pressure(dof)
                             : pressure(dof + pressure_active_dof);
    }
    std::atomic<bool> update_status{true};
    std::atomic<std::size_t> phase_limited_particles{0};
    std::atomic<std::size_t> capillary_limited_particles{0};
    mesh_->iterate_over_particles(
        [&](const std::shared_ptr<mpm::ParticleBase<Tdim>>& base_particle) {
          const auto particle = std::dynamic_pointer_cast<
              mpm::ThreePhasePressureParticle<Tdim>>(base_particle);
          if (!particle) return;
          const int limiter_flags = particle->update_semi_implicit_pressure(
              liquid_pressure_increment, gas_pressure_increment,
              nodal_liquid_pressure, nodal_gas_pressure, dt);
          if (limiter_flags < 0) {
            update_status.store(false);
            return;
          }
          if (limiter_flags & 1) ++phase_limited_particles;
          if (limiter_flags & 2) ++capillary_limited_particles;
        });
    if (!update_status.load())
      throw std::runtime_error("Particle pressure update failed");

    if (log_pressure_solver_ && step_ % output_steps_ == 0)
      console_->info(
          "Semi-implicit pressure: dofs={}, mesh_active_nodes={}, "
          "fixed_gas={}, phase_limited={}, capillary_limited={}, residual={}",
          system_size, mesh_active_dof, fixed_gas_pressure,
          phase_limited_particles.load(), capillary_limited_particles.load(),
          relative_residual);
    return true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__,
                    __func__, exception.what());
    return false;
  }
}

//! MPM Explicit two-phase pressure smoothing
template <unsigned Tdim>
void mpm::ThermoMPMExplicitThreePhaseNew<Tdim>::pressure_smoothing(unsigned phase) {

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
void mpm::ThermoMPMExplicitThreePhaseNew<Tdim>::compute_critical_timestep_size(double dt) {
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
