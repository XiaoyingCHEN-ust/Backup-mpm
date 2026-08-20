//==============================================================================
// CONSTRUCT AND DESTRUCT A NODE
//==============================================================================

// Constructor with id, coordinates and dof
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
mpm::Node<Tdim, Tdof, Tnphases>::Node(
    Index id, const Eigen::Matrix<double, Tdim, 1>& coord)
    : NodeBase<Tdim>(id, coord) {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
  id_ = id;
  coordinates_ = coord;
  dof_ = Tdof;

  // Logger
  std::string logger =
      "node" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);

  // Clear any velocity constraints
  velocity_constraints_.clear();
  temperature_constraints_.clear();
  concentrated_force_.setZero();
  this->initialise();
}

//==============================================================================
// INITIALISE NODE
//==============================================================================

// Initialise nodal properties
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::initialise() noexcept {
  mass_.setZero();
  volume_.setZero();
  density_.setZero();   
  mass_source_.setZero(); 
  external_force_.setZero();
  internal_force_.setZero();
  reaction_force_.setZero();
  reaction_dir_.setZero();
  pressure_.setZero();
  velocity_.setZero();
  momentum_.setZero();
  acceleration_.setZero();
  status_ = false;
  contact_ = false;
  rigid_velocity_.setZero();
  material_ids_.clear();
  contact_normal_.setZero();
  contact_tangential_.setZero(); 
  rigid_acceleration_.setZero(); 
  rigid_mass_ = 0;
  rigid_mass2_ = 0;  
  mass_contact_normal_.setZero();
  mass_rigid_velocity_.setZero();
  mass_pressure_.setZero();
  gas_saturation_ = 0;
  pressure_.setZero();
  pressure_acceleration_.setZero();
  hydraulic_conduction_.setZero();
  K_coeff_.setZero();
  volumetric_strain_.setZero();
  mean_length_ = 0.; 
  scaler_.setZero();

  // Specific variables for two phase
  free_surface_ = false;
  phase_kinematic_boundary_ = false;
  set_velocity_constraints_ = false;
  drag_force_coefficient_.setZero();
  drag_force_.setZero();
  force_total_inter_.setZero();
  force_fluid_inter_.setZero();
  velocity_inter_.setZero();
  acceleration_inter_.setZero();
  pore_pressure_increment_ = 0.;
  force_cor_.setZero();

  // Specific variables thermal part
  temperature_.setZero();
  temperature_acceleration_.setZero();
  heat_capacity_.setZero();
  heat_.setZero();
  heat_conduction_.setZero();
  hydraulic_conduction_.setZero();
  mass_convection_.setZero();  
  heat_convection_.setZero();
  heat_source_.setZero();
  latent_heat_.setZero();
  convective_heat_flux_.setZero();
  plastic_work_.setZero();

  // Link data with NAME
  this->nodal_scalar_property_ = {
    {"temperature",         [&]() {return this->temperature_(0);}},
    {"temperature_acc",     [&]() {return this->temperature_acceleration_(0);}},
    {"mixture_mass",        [&]() {return this->mass_(0);}},
    {"free_surface",        [&]() {return this->free_surface();}}
  };
  this->nodal_vector_property_ = {
    {"solid_velocity",      [&]() {return this->velocity_.col(0);}},
    {"solid_acceleration",  [&]() {return this->acceleration_.col(0);}},
    {"mix_ext_force",       [&]() {return this->external_force_.col(0);}},
    {"mix_int_force",       [&]() {return this->internal_force_.col(0);}},
    {"reaction_force",      [&]() {return this->reaction_force_;}}
  };

  // Liquid phase
  if (Tnphases > 1) {
    this->nodal_scalar_property_.insert({
      {"liquid_mass",         [&]() {return this->mass_(1);}},
      {"liquid_pressure",     [&]() {return this->pressure_(1);}},
      {"liquid_pressure_acc", [&]() {return this->pressure_acceleration_(1);}},
      {"liquid_mass_source",  [&]() {return this->mass_source_(1);}}, 
      {"liquid_hydraulic_conduction",    [&]() {return this->hydraulic_conduction_(1);}}
    });
    this->nodal_vector_property_.insert({
      {"liquid_velocity",     [&]() {return this->velocity_.col(1);}},
      {"liquid_acceleration", [&]() {return this->acceleration_.col(1);}},
      {"liquid_ext_force",    [&]() {return this->external_force_.col(1);}},
      {"liquid_int_force",    [&]() {return this->internal_force_.col(1);}},
      {"drag_force_liquid",   [&]() {return this->drag_force_liquid_;}}
    });
  }

  // Gas phase
  if (Tnphases > 2) {
    this->nodal_scalar_property_.insert({
      {"gas_mass",            [&]() {return this->mass_(2);}},
      {"gas_pressure",        [&]() {return this->pressure_(2);}},
      {"gas_pressure_acc",    [&]() {return this->pressure_acceleration_(2);}},
      {"gas_mass_source",     [&]() {return this->mass_source_(2);}},
      {"gas_hydraulic_conduction",    [&]() {return this->hydraulic_conduction_(2);}}
    });
    this->nodal_vector_property_.insert({
      {"gas_velocity",        [&]() {return this->velocity_.col(2);}},
      {"gas_acceleration",    [&]() {return this->acceleration_.col(2);}},
      {"gas_ext_force",       [&]() {return this->external_force_.col(2);}},
      {"gas_int_force",       [&]() {return this->internal_force_.col(2);}},
      {"drag_force_gas",      [&]() {return this->drag_force_gas_;}}
    });
  }

}

// Add material id from material points to material_ids_
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::append_material_id(unsigned id) {
  std::lock_guard<std::mutex> guard(node_mutex_);
  material_ids_.emplace(id);
}

// Assign MPI rank to node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::mpi_rank(unsigned rank) {
  std::lock_guard<std::mutex> guard(node_mutex_);
  auto status = this->mpi_ranks_.insert(rank);
  return status.second;
}

//==============================================================================
// UPDATE NODAL INFORMATION FROM PARTICLES
//==============================================================================

// MPI// Update mass at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_mass_momentum(
                    bool update, unsigned phase, double mass,
                    const Eigen::Matrix<double, Tdim, 1>& momentum) noexcept {
  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign mass
  std::lock_guard<std::mutex> guard(node_mutex_);
  mass_(phase) = mass_(phase) * factor + mass;
  momentum_.col(phase) = momentum_.col(phase) * factor + momentum;
}

// MPI// Update mass at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_mass(bool update, unsigned phase,
                                                  double mass) noexcept {
  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign mass
  std::lock_guard<std::mutex> guard(node_mutex_);
  mass_(phase) = mass_(phase) * factor + mass;
}

// MPI// Assign nodal momentum
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_momentum(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& momentum) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign momentum
  std::lock_guard<std::mutex> guard(node_mutex_);
  momentum_.col(phase) = momentum_.col(phase) * factor + momentum;
}

// MPI// Update external force (body force / traction force)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_external_force(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& force) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign external force
  std::lock_guard<std::mutex> guard(node_mutex_);
  external_force_.col(phase) = external_force_.col(phase) * factor + force;
}

// MPI// Update internal force
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_internal_force(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& force) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign internal force
  std::lock_guard<std::mutex> guard(node_mutex_);
  internal_force_.col(phase) = internal_force_.col(phase) * factor + force;
}

///! MPI// Assign/update heat capacity at nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_heat_capacity(
    bool update, unsigned phase, double heat_capacity) noexcept {
  
  // Assert
  assert(phase < Tnphases);
  
  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;
  
  // Update/assign mass
  std::lock_guard<std::mutex> guard(node_mutex_);
  heat_capacity_(phase) = heat_capacity_(phase) * factor + heat_capacity;
}

///! MPI// Assign/update heat
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_heat(
    bool update, unsigned phase, double heat) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign heat
  std::lock_guard<std::mutex> guard(node_mutex_);
  heat_(phase) = heat_(phase) * factor + heat;
}

///! MPI// Assign/update heat
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_energy(
    bool update, unsigned phase, double heat_capacity, double heat) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign heat
  std::lock_guard<std::mutex> guard(node_mutex_);
  heat_capacity_(phase) = heat_capacity_(phase) * factor + heat_capacity;  
  heat_(phase) = heat_(phase) * factor + heat;
}

///! Assign/update heat conduction
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_heat_conduction(
    bool update, unsigned phase, const double heat_conduction) noexcept {

  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign heat conduction
  std::lock_guard<std::mutex> guard(node_mutex_);
  heat_conduction_(phase) = heat_conduction_(phase) * factor + heat_conduction;  
}

///! Assign/update hydraulic conduction
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_hydraulic_conduction(
    bool update, unsigned phase, const double hydraulic_conduction) noexcept {

  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign hydraulic conduction
  std::lock_guard<std::mutex> guard(node_mutex_);
  hydraulic_conduction_(phase) = hydraulic_conduction_(phase) * factor + hydraulic_conduction;  
}

///! Assign/update heat convection (only for twophase)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_heat_convection(
    bool update, unsigned phase, double heat_convection) noexcept {

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign heat convection
  std::lock_guard<std::mutex> guard(node_mutex_);
  heat_convection_(phase) = heat_convection_(phase) * factor + heat_convection;
}

///! Assign/update mass convection (only for twophase)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_mass_convection(
    bool update, unsigned phase, double mass_convection) noexcept {

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign mass convection
  std::lock_guard<std::mutex> guard(node_mutex_);
  mass_convection_(phase) = mass_convection_(phase) * factor + mass_convection;
}

///! Assign/update heat convection (only for twophase)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_covective_heat_flux(
    bool update, unsigned phase, double covective_heat_flux) noexcept {

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign heat convection
  std::lock_guard<std::mutex> guard(node_mutex_);
  convective_heat_flux_(phase) = convective_heat_flux_(phase) * factor + covective_heat_flux;
}

///! Assign/update latent heat
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_latent_heat(
    bool update, unsigned phase, double latent_heat) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign heat conduction
  std::lock_guard<std::mutex> guard(node_mutex_);
  latent_heat_(phase) = latent_heat_(phase) * factor + latent_heat;
}

///! Assign/update heat source
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_heat_source(
    bool update, unsigned phase, double heat_source) noexcept {

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign heat convection
  std::lock_guard<std::mutex> guard(node_mutex_);
  heat_source_(phase) = heat_source_(phase) * factor + heat_source;
}

///! Assign/update heat source
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_mass_source(
    bool update, unsigned phase, double mass_source) noexcept {

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign heat convection
  std::lock_guard<std::mutex> guard(node_mutex_);
  mass_source_(phase) = mass_source_(phase) * factor + mass_source;
}

///! Assign/update heat source
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_volumetric_strain(
    bool update, unsigned phase, double volumetric_strain) noexcept {

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign heat convection
  std::lock_guard<std::mutex> guard(node_mutex_);
  volumetric_strain_(phase) = volumetric_strain_(phase) * factor + volumetric_strain;
}

///! Assign/update plastic work
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_plastic_work(
    bool update, unsigned phase, double plastic_work) noexcept {

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign heat convection
  std::lock_guard<std::mutex> guard(node_mutex_);
  plastic_work_(phase) = plastic_work_(phase) * factor + plastic_work;
}

//! Update drag force
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_drag_force_coefficient(
    bool update, unsigned phase, 
    const double drag_force_coefficient) {
  bool status = false;
  try {
    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    // Update/assign drag force coefficient
    std::lock_guard<std::mutex> guard(node_mutex_);
    drag_force_coefficient_(phase) =
        drag_force_coefficient_(phase) * factor + drag_force_coefficient;
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Update pressure at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_mass_pressure(
    unsigned phase, double mass_pressure, double current_time) noexcept {
  // Assert
  assert(phase < Tnphases);

  const double tolerance = 1.E-15;
  // Compute pressure from mass*pressure
  if (mass_(phase) > tolerance) {
    std::lock_guard<std::mutex> guard(node_mutex_);
    pressure_(phase) += mass_pressure / mass_(phase);

    if (pressure_constraints_.find(phase) != pressure_constraints_.end()) {
      const double scalar =
          (pressure_function_.find(phase) != pressure_function_.end())
              ? pressure_function_[phase]->value(current_time)
              : 1.0;
      this->pressure_(phase) = scalar * pressure_constraints_[phase];
    }
  }
  // if(this->free_surface()) this->pressure_.setZero();
}

// Update pressure at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_gas_saturation(unsigned phase, 
    double gas_saturation) noexcept {
  // Assert
  assert(phase < Tnphases);
  const double tolerance = 1.E-15;
  if (mass_(phase) > tolerance) {
    std::lock_guard<std::mutex> guard(node_mutex_);
    gas_saturation_ += gas_saturation / mass_(phase);
  }
}

// Update pressure at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_pressure(bool update, 
    unsigned phase, double mass_pressure) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign pressure
  std::lock_guard<std::mutex> guard(node_mutex_);
  mass_pressure_(phase) = mass_pressure_(phase) * factor + mass_pressure;
}

// Update pressure at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_density(bool update, 
    unsigned phase, double density) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign pressure
  std::lock_guard<std::mutex> guard(node_mutex_);
  density_(phase) =density_(phase) * factor + density;
}

// Update pressure at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_K_coeff(bool update, 
    unsigned phase, double K_coeff) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign pressure
  std::lock_guard<std::mutex> guard(node_mutex_);
  K_coeff_(phase) = K_coeff_(phase) * factor + K_coeff;
}

// Update volume at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_volume(bool update, unsigned phase,
                                    double volume, double mean_length) noexcept {
  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign volume
  std::lock_guard<std::mutex> guard(node_mutex_);
  volume_(phase) = volume_(phase) * factor + volume;
  mean_length_ = mean_length_ * factor + mean_length;

}

// Update smoothed scaler properties at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_scalers(
    bool update, unsigned scaler_id, double scaler) noexcept {
  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Compute pressure from mass*pressure
  std::lock_guard<std::mutex> guard(node_mutex_);
  scaler_(scaler_id) = scaler_(scaler_id) + scaler * factor;
}

//==============================================================================
// COMPUTE NODAL VARIABLES
//==============================================================================

// Compute velocity from momentum
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_velocity(double dt) {
  const double tolerance = 1.E-15;
  if (minimum_nodal_density_ > 0. && volume_(0) > tolerance &&
      mass_(0) / volume_(0) < minimum_nodal_density_) {
    velocity_.setZero();
    acceleration_.setZero();
    this->apply_velocity_constraints();
    return;
  }
  for (unsigned phase = 0; phase < Tnphases; ++phase) {
    if (mass_(phase) > tolerance) {
      velocity_.col(phase) = momentum_.col(phase) / mass_(phase);

      // Check to see if value is below threshold
      for (unsigned i = 0; i < velocity_.rows(); ++i)
        if (std::abs(velocity_.col(phase)(i)) < 1.E-15)
          velocity_.col(phase)(i) = 0.;
    }
  }

  // At a three-phase free-surface node the pore phases share the skeleton
  // velocity. Apply the condition only after every phase velocity has been
  // reconstructed; doing it inside the phase loop can be overwritten by a
  // later momentum-to-velocity calculation.
  if (this->phase_kinematic_boundary_ && Tnphases == 3) {
    velocity_.col(1) = velocity_.col(0);
    velocity_.col(2) = velocity_.col(0);
  }

  // Contact is authoritative over the unconstrained free-surface update. A
  // prescribed velocity constraint is applied below and remains the final
  // boundary condition in case both are present at the same node.
  if (this->contact_) {
    rigid_acceleration_ = (rigid_velocity_ - velocity_) / dt;
    velocity_ = rigid_velocity_;
    acceleration_ = rigid_acceleration_;
  }

  // Apply velocity constraints last; this also sets acceleration to zero.
  this->apply_velocity_constraints();
  // if (this->contact_) velocity_ = rigid_velocity_;   

  // // Debug
  // if (id_ <= 100 && id_ >= 90) {
  //   std::cout << "nodes_.id_: " << id_ << std::endl;
  //   std::cout << "    velocity_.col(liquid): " << velocity_.col(1) << std::endl;
  //   std::cout << "    velocity_.col(soil): " << velocity_.col(0) << std::endl;
  //   //std::cout << "    velocity_.col(gas): " << velocity_(mpm::ParticlePhase::Gas) << std::endl;
  // }

}

// Compute nodal temperature from heat and heat capacity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_temperature(unsigned phase) {
  const double tolerance = 1.E-15;
  if (this->heat_capacity_(phase) > tolerance) {
    this->temperature_(phase) = this->heat_(phase) / this->heat_capacity_(phase);
  }

  if ((this->free_surface()) || (this->contact_)) {
      if (this->current_time_ >= 1.0) {
        this->temperature_(phase) = 10;
        this->temperature_acceleration_(phase) = 0;
      }
  }

  // // Debug
  // if (id_ <= 100 && id_ >= 90) {
  //   std::cout << "nodes_.id_: " << id_ << std::endl;
  //   std::cout << "    temperature_(phase): " << temperature_(phase) << std::endl;
  //   std::cout << "    heat_(phase): " << this->heat_(phase) << std::endl;
  //   std::cout << "    heat_capacity_(phase): " << this->heat_capacity_(phase) << std::endl;
  //   std::cout << "    temperature_acceleration_(phase): " << this->temperature_acceleration_(phase) << std::endl;
  // }

}

// Compute nodal scalers
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_scalers() {
  const double tolerance = 1.E-15;
  if (scaler_(0) > tolerance) {
    this->scaler_(1) /= this->scaler_(0);
    this->scaler_(2) /= this->scaler_(0);
    this->scaler_(3) /= this->scaler_(0);
    this->scaler_(4) /= this->scaler_(0);
    this->scaler_(5) /= this->scaler_(0);
    this->scaler_(6) /= this->scaler_(0);
    this->scaler_(7) /= this->scaler_(0);
  } 
}

// Compute nodal pressure 
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_pressure(unsigned phase) {
  const double tolerance = 1.E-15;
  if (scaler_(0) > tolerance) {
    this->pressure_(phase) = this->mass_pressure_(phase) / this->scaler_(0);
  } 
}

// Compute acceleration and velocity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_acceleration_velocity(
    unsigned phase, double dt) noexcept {
  bool status = false;
  const double tolerance = 1.0E-15;
  if (mass_(phase) > tolerance) {

    // set zero total force for rigid particle influence node
    // get reaction force
    reaction_force_.setZero();

    for (unsigned i = 0; i < Tdim; ++i) {
      if (reaction_dir_(i) == 1) {
        reaction_force_(i) = -(this->external_force_.col(0)(i) +
                               this->internal_force_.col(0)(i));
      }
    }
  
    // acceleration = (unbalaced force / mass)
    this->acceleration_.col(phase) =
        (this->external_force_.col(phase) + this->internal_force_.col(phase) +
         this->reaction_force_) /
        this->mass_(phase);

    // Apply friction constraints
    this->apply_friction_constraints(dt);

    // Velocity += acceleration * dt
    this->velocity_.col(phase) += this->acceleration_.col(phase) * dt;
    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints();

    // Set a threshold
    for (unsigned i = 0; i < Tdim; ++i) {
      if (std::abs(velocity_.col(phase)(i)) < tolerance)
        velocity_.col(phase)(i) = 0.;
      if (std::abs(acceleration_.col(phase)(i)) < tolerance)
        acceleration_.col(phase)(i) = 0.;
    }
    status = true;
  }
  return status;
}

// Compute nodal temperature acceleration and update nodal temperature
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_acceleration_temperature(
    unsigned phase, double dt) noexcept {
  bool status = false;
  const double tolerance = 1.0E-15;

  // rock freezing
  // this->current_time_ += dt * 1E5;
  // double ambient_T = -20;
  // if (current_time_ > 259200) ambient_T = 20;
  // if (id_ < 51){
  //   this->convective_heat_flux_(phase) = (ambient_T - this->temperature_(phase)) * 100 * 0.002258 * 1E5;    
  // }
  // if (id_ == 50){
  //   this->convective_heat_flux_(phase) = (ambient_T - this->temperature_(phase)) * 100 * 0.002258 * 1E5 / 2;    
  // }
  // if (id_ == 34){
  //   this->convective_heat_flux_(phase) = (ambient_T - this->temperature_(phase)) * 100 * 0.002258 * 1E5 / 2;    
  // }

  // rock freezing
  // this->current_time_ += dt * 1E5;
  // double ambient_T = -20;
  // if (current_time_ > 259200) ambient_T = 20;
  // if (id_ < 51){
  //   this->convective_heat_flux_(phase) = (ambient_T - this->temperature_(phase)) * 100 * 0.002258 * 1E5;    
  // }
  // if (id_ == 50){
  //   this->convective_heat_flux_(phase) = (ambient_T - this->temperature_(phase)) * 100 * 0.002258 * 1E5 / 2;    
  // }
  // if (id_ == 34){
  //   this->convective_heat_flux_(phase) = (ambient_T - this->temperature_(phase)) * 100 * 0.002258 * 1E5 / 2;    
  // }


  if ((this->heat_capacity_(phase) + this->latent_heat_(phase)) > tolerance) {
    // temperature acceleration = (total_heat / heat_capacity)

    this->temperature_acceleration_(phase) =
        (this->heat_conduction_(phase) + this->heat_source_(phase) 
        + this->convective_heat_flux_(phase) + this->plastic_work_(phase))  / 
        this->heat_capacity_(phase);

    // this->temperature_acceleration_(phase) =
    //     (this->heat_conduction_(phase) 
    //     + this->convective_heat_flux_(phase) + 
    //     this->heat_source_(phase) + this->plastic_work_(phase)) / 
    //     (this->heat_capacity_(phase) + this->latent_heat_(phase));

    this->temperature_ += this->temperature_acceleration_ * dt;

    this->current_time_ += dt;
    if ((this->free_surface()) || (this->contact_)) {
      if (this->current_time_ >= 1.0) {
        this->temperature_(phase) = 10;
        this->temperature_acceleration_(phase) = 0;
      }
    }

    status = true;
  }
  return status;
}

// Compute nodal pressure acceleration and update nodal pressure
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_acceleration_pressure(
    unsigned phase, double dt) noexcept {
  bool status = false;
  const double tolerance = 1.0E-15;

  if (this->mass_(phase) > tolerance) {

    this->pressure_acceleration_(phase) =
        this->hydraulic_conduction_(phase) / this->mass_(phase);

    this->pressure_ += this->pressure_acceleration_ * dt;

    if (this->free_surface_) {
      this->pressure_acceleration_.setZero();
      this->pressure_.setZero();
    }

    status = true;
  }
  return status;
}

// Compute mass density (Z. Wiezckowski, 2004)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_density() {
  try {
    const double tolerance = 1.E-15;  // std::numeric_limits<double>::lowest();

    for (unsigned phase = 0; phase < Tnphases; ++phase) {
      if (mass_(phase) > tolerance) {
        if (volume_(phase) > tolerance)
          density_(phase) = mass_(phase) / volume_(phase);

        // Check to see if value is below threshold
        if (std::abs(density_(phase)) < 1.E-15) density_(phase) = 0.;
      } else {
        density_(phase) = 0.;
      }
    }

    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints();
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
}

//==============================================================================
// ASSIGN AND APPLY NODAL BOUDARY CONDITIONS
//==============================================================================

// Assign velocity constraint
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_velocity_constraint(
    unsigned dir, double velocity) {
  bool status = true;
  try {
    // Constrain directions can take values between 0 and Dim * Nphases
    if (dir < (Tdim * Tnphases)) {
      // this->velocity_constraints_.insert(std::make_pair<unsigned, double>(
      //     static_cast<unsigned>(dir), static_cast<double>(velocity)));
      this->velocity_constraints_[dir] = velocity;
      this->set_velocity_constraints_ = true; 
    }
    else
      throw std::runtime_error("Constraint direction is out of bounds");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Apply velocity constraints
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::apply_velocity_constraints() {
  // Set velocity constraint
  for (const auto& constraint : this->velocity_constraints_) {
    // Direction value in the constraint (0, Dim * Nphases)
    const unsigned dir = constraint.first;
    // Direction: dir % Tdim (modulus)
    const auto direction = static_cast<unsigned>(dir % Tdim);
    // Phase: Integer value of division (dir / Tdim)
    const auto phase = static_cast<unsigned>(dir / Tdim);

    if (!generic_boundary_constraints_) {
      // Velocity constraints are applied on Cartesian boundaries
      this->velocity_(direction, phase) = constraint.second;
      this->velocity_inter_(direction, phase) = constraint.second;
      // Set acceleration to 0 in direction of velocity constraint
      this->acceleration_(direction, phase) = 0.;
      this->acceleration_inter_(direction, phase) = 0.;

    } else {
      // Velocity constraints on general boundaries
      // Compute inverse rotation matrix
      const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
          rotation_matrix_.inverse();
      // Transform to local coordinate
      Eigen::Matrix<double, Tdim, Tnphases> local_velocity =
          inverse_rotation_matrix * this->velocity_;
      Eigen::Matrix<double, Tdim, Tnphases> local_acceleration =
          inverse_rotation_matrix * this->acceleration_;
      // Apply boundary condition in local coordinate
      local_velocity(direction, phase) = constraint.second;
      local_acceleration(direction, phase) = 0.;
      // Transform back to global coordinate
      this->velocity_ = rotation_matrix_ * local_velocity;
      this->acceleration_ = rotation_matrix_ * local_acceleration;
    }
  }
}

// Assign temperature constraint
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_temperature_constraint(
    const unsigned phase, const double temperature,
    const std::shared_ptr<FunctionBase>& function) {
  bool status = true;
  try {
    // Constrain directions can take values between 0 and Tnphases
    if (phase < Tnphases * 2) {
      this->temperature_constraints_.insert(std::make_pair<unsigned, double>(
           static_cast<unsigned>(phase), static_cast<double>(temperature)));
      // Assign temperature function
      if (function != nullptr)
        this->temperature_function_.insert(
            std::make_pair<unsigned, std::shared_ptr<FunctionBase>>(
                static_cast<unsigned>(phase),
                static_cast<std::shared_ptr<FunctionBase>>(function)));
    } else
      throw std::runtime_error("Temperature constraint phase is out of bounds");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign temperature constraint
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_convective_heat_constraint(
    const unsigned phase, const double temperature,
    const std::shared_ptr<FunctionBase>& function, 
    const double coeff, const int set_id) {
  bool status = true;
  try {
    // Constrain directions can take values between 0 and Tnphases
    if (phase < Tnphases * 2) {
      this->convective_heat_constraints_.insert(std::make_pair<unsigned, double>(
          static_cast<unsigned>(coeff), static_cast<double>(temperature)));
      if (set_id == -1) this->convective_heat_boundary_ = true;
      // Assign temperature function
      if (function != nullptr)
        this->temperature_function_.insert(
            std::make_pair<unsigned, std::shared_ptr<FunctionBase>>(
                static_cast<unsigned>(phase),
                static_cast<std::shared_ptr<FunctionBase>>(function)));
    } else
      throw std::runtime_error("Temperature constraint phase is out of bounds");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign pressure constraint
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_pressure_constraint(
    const unsigned phase, const double pressure,
    const std::shared_ptr<FunctionBase>& function) {
  bool status = true;     
  try {
    // Constrain directions can take values between 0 and Tnphases
    if (phase < Tnphases * 2) {
      this->pressure_constraints_.insert(std::make_pair<unsigned, double>(
          static_cast<unsigned>(phase), static_cast<double>(pressure)));
      // Assign pressure function
      if (function != nullptr)
        this->pressure_function_.insert(
            std::make_pair<unsigned, std::shared_ptr<FunctionBase>>(
                static_cast<unsigned>(phase),
                static_cast<std::shared_ptr<FunctionBase>>(function)));
    } else
      throw std::runtime_error("Pressure constraint phase is out of bounds");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Apply temperature constraint
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::apply_temperature_constraints(
    const unsigned phase, const double current_time) {
  bool status = true;
  try {
    // Set temperature constraint
    for (const auto& constraint : this->temperature_constraints_) {
      const double scalar =
          (temperature_function_.find(phase) != temperature_function_.end())
              ? temperature_function_[phase]->value(current_time)
              : 1.0;     
      // Phase
      const auto phase = static_cast<unsigned>(constraint.first);
      // temperature constraints are applied
      this->temperature_(phase) = scalar * constraint.second;
      this->temperature_acceleration_(phase) = 0;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Apply convective_heat constraint
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::apply_convective_heat_constraints(
    const unsigned phase, const double current_time) 
{
  bool status = true;
  try 
  {
    // Set temperature constraint
    for (const auto& constraint : this->convective_heat_constraints_) 
    {
      const double scalar =
          (temperature_function_.find(phase) != temperature_function_.end())
              ? temperature_function_[phase]->value(current_time)
              : 1.0;     
      // Phase
      const double conv_coeff = constraint.first;
      const double temperature = scalar * constraint.second;
      // temperature constraints are applied
      if (!this->convective_heat_boundary_ ||
          (this->convective_heat_boundary_ & this->free_surface_))
      {
        this->convective_heat_flux_(phase) = 
          (temperature - this->temperature_(phase)) * conv_coeff * mean_length_;
        // this->convective_heat_flux_(phase) = conv_coeff * mean_length_; 
      }
    }
  } 
  catch (std::exception& exception)
  {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Apply pressure constraints
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::apply_pressure_constraints(
    const unsigned phase, const double current_time) {
  bool status = true;
  try {
    // Set pressure constraint
    for (const auto& constraint : this->pressure_constraints_) {
      const double scalar =
          (pressure_function_.find(phase) != pressure_function_.end())
              ? pressure_function_[phase]->value(current_time)
              : 1.0;
      // Phase
      const auto phase = static_cast<unsigned>(constraint.first);
      // Pressure constraints are applied
      this->pressure_(phase) = scalar * constraint.second;
      this->pressure_acceleration_(phase) = 0;

      if (Tnphases > 2) {
        this->pressure_(phase + 1) = scalar * constraint.second;
        this->pressure_acceleration_(phase + 1) = 0;
    }

    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign concentrated force to the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_concentrated_force(
    unsigned phase, unsigned direction, double concentrated_force,
    const std::shared_ptr<FunctionBase>& function) {
  bool status = false;
  try {
    if (phase >= Tnphases || direction >= Tdim) {
      throw std::runtime_error(
          "Cannot assign nodal concentrated forcey: Direction / phase is "
          "invalid");
    }
    // Assign concentrated force
    concentrated_force_(direction, phase) = concentrated_force;
    status = true;
    this->force_function_ = function;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Apply concentrated force to the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::apply_concentrated_force(
    unsigned phase, double current_time) {
  const double scalar =
      (force_function_ != nullptr) ? force_function_->value(current_time) : 1.0;
  this->update_external_force(true, phase,
                              scalar * concentrated_force_.col(phase));
}

// Assign friction constraint
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_friction_constraint(
    unsigned dir, int sign_n, double friction) {
  bool status = true;
  try {
    // Constrain directions can take values between 0 and Dim * Nphases
    if (dir < Tdim) {
      this->friction_constraint_ =
          std::make_tuple(static_cast<unsigned>(dir), static_cast<int>(sign_n),
                          static_cast<double>(friction));
      this->friction_ = true;
    } else
      throw std::runtime_error("Constraint direction is out of bounds");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Apply friction constraints
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::apply_friction_constraints(double dt) {
  if (friction_) {
    auto sign = [](double value) { return (value > 0.) ? 1. : -1.; };

    // Set friction constraint
    // Direction value in the constraint (0, Dim)
    const unsigned dir_n = std::get<0>(this->friction_constraint_);

    // Normal direction of friction: outward direction
    const double sign_dir_n = sign(std::get<1>(this->friction_constraint_));

    // Friction co-efficient
    const double mu = std::get<2>(this->friction_constraint_);

    for (unsigned phase = 0; phase < Tnphases; ++phase) {

      // Acceleration and velocity
      double acc_n, acc_t, vel_t;

      if (Tdim == 2) {
        // tangential direction to boundary
        const unsigned dir_t = (Tdim - 1) - dir_n;  // TODO

        if (!generic_boundary_constraints_) {
          // Cartesian case
          // Normal and tangential acceleration
          acc_n = this->acceleration_(dir_n, phase);
          acc_t = this->acceleration_(dir_t, phase);
          // Velocity tangential
          vel_t = this->velocity_(dir_t, phase);
        } else {
          // General case, transform to local coordinate
          // Compute inverse rotation matrix
          const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
              rotation_matrix_.inverse();
          // Transform to local coordinate
          Eigen::Matrix<double, Tdim, Tnphases> local_acceleration =
              inverse_rotation_matrix * this->acceleration_;
          Eigen::Matrix<double, Tdim, Tnphases> local_velocity =
              inverse_rotation_matrix * this->velocity_;
          // Normal and tangential acceleration
          acc_n = local_acceleration(dir_n, phase);
          acc_t = local_acceleration(dir_t, phase);
          // Velocity tangential
          vel_t = local_velocity(dir_t, phase);
        }

        if ((acc_n * sign_dir_n) > 0.0) {
          if (vel_t != 0.0) {  // kinetic friction
            const double vel_net = dt * acc_t + vel_t;
            const double vel_frictional = dt * mu * std::abs(acc_n);
            if (std::abs(vel_net) <= vel_frictional)
              acc_t = -vel_t / dt;
            else
              acc_t -= sign(vel_net) * mu * std::abs(acc_n);
          } else {  // static friction
            if (std::abs(acc_t) <= mu * std::abs(acc_n))
              acc_t = 0.0;
            else
              acc_t -= sign(acc_t) * mu * std::abs(acc_n);
          }

          if (!generic_boundary_constraints_) {
            // Cartesian case
            this->acceleration_(dir_t, phase) = acc_t;
          } else {
            // Local acceleration in terms of tangential and normal
            Eigen::Matrix<double, Tdim, Tnphases> acc;
            acc(dir_t, phase) = acc_t;
            acc(dir_n, phase) = acc_n;

            // General case, transform to global coordinate
            this->acceleration_.col(phase) = rotation_matrix_ * acc.col(phase);
          }
        }
      
      } else if (Tdim == 3) {
        Eigen::Matrix<int, 3, 2> dir;
        dir(0, 0) = 1;
        dir(0, 1) = 2;  // tangential directions for dir_n = 0
        dir(1, 0) = 0;
        dir(1, 1) = 2;  // tangential directions for dir_n = 1
        dir(2, 0) = 0;
        dir(2, 1) = 1;  // tangential directions for dir_n = 2

        const unsigned dir_t0 = dir(dir_n, 0);
        const unsigned dir_t1 = dir(dir_n, 1);

        Eigen::Matrix<double, Tdim, 1> acc, vel;
        if (!generic_boundary_constraints_) {
          // Cartesian case
          acc = this->acceleration_.col(phase);
          vel = this->velocity_.col(phase);
        } else {
          // General case, transform to local coordinate
          // Compute inverse rotation matrix
          const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
              rotation_matrix_.inverse();
          // Transform to local coordinate
          acc = inverse_rotation_matrix * this->acceleration_.col(phase);
          vel = inverse_rotation_matrix * this->velocity_.col(phase);
        }

        const auto acc_n = acc(dir_n);
        auto acc_t =
            std::sqrt(acc(dir_t0) * acc(dir_t0) + acc(dir_t1) * acc(dir_t1));
        const auto vel_t =
            std::sqrt(vel(dir_t0) * vel(dir_t0) + vel(dir_t1) * vel(dir_t1));

        if (acc_n * sign_dir_n > 0.0) {
          // kinetic friction
          if (std::abs(vel_t) > 1.E-15) {
            Eigen::Matrix<double, 2, 1> vel_net;
            // friction is applied opposite to the vel_net
            vel_net(0) = vel(dir_t0) + acc(dir_t0) * dt;
            vel_net(1) = vel(dir_t1) + acc(dir_t1) * dt;
            const double vel_net_t =
                sqrt(vel_net(0) * vel_net(0) + vel_net(1) * vel_net(1));
            const double vel_fricion = mu * std::abs(acc_n) * dt;

            if (vel_net_t <= vel_fricion) {
              acc(dir_t0) =
                  -vel(dir_t0) / dt;  // To set particle velocity to zero
              acc(dir_t1) = -vel(dir_t1) / dt;
            } else {
              acc(dir_t0) -= mu * std::abs(acc_n) * (vel_net(0) / vel_net_t);
              acc(dir_t1) -= mu * std::abs(acc_n) * (vel_net(1) / vel_net_t);
            }
          } else {                                // static friction
            if (acc_t <= mu * std::abs(acc_n)) {  // since acc_t is positive
              acc(dir_t0) = 0;
              acc(dir_t1) = 0;
            } else {
              acc_t -= mu * std::abs(acc_n);
              acc(dir_t0) -= mu * std::abs(acc_n) * (acc(dir_t0) / acc_t);
              acc(dir_t1) -= mu * std::abs(acc_n) * (acc(dir_t1) / acc_t);
            }
          }

          if (!generic_boundary_constraints_) {
            // Cartesian case
            this->acceleration_.col(phase) = acc;
          } else {
            // General case, transform to global coordinate
            this->acceleration_.col(phase) = rotation_matrix_ * acc;
          }
        }
      }
    }
  }
}

// Assign heat source
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_heat_source(
    const unsigned phase, const double heat_source,
    const std::shared_ptr<FunctionBase>& function) {
  bool status = true;
  try {
    // Constrain directions can take values between 0 and Tnphases
    if (phase < Tnphases * 2) {
      this->heat_source_constraints_.insert(std::make_pair<unsigned, double>(
          static_cast<unsigned>(phase), static_cast<double>(heat_source)));
      // Assign heat source function
      if (function != nullptr)
        this->heat_source_function_.insert(
            std::make_pair<unsigned, std::shared_ptr<FunctionBase>>(
                static_cast<unsigned>(phase),
                static_cast<std::shared_ptr<FunctionBase>>(function)));
    } else
      throw std::runtime_error("Heat_source phase is out of bounds");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Apply heat source
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::apply_heat_source(
    const unsigned phase, const double current_time) {
  bool status = true;
  try {
    // Set heat_source
    for (const auto& constraint : this->heat_source_constraints_) {
      const double scalar =
          (heat_source_function_.find(phase) != heat_source_function_.end())
              ? heat_source_function_[phase]->value(current_time)
              : 1.0;
      // Phase
      const auto phase = static_cast<unsigned>(constraint.first);
      // heat_source are applied
      this->heat_source_(phase) = scalar * constraint.second; 
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign water table
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_water_table(
    const std::shared_ptr<FunctionBase>& function, const unsigned dir,
    const double h0) {
  bool status = true;
  const unsigned pore_liquid = 1;
  try {
    // Assign zero pore pressure level
    this->h0_ = std::make_pair<unsigned, double>(static_cast<unsigned>(dir),
                                                 static_cast<double>(h0));
    // Assign math function for dynamic water table
    if (function != nullptr) {
      this->water_table_function_ = function;
    } else {
      if (this->coordinates_(dir) > h0) {
        // Assign 0 for node over zero pore pressure level
        pressure_constraints_[pore_liquid] = 0.;
      } else {
        // Compute pore pressure at node
        pressure_constraints_[pore_liquid] =
            (h0 - this->coordinates_(dir)) * 1000 * 9.81;
      }
    }
    // Activate water table
    this->water_table_ = true;

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Apply water table
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::apply_water_table(
    const double current_time) {
  bool status = true;
  const unsigned pore_liquid = 1;
  try {
    if (this->water_table_ && water_table_function_ != nullptr) {
      // Compute current h0
      const double h0 =
          this->h0_.second * water_table_function_->value(current_time);
      // Assign pore pressure constraint
      if (this->coordinates_(this->h0_.first) > h0) {
        // Assign 0 for node over zero pore pressure level
        pressure_constraints_[pore_liquid] = 0.;
      } else {
        // Compute pore pressure at node
        pressure_constraints_[pore_liquid] =
            (h0 - this->coordinates_(this->h0_.first)) * 1000 * 9.81;
      }
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign pressure at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_pressure(unsigned phase,
                                                      double pressure) {
  const double tolerance = 1.E-15;

  // Compute pressure from mass*pressure
  std::lock_guard<std::mutex> guard(node_mutex_);
  pressure_(phase) = pressure;
}

template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_nodal_pressure_increment_constraints(
                                        const Index step, const double dt) {
  bool status = true;
  try {
    const unsigned pore_liquid = 1;
    if (step == ref_step_) {

      const double scalar =
          (pressure_function_.find(pore_liquid + Tnphases) !=
           pressure_function_.end())
              ? pressure_function_[pore_liquid + Tnphases]->value(step * dt)
              : 1.0;

      pressure_constraints_[pore_liquid] =
          pressure_(pore_liquid) +
          scalar * pressure_constraints_.at(pore_liquid + Tnphases);

    } else if (step > ref_step_) {
      const double scalar =
          (pressure_function_.find(pore_liquid + Tnphases) !=
           pressure_function_.end())
              ? pressure_function_[pore_liquid + Tnphases]->value(step * dt)
              : 1.0;

      pressure_constraints_[pore_liquid] +=
          scalar * pressure_constraints_.at(pore_liquid + Tnphases);
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//============================================================================
// CONTACT： TO BE OPTIMIZED

template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_velocity_from_rigid(
    unsigned dir, double velocity, double dt) {
  std::lock_guard<std::mutex> guard(node_mutex_);

  switch (Tnphases) {
    // One phase
    case 1: {
      rigid_velocity_.col(0)[dir] = velocity;
      break;
    }
    // Two phase
    case 2: {
      if (dir < Tdim) {
        rigid_velocity_.col(0)[dir] = velocity;
        reaction_dir_[dir] = 1;
      } 
      else {
        rigid_velocity_.col(1)[dir - Tdim] = velocity;
      }
      break;
    }
    // Three phase
    case 3: {
      if (dir < Tdim) {
        rigid_velocity_.col(0)[dir] = velocity;
        reaction_dir_[dir] = 1;
      } 
      else if (dir < (Tdim * 2)) {
        rigid_velocity_.col(1)[dir - Tdim] = velocity;
      }
      else {
        rigid_velocity_.col(2)[dir - Tdim * 2] = velocity;
      }
      break;
    }
  }
  this->contact_ = true;
}

template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_contact_normal(
      double mass, const Eigen::Matrix<double, Tdim, 1>& mass_contact_normal) {
  // contact_normal_ = contact_normal;
  rigid_mass2_ += mass;
  mass_contact_normal_ += mass_contact_normal;  
}

// To be deleted
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_rigid_mass_momentum(
        double mass, const Eigen::MatrixXd& momentum) {
  std::lock_guard<std::mutex> guard(node_mutex_);
  rigid_mass_ += mass;
  rigid_momentum_ += momentum;
}

// Compute rigid velocity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_rigid_velocity(double dt) {
  const double tolerance = 1.E-15;
  if (this->contact_) {
    // Compute acceleration
    rigid_acceleration_ = (rigid_velocity_ - velocity_) / dt;
    // velocity_ = rigid_velocity_;
    // acceleration_ = rigid_acceleration_;
  }
}

template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_corrected_velocity_from_rigid(
                                          double dt) {
  if (this->contact_) {
    velocity_ = rigid_velocity_;
    acceleration_ = rigid_acceleration_;
  }
}

template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_intermediate_velocity_from_rigid(
                                          double dt) {
  if (this->contact_) {
    velocity_inter_ = rigid_velocity_;
    acceleration_inter_ = rigid_acceleration_;
  }
}
