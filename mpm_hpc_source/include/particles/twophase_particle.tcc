// Construct a two phase particle with id and coordinates
template <unsigned Tdim>
mpm::TwoPhaseParticle<Tdim>::TwoPhaseParticle(Index id, const VectorDim& coord)
    : mpm::Particle<Tdim>(id, coord) {
  this->initialise_liquid_phase();

  // Set material pointer to null
  liquid_material_ = nullptr;
  // Logger
  std::string logger =
      "twophaseparticle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//==============================================================================
// ASSIGN INITIAL CONDITIONS
//==============================================================================

// Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::initialise_particle(
    const HDF5Particle& particle) {
  // Derive from particle
  mpm::Particle<Tdim>::initialise_particle(particle);

  // // Liquid density
  // this->liquid_density_ = particle.liquid_density;
  // Liquid mass
  this->liquid_mass_ = particle.liquid_mass;
  // Liquid mass Density
  this->liquid_mass_density_ = particle.liquid_mass / particle.volume;

  // Pore pressure
  this->pore_pressure_ = particle.pore_pressure;

  // Liquid velocity
  Eigen::Vector3d liquid_velocity;
  liquid_velocity << particle.liquid_velocity_x, particle.liquid_velocity_y,
      particle.liquid_velocity_z;
  // Initialise velocity
  for (unsigned i = 0; i < Tdim; ++i)
    this->liquid_velocity_(i) = liquid_velocity(i);

  // Liquid strain
  this->liquid_strain_[0] = particle.liquid_strain_xx;
  this->liquid_strain_[1] = particle.liquid_strain_yy;
  this->liquid_strain_[2] = particle.liquid_strain_zz;
  this->liquid_strain_[3] = particle.liquid_gamma_xy;
  this->liquid_strain_[4] = particle.liquid_gamma_yz;
  this->liquid_strain_[5] = particle.liquid_gamma_xz;

  // Liquid material id
  this->liquid_material_id_ = particle.liquid_material_id;

  return true;
}

// Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::initialise_particle(
    const HDF5Particle& particle,
    const std::shared_ptr<mpm::Material<Tdim>>& material) {

  bool status = this->initialise_particle(particle);
  if (material != nullptr) {
    if (this->material_id_ == material->id() ||
        this->material_id_ == std::numeric_limits<unsigned>::max()) {
      material_ = material;
      // Reinitialize state variables
      auto mat_state_vars = material_->initialise_state_variables();
      if (mat_state_vars.size() == particle.nstate_vars) {
        unsigned i = 0;
        for (const auto& mat_state_var : mat_state_vars) {
          this->state_variables_[mat_state_var.first] = particle.svars[i];
          ++i;
        }
      }
    } else {
      status = false;
      throw std::runtime_error("Material is invalid to assign to particle!");
    }
  }
  return status;
}

// Initialise liquid phase particle properties
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::initialise_liquid_phase() {

  // Scalar data
  liquid_mass_ = 0.;
  liquid_density_ = 0.;
  liquid_fraction_ = 0.;
  pore_pressure_ = 0.;
  permeability_ = std::numeric_limits<double>::max();
  liquid_volumetric_strain_ = 0.;
  liquid_dvolumetric_strain_ = 0.;

  // Vctor data
  liquid_velocity_.setZero();
  liquid_traction_.setZero(); 

  // Tensor data
  liquid_strain_rate_.setZero();
  liquid_strain_.setZero();
  liquid_flux_.setZero();  

  // Bool data
  set_liquid_traction_ = false;
  set_mixture_traction_ = false;
  set_contact_ = false;  

  // Link data with NAME
  this->scalar_property_.insert({
    {"pore_pressures",    [&]() {return this->pore_pressure_;}},
    {"liquid_densities",  [&]() {return this->liquid_density_;}},
    {"liquid_fractions",  [&]() {return this->porosity_;}},
    {"permeabilities",    [&]() {return this->permeability_;}},
    {"PIC_pore_pressures",[&]() {return this->PIC_pore_pressure_;}},
    {"liquid_volumetric_strains",     [&]() {return this->liquid_volumetric_strain_;}  },
    {"liquid_dvolumetric_strains",     [&]() {return this->liquid_dvolumetric_strain_;}  },
  });

  this->vector_property_.insert({
    {"liquid_velocities", [&]() {return this->liquid_velocity_;}},
    {"liquid_strains",    [&]() {return this->liquid_strain_;}},
    {"liquid_fluxes",     [&]() {return this->liquid_flux_;}}
  }); 
}

// Return particle data in HDF5 format
template <unsigned Tdim>
// cppcheck-suppress *
mpm::HDF5Particle mpm::TwoPhaseParticle<Tdim>::hdf5() {
  // Derive from particle
  auto particle_data = mpm::Particle<Tdim>::hdf5();
  // Particle liquid velocity
  Eigen::Vector3d liquid_velocity;
  liquid_velocity.setZero();
  for (unsigned j = 0; j < Tdim; ++j)
    liquid_velocity[j] = this->liquid_velocity_[j];
  // Particle liquid strain
  Eigen::Matrix<double, 6, 1> liquid_strain = this->liquid_strain_;
  // // Particle liquid density
  // particle_data.liquid_density = this->liquid_density_;
  // Particle liquid mass
  particle_data.liquid_mass = this->liquid_mass_;
  // Particle pore pressure
  particle_data.pore_pressure = this->pore_pressure_;
  particle_data.liquid_velocity_x = liquid_velocity[0];
  particle_data.liquid_velocity_y = liquid_velocity[1];
  particle_data.liquid_velocity_z = liquid_velocity[2];
  // Particle liquid strain
  particle_data.liquid_strain_xx = liquid_strain[0];
  particle_data.liquid_strain_yy = liquid_strain[1];
  particle_data.liquid_strain_zz = liquid_strain[2];
  particle_data.liquid_gamma_xy = liquid_strain[3];
  particle_data.liquid_gamma_yz = liquid_strain[4];
  particle_data.liquid_gamma_xz = liquid_strain[5];
  // Particle liquid material id
  particle_data.liquid_material_id = this->liquid_material_id_;

  return particle_data;
}

// Assign a liquid material to particle
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_liquid_material(
    const std::shared_ptr<Material<Tdim>>& material) {
  bool status = false;
  try {
    // Check if material is valid and properties are set
    if (material != nullptr) {
      liquid_material_ = material;
      liquid_material_id_ = liquid_material_->id();
      status = true;
    } else {
      throw std::runtime_error("Material is undefined!");
    }
   
    // Assign permeability   
    this->assign_permeability();

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

// Compute mass of particle (both solid and fluid)
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_mass() {
  mpm::Particle<Tdim>::compute_mass();
  this->compute_liquid_mass();
}

// Compute fluid mass of particle
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_liquid_mass() noexcept {
  // Check if particle volume is set and liquid material ptr is valid
  assert(volume_ != std::numeric_limits<double>::max() &&
          liquid_material_ != nullptr);

  this->liquid_density_ = 
    liquid_material_->template property<double>(std::string("density"));

  // Mass = volume of particle * bulk_density
  this->liquid_mass_density_ =
      porosity_ * this->liquid_density_;
  this->liquid_mass_ = volume_ * liquid_mass_density_;
}

// Assign particle permeability
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_permeability() {
  bool status = true;
  try {
    // Check if material ptr is valid
    if (material_ != nullptr) {
      // Porosity parameter
      permeability_ = material_->template property<double>("k_x");
    } else {
      throw std::runtime_error("Material is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Initial pore pressure
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::initialise_pore_pressure_watertable(
    const unsigned dir_v, const unsigned dir_h,
    std::map<double, double>& refernece_points) {
  bool status = true;
  try {
    // Initialise left boundary position (coordinate) and h0
    double left_boundary = std::numeric_limits<double>::lowest();
    double h0_left = 0.;
    // Initialise right boundary position (coordinate) and h0
    double right_boundary = std::numeric_limits<double>::max();
    double h0_right = 0.;
    // Position and h0 of particle (coordinate)
    const double position = this->coordinates_(dir_h);
    // Iterate over each refernece_points
    for (const auto& refernece_point : refernece_points) {
      // Find boundary
      if (refernece_point.first > left_boundary &&
          refernece_point.first <= position) {
        // Left boundary position and h0
        left_boundary = refernece_point.first;
        h0_left = refernece_point.second;
      } else if (refernece_point.first > position &&
                 refernece_point.first <= right_boundary) {
        // Right boundary position and h0
        right_boundary = refernece_point.first;
        h0_right = refernece_point.second;
      }
    }

    if (left_boundary != std::numeric_limits<double>::lowest()) {
      // Particle with left and right boundary
      if (right_boundary != std::numeric_limits<double>::max()) {
        this->pore_pressure_ =
            ((h0_right - h0_left) / (right_boundary - left_boundary) *
                 (position - left_boundary) +
             h0_left - this->coordinates_(dir_v)) *
            1000 * 9.81;
      } else
        // Particle with only left boundary
        this->pore_pressure_ =
            (h0_left - this->coordinates_(dir_v)) * 1000 * 9.81;
    }
    // Particle with only right boundary
    else if (right_boundary != std::numeric_limits<double>::max())
      this->pore_pressure_ =
          (h0_right - this->coordinates_(dir_v)) * 1000 * 9.81;

    else
      throw std::runtime_error(
          "Particle pore pressure can not be initialised by water table");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//==============================================================================
// APPLY BOUNDARY CONDITIONS
//==============================================================================

// Assign traction
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_particle_traction(unsigned direction,
                                                  double traction) {
  bool status = true;
  this->assign_mixture_traction(direction, traction);
  return status;
}

// Assign traction to the liquid phase
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_liquid_traction(unsigned direction,
                                                          double traction) {
  bool status = false;
  try {
    if (direction >= Tdim ||
        this->volume_ == std::numeric_limits<double>::max()) {
      throw std::runtime_error(
          "Particle liquid traction property: volume / direction is invalid");
    }
    // Assign liquid traction
    liquid_traction_(direction) =
        traction * this->volume_ / this->size_(direction);
    status = true;
    this->set_liquid_traction_ = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign traction to the mixture
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_mixture_traction(unsigned direction,
                                                          double traction) {
  bool status = false;
  try {
    if (direction >= Tdim ||
        this->volume_ == std::numeric_limits<double>::max()) {
      throw std::runtime_error(
          "Particle mixture traction property: volume / direction is invalid");
    }
    // Assign mixture traction
    mixture_traction_(direction) =
        traction * this->volume_ / this->size_(direction);
    status = true;
    this->set_mixture_traction_ = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign particle liquid phase velocity constraint
// Constrain directions can take values between 0 and Dim
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_particle_liquid_velocity_constraint(
    unsigned dir, double velocity) {
  bool status = true;
  try {
    // Constrain directions can take values between 0 and Dim
    if (dir < Tdim)
      this->liquid_velocity_constraints_.insert(
          std::make_pair<unsigned, double>(static_cast<unsigned>(dir),
                                           static_cast<double>(velocity)));
    else
      throw std::runtime_error(
          "Particle liquid velocity constraint direction is out of bounds");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Apply particle velocity constraints
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::apply_particle_liquid_velocity_constraints() {
  // Set particle velocity constraint
  for (const auto& constraint : this->liquid_velocity_constraints_) {
    // Direction value in the constraint (0, Dim)
    const unsigned dir = constraint.first;
    // Direction: dir % Tdim (modulus)
    const auto direction = static_cast<unsigned>(dir % Tdim);
    this->liquid_velocity_(direction) = constraint.second;
  }
}

// Assign particle pressure constraints
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_particle_pore_pressure_constraint(
    double pressure) {
  bool status = true;
  try {
    this->pore_pressure_constraint_ = pressure;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Apply particle pore pressure constraints
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::apply_particle_pore_pressure_constraints(double pore_pressure) {
  // Set particle temperature constraint
  this->pore_pressure_ = pore_pressure;
  this->set_pressure_constraint_ = true;
}

//==============================================================================
// MAP PARTICLE INFORMATION TO NODES
//==============================================================================

// Map particle mass and momentum to nodes
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_mass_momentum_to_nodes() noexcept {
  mpm::Particle<Tdim>::map_mass_momentum_to_nodes();
  this->map_liquid_mass_momentum_to_nodes();
}

// Map liquid mass and momentum to nodes
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_liquid_mass_momentum_to_nodes() noexcept {
  // Check if liquid mass is set and positive
  assert(liquid_mass_ != std::numeric_limits<double>::max());
  if (this->material_id_ != 999){
    // Map liquid mass and momentum to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      nodes_[i]->update_mass(true, mpm::ParticlePhase::Liquid,
                            liquid_mass_ * shapefn_[i]);
      nodes_[i]->update_momentum(true, mpm::ParticlePhase::Liquid,
                                liquid_mass_ * shapefn_[i] * liquid_velocity_);
    }

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   std::cout << "    liquid_velocity_: " << liquid_velocity_ << std::endl;
  //   std::cout << "    liquid_mass_: " << liquid_mass_ << std::endl;
  // }

  }
}

// Map body force for both mixture and liquid
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_external_force(
    const VectorDim& pgravity) {
  this->map_mixture_body_force(mpm::ParticlePhase::Mixture, pgravity);
  this->map_liquid_body_force(pgravity);
}

// Map liquid phase body force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_liquid_body_force(
    const VectorDim& pgravity) noexcept {
  if (this->material_id_ != 999){        
    // Compute nodal liquid body forces
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_external_force(
          true, mpm::ParticlePhase::Liquid,
          (pgravity * this->liquid_mass_ * shapefn_(i)));   
          
  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   std::cout << "    liquid_force: " << pgravity * this->liquid_mass_ * shapefn_(0) << std::endl;
  // }              

  }
}

// Map mixture body force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_mixture_body_force(
    unsigned mixture, const VectorDim& pgravity) noexcept {
  if (this->material_id_ != 999){        
    // Compute nodal mixture body forces
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_external_force(
          true, mixture,
          (pgravity * (this->liquid_mass_ + this->mass_) * shapefn_(i)));      
          
  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   std::cout << "    mixture_force: " << pgravity * (this->liquid_mass_ + this->mass_) * shapefn_(0) << std::endl;
  // }      

  }
}

// Map traction force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_traction_force() noexcept {
  if (this->material_id_ != 999)  
    this->map_mixture_traction_force(mpm::ParticlePhase::Mixture);
}

// Map liquid traction force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_liquid_traction_force() noexcept {
  if (this->set_liquid_traction_) {
    // Map particle liquid traction forces to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_external_force(
          true, mpm::ParticlePhase::Liquid,
          (-1. * shapefn_[i] * porosity_ * this->liquid_traction_));
  }
}

// Map mixture traction force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_mixture_traction_force(
    unsigned mixture) noexcept {
  if (this->set_mixture_traction_) {
    // Map particle mixture traction forces to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      nodes_[i]->update_external_force(true, mixture,
                                       (shapefn_[i] * this->mixture_traction_)); 
    }                                  
  }
}

// Map both mixture and liquid internal force
template <unsigned Tdim>
inline void mpm::TwoPhaseParticle<Tdim>::map_internal_force() {
  if (this->material_id_ != 999){    
    mpm::TwoPhaseParticle<Tdim>::map_mixture_internal_force(
        mpm::ParticlePhase::Mixture);
    mpm::TwoPhaseParticle<Tdim>::map_liquid_internal_force();
  }
}

// Map both mixture and liquid internal force
template <unsigned Tdim>
inline void mpm::TwoPhaseParticle<Tdim>::map_internal_force_semi(double beta) {
  if (this->material_id_ != 999){    
    mpm::TwoPhaseParticle<Tdim>::map_mixture_internal_force(
        mpm::ParticlePhase::Mixture, beta);
    mpm::TwoPhaseParticle<Tdim>::map_liquid_internal_force(beta);
  }
}

// Map liquid phase internal force
template <>
void mpm::TwoPhaseParticle<2>::map_liquid_internal_force(double beta = 1) {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> pressure;
  pressure.setZero();
  pressure[0] = -this->pore_pressure_ * beta;
  pressure[1] = -this->pore_pressure_ * beta;

  if (is_axisymmetric_) pressure[2] = -this->pore_pressure_ * beta;

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 2, 1> force;
    force[0] = dn_dx_(i, 0) * pressure[0] * this->porosity_;
    force[1] = dn_dx_(i, 1) * pressure[1] * this->porosity_;    

    if (is_axisymmetric_) force[0] += shapefn_[i] / this->coordinates_(0) * 
                                      pressure[2] * this->porosity_;

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, force);
  }

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   // std::cout << "    liquid_force/v: " << -this->pore_pressure_ * beta * this->porosity_<< std::endl;
  //   std::cout << "    liquid_pressure_: " << pore_pressure_ << std::endl;
  //   // std::cout << "    dn_dx_(0,0): " << dn_dx_(0,0) << std::endl;
  //   // std::cout << "    dn_dx_(0,1): " << dn_dx_(0,1) << std::endl;
  //   // std::cout << "    dn_dx_(1,0): " << dn_dx_(1,0) << std::endl;
  //   // std::cout << "    dn_dx_(1,1): " << dn_dx_(1,1) << std::endl;
  // }

}

template <>
void mpm::TwoPhaseParticle<3>::map_liquid_internal_force(double beta = 1) {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> pressure;
  pressure.setZero();
  pressure[0] = -this->pore_pressure_ * beta;
  pressure[1] = -this->pore_pressure_ * beta;
  pressure[2] = -this->pore_pressure_ * beta;

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 3, 1> force;
    force[0] = dn_dx_(i, 0) * pressure[0] * this->porosity_;
    force[1] = dn_dx_(i, 1) * pressure[1] * this->porosity_;
    force[2] = dn_dx_(i, 2) * pressure[2] * this->porosity_;

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, force);
  }
}

// Map mixture internal force
template <>
void mpm::TwoPhaseParticle<2>::map_mixture_internal_force(unsigned mixture,
                                                          double beta = 1) {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress[0] -= this->pore_pressure_ * beta;
  total_stress[1] -= this->pore_pressure_ * beta; 

  if (is_axisymmetric_) total_stress[2] -= this->pore_pressure_ * beta;

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 2, 1> force;
    force[0] = dn_dx_(i, 0) * total_stress[0] + dn_dx_(i, 1) * total_stress[3];
    force[1] = dn_dx_(i, 1) * total_stress[1] + dn_dx_(i, 0) * total_stress[3];

    // force[0] = -dn_dx_centroid_(i, 0) * this->pore_pressure_ * beta;
    // force[1] = -dn_dx_centroid_(i, 1) * this->pore_pressure_ * beta;

    // force[0] += (dn_dx_(i, 0) + (dn_dx_centroid_(i, 0) - dn_dx_(i, 0)) / 2.) * stress_[0] +
    //           (dn_dx_centroid_(i, 0) - dn_dx_(i, 0)) / 2. * stress_[1] +
    //           dn_dx_(i, 1) * stress_[3];
    // force[1] += (dn_dx_centroid_(i, 1) - dn_dx_(i, 1)) / 2. * stress_[0] +
    //           (dn_dx_(i, 1) + (dn_dx_centroid_(i, 1) - dn_dx_(i, 1)) / 2.) * stress_[1] +
    //           dn_dx_(i, 0) * stress_[3];

    if (is_axisymmetric_) force[0] += shapefn_[i]/this->coordinates_(0) * total_stress[2];     

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mixture, force);
  }
}

template <>
void mpm::TwoPhaseParticle<3>::map_mixture_internal_force(unsigned mixture,
                                                          double beta = 1) {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -= this->pore_pressure_ * beta;
  total_stress(1) -= this->pore_pressure_ * beta;
  total_stress(2) -= this->pore_pressure_ * beta;

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 3, 1> force;
    force[0] = dn_dx_(i, 0) * total_stress[0] + dn_dx_(i, 1) * total_stress[3] +
               dn_dx_(i, 2) * total_stress[5];

    force[1] = dn_dx_(i, 1) * total_stress[1] + dn_dx_(i, 0) * total_stress[3] +
               dn_dx_(i, 2) * total_stress[4];

    force[2] = dn_dx_(i, 2) * total_stress[2] + dn_dx_(i, 1) * total_stress[4] +
               dn_dx_(i, 0) * total_stress[5];

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mixture, force);
  }
}

// Map particle pore liquid pressure to nodes
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_pore_pressure_to_nodes(
    double current_time) noexcept {

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   std::cout << "    update_pressure: " << liquid_mass_ * pore_pressure_<< std::endl;
  // }

  // Check if particle mass is set
  assert(liquid_mass_ != std::numeric_limits<double>::max());
  bool status = true;
  if (this->material_id_ != 999){
    // Map particle liquid mass and pore pressure to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_mass_pressure(mpm::ParticlePhase::Liquid,
                                      shapefn_[i] * liquid_mass_ * pore_pressure_,
                                      current_time);

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   //std::cout << "    update_pressure: " << liquid_mass_ * pore_pressure_<< std::endl;
  //   std::cout << "    nodes_[0]->pressure: " << nodes_[0]->pressure(mpm::ParticlePhase::Liquid) << std::endl;
  // }

  }
  return status;
}

// Map particle pore liquid pressure to nodes
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_mass_pressure_to_nodes() {
  // Check if particle mass is set
  assert(liquid_mass_ != std::numeric_limits<double>::max());
  if (this->material_id_ != 999){
    // Map particle liquid mass and pore pressure to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_pressure(true, mpm::ParticlePhase::Liquid,
                                      shapefn_[i] * liquid_mass_ * pore_pressure_);
  }
}

// Map drag force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_drag_force_coefficient() {
  if (this->material_id_ != 999){  
    try {
      // Update permeability
      this->update_permeability();
      // Initialise drag force coefficient
      double drag_force_coefficient;
      drag_force_coefficient = 0;

      // Check if permeability coefficient is valid
      if (permeability_ > 0.)
        drag_force_coefficient = porosity_ * porosity_ * 9.81 *
                                    this->liquid_density_ /
                                    permeability_;
      else throw std::runtime_error("Permeability coefficient is invalid");

      // Map drag forces from particle to nodes
      for (unsigned j = 0; j < nodes_.size(); ++j)
        nodes_[j]->update_drag_force_coefficient(
            true, mpm::ParticlePhase::Liquid, 
            drag_force_coefficient * this->volume_ * shapefn_(j));

    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    }
  }
}

// Map particle heat capacity and heat to nodes
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_heat_to_nodes() {
  if (this->material_id_ != 999){    
    mpm::Particle<Tdim>::map_heat_to_nodes();
    this->map_liquid_heat_to_nodes();
  }
}

// Map liquid heat capacity and heat to nodes
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_liquid_heat_to_nodes() noexcept {
  // Check if liquid mass is set and positive
  assert(liquid_mass_ != std::numeric_limits<double>::max());
  // get the specific_heat of liquid       
  double liquid_specific_heat_ = 
         liquid_material_->template property<double>(std::string("liquid_specific_heat"));
       
  // Map liquid heat capacity and heat to nodes
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    nodes_[i]->update_heat_capacity(true, mpm::ParticlePhase::Mixture,
                   liquid_mass_ * liquid_specific_heat_ * shapefn_[i]);
    nodes_[i]->update_heat(true, mpm::ParticlePhase::Mixture,
                   liquid_mass_ * liquid_specific_heat_ * shapefn_[i] * this->temperature_);      
  }
  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   std::cout << "    shapefn_[0]: " << shapefn_[0] << std::endl;
  //   std::cout << "    liquid_heat_capacity_: " << liquid_mass_ * liquid_specific_heat_ << std::endl;
  // }
}

// Map particle heat conduction to node
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_heat_conduction() {
  if (this->material_id_ != 999){    
    mpm::Particle<Tdim>::map_heat_conduction();
    mpm::TwoPhaseParticle<Tdim>::map_liquid_heat_conduction();
  }
}

// Map liquid phase heat conduction 
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_liquid_heat_conduction() noexcept {
  if (liquid_material_id_ != 999) {
    // Assign the liquid thermal conductivity
    const double liquid_k_conductivity = 
           liquid_material_->template property<double>(std::string("liquid_thermal_conductivity"));    
    
    // Compute nodal heat conduction
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      double heat_conduction = 0;
      for (unsigned j = 0; j < Tdim; ++j){
        heat_conduction += dn_dx_(i, j)  * this->temperature_gradient_[j]; 
      }
      heat_conduction *= -1 * this->volume_ * porosity_ * liquid_k_conductivity;
      nodes_[i]->update_heat_conduction(true, mpm::ParticlePhase::Mixture, heat_conduction);     
    }

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   //std::cout << "particle id: " << this->id_ << std::endl;
  //   //std::cout << "    volume_: " <<this->volume_;
  //   //std::cout << "    Lheat_conduction: " << -1 * this->volume_ * porosity_ * liquid_k_conductivity << std::endl;
  //   std::cout << "    temperature_gradient_y: " << this->temperature_gradient_[1] << std::endl;
  // }

  }
}

// Map heat convection of mixture
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_heat_convection() {
  if (liquid_material_id_ != 999) {
    // get the specific heat of liquid
    double liquid_specific_heat = 
          liquid_material_->template property<double>(std::string("liquid_specific_heat"));
    // Compute nodal heat convection
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      double heat_convection = 0;
      for (unsigned j = 0; j < Tdim; ++j){
        heat_convection += 
              shapefn_[i] * temperature_gradient_[j] * (liquid_velocity_[j] - this->velocity_[j]);
      }
      heat_convection *= -1 * liquid_mass_ * liquid_specific_heat;
      nodes_[i]->update_heat_convection(true, mpm::ParticlePhase::Mixture, heat_convection);
    }
  }
}

// Map hydraulic conduction
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_hydraulic_conduction() {
  // Assign the liquid hydraulic conductivity
  const double permeability = 
          material_->template property<double>(std::string("intrinsic_permeability"));
  const double viscosity = 
          liquid_material_->template property<double>(std::string("liquid_viscosity"));
  const double liquid_modulus = 
          liquid_material_->template property<double>(std::string("bulk_modulus"));

  const double hydraulic_conductivity = permeability / viscosity * liquid_modulus * liquid_density_;  

  // Assign liquid gradient
  this->compute_pressure_gradient(mpm::ParticlePhase::Liquid);

  // Compute nodal hydraulic conduction
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    double hydraulic_conduction = 0;
    for (unsigned j = 0; j < Tdim; ++j){
      hydraulic_conduction += dn_dx_(i, j) * this->pressure_gradient_[j]; 
    }
    hydraulic_conduction *= -1 * this->volume_ * hydraulic_conductivity;
    nodes_[i]->update_hydraulic_conduction(true, mpm::ParticlePhase::Liquid, hydraulic_conduction);     
  }
}

// Compute pressure gradient of the particle
template <unsigned Tdim>
inline Eigen::Matrix<double, Tdim, 1> mpm::TwoPhaseParticle<Tdim>::
                      compute_pressure_gradient(unsigned phase) noexcept {

  Eigen::Matrix<double, Tdim, 1> pressure_gradient;
  pressure_gradient.setZero();
  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    double pressure = nodes_[i]->pressure(phase);
    for (unsigned j = 0; j < Tdim; ++j) {
      // pressure_gradient = partial p / partial X = p_{i,j}
      pressure_gradient[j] += dn_dx_(i, j) * pressure;
      if (std::fabs(pressure_gradient[j]) < 1.E-15)
        pressure_gradient[j] = 0.;
    }
  }
  this->pressure_gradient_ = pressure_gradient;
  return pressure_gradient_;
}

//------------------------------------------------------------
// Implict mpm

// Map heat laplacian matrix to cell
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_MTT_to_cell() {
  bool status = true;
  if (this->material_id_ != 999){  
    try {
      // get the specific_heat of solid  
      double solid_specific_heat = 
            material_->template property<double>(std::string("specific_heat"));
      // get the specific_heat of liquid       
      double liquid_specific_heat = 
            liquid_material_->template property<double>(std::string("liquid_specific_heat"));

      double specific_heat_m = liquid_mass_density_ * liquid_specific_heat + 
                              mass_density_ * solid_specific_heat;
      //cell_->compute_MTT_element(shapefn_, dn_dx_, volume_, specific_heat_m); 
    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

// Map heat laplacian matrix to cell
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_KTT_to_cell() {
  bool status = true;
  if (this->material_id_ != 999){  
    try {
      const double solid_k_conductivity = 
           material_->template property<double>(std::string("thermal_conductivity"));  
      // get the thermal conductivity coefficient of solid
      const double liquid_k_conductivity = 
           liquid_material_->template property<double>(std::string("liquid_thermal_conductivity"));  

      double beta_m = porosity_ * liquid_k_conductivity + (1 - porosity_) * solid_k_conductivity;
      cell_->compute_KTT_element(shapefn_, dn_dx_, volume_, beta_m); 
    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

//------------------------------------------------------------
// Semi-implict mpm
// Compute local matrix of K_inter
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_K_inter_to_cell() {
  bool status = true;
  if (this->material_id_ != 999){  
    try {
      // Initialise multiplier
      VectorDim multiplier;
      multiplier.setZero();
      // Compute multiplier
      for (unsigned i = 0; i < Tdim; ++i)
        multiplier(i) =
            porosity_ * porosity_ * 9.81 *
            liquid_material_->template property<double>(std::string("density")) /
            permeability_;
      // Compute local matrix of K_inter
      cell_->compute_K_inter_element(shapefn_, volume_, multiplier);
      cell_->compute_average_shapefn(shapefn_, volume_);  
    }

    catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

// Compute laplacian matrix element
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_L_to_cell(double dt, double alpha) {
  bool status = true;
  if (this->material_id_ != 999){  
    try {
      if (alpha) {
          double drag = porosity_ * porosity_ * 9.81 * liquid_material_->template 
                        property<double>(std::string("density")) / 
                        material_->template property<double>("k_x");
          double ksi_s = (1 + dt * drag / liquid_mass_density_ / (1 - porosity_)) / 
                      (1 + dt * drag * (1 / mass_density_ + 1 / liquid_mass_density_));
          double ksi_l = (1 + dt * drag / mass_density_ / porosity_) / 
                      (1 + dt * drag * (1 / mass_density_ + 1 / liquid_mass_density_));
          const double multiplier_s = ksi_s * (1 - porosity_) /
                        material_->template property<double>(std::string("density"));;
          const double multiplier_l = ksi_l * porosity_ /
                        liquid_material_->template property<double>(std::string("density"));;
          const double alpha_liquid = porosity_ *
                liquid_material_->template property<double>(std::string("liquid_compressibility"));

          const double youngs_modulus = material_->template property<double>(std::string("youngs_modulus"));      
          const double poisson_ratio = material_->template property<double>(std::string( "poisson_ratio"));
          const double stab_para = (1 + poisson_ratio) / youngs_modulus * liquid_material_->template property<double>(std::string("stab_para"));    
          // Compute local matrix of Laplacian
          cell_->compute_L_element(shapefn_, dn_dx_, volume_, multiplier_s, multiplier_l, alpha_liquid);
          cell_->compute_stab_element(shapefn_, shapefn_centroid_, stab_para);                      
        } else {
          const double multiplier_s = (1 - porosity_) /
                  material_->template property<double>(std::string("density"));
          const double multiplier_l = porosity_ /
                  liquid_material_->template property<double>(std::string("density"));
          const double alpha_liquid = porosity_ *
                  liquid_material_->template property<double>(std::string("liquid_compressibility"));
          const double youngs_modulus = material_->template property<double>(std::string("youngs_modulus"));      
          const double poisson_ratio = material_->template property<double>(std::string( "poisson_ratio"));
          const double stab_para = (1 + poisson_ratio) / youngs_modulus * liquid_material_->template property<double>(std::string("stab_para"));     
          // Compute local matrix of Laplacian
          cell_->compute_L_element(shapefn_, dn_dx_, volume_, multiplier_s, multiplier_l, alpha_liquid); 
          cell_->compute_stab_element(shapefn_, shapefn_centroid_, stab_para);

        }   
    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

// Compute element F_s_element and F_m_element
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_F_to_cell() {
  bool status = true;
  if (this->material_id_ != 999){  
    try {
      // Compute local matrix of F
      cell_->compute_F_element(shapefn_, dn_dx_, volume_, 1, porosity_);

    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

// Map element T_element
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_P_to_cell(double beta) {
  bool status = true;
  if (this->material_id_ != 999){  
    try {
      double multiplier = 0;    
      if (!beta) {                    
        double alpha_liquid = porosity_ * liquid_material_->template 
                                  property<double>(std::string("liquid_compressibility"));
        multiplier = alpha_liquid * pore_pressure_; 
      }
      cell_->compute_P_element(shapefn_, volume_, multiplier);
    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

// Compute element T_element
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_T_to_cell() {
  bool status = true;
  if (this->material_id_ != 999){  
    try {
      // get the thermal conductivity coefficient of liquid   
      double beta_liquid =
          liquid_material_->template property<double>(std::string("liquid_thermal_expansivity"));
      // get the thermal conductivity coefficient of solid
      double beta_solid =
          material_->template property<double>(std::string("thermal_expansivity"));
      double beta =  beta_liquid * porosity_ + beta_solid * (1 - porosity_);
      // Compute local matrix of T
      cell_->compute_T_element(shapefn_, volume_, beta);   
    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

// Compute element K_cor_w_element_
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_K_cor_to_cell(double dt, double alpha) {
  bool status = true;
  if (this->material_id_ != 999){  
    try {
      if (alpha) {
        double drag = porosity_ * porosity_ * 9.81 * liquid_material_->template 
                      property<double>(std::string("density")) / 
                      material_->template property<double>("k_x");
        double ksi_s = (1 + dt * drag / liquid_mass_density_ / (1 - porosity_)) / 
                      (1 + dt * drag * (1 / mass_density_ + 1 / liquid_mass_density_));
        double ksi_l = (1 + dt * drag / mass_density_ / porosity_) / 
                      (1 + dt * drag * (1 / mass_density_ + 1 / liquid_mass_density_));
        const double multiplier_s = ksi_s * (1 - porosity_);
        const double multiplier_l = ksi_l * porosity_;
        // const double multiplier_s = ((1 - porosity_) + dt * drag / liquid_mass_density_) / 
        //             (1 + dt * drag * (1 / mass_density_ + 1 / liquid_mass_density_));
        // const double multiplier_l = (porosity_ + dt * drag / mass_density_) / 
        //             (1 + dt * drag * (1 / mass_density_ + 1 / liquid_mass_density_)); 
        cell_->compute_K_cor_element(shapefn_, dn_dx_, volume_, multiplier_s, multiplier_l, 
                                    this->coordinates_(0), is_axisymmetric_);                                                                      
      } else {
        const double multiplier_s = 1 - porosity_;
        const double multiplier_l = porosity_;
        cell_->compute_K_cor_element(shapefn_, dn_dx_, volume_, multiplier_s, multiplier_l, 
                                    this->coordinates_(0), is_axisymmetric_);    
      }
    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    }
  }
  return status;
}




//==============================================================================
// UPDATE PARTICLE INFORMATION
//==============================================================================

// Compute updated position of the particle and kinematics of both solid and
// liquid phase
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_updated_velocity(
    double dt, double pic, double damping_factor) {
  mpm::Particle<Tdim>::compute_updated_velocity(dt, pic, damping_factor);
  this->compute_updated_liquid_velocity(dt, pic, damping_factor);

}

// Compute updated velocity of the liquid
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_updated_liquid_velocity(
    double dt, double pic, double damping_factor) {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  assert((-1.E-15) <= pic && pic <= (1 + 1.E-15));

  // Get PIC velocity
  Eigen::Matrix<double, Tdim, 1> pic_velocity =
      Eigen::Matrix<double, Tdim, 1>::Zero();
  for (unsigned i = 0; i < nodes_.size(); ++i)
    pic_velocity +=
        shapefn_[i] * nodes_[i]->velocity(mpm::ParticlePhase::Liquid);

  // Get interpolated nodal acceleration
  Eigen::Matrix<double, Tdim, 1> nodal_acceleration =
      Eigen::Matrix<double, Tdim, 1>::Zero();
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodal_acceleration +=
        shapefn_[i] * nodes_[i]->acceleration(mpm::ParticlePhase::Liquid);

  // Applying particle damping
  nodal_acceleration -= damping_factor * this->liquid_velocity_;

  // Get FLIP velocity
  Eigen::Matrix<double, Tdim, 1> flip_velocity =
      this->liquid_velocity_ + nodal_acceleration * dt;

  // Update particle velocity based on PIC value
  this->liquid_velocity_ = pic * pic_velocity + (1. - pic) * flip_velocity;

  // Apply particle velocity constraints
  this->apply_particle_liquid_velocity_constraints();

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   std::cout << "    liquid_velocity_: " << this->liquid_velocity_ * dt<< std::endl;
  // }

}

// Compute updated pore_pressure of the particle
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::update_particle_pore_pressure(
    double dt, double pic_t) noexcept {
  if (this->material_id_ != 999) {      
    // Check if particle has a valid cell ptr and pic_t value
    assert(cell_ != nullptr);
    assert((-1.E-15) <= pic_t && pic_t <= (1 + 1.E-15));

    // Get PIC pressure
    double PIC_pore_pressure = 0;
    for (unsigned i = 0; i < nodes_.size(); ++i)
      PIC_pore_pressure +=
          shapefn_[i] * nodes_[i]->pressure(mpm::ParticlePhase::Liquid);
    // pressure increment
    this->pore_pressure_increment_ = PIC_pore_pressure - this->PIC_pore_pressure_; 

    this->PIC_pore_pressure_= PIC_pore_pressure;
    // Get interpolated nodal pressure acceleration
    this->pore_pressure_acceleration_ = 0.;
    for (unsigned i = 0; i < nodes_.size(); ++i)
      this->pore_pressure_acceleration_ +=
          shapefn_[i] * nodes_[i]->pressure_acceleration(mpm::ParticlePhase::Liquid);

    // Get FLIP pressure
    double FLIP_pore_pressure =
        this->pore_pressure_ + this->pore_pressure_acceleration_ * dt;
    this->FLIP_pore_pressure_ = FLIP_pore_pressure;

    // Update particle pressure based on PIC value
    this->pore_pressure_ = pic_t * PIC_pore_pressure + (1. - pic_t) * FLIP_pore_pressure;

  }
}

// Compute updated pore pressure of the particle based on nodal pressure
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::compute_updated_pore_pressure(double beta) {
  bool status = true;
  if (this->material_id_ != 999){  
    try {
      double pore_pressure_increment = 0.0;
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        pore_pressure_increment +=
            shapefn_(i) * nodes_[i]->pore_pressure_increment();
            // shapefn_centroid_(i) * nodes_[i]->pore_pressure_increment();
      }
      // Get interpolated nodal pore pressure
      this->pore_pressure_ = pore_pressure_ * beta + pore_pressure_increment;

      // this->pore_pressure_ = 0.0;
      // for (unsigned i = 0; i < nodes_.size(); ++i) {
      //   this->pore_pressure_ +=
      //       shapefn_(i) * nodes_[i]->pressure(mpm::ParticlePhase::Liquid);
      // }

      // Apply free surface
      if (this->free_surface()) {
        this->pore_pressure_ = 0.0;
      }

    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

// Compute pore pressure
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_pore_pressure(double dt) {
  // Check if liquid material and cell pointer are set and positive
  assert(liquid_material_ != nullptr && cell_ != nullptr);
  if (this->material_id_ != 999){
    // get the bulk modulus of liquid
    double K =
        liquid_material_->template property<double>(std::string("bulk_modulus"));

    // Compute at centroid
    // get liquid phase strain rate at cell centre
    auto liquid_strain_rate_centroid =
        this->compute_strain_rate(dn_dx_centroid_, shapefn_centroid_, mpm::ParticlePhase::Liquid);

    // update pressure
    if (!is_axisymmetric_) {
      this->pore_pressure_ +=
          -dt * (K / porosity_) *
          ((1 - porosity_) * strain_rate_.head(Tdim).sum() +
          porosity_ * liquid_strain_rate_centroid.head(Tdim).sum());
    } else {
        this->pore_pressure_ +=
        -dt * (K / porosity_) *
        ((1 - porosity_) * strain_rate_.head(3).sum() +
        porosity_ * liquid_strain_rate_centroid.head(3).sum());
    }
  }     
}

// Compute pore liquid pressure smoothing based on nodal pressure
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::compute_pore_pressure_smoothing() noexcept {

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   //std::cout << "    update_pressure: " << liquid_mass_ * pore_pressure_<< std::endl;
  //   std::cout << "    before:nodes_[0]->pressure: " << nodes_[0]->pressure(mpm::ParticlePhase::Liquid) << std::endl;
  //   std::cout << "    liquid_pressure_: " << this->pore_pressure_ << std::endl;
  // }

  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  bool status = true;
  if (this->material_id_ != 999){
    double pressure = 0;
    for (unsigned i = 0; i < nodes_.size(); ++i)
      pressure += shapefn_(i) * nodes_[i]->pressure(mpm::ParticlePhase::Liquid);

    // Update pore liquid pressure to interpolated nodal pressure
    this->pore_pressure_ = pressure;
    // Apply free surface
      if (this->free_surface()) {
        this->pore_pressure_ = 0.0;
      }

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   //std::cout << "    update_pressure: " << liquid_mass_ * pore_pressure_<< std::endl;
  //   std::cout << "    after:nodes_[0]->pressure: " << nodes_[0]->pressure(mpm::ParticlePhase::Liquid) << std::endl;
  //   std::cout << "    liquid_pressure_: " << this->pore_pressure_ << std::endl;
  // }

  }
  return status;
}

// Update particle permeability
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::update_permeability() {
  if (this->material_id_ != 999){  
    try {
      double porosity_0 = material_->template property<double>("porosity");
      // Porosity parameter
      double k_para = std::pow(this->porosity_, 3) / std::pow((1. - this->porosity_), 2) /
                (std::pow(porosity_0, 3) / std::pow((1. - porosity_0), 2));
      // Update permeability
      permeability_ *= k_para;
    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    }
  }
}

// update density of the particle
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::update_particle_density(double dt) {
  if (this->material_id_ != 999){    
    mpm::Particle<Tdim>::update_particle_density(dt);
    this->update_liquid_density();
  }
}

// update liquid density of the particle
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::update_liquid_density() noexcept {
  // Check if liquid material ptr is valid
  assert(liquid_material_ != nullptr);
  // get liquid initial density
  double liquid_density_0 = 
      liquid_material_->template property<double>(std::string("density"));
  // get the compressibility of liquid
  double alpha_liquid =
      liquid_material_->template property<double>(std::string("liquid_compressibility"));
  // get the thermal conductivity coefficient of liquid
  double beta_liquid =
      liquid_material_->template property<double>(std::string("liquid_thermal_expansivity"));
  this->liquid_density_ = liquid_density_0 *
                          std::exp(pore_pressure_ * alpha_liquid - 
                          3 * beta_liquid * PIC_temperature_);    
  // this->liquid_density_ *= std::exp(pore_pressure_increments_ * alpha_liquid - 
  //                         3 * beta_liquid * temperature_increment_);         

}

// update mass of the particle
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::update_particle_volume()  {
  if (this->material_id_ != 999) {
    mpm::Particle<Tdim>::update_particle_volume();  
    this->update_liquid_mass();
  }  
}

// update liquid mass of particle
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::update_liquid_mass() noexcept {
  // Check if particle volume is set and liquid material ptr is valid
  assert(volume_ != std::numeric_limits<double>::max() &&
         liquid_material_ != nullptr);

  // Mass = volume of particle * bulk_density
  this->liquid_mass_density_ = porosity_ * liquid_density_;
  this->liquid_mass_ = volume_ * liquid_mass_density_;

}

// Compute pore pressure
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_thermal_pore_pressure(double dt) noexcept {
  // Check if material and cell pointer are set and positive
  assert(liquid_material_ != nullptr && cell_ != nullptr);
  // Check if liquid material and cell pointer are set and positive
  assert(liquid_material_ != nullptr && cell_ != nullptr);
  if (this->material_id_ != 999){  

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   std::cout << "    pore_pressure_: " << pore_pressure_ << std::endl;
  // }

    // get the bulk modulus of liquid
    auto K =
        liquid_material_->template property<double>(std::string("bulk_modulus"));
    // get the thermal conductivity coefficient of liquid
    auto beta_liquid =
        liquid_material_->template property<double>(std::string("liquid_thermal_expansivity"));
    // get the thermal conductivity coefficient of solid
    auto beta_solid =
        material_->template property<double>(std::string("thermal_expansivity"));
    
    // Compute at centroid
    // get liquid phase strain rate at cell centre
    auto liquid_strain_rate_centroid =
        this->compute_strain_rate(dn_dx_, shapefn_, mpm::ParticlePhase::Liquid);

    auto strain_rate_centroid =
        this->compute_strain_rate(dn_dx_centroid_, shapefn_centroid_, mpm::ParticlePhase::Solid);

    // Compute mass-convection-induced pressure
    double convective_term = 0;
    // for (unsigned j = 0; j < Tdim; ++j){
    //     double convective_term += 
    //           (velocity_[i] - liquid_velocity_[i]) * porosity_gradient_[i];
    //}

    // Compute thermal-expansion-induced pressure
    double thermal_term = (3 * porosity_ * beta_liquid + 
                            3 * (1 - porosity_) * beta_solid) 
                          * this->temperature_increment_cent_;

    liquid_dvolumetric_strain_ = dt * liquid_strain_rate_centroid.head(Tdim).sum();  
                                  //+ 3 * beta_liquid * this->temperature_increment_;
    liquid_volumetric_strain_ += liquid_dvolumetric_strain_;

    double strain_term = 0;
    // Compute strain-induced pressure
    if (is_axisymmetric_) {
      strain_term = dt * ((1 - porosity_) * strain_rate_.head(3).sum() +
                  porosity_ * liquid_strain_rate_centroid.head(3).sum());  
    } else {           
      strain_term = dt * ((1 - porosity_) * strain_rate_.head(Tdim).sum() +
                  porosity_ * liquid_strain_rate_centroid.head(Tdim).sum());
    } 
    // update pressure
    if (porosity_ > 1E-3) {
      double pore_pressure_increments_ = 
        (K / porosity_) * (thermal_term - strain_term);
      this->pore_pressure_ += pore_pressure_increments_; 
      
      // this->pore_pressure_ = (K / porosity_) * (-porosity_ * liquid_volumetric_strain_ - (1 - porosity_) * volumetric_strain_);

    }

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   //std::cout << "    K_ww: " << porosity_ / K << std::endl;
  //   //std::cout << "    thermal_term: " << thermal_term << std::endl;
  //   //std::cout << "    strain_term: " << strain_term << std::endl;
  //   std::cout << "    strain_rate_.head(Tdim).sum(): " << strain_rate_.head(Tdim).sum() << std::endl;
  //   std::cout << "    liquid_strain_rate_centroid.head(Tdim).sum(): " << liquid_strain_rate_centroid.head(Tdim).sum() << std::endl;
  //   std::cout << "    after:pore_pressure_: " << pore_pressure_ << std::endl;
  //   std::cout << "    pore_pressure_increments_: " << (K / porosity_) * (thermal_term - strain_term) << std::endl;
  //   //std::cout << "    pore:temperature_increment_cent_: " << this->temperature_increment_cent_ << std::endl;
  // }

  }
}

//==============================================================================
// CONTACT
//==============================================================================

// Assign contact to the particle
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_particle_contact(unsigned dir, double normal) {
  bool status = false;
  try {
    // Assign contact
    contact_normal_(dir) = normal;
    status = true;
    this->set_contact_ = true;
    for (unsigned i = 0; i < std::pow(2, Tdim); ++i) {       
      nodes_[i]->update_contact_normal(shapefn_[i] * this->mass_,
                             shapefn_[i] * this->mass_ * contact_normal_);
    }    
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Overwrite node velocity to get strain correct
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_moving_rigid_velocity_to_nodes(
    unsigned dir, double velocity, double dt) noexcept {   
  if (this->material_id_ == 999){  
    for (unsigned i = 0; i < std::pow(2, Tdim); ++i) {
    // for (unsigned i = 0; i < nodes_.size(); ++i) {
      nodes_[i]->assign_velocity_from_rigid(dir, velocity, dt);
    }
  }
}

// Overwrite node velocity to get strain correct
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_rigid_mass_momentum_to_nodes() noexcept {
  if (this->material_id_ == 999 && this->set_contact_){  
    for (unsigned i = 0; i < std::pow(2, Tdim); ++i) {       
      nodes_[i]->update_rigid_mass_momentum(shapefn_[i] * this->mass_, 
                                  shapefn_[i] * this->mass_ * this->velocity_);
    }       
  }
}