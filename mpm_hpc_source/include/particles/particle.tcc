// Construct a particle with id and coordinates
template <unsigned Tdim>
mpm::Particle<Tdim>::Particle(Index id, const VectorDim& coord)
    : mpm::ParticleBase<Tdim>(id, coord) {
  this->initialise();
  // Clear cell ptr
  cell_ = nullptr;
  // Nodes
  nodes_.clear();
  // Set material pointer to null
  material_ = nullptr;
  // Logger
  std::string logger =
      "particle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

// Construct a particle with id, coordinates and status
template <unsigned Tdim>
mpm::Particle<Tdim>::Particle(Index id, const VectorDim& coord, bool status)
    : mpm::ParticleBase<Tdim>(id, coord, status) {
  this->initialise();
  cell_ = nullptr;
  nodes_.clear();
  material_ = nullptr;
  // Logger
  std::string logger =
      "particle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

// Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::Particle<Tdim>::initialise_particle(const HDF5Particle& particle) {
  
  // Particle ID
  this->id_ = particle.id;
  this->cell_id_ = particle.cell_id;
  this->material_id_ = particle.material_id;
  this->current_time_ = particle.current_time;
  this->status_ = particle.status;
  this->cell_ = nullptr;
  this->nodes_.clear();
  this->free_surface_ = 0;

  // Scaler properties
  this->mass_ = particle.mass;
  this->volume_ = particle.volume;
  this->density_ = particle.density;
  this->porosity_ = particle.porosity;
  this->solid_fraction_ = 1. - particle.porosity; 
  this->mass_density_ = particle.mass / particle.volume;
  this->volumetric_strain_ = particle.epsilon_v;
  this->deviatoric_strain_ = 0.;
  this->thermal_volumetric_strain_ = particle.thermal_epsilon_v;
  this->temperature_ = particle.temperature;
  this->PIC_temperature_ = particle.PIC_temperature;

  this->assign_initial_volume(particle.volume); 

  // Vector properties
  Eigen::Vector3d psize, coordinates, velocity, displacement;
  psize        << particle.nsize_x, particle.nsize_y, particle.nsize_z;
  coordinates  << particle.coord_x, particle.coord_y, particle.coord_z;
  velocity     << particle.velocity_x, particle.velocity_y, particle.velocity_z;
  displacement << particle.displacement_x, particle.displacement_y, particle.displacement_z;
  for (unsigned i = 0; i < Tdim; ++i) {
    this->natural_size_[i] = psize[i];
    this->coordinates_[i] = coordinates[i];
    this->displacement_[i] = displacement[i];
    this->velocity_[i] = velocity[i];
  } 

  // Tensor properties
  // Stress
  this->stress_[0] = particle.stress_xx;
  this->stress_[1] = particle.stress_yy;
  this->stress_[2] = particle.stress_zz;
  this->stress_[3] = particle.tau_xy;
  this->stress_[4] = particle.tau_yz;
  this->stress_[5] = particle.tau_xz;
  // Mechanical strain
  this->strain_[0] = particle.strain_xx;
  this->strain_[1] = particle.strain_yy;
  this->strain_[2] = particle.strain_zz;
  this->strain_[3] = particle.gamma_xy;
  this->strain_[4] = particle.gamma_yz;
  this->strain_[5] = particle.gamma_xz;
  // Thermal strain
  this->thermal_strain_[0] = particle.thermal_strain_xx;
  this->thermal_strain_[1] = particle.thermal_strain_yy;
  this->thermal_strain_[2] = particle.thermal_strain_zz;
  this->thermal_strain_[3] = particle.thermal_gamma_xy;
  this->thermal_strain_[4] = particle.thermal_gamma_yz;
  this->thermal_strain_[5] = particle.thermal_gamma_xz;
  // Deformation gradient
  this->deformation_gradient_(0, 0) = particle.fxx;
  if (Tdim != 1) {
    this->deformation_gradient_(0, 1) = particle.fxy;
    this->deformation_gradient_(1, 0) = particle.fyx;
    this->deformation_gradient_(1, 1) = particle.fyy;
  }
  if (Tdim == 3) {
    this->deformation_gradient_(0, 2) = particle.fxz;
    this->deformation_gradient_(1, 3) = particle.fyz;
    this->deformation_gradient_(2, 0) = particle.fzx;
    this->deformation_gradient_(2, 1) = particle.fzy;
    this->deformation_gradient_(2, 2) = particle.fzz;
  }
  return true;
}

// Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::Particle<Tdim>::initialise_particle(const HDF5Particle& particle,
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

// Return particle data in HDF5 format
template <unsigned Tdim>
mpm::HDF5Particle mpm::Particle<Tdim>::hdf5() {

  mpm::HDF5Particle particle_data;

  // Particle information
  particle_data.id = this->id_;
  particle_data.status = this->status_;  
  particle_data.cell_id = this->cell_id_;
  particle_data.material_id = this->material_id_;
  particle_data.current_time = this->current_time_;

  // Scaler properties
  particle_data.mass = this->mass_;
  particle_data.volume = this->volume_;
  particle_data.density = this->density_;
  particle_data.porosity = this->porosity_;
  particle_data.temperature = this->temperature_;
  particle_data.epsilon_v = this->volumetric_strain_;  
  particle_data.PIC_temperature = this->PIC_temperature_;
  particle_data.thermal_epsilon_v = this->thermal_volumetric_strain_;

  // Vector properties
  Eigen::Vector3d nsize, coordinates, velocity, displacement;
  nsize.setZero();
  coordinates.setZero();
  velocity.setZero();
  displacement.setZero();
  for (unsigned i = 0; i < Tdim; ++i) {
    nsize[i] = this->natural_size_[i];
    coordinates[i] = this->coordinates_[i];
    displacement[i] = this->displacement_[i];
    velocity[i] = this->velocity_[i];
  } 

  particle_data.coord_x = coordinates[0];
  particle_data.coord_y = coordinates[1];
  particle_data.coord_z = coordinates[2];

  particle_data.displacement_x = displacement[0];
  particle_data.displacement_y = displacement[1];
  particle_data.displacement_z = displacement[2];

  particle_data.nsize_x = nsize[0];
  particle_data.nsize_y = nsize[1];
  particle_data.nsize_z = nsize[2];

  particle_data.velocity_x = velocity[0];
  particle_data.velocity_y = velocity[1];
  particle_data.velocity_z = velocity[2];

  particle_data.rotation_xy = this->rotation_[0];
  particle_data.rotation_yz = this->rotation_[1];
  particle_data.rotation_xz = this->rotation_[2];  

  // Tensor properties
  particle_data.stress_xx = this->stress_[0];
  particle_data.stress_yy = this->stress_[1];
  particle_data.stress_zz = this->stress_[2];
  particle_data.tau_xy = this->stress_[3];
  particle_data.tau_yz = this->stress_[4];
  particle_data.tau_xz = this->stress_[5];

  particle_data.strain_xx = this->strain_[0];
  particle_data.strain_yy = this->strain_[1];
  particle_data.strain_zz = this->strain_[2];
  particle_data.gamma_xy = this->strain_[3];
  particle_data.gamma_yz = this->strain_[4];
  particle_data.gamma_xz = this->strain_[5];

  particle_data.thermal_strain_xx = this->thermal_strain_[0];
  particle_data.thermal_strain_yy = this->thermal_strain_[1];
  particle_data.thermal_strain_zz = this->thermal_strain_[2];
  particle_data.thermal_gamma_xy = this->thermal_strain_[3];
  particle_data.thermal_gamma_yz = this->thermal_strain_[4];
  particle_data.thermal_gamma_xz = this->thermal_strain_[5];

  Eigen::Matrix<double, 9, 1> deformation_gradient;
  deformation_gradient.setZero();
  for (unsigned j = 0; j < Tdim; ++j) {
    for (unsigned k = 0; k < Tdim; ++k) {
      deformation_gradient[3 * j + k] = this->deformation_gradient_(j, k);
    }
  }

  particle_data.fxx = deformation_gradient[0];
  particle_data.fxy = deformation_gradient[1];
  particle_data.fxz = deformation_gradient[2];
  particle_data.fyx = deformation_gradient[3];
  particle_data.fyy = deformation_gradient[4];
  particle_data.fyz = deformation_gradient[5];
  particle_data.fzx = deformation_gradient[6];
  particle_data.fzy = deformation_gradient[7];
  particle_data.fzz = deformation_gradient[8];

  particle_data.pressure =
      (state_variables_.find("pressure") != state_variables_.end())
          ? state_variables_.at("pressure")
          : 0.;

  // Write state variables
  if (material_ != nullptr) {
    particle_data.nstate_vars = state_variables_.size();
    if (state_variables_.size() > 20)
      throw std::runtime_error("# of state variables cannot be more than 20");
    unsigned i = 0;
    for (const auto& state_var : this->state_variables_) {
      particle_data.svars[i] = state_var.second;
      ++i;
    }
  }
  return particle_data;
}

// Initialise particle properties
template <unsigned Tdim>
void mpm::Particle<Tdim>::initialise() {

  // Scalar properties
  this->current_time_ = 0.; 
  this->volume_ = 0.;
  this->density_ = 1.;
  this->porosity_ = 0.;
  this->temperature_ = 0.;
  this->PIC_temperature_ = 0.;
  this->heat_capacity_= 0.; 
  this->mass_ = 1.;
  this->mass_density_ = 1.;
  this->solid_fraction_ = 1.;
  this->volumetric_strain_ = 0.;
  this->deviatoric_strain_ = 0.;
  this->thermal_volumetric_strain_ = 0.;
  this->dthermal_volumetric_strain_ = 0.;
  this->temperature_increment_ = 0.;
  this->heat_source_ = 0;
  this->temperature_acceleration_ = 0.;

  // Vector properties
  this->size_.setZero();  
  this->natural_size_.setZero(); 
  this->velocity_.setZero();
  // this->displacement_.setZero();
  this->displacement_.setZero();
  this->rotation_.setZero();
  this->temperature_gradient_.setZero();
  this->mass_gradient_.setZero();  
  this->contact_normal_.setZero();
  this->contact_tangential_.setZero();

  // Tensor properties
  this->strain_rate_.setZero();
  this->strain_.setZero();
  this->stress_.setZero();  
  this->thermal_strain_.setZero();
  this->heat_flux_.setZero();  
  this->fabric_CN_.setZero();
  this->fabric_PO_.setZero();
  this->C_matrix_.setZero();
  this->displacement_gradient_.setZero();   
  this->deformation_gradient_ = Eigen::Matrix<double, Tdim, Tdim>::Identity();

  // Bool properties
  this->set_traction_ = false;
  this->set_heat_source_ = false;
  this->set_contact_ = false;
  this->affine_mpm_ = false;

  // Scalar properties
  this->scalar_property_ = {
      {"ids",                    [&]() {return static_cast<double>(this->id_);} },
      {"current_time",           [&]() {return this->current_time_;}       },
      {"volumes",                [&]() {return this->volume_;}           },
      {"masses",                 [&]() {return this->mass_;}           },      
      {"porosities",             [&]() {return this->porosity_;}           },
      {"solid_fractions",        [&]() {return this->solid_fraction_;}     },      
      {"densities",              [&]() {return this->density_;}            },
      {"temperatures",           [&]() {return this->temperature_;}        },
      {"PIC_temperatures",       [&]() {return this->PIC_temperature_;}    },
      {"temperature_accelerations", [&]() {return this->temperature_acceleration_;}   },
      {"volumetric_strains",     [&]() {return this->volumetric_strain_;}  },
      {"dvolumetric_strains",     [&]() {return this->dvolumetric_strain_;}  },
      {"thermal_volumetric_strains",     [&]() {return this->thermal_volumetric_strain_;}  },
      {"dthermal_volumetric_strains",     [&]() {return this->dthermal_volumetric_strain_;}  },
      {"free_surfaces",          [&]() {return this->free_surface();}      } 
  };    
    
  this->vector_property_ = {    
      {"displacements",          [&]() {return this->displacement_;}       },     
      {"rotations",              [&]() {return this->rotation_;}           },
      {"velocities",             [&]() {return this->velocity_;}           },
      {"accelerations",          [&]() {return this->acceleration_;}           },
      {"stresses",               [&]() {return this->stress_;}             },
      {"strains",                [&]() {return this->strain_;}             },
      {"thermal_strains",        [&]() {return this->thermal_strain_;}     },
      {"temperature_gradients",  [&]() {return this->temperature_gradient_;}},
      {"outward_normals",         [&]() {return this->outward_normal_;}},             
      {"mass_gradients",         [&]() {return this->mass_gradient_;}},             
      {"heat_fluxes",            [&]() {return this->heat_flux_;}       },      
      {"fabric_CNs",             [&]() {return reshape_tensor(this->fabric_CN_);}},
      {"fabric_POs",             [&]() {return reshape_tensor(this->fabric_PO_);}},
      {"velocity_gradients",     [&]() {return reshape_tensor(this->velocity_gradient_);}},      
      {"displacement_gradients", [&]() {return reshape_tensor(this->displacement_gradient_);}},
      {"deformation_gradients",  [&]() {return reshape_tensor(this->deformation_gradient_);}}
  };
}

//==============================================================================
//  PART 0: ASSIGN INITIAL CONDITIONS
//==============================================================================

// Assign a material to particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_material(
    const std::shared_ptr<Material<Tdim>>& material) {
  bool status = false;
  try {
    // Check if material is valid and properties are set
    if (material != nullptr) {
      material_ = material;
      material_id_ = material_->id();
      state_variables_ = material_->initialise_state_variables();
      status = true;
    } else {
      throw std::runtime_error("Material is undefined!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

// Assign initial volume to the particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_initial_volume(double volume) {
  bool status = true;
  try {
    if (volume <= 0.)
      throw std::runtime_error("Particle volume cannot be negative");

    this->volume_ = volume;
    // Compute size of particle in each direction
    const double length =
        std::pow(this->volume_, static_cast<double>(1. / Tdim));
    // Set particle size as length on each side
    this->size_.fill(length);

    if (cell_ != nullptr) {
      // Get element ptr of a cell
      const auto element = cell_->element_ptr();

      // Set local particle size based on length of element in natural
      // coordinates (cpGIMP Bardenhagen 2008 (pp485))
      this->natural_size_.fill(
          element->unit_element_length() /
          std::pow(cell_->nparticles(), static_cast<double>(1. / Tdim)));
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute initial volume of the particle
template <unsigned Tdim>
void mpm::Particle<Tdim>::compute_volume(bool is_axisymmetric) noexcept {
  try {  
    // Check if particle has a valid cell ptr
    assert(cell_ != nullptr);
    // Volume of the cell / # of particles
    if (this->volume_ < 1E-16)
      this->assign_initial_volume(cell_->volume() / cell_->nparticles());

    // If axisymmetric condition, modify the volume according to ridial dimension
    if (is_axisymmetric) {
      this->is_axisymmetric_ = is_axisymmetric;
      volume_ *= coordinates_(0);
    };
  } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__, 
                      exception.what());
  }
}

// Compute mass of particle
template <unsigned Tdim>
void mpm::Particle<Tdim>::compute_mass() {
  // Check if particle volume is set and material ptr is valid
  assert(volume_ != std::numeric_limits<double>::max() && material_ != nullptr);

  // Assign intial porosity and density
  this->porosity_ = material_->template 
                          property<double>(std::string("porosity"));
  this->density_ = material_->template 
                          property<double>(std::string("density"));

  // Throw error if porosity is negative or larger than one 
  if (porosity_ < 0. || porosity_ > 1. || std::isnan(porosity_)) {
    throw std::runtime_error(
          "Particle porosity cannot be negative or larger than one");
    std::exit(1);
  }

  // Calculate mass density and mass
  this->solid_fraction_ = 1 - this->porosity_;
  this->mass_density_ = this->density_ * this->solid_fraction_;
  this->mass_ = this->volume_ * this->mass_density_; 
}

//==============================================================================
//  PART 1: LOCATE PARTICLES & COMPUTE SHAPE FUNCTIONS
//==============================================================================

// Assign a cell to particle
// If point is in new cell, assign new cell and remove particle id from old
// cell. If point can't be found in the new cell, check if particle is still
// valid in the old cell, if it is leave it as is. If not, set cell as null
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_cell(
    const std::shared_ptr<Cell<Tdim>>& cellptr) {
  bool status = true;
  try {
    Eigen::Matrix<double, Tdim, 1> xi;
    // Assign cell to the new cell ptr, if point can be found in new cell
    if (cellptr->is_point_in_cell(this->coordinates_, &xi)) {
      // if a cell already exists remove particle from that cell
      if (cell_ != nullptr) cell_->remove_particle_id(this->id_);

      cell_ = cellptr;
      cell_id_ = cellptr->id();
      //dn_dx centroid
      dn_dx_centroid_ = cell_->dn_dx_centroid();
      // shapefn centroid
      shapefn_centroid_ = cell_->shapefn_centroid();
      // Copy nodal pointer to cell
      nodes_.clear();
      nodes_ = cell_->nodes();

      // Compute reference location of particle
      bool xi_status = this->compute_reference_location();
      if (!xi_status) return false;
      status = cell_->add_particle_id(this->id());
    } else {
      throw std::runtime_error("Point cannot be found in cell!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign a cell to particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_cell_xi(
    const std::shared_ptr<Cell<Tdim>>& cellptr,
    const Eigen::Matrix<double, Tdim, 1>& xi) {
  bool status = true;
  try {
    // Assign cell to the new cell ptr, if point can be found in new cell
    if (cellptr != nullptr) {
      // if a cell already exists remove particle from that cell
      if (cell_ != nullptr) cell_->remove_particle_id(this->id_);

      cell_ = cellptr;
      cell_id_ = cellptr->id(); 
      // dn_dx centroid
      dn_dx_centroid_ = cell_->dn_dx_centroid();
      // shapefn centroid
      shapefn_centroid_ = cell_->shapefn_centroid();      
      // Copy nodal pointer to cell
      nodes_.clear();
      nodes_ = cell_->nodes();

      // Assign the reference location of particle
      bool xi_nan = false;

      // Check if point is within the cell
      for (unsigned i = 0; i < xi.size(); ++i)
        if (xi[i] < -1. || xi[i] > 1. || std::isnan(xi[i])) xi_nan = true;

      if (xi_nan == false)
        this->xi_ = xi;
      else
        return false;

      status = cell_->add_particle_id(this->id());
    } else {
      throw std::runtime_error("Point cannot be found in cell!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign a cell id to particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_cell_id(mpm::Index id) {
  bool status = false;
  try {
    // if a cell ptr is null
    if (cell_ == nullptr && id != std::numeric_limits<Index>::max()) {
      cell_id_ = id;
      status = true;
    } else {
      throw std::runtime_error("Invalid cell id or cell is already assigned!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Remove cell for the particle
template <unsigned Tdim>
void mpm::Particle<Tdim>::remove_cell() {
  // if a cell is not nullptr
  if (cell_ != nullptr) cell_->remove_particle_id(this->id_);
  cell_id_ = std::numeric_limits<Index>::max();
  // Clear all the nodes
  nodes_.clear();
}

// Compute reference location cell to particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::compute_reference_location() noexcept {
  // Set status of compute reference location
  bool status = false;
  // Compute local coordinates
  Eigen::Matrix<double, Tdim, 1> xi;
  // Check if the point is in cell
  if (cell_ != nullptr && cell_->is_point_in_cell(this->coordinates_, &xi)) {
    this->xi_ = xi;
    status = true;
  }

  return status;
}

///! Compute shape functions and gradients
template <unsigned Tdim>
void mpm::Particle<Tdim>::compute_shapefn() noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  // Get element ptr of a cell
  const auto element = cell_->element_ptr();

  // Zero matrix
  Eigen::Matrix<double, Tdim, 1> zero = Eigen::Matrix<double, Tdim, 1>::Zero();

  // Compute shape function of the particle
  shapefn_ = element->shapefn(this->xi_, this->natural_size_, zero);

  // Compute dN/dx
  dn_dx_ = element->dn_dx(this->xi_, cell_->nodal_coordinates(),
                          this->natural_size_, zero);

}

//==============================================================================
//  PART 2: MAP PARTICLE INFORMATION TO NODES
//==============================================================================

// Map particle mass and momentum to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_mass_momentum_to_nodes() noexcept {
  // Check if particle mass is set
  assert(mass_ != std::numeric_limits<double>::max());
  if (this->material_id_ != 999) {
    // Map mass and momentum to nodes
    if (!affine_mpm_) {
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        nodes_[i]->update_mass(true, mpm::ParticlePhase::Solid,
                              mass_ * shapefn_[i]);
        nodes_[i]->update_momentum(true, mpm::ParticlePhase::Solid,
                                  mass_ * shapefn_[i] * velocity_);
      }
    } else {
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        Eigen::Matrix<double, Tdim, 1> affine_vel = velocity_ + 
                  C_matrix_ * (nodes_[i]->coordinates() - this->coordinates_);
        nodes_[i]->update_mass(true, mpm::ParticlePhase::Solid,
                              mass_ * shapefn_[i]);
        nodes_[i]->update_momentum(true, mpm::ParticlePhase::Solid,
                                  mass_ * shapefn_[i] * affine_vel);
      }
    }

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   std::cout << "    solid_velocity_: " << velocity_ << std::endl;
  //   std::cout << "    solid_mass_: " << mass_ << std::endl;
  // }

  }  
}

// Map body force to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_external_force(const VectorDim& pgravity) {
  if (material_id_ != 999) {
    // Compute nodal body forces
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_external_force(true, mpm::ParticlePhase::Solid,
                                       (pgravity * mass_ * shapefn_[i]));
  }
}

// Map internal force 1D to nodes
template <>
inline void mpm::Particle<1>::map_internal_force() {
  if (material_id_ != 999) {
    // Compute nodal internal forces
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      // Compute force: -pstress * volume
      Eigen::Matrix<double, 1, 1> force;
      force[0] = -1. * dn_dx_(i, 0) * volume_ * stress_[0];

      nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Solid, force);
    }
  }
}

//! Map internal force 2D
template <>
inline void mpm::Particle<2>::map_internal_force() {
  if (material_id_ != 999) {
    // Compute nodal internal forces
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      // Compute force: -pstress * volume
      Eigen::Matrix<double, 2, 1> force;
      force[0] = dn_dx_(i, 0) * stress_[0] + dn_dx_(i, 1) * stress_[3];
      force[1] = dn_dx_(i, 1) * stress_[1] + dn_dx_(i, 0) * stress_[3];

      // force[0] = (dn_dx_(i, 0) + (dn_dx_centroid_(i, 0) - dn_dx_(i, 0)) / 2.) * stress_[0] +
      //           (dn_dx_centroid_(i, 0) - dn_dx_(i, 0)) / 2. * stress_[1] +
      //           dn_dx_(i, 1) * stress_[3];
      // force[1] = (dn_dx_centroid_(i, 1) - dn_dx_(i, 1)) / 2. * stress_[0] +
      //           (dn_dx_(i, 1) + (dn_dx_centroid_(i, 1) - dn_dx_(i, 1)) / 2.) * stress_[1] +
      //           dn_dx_(i, 0) * stress_[3];

      if (is_axisymmetric_) force[0] += shapefn_[i]/this->coordinates_(0) * stress_[2]; 

      force *= -1. * this->volume_;

      nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Solid, force);

    }    
  }
}

// Map internal force 3D to nodes
template <>
inline void mpm::Particle<3>::map_internal_force() {
  if (material_id_ != 999) {
    // Compute nodal internal forces
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      // Compute force: -pstress * volume
      Eigen::Matrix<double, 3, 1> force;
      force[0] = dn_dx_(i, 0) * stress_[0] + dn_dx_(i, 1) * stress_[3] +
                 dn_dx_(i, 2) * stress_[5];

      force[1] = dn_dx_(i, 1) * stress_[1] + dn_dx_(i, 0) * stress_[3] +
                 dn_dx_(i, 2) * stress_[4];

      force[2] = dn_dx_(i, 2) * stress_[2] + dn_dx_(i, 1) * stress_[4] +
                 dn_dx_(i, 0) * stress_[5];

    // force[0] = (dn_dx_(i, 0) + (dn_dx_centroid_(i, 0) - dn_dx_(i, 0)) / 3.) * stress_[0] +
    //            (dn_dx_centroid_(i, 0) - dn_dx_(i, 0)) / 3. * stress_[1] +
    //            (dn_dx_centroid_(i, 0) - dn_dx_(i, 0)) / 3. * stress_[2] +
    //            dn_dx_(i, 1) * stress_[3] + dn_dx_(i, 2) * stress_[5];

    // force[1] = (dn_dx_centroid_(i, 1) - dn_dx_(i, 1)) / 3. * stress_[0] +
    //            (dn_dx_(i, 1) + (dn_dx_centroid_(i, 1) - dn_dx_(i, 1)) / 3.) * stress_[1] +
    //            (dn_dx_centroid_(i, 1) - dn_dx_(i, 1)) / 3. * stress_[2] +
    //            dn_dx_(i, 0) * stress_[3] + dn_dx_(i, 2) * stress_[4];

    // force[2] = (dn_dx_centroid_(i, 2) - dn_dx_(i, 2)) / 3. * stress_[0] +
    //            (dn_dx_centroid_(i, 2) - dn_dx_(i, 2)) / 3. * stress_[1] +
    //            (dn_dx_(i, 2) + (dn_dx_centroid_(i, 2) - dn_dx_(i, 2)) / 3.) * stress_[2] +
    //            dn_dx_(i, 1) * stress_[4] + dn_dx_(i, 0) * stress_[5];

      force *= -1. * this->volume_;

      nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Solid, force);
    }
  }
}

// Map traction force to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_traction_force() noexcept {
  if (this->set_traction_) {
    // Map particle traction forces to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_external_force(true, mpm::ParticlePhase::Solid,
                                       (shapefn_[i] * traction_));
  }
}

// Map particle heat capacity and heat to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_heat_to_nodes() {
  if (material_id_ != 999) {  
    // Check if mass is set and positive
    assert(mass_ != std::numeric_limits<double>::max());
    // get the specific_heat 
    const double specific_heat_ = 
          material_->template property<double>(std::string("specific_heat"));

    if (this->material_id_ != 999) {
      // Map heat capacity and heat to nodes
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        nodes_[i]->update_heat_capacity(true, mpm::ParticlePhase::Solid,
                                    mass_ * specific_heat_ * shapefn_[i]);
        nodes_[i]->update_heat(true, mpm::ParticlePhase::Solid,
                      mass_ * specific_heat_ * shapefn_[i] * temperature_);
      }
  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   //std::cout << "particle id: " << this->id_ << std::endl;
  //   std::cout << "    solid_heat_capacity_: " << mass_ * specific_heat_ << std::endl;
  //   //std::cout << "    temperature_: " << this->temperature_ << std::endl;
  // }
    }
  }
}

// Map heat conduction to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_heat_conduction() {
  if (material_id_ != 999) {
    
    // Assign the thermal conductivity
    const double k_conductivity = 
            material_->template property<double>(std::string("thermal_conductivity"));
    
    // Assign temperature gradient
    this->compute_temperature_gradient(mpm::ParticlePhase::Solid);

    // Compute nodal heat conduction
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      double heat_conduction = 0;
      for (unsigned j = 0; j < Tdim; ++j){
        heat_conduction += dn_dx_(i, j) * this->temperature_gradient_[j]; 
      }

      heat_conduction *= -1 * this->volume_ * (1 - porosity_) * k_conductivity;
      nodes_[i]->update_heat_conduction(true, mpm::ParticlePhase::Solid, heat_conduction);
    }

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   std::cout << "    Sheat_conduction: " << -1 * this->volume_ * (1 - porosity_) * k_conductivity << std::endl;
  //   //std::cout << "    volume_: " <<this->volume_;
  //   std::cout << "    solid_fraction_: " << (1 - porosity_) << std::endl;
  // }

  }
}

// Map heat conduction to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_virtual_heat_flux(bool convective, 
                                          const double para_1,
                                          const double para_2) {
  if (material_id_ != 999) {
    
    this->compute_mass_gradient(mpm::ParticlePhase::Solid);
    this->outward_normal_.setZero();
    this->free_surface_ = false;    
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      if (nodes_[i]->free_surface()) {
        this->outward_normal_ = -mass_gradient_/mass_gradient_.norm();
        this->free_surface_ = true;

        // this->outward_normal_[0] = -1;
        // this->outward_normal_[1] = 0; 
      } 
    }

    for (unsigned i = 0; i < nodes_.size(); ++i) {
      // Judge the heat flux type
      // Type 1: convective heat flux 
      if (convective) {
        const double heat_transfer_coeff = para_1;
        const double ambient_temperature = para_2;
        // Only the particles related surface node matter
        if (nodes_[i]->free_surface()) {
          double heat_conduction_bc = 0;
          for (unsigned j = 0; j < Tdim; ++j){
            heat_conduction_bc += dn_dx_(i, j) * heat_transfer_coeff * (ambient_temperature - 
                                    this->temperature_) * outward_normal_[j];
          }
            heat_conduction_bc *= this->volume_ * (1. - porosity_);
            nodes_[i]->update_heat_conduction(true, mpm::ParticlePhase::Solid, heat_conduction_bc);
        }
      } 
      // Type 2: conductive heat flux
      else {
        const double const_heat_flux = para_1;
        // Only the particles related surface node matter
        if (nodes_[i]->free_surface()) {
          double heat_conduction_bc = 0;
          for (unsigned j = 0; j < Tdim; ++j){
            heat_conduction_bc += dn_dx_(i, j) * const_heat_flux * outward_normal_[j];
          }
          heat_conduction_bc *= this->volume_ * (1. - porosity_);
          nodes_[i]->update_heat_conduction(true, mpm::ParticlePhase::Solid, heat_conduction_bc);
        }
      }

      // if (nodes_[i]->free_surface()) {
      //   double heat_conduction_bc = 0;
      //   double heat_flux_bc = material_->template property<double>(std::string("heat_flux_bc"));
      //   heat_conduction_bc = shapefn_[i] * heat_flux_bc * (50 - this->temperature_) / std::pow(volume_, 1/Tdim);              
      //   heat_conduction_bc *= -1 * this->volume_ * (1 - porosity_);
      //   nodes_[i]->update_heat_conduction(true, mpm::ParticlePhase::Solid, heat_conduction_bc);
      // }

      // // Particle boundary
      // if (this->free_surface_) {
      //   double heat_conduction_bc = 0;
      //   double heat_flux_bc = material_->template property<double>(std::string("heat_flux_bc"));
      //   heat_conduction_bc = shapefn_[i] * heat_flux_bc / this->size_(0);         
      //   // heat_conduction_bc = shapefn_[i] * heat_flux_bc * (50 - this->temperature_) / std::pow(volume_, 1/Tdim);              
      //   heat_conduction_bc *= -1 * this->volume_ * (1 - porosity_);
      //   nodes_[i]->update_heat_conduction(true, mpm::ParticlePhase::Solid, heat_conduction_bc);
      // }

    }
  }
}

// Map heat source to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_heat_source() {
  if (material_id_ != 999) {  
    if (this->set_heat_source_) {
      // Map particle heat source forces to nodes
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        nodes_[i]->update_heat_source(true, mpm::ParticlePhase::Solid,
                                        (this->volume_ * shapefn_[i] * heat_source_));                                 
      }
    }
  }
}

// Map plastic work to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_plastic_work(double dt) noexcept {
  if (material_id_ != 999) {  
      // Map particle heat source forces to nodes
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        nodes_[i]->update_plastic_work(true, mpm::ParticlePhase::Solid,
                                        (this->volume_ * shapefn_[i] * plastic_work_));                              
    }
  }
}

//==============================================================================
//  PART 3: UPDATE PARTICLE INFORMATION
//==============================================================================

// Compute strain rate of the particle 1D
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Particle<1>::compute_strain_rate(
    const Eigen::MatrixXd& dn_dx, Eigen::VectorXd& shapefn, 
    unsigned phase) noexcept {
  // Define strain rate
  Eigen::Matrix<double, 6, 1> strain_rate = Eigen::Matrix<double, 6, 1>::Zero();

  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    Eigen::Matrix<double, 1, 1> vel = nodes_[i]->velocity(phase);
    strain_rate[0] += dn_dx(i, 0) * vel[0];
  }

  if (std::fabs(strain_rate(0)) < 1.E-15) strain_rate[0] = 0.;
  return strain_rate;
}

// Compute strain rate of the particle 2D
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Particle<2>::compute_strain_rate(
    const Eigen::MatrixXd& dn_dx, Eigen::VectorXd& shapefn, 
    unsigned phase) noexcept {
  // Define strain rate
  Eigen::Matrix<double, 6, 1> strain_rate = Eigen::Matrix<double, 6, 1>::Zero();

  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    Eigen::Matrix<double, 2, 1> vel = nodes_[i]->velocity(phase);
    strain_rate[0] += dn_dx(i, 0) * vel[0];
    strain_rate[1] += dn_dx(i, 1) * vel[1];
    strain_rate[3] += dn_dx(i, 1) * vel[0] + dn_dx(i, 0) * vel[1];
    if (is_axisymmetric_) strain_rate[2] += shapefn[i] * vel[0] / this->coordinates_(0);
  } 

  if (std::fabs(strain_rate[0]) < 1.E-15) strain_rate[0] = 0.;
  if (std::fabs(strain_rate[1]) < 1.E-15) strain_rate[1] = 0.;
  if (std::fabs(strain_rate[2]) < 1.E-15) strain_rate[2] = 0.;
  if (std::fabs(strain_rate[3]) < 1.E-15) strain_rate[3] = 0.;  

  return strain_rate;
}

// Compute strain rate of the particle 3D
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Particle<3>::compute_strain_rate(
    const Eigen::MatrixXd& dn_dx, Eigen::VectorXd& shapefn, 
    unsigned phase) noexcept {
  // Define strain rate
  Eigen::Matrix<double, 6, 1> strain_rate = Eigen::Matrix<double, 6, 1>::Zero();

  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    Eigen::Matrix<double, 3, 1> vel = nodes_[i]->velocity(phase);
    strain_rate[0] += dn_dx(i, 0) * vel[0];
    strain_rate[1] += dn_dx(i, 1) * vel[1];
    strain_rate[2] += dn_dx(i, 2) * vel[2];
    strain_rate[3] += dn_dx(i, 1) * vel[0] + dn_dx(i, 0) * vel[1];
    strain_rate[4] += dn_dx(i, 2) * vel[1] + dn_dx(i, 1) * vel[2];
    strain_rate[5] += dn_dx(i, 2) * vel[0] + dn_dx(i, 0) * vel[2];
  }

  for (unsigned i = 0; i < strain_rate.size(); ++i)
    if (std::fabs(strain_rate[i]) < 1.E-15) strain_rate[i] = 0.;
  return strain_rate;
}

// Compute jaumann stress
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Particle<2>::compute_jaumann_stress() noexcept {
  if (material_id_ != 999) {
    // Check if material ptr is valid
    assert(material_ != nullptr);

    // stress tensor
    Eigen::Matrix<double, 3, 3> stress_tensor = Eigen::Matrix<double, 3, 3>::Zero();    
    stress_tensor(0,0) = this->stress_[0];
    stress_tensor(0,1) = this->stress_[3];
    stress_tensor(0,2) = this->stress_[4];
    stress_tensor(1,0) = this->stress_[3];
    stress_tensor(1,1) = this->stress_[1]; 
    stress_tensor(1,2) = this->stress_[5];
    stress_tensor(2,0) = this->stress_[4];
    stress_tensor(2,1) = this->stress_[5];
    stress_tensor(2,2) = this->stress_[2];
      
    // Define spin tensor
    Eigen::Matrix<double, 3, 3> spin_tensor = Eigen::Matrix<double, 3, 3>::Zero();
    if (!affine_mpm_) {
      for (unsigned i = 0; i < this->nodes_.size(); ++i) {
        Eigen::Matrix<double, 2, 1> vel = nodes_[i]->velocity(0);
        spin_tensor(0,0) += 0;
        spin_tensor(1,1) += 0;
        spin_tensor(2,2) += 0;
        spin_tensor(0,1) += (dn_dx_(i, 1) * vel[0] - dn_dx_(i, 0) * vel[1])/2;
        spin_tensor(1,0) += (dn_dx_(i, 0) * vel[1] - dn_dx_(i, 1) * vel[0])/2;
      }      
    } else {
      Eigen::Matrix<double, 2, 2> spin_rate_tensor = Eigen::Matrix<double, 2, 2>::Zero();
      spin_rate_tensor = (C_matrix_ - C_matrix_.transpose()) / 2;
      spin_tensor(0,0) = 0;
      spin_tensor(1,1) = 0;
      spin_tensor(2,2) = 0;
      spin_tensor(0,1) = spin_rate_tensor(0,1);
      spin_tensor(1,0) = spin_rate_tensor(1,0);
    }

    // Define Jaumann stress
    Eigen::Matrix<double, 3, 3> jaumann_stress_tensor = Eigen::Matrix<double, 3, 3>::Zero();
    jaumann_stress_tensor = stress_tensor * spin_tensor - spin_tensor * stress_tensor;

    Eigen::Matrix<double, 6, 1> jaumann_stress = Eigen::Matrix<double, 6, 1>::Zero(); 
    jaumann_stress[0] = jaumann_stress_tensor(0,0);
    jaumann_stress[1] = jaumann_stress_tensor(1,1); 
    jaumann_stress[2] = jaumann_stress_tensor(2,2); 
    jaumann_stress[3] = jaumann_stress_tensor(0,1);  
    jaumann_stress[4] = jaumann_stress_tensor(0,2);
    jaumann_stress[5] = jaumann_stress_tensor(1,2);

    return jaumann_stress; 
  }
}

// Compute spin tensor
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Particle<3>::compute_jaumann_stress() noexcept {
  if (material_id_ != 999) {
    // Check if material ptr is valid
    assert(material_ != nullptr);

    // stress tensor
    Eigen::Matrix<double, 3, 3> stress_tensor = Eigen::Matrix<double, 3, 3>::Zero();    
    stress_tensor(0,0) = this->stress_[0];
    stress_tensor(0,1) = this->stress_[3];
    stress_tensor(0,2) = this->stress_[4];
    stress_tensor(1,0) = this->stress_[3];
    stress_tensor(1,1) = this->stress_[1]; 
    stress_tensor(1,2) = this->stress_[5];
    stress_tensor(2,0) = this->stress_[4];
    stress_tensor(2,1) = this->stress_[5];
    stress_tensor(2,2) = this->stress_[2];
      
    // Define spin tensor
    Eigen::Matrix<double, 3, 3> spin_tensor = Eigen::Matrix<double, 3, 3>::Zero();
    for (unsigned i = 0; i < this->nodes_.size(); ++i) {
      Eigen::Matrix<double, 3, 1> vel = nodes_[i]->velocity(0);
      spin_tensor(0,0) += 0;
      spin_tensor(1,1) += 0;
      spin_tensor(2,2) += 0;
      spin_tensor(0,1) += (dn_dx_(i, 1) * vel[0] - dn_dx_(i, 0) * vel[1])/2;
      spin_tensor(1,0) += (dn_dx_(i, 0) * vel[1] - dn_dx_(i, 1) * vel[0])/2;
      spin_tensor(0,2) += (dn_dx_(i, 2) * vel[0] - dn_dx_(i, 0) * vel[2])/2;    
      spin_tensor(1,2) += (dn_dx_(i, 2) * vel[1] - dn_dx_(i, 1) * vel[2])/2;    
      spin_tensor(2,0) += (dn_dx_(i, 0) * vel[2] - dn_dx_(i, 2) * vel[0])/2;
      spin_tensor(2,1) += (dn_dx_(i, 1) * vel[2] - dn_dx_(i, 2) * vel[1])/2;         
    }  
    // Define Jaumann stress
    Eigen::Matrix<double, 3, 3> jaumann_stress_tensor = Eigen::Matrix<double, 3, 3>::Zero();
    jaumann_stress_tensor = stress_tensor * spin_tensor - spin_tensor * stress_tensor;

    Eigen::Matrix<double, 6, 1> jaumann_stress = Eigen::Matrix<double, 6, 1>::Zero(); 
    jaumann_stress[0] = jaumann_stress_tensor(0,0);
    jaumann_stress[1] = jaumann_stress_tensor(1,1);
    jaumann_stress[2] = jaumann_stress_tensor(2,2);
    jaumann_stress[3] = jaumann_stress_tensor(0,1);
    jaumann_stress[4] = jaumann_stress_tensor(0,2);
    jaumann_stress[5] = jaumann_stress_tensor(1,2);

    return jaumann_stress; 
  }
}

// Compute temperature gradient of the particle
template <unsigned Tdim>
inline Eigen::Matrix<double, Tdim, 1> mpm::Particle<Tdim>::compute_temperature_gradient(unsigned phase) noexcept {

  Eigen::Matrix<double, Tdim, 1> temperature_gradient;
  temperature_gradient.setZero();
  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    double temperature = nodes_[i]->temperature(phase);
    for (unsigned j = 0; j < Tdim; ++j) {
      // temperature_gradient = partial T / partial X = T_{i,j}
      temperature_gradient[j] += dn_dx_(i, j) * temperature;
      if (std::fabs(temperature_gradient[j]) < 1.E-15)
        temperature_gradient[j] = 0.;
    }
  }
  this->temperature_gradient_ = temperature_gradient;

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   std::cout << "    phase: " << phase << std::endl;
  //   std::cout << "    temperature: " << nodes_[0]->temperature(phase) << std::endl;
  //   std::cout << "    temperature_gradient_y: " << temperature_gradient[1] << std::endl;
  // }

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   std::cout << "    temperature: " << nodes_[0]->temperature(mpm::ParticlePhase::Solid) << std::endl;
  // }

  return temperature_gradient_;
}

// Compute mass gradient of the particle
template <unsigned Tdim>
inline Eigen::Matrix<double, Tdim, 1> mpm::Particle<Tdim>::compute_mass_gradient(unsigned phase) noexcept {

  Eigen::Matrix<double, Tdim, 1> mass_gradient;
  mass_gradient.setZero();
  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    double mass = nodes_[i]->mass(phase);
    for (unsigned j = 0; j < Tdim; ++j) {
      // mass_gradient = partial T / partial X = m_{i,j}
      mass_gradient[j] += dn_dx_(i, j) * mass;
      if (std::fabs(mass_gradient[j]) < 1.E-15)
        mass_gradient[j] = 0.;
    }
  }
  this->mass_gradient_ = mass_gradient;
  return mass_gradient_;
}

// Compute updated velocity and position of the particle
template <unsigned Tdim>
void mpm::Particle<Tdim>::compute_updated_velocity(
    double dt, double pic, double damping_factor) {

  if (material_id_ == 999) {

    auto displacement_increment = this->velocity_ * dt;
    this->coordinates_ += displacement_increment;
    this->displacement_ += displacement_increment;

  } else {
    // Check if particle has a valid cell ptr and pic value
    assert(cell_ != nullptr);
    assert((-1.E-15) <= pic && pic <= (1 + 1.E-15));

    // Get PIC velocity
    Eigen::Matrix<double, Tdim, 1> pic_velocity =
        Eigen::Matrix<double, Tdim, 1>::Zero();
    for (unsigned i = 0; i < nodes_.size(); ++i)
      pic_velocity +=
          shapefn_[i] * nodes_[i]->velocity(mpm::ParticlePhase::Solid);

    // Get interpolated nodal acceleration
    Eigen::Matrix<double, Tdim, 1> nodal_acceleration =
        Eigen::Matrix<double, Tdim, 1>::Zero();
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodal_acceleration +=
          shapefn_[i] * nodes_[i]->acceleration(mpm::ParticlePhase::Solid);

    // Applying particle damping
    nodal_acceleration -= damping_factor * this->velocity_;

    // The damping acceleration above contributes to the FLIP update, but it
    // would otherwise be discarded by a pure PIC update (pic == 1). Apply the
    // same first-order damping increment to the PIC velocity so that Cundall
    // damping remains effective during PIC equilibrium initialization.
    pic_velocity -= damping_factor * this->velocity_ * dt;

    this->acceleration_ = nodal_acceleration;

    // Get FLIP velocity
    Eigen::Matrix<double, Tdim, 1> flip_velocity =
        this->velocity_ + nodal_acceleration * dt;

    // Update particle velocity based on PIC value
    this->velocity_ = pic * pic_velocity + (1. - pic) * flip_velocity;

    // Displacement increment
    // Use mid-point scheme
    Eigen::Matrix<double, Tdim, 1> displacement_increment =
        (pic_velocity - 0.5 * nodal_acceleration * dt) * dt;
    // auto displacement_increment = pic_velocity * dt;

    if (affine_mpm_) {
      C_matrix_ = this->compute_affine_matrix(shapefn_, mpm::ParticlePhase::Solid);
      displacement_increment = this->velocity_ * dt;
    }
    
    this->coordinates_ += displacement_increment;
    this->displacement_ += displacement_increment;

    // // if (this->current_time_ < (2 * dt)) {
    //   if (Tdim == 2) {
    //     Eigen::Matrix<double, 2, 2> rotation;
    //     rotation.setZero();
    //     // const double PI = std::atan(1.0) * 4;

    //     const double angle = material_->template property<double>(std::string("rotation_angle")) / 180. * 3.1415926535 * dt;
    //     rotation(0,0) = cos(angle);
    //     rotation(0,1) = -sin(angle); 
    //     rotation(1,0) = sin(angle);
    //     rotation(1,1) = cos(angle);
    //     Eigen::Matrix<double, 2, 1> coordinates;
    //     coordinates[0] = rotation(0,0) * this->coordinates_[0] + rotation(0,1) * this->coordinates_[1];
    //     coordinates[1] = rotation(1,0) * this->coordinates_[0] + rotation(1,1) * this->coordinates_[1];
    //     if (id_ == 0) 
    //       std::cout <<  this->coordinates_ << "\n";           
    //     this->coordinates_[0] = coordinates[0];
    //     this->coordinates_[1] = coordinates[1];        
    //     if (id_ == 0) 
    //       std::cout << rotation << "\t" 
    //                 << this->coordinates_ << "\n";                  
    //     // std::cout <<  cos(0 / 180 * 3.1415926535) << "\n";             
    //   // }
    // }
  }
}

// Compute updated temperature of the particle
template <unsigned Tdim>
void mpm::Particle<Tdim>::update_particle_temperature(
    double dt, double pic_t) noexcept {
  if (material_id_ != 999) {  
    
  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   //std::cout << "    temperature_acceleration_: " << this->temperature_acceleration_ << std::endl;
  //   //std::cout << "    temperature_gradient_y: " << this->temperature_gradient_[1] << std::endl;
  //   //std::cout << "    temperature_acceleration: " << nodes_[0]->temperature_acceleration(mpm::ParticlePhase::Solid) << std::endl;
  //   std::cout << "    temperature_: " << this->temperature_ << std::endl;
  //   std::cout << "    before:temperature_increment_cent_: " << this->temperature_increment_cent_ << std::endl;
  // }    

    // Check if particle has a valid cell ptr and pic_t value
    assert(cell_ != nullptr);
    assert((-1.E-15) <= pic_t && pic_t <= (1 + 1.E-15));

    // Get PIC temperature
    double PIC_temperature = 0;
    double temperature_acceleration_cent = 0;
    for (unsigned i = 0; i < nodes_.size(); ++i)
      PIC_temperature +=
          shapefn_[i] * nodes_[i]->temperature(mpm::ParticlePhase::Solid);
    // temperature increment
    this->temperature_increment_ = PIC_temperature - this->PIC_temperature_; 

    this->PIC_temperature_= PIC_temperature;
    // Get interpolated nodal temperature acceleration
    this->temperature_acceleration_ = 0.;
    for (unsigned i = 0; i < nodes_.size(); ++i){
      this->temperature_acceleration_ +=
          shapefn_[i] * nodes_[i]->temperature_acceleration(mpm::ParticlePhase::Solid);
      temperature_acceleration_cent +=
          shapefn_centroid_[i] * nodes_[i]->temperature_acceleration(mpm::ParticlePhase::Solid);
    }

    // Get FLIP temperature
    double FLIP_temperature =
        this->temperature_ + this->temperature_acceleration_ * dt;
    this->FLIP_temperature_ = FLIP_temperature;

    // Update particle temperature based on PIC value
    this->temperature_ = pic_t * PIC_temperature + (1. - pic_t) * FLIP_temperature;

    // temperature increment
    this->temperature_increment_cent_ = temperature_acceleration_cent * dt; 
  }
  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   //std::cout << "particle id: " << this->id_ << std::endl;
  //   //std::cout << "    temperature_acceleration_: " << this->temperature_acceleration_ << std::endl;
  //   //std::cout << "    temperature_gradient_y: " << this->temperature_gradient_[1] << std::endl;
  //   //std::cout << "    temperature_acceleration: " << nodes_[0]->temperature_acceleration(mpm::ParticlePhase::Solid) << std::endl;
  //   std::cout << "    temperature_: " << this->temperature_ << std::endl;
  //   std::cout << "    after:temperature_increment_cent_: " << this->temperature_increment_cent_ << std::endl;
  // }
}

// Update particle strain
template <unsigned Tdim>
void mpm::Particle<Tdim>::update_particle_strain(double dt) noexcept {
  if (material_id_ != 999) {

    // Assign strain rate
    strain_rate_ = this->compute_strain_rate(dn_dx_, shapefn_, mpm::ParticlePhase::Solid);
    // Update dstrain
    dstrain_ = strain_rate_ * dt;
    // Update strain
    strain_ += dstrain_;
    // // Update displacement gradient
    // displacement_gradient_ = velocity_gradient_ * dt;
    // // Update deformation gradient
    // deformation_gradient_ *= (Eigen::Matrix<double, Tdim, Tdim>::Identity() +
    //                           displacement_gradient_);

    // Compute at centroid
    // Strain rate for reduced integration
    const Eigen::Matrix<double, 6, 1> strain_rate_centroid =
        this->compute_strain_rate(dn_dx_centroid_, shapefn_centroid_, mpm::ParticlePhase::Solid);

    // Assign volumetric strain at centroid (infinitesimal assumption)
    dvolumetric_strain_ = dt * strain_rate_.head(Tdim).sum();
    volumetric_strain_ += dvolumetric_strain_;

    if (is_axisymmetric_) dvolumetric_strain_ = dt * strain_rate_centroid.head(3).sum();

    // jaumann_stress_ = this->compute_jaumann_stress() * dt;  
  }
}

// Compute thermal strain of the particle
template <unsigned Tdim>
void mpm::Particle<Tdim>::update_particle_thermal_strain() noexcept {
  if (material_id_ != 999) {
    // get the thermal conductivity coefficient
    const double beta_solid =
      material_->template property<double>(std::string("thermal_expansivity"));

    // compute thermal strain increment
    for (unsigned i = 0; i < 3; i++) {
      dthermal_strain_[i] = -1 * beta_solid * this->temperature_increment_;
    }

    // Compute volumetric thermal strain
    dthermal_volumetric_strain_ = dthermal_strain_.head(Tdim).sum();
    if (is_axisymmetric_) dthermal_volumetric_strain_ = dthermal_strain_.head(3).sum();    
    
    // update thermal strain 
    thermal_strain_ += dthermal_strain_;
    thermal_volumetric_strain_ += dthermal_volumetric_strain_;

    // compute total strain increment
    dstrain_ += dthermal_strain_;
    // compute total strain
    strain_ += dthermal_strain_;
    // compute total volumetric strain increment
    dvolumetric_strain_ += dthermal_volumetric_strain_;
    // compute total volumetric strain
    volumetric_strain_ += dthermal_volumetric_strain_;

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   //std::cout << "    dstrain_: " << dstrain_ << std::endl;
  //   std::cout << "    temperature_increment_: " << this->temperature_increment_ << std::endl;
  //   std::cout << "    dvolumetric_strain_: " << dvolumetric_strain_ << std::endl;
  // }

  }

}

// Update particle stress
template <unsigned Tdim>
void mpm::Particle<Tdim>::update_particle_stress() noexcept {
  if (material_id_ != 999) {
    // Check if material ptr is valid
    assert(material_ != nullptr);
    // Calculate stress    
    Eigen::Matrix<double, 6, 1> updated_stress = 
          material_->compute_stress(stress_, dstrain_, this, &state_variables_);
    this->stress_ = updated_stress;

    // get theta
    const double theta = 
        material_->template property<double>(std::string("theta"));  

    // // this->stress_ += jaumann_stress_;
    this->plastic_work_ = strain_rate_(0) * stress_(0) +
                          strain_rate_(1) * stress_(1) +
                          strain_rate_(2) * stress_(2) +
                          2 * strain_rate_(3) * stress_(3) +
                          2 * strain_rate_(4) * stress_(4) +
                          2 * strain_rate_(5) * stress_(5);

    this->plastic_work_ *= theta;
  }
}

// Compute updated update porosity of the particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::update_particle_porosity(double dt) {
  bool status = true;
  if (material_id_ != 999) {
    try {
      // Update particle porosity
      const double porosity = 
          1 - (1 - this->porosity_) / (1 + dvolumetric_strain_);
      // Check if the value is valid
      if (porosity > 0 && porosity < 1) {
        this->porosity_ = porosity;
        this->solid_fraction_ = 1.0 - this->porosity_;
      }
      // Throw error if porosity is negative or larger than one 
      else {
        throw std::runtime_error(
              "Particle porosity cannot be negative or larger than one");
        std::exit(1);
      }

    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

// Update particle volume of the particle 
template <unsigned Tdim>
void mpm::Particle<Tdim>::update_particle_volume()  {
  if (material_id_ != 999) {
    // Check if particle has a valid cell ptr and a valid volume
    assert(cell_ != nullptr && volume_ != std::numeric_limits<double>::max());
    // Compute at centroid
    // Strain rate for reduced integration
    this->volume_ *= (1. + dvolumetric_strain_);
    this->mass_density_ = this->mass_density_ / (1. + dvolumetric_strain_);
  } 
}

// update material density of particle (considering thermal expansion)
template <unsigned Tdim>
void mpm::Particle<Tdim>::update_particle_density(double dt) {
  // Check if material ptr is valid
  assert(material_ != nullptr);
  if (material_id_ != 999) {
    // get initial density
    double density_0 = 
        material_->template property<double>(std::string("density"));
    // get the thermal conductivity coefficient
    double beta_solid =
        material_->template property<double>(std::string("thermal_expansivity"));
    // material density is a function of temperature
    this->density_ = density_0 * std::exp(-3 * beta_solid * PIC_temperature_);                
    // this->density_ *= std::exp(-1 * 3 * beta_solid * this->temperature_increment_);                 
  }
}

//==============================================================================
//  APPENDIX 2: APPLY BOUNDARY CONDITIONS
//==============================================================================

// Assign traction to the particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_particle_traction(unsigned direction, double traction) {
  bool status = false;
  try {
    if (direction >= Tdim ||
        this->volume_ == std::numeric_limits<double>::max()) {
      throw std::runtime_error(
          "Particle traction property: volume / direction is invalid");
    }
    // Assign traction
    traction_(direction) = traction * this->volume_ / this->size_(direction);
    status = true;
    this->set_traction_ = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute free surface by density method
template <unsigned Tdim>
bool mpm::Particle<Tdim>::compute_particle_free_surface() {
  bool status = true;
  try {
    this->free_surface_ = false;
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr) {
      // 1. Simple approach of density comparison (Hamad, 2015)
      // Get interpolated nodal density
      double nodal_mass_density = 0;

      for (unsigned i = 0; i < nodes_.size(); ++i)
        nodal_mass_density +=
            shapefn_[i] * nodes_[i]->density(mpm::ParticlePhase::Solid);

    // get the density_ratio
    const double density_ratio =
      material_->template property<double>(std::string("density_ratio"));
      // Compare smoothen density to actual particle mass density
      if ((nodal_mass_density / mass_density_ ) <= density_ratio){
        if (cell_->free_surface()) this->free_surface_ = true;
      }
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
};

//! Assign free surface manually
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_particle_free_surfaces() {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr) {
      this->free_surface_ = false;
      if (this->initial_free_surface() && !this->initial_nonfree_surface_)
        this->free_surface_ = true;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign contact to the particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_particle_contact(unsigned dir, double normal) {
  bool status = false;
  try {
    // Assign contact
    contact_normal_(dir) = normal;
    status = true;
    this->set_contact_ = true;
    for (unsigned i = 0; i < nodes_.size(); ++i) {       
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
void mpm::Particle<Tdim>::map_moving_rigid_velocity_to_nodes(
    unsigned dir, double velocity, double dt) noexcept {   
  if (material_id_ == 999) {
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      nodes_[i]->assign_velocity_from_rigid(dir, velocity, dt);
    }
  }
}

// Overwrite node velocity to get strain correct
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_rigid_mass_momentum_to_nodes() noexcept {
  if (this->material_id_ == 999){  
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      nodes_[i]->update_rigid_mass_momentum(shapefn_[i] * this->mass_, 
                                  shapefn_[i] * this->mass_ * velocity_);
    }
  }
}

// Assign heat source to the particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_particle_heat_source(double heat_source, double dt) {
  bool status = false;
  try {
    // Assign heat source
    heat_source_ = heat_source;
    status = true;
    this->set_heat_source_ = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//==============================================================================
//  APPENDIX 3: SMOOTHING METHOD
//==============================================================================

// Map particle pressure to nodes 
template <unsigned Tdim>
bool mpm::Particle<Tdim>::map_pressure_to_nodes(double current_time) noexcept {
  // Mass is initialized
  assert(mass_ != std::numeric_limits<double>::max());

  bool status = false;
  // Check if particle mass is set and state variable pressure is found
  if (mass_ != std::numeric_limits<double>::max() &&
      (state_variables_.find("pressure") != state_variables_.end())) {
    // Map particle pressure to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_mass_pressure(
          mpm::ParticlePhase::Solid,
          shapefn_[i] * mass_ * state_variables_.at("pressure"), current_time);

    status = true;
  }
  return status;
}

// Compute pressure smoothing of the particle based on nodal pressure 
template <unsigned Tdim>
bool mpm::Particle<Tdim>::compute_pressure_smoothing() noexcept {
  // Assert
  assert(cell_ != nullptr);

  bool status = false;
  // Check if particle has a valid cell ptr
  if (cell_ != nullptr &&
      (state_variables_.find("pressure") != state_variables_.end())) {

    double pressure = 0.;
    // Update particle pressure to interpolated nodal pressure
    for (unsigned i = 0; i < this->nodes_.size(); ++i)
      pressure += shapefn_[i] * nodes_[i]->pressure(mpm::ParticlePhase::Solid);

    state_variables_.at("pressure") = pressure;
    status = true;
  }
  return status;

}

//==============================================================================
//  APPENDIX 4: MULTISCALE FUNCTIONS
//==============================================================================

// // Compute velocity gradient of the particle 1D
// template <>
// inline Eigen::Matrix<double, 1, 1> mpm::Particle<1>::compute_velocity_gradient(
//     const Eigen::MatrixXd& dn_dx, unsigned phase) noexcept {
//   // Define strain rate
//   Eigen::Matrix<double, 1, 1> velocity_gradient =
//       Eigen::Matrix<double, 1, 1>::Zero();

//   for (unsigned i = 0; i < this->nodes_.size(); ++i) {
//     Eigen::Matrix<double, 1, 1> vel = nodes_[i]->velocity(phase);
//     // vel_gradient = partial x / partial X
//     velocity_gradient[0] +=dn_dx(i, 0) * vel[0];
//   }

//   if (std::fabs(velocity_gradient[0]) < 1.E-15) velocity_gradient[0] = 0.;
//   return velocity_gradient;
// }

// // Compute velocity gradient of the particle 2D
// template <>
// inline Eigen::Matrix<double, 2, 2> mpm::Particle<2>::compute_velocity_gradient(
//     const Eigen::MatrixXd& dn_dx, unsigned phase) noexcept {
//   Eigen::Matrix<double, 2, 2> velocity_gradient =
//       Eigen::Matrix<double, 2, 2>::Zero();

//   for (unsigned i = 0; i < this->nodes_.size(); ++i) {
//     Eigen::Matrix<double, 2, 1> vel = nodes_[i]->velocity(phase);
//     // vel_gradient = partial v / partial X = v_{i,j}
//     velocity_gradient(0, 0) +=dn_dx(i, 0) * vel[0];
//     velocity_gradient(0, 1) +=dn_dx(i, 1) * vel[0];
//     velocity_gradient(1, 0) +=dn_dx(i, 0) * vel[1];
//     velocity_gradient(1, 1) +=dn_dx(i, 1) * vel[1];
//   }

//   if (std::fabs(velocity_gradient(0, 0)) < 1.E-15) velocity_gradient(0, 0) = 0.;
//   if (std::fabs(velocity_gradient(0, 1)) < 1.E-15) velocity_gradient(0, 1) = 0.;
//   if (std::fabs(velocity_gradient(1, 0)) < 1.E-15) velocity_gradient(1, 0) = 0.;
//   if (std::fabs(velocity_gradient(1, 1)) < 1.E-15) velocity_gradient(1, 0) = 0.;

//   return velocity_gradient;
// }

// // Compute velocity gradient of the particle 3D
// template <>
// inline Eigen::Matrix<double, 3, 3> mpm::Particle<3>::compute_velocity_gradient(
//     const Eigen::MatrixXd& dn_dx, unsigned phase) noexcept {
//   // Define strain rate
//   Eigen::Matrix<double, 3, 3> velocity_gradient =
//       Eigen::Matrix<double, 3, 3>::Zero();

//   for (unsigned i = 0; i < this->nodes_.size(); ++i) {
//     Eigen::Matrix<double, 3, 1> vel = nodes_[i]->velocity(phase);
//     // vel_gradient = partial v / partial X = v_{i,j}
//     velocity_gradient(0, 0) +=dn_dx(i, 0) * vel[0];
//     velocity_gradient(0, 1) +=dn_dx(i, 1) * vel[0];
//     velocity_gradient(0, 2) +=dn_dx(i, 2) * vel[0];
//     velocity_gradient(1, 0) +=dn_dx(i, 0) * vel[1];
//     velocity_gradient(1, 1) +=dn_dx(i, 1) * vel[1];
//     velocity_gradient(1, 2) +=dn_dx(i, 2) * vel[1];
//     velocity_gradient(2, 0) +=dn_dx(i, 0) * vel[2];
//     velocity_gradient(2, 1) +=dn_dx(i, 1) * vel[2];
//     velocity_gradient(2, 2) +=dn_dx(i, 2) * vel[2];
//   }

//   for (unsigned i = 0; i < 3; ++i) {
//     for (unsigned j = 0; j < 3; ++j) {
//       if (std::fabs(velocity_gradient(i, j)) < 1.E-15)
//         velocity_gradient(i, j) = 0.;
//     }
//   }
//   return velocity_gradient;
// }

// Compute strain of the particle
template <unsigned Tdim>
void mpm::Particle<Tdim>::compute_displacement_gradient(double dt, bool thermal) noexcept {

  if (material_id_ != 999) {
    // Assign strain rate
    auto velocity_gradient =
        this->compute_velocity_gradient(dn_dx_, mpm::ParticlePhase::Solid);

    // Update displacement gradient
    displacement_gradient_ = velocity_gradient * dt;

    if (thermal) {
      // get the thermal conductivity coefficient
      const double beta_solid =
      material_->template property<double>(std::string("thermal_expansivity"));
      // compute thermal strain increment
      for (unsigned i = 0; i < Tdim; i++) {
        displacement_gradient_(i,i) += -1 * beta_solid * this->temperature_increment_;   
      }
    }
    
    // Update deformation gradient
    deformation_gradient_ *= (Eigen::Matrix<double, Tdim, Tdim>::Identity() +
                              displacement_gradient_);

    // auto strain_tensor =
    //     0.5 * (displacement_gradient_ + displacement_gradient_.transpose());

    // strain_[0] = strain_tensor(0, 0);
    // if (Tdim == 2) {
    //   strain_[1] = strain_tensor(1, 1);
    //   strain_[3] = 2 * strain_tensor(0, 1);
    // } else if (Tdim == 3) {
    //   strain_[1] = strain_tensor(1, 1);
    //   strain_[2] = strain_tensor(2, 2);
    //   strain_[3] = 2 * strain_tensor(0, 1);
    //   strain_[4] = 2 * strain_tensor(1, 2);
    //   strain_[5] = 2 * strain_tensor(0, 2);
    // }
  }
}

// Compute stress
template <unsigned Tdim>
void mpm::Particle<Tdim>::set_stress(const Eigen::MatrixXd& stresses,
                                     bool increment) noexcept {
  if (material_id_ != 999) {
    // set stress
    if (increment) {
      this->stress_ += stresses.col(id_);
    } else {
      this->stress_ = stresses.col(id_);
    }
    // get theta
    const double theta = 
        material_->template property<double>(std::string("theta"));  

    // // this->stress_ += jaumann_stress_;
    this->plastic_work_ = strain_rate_(0) * stress_(0) +
                          strain_rate_(1) * stress_(1) +
                          strain_rate_(2) * stress_(2) +
                          2 * strain_rate_(3) * stress_(3) +
                          2 * strain_rate_(4) * stress_(4) +
                          2 * strain_rate_(5) * stress_(5);

    this->plastic_work_ *= theta;;
  }
}

// Set fabric
template <unsigned Tdim>
void mpm::Particle<Tdim>::set_fabric(std::string fabric_type,
                                     const Eigen::MatrixXd& fabrics) {
  if (material_id_ != 999) {
    Eigen::MatrixXd temp_fabric;
    temp_fabric = fabrics.col(id_);
    temp_fabric.resize(Tdim, Tdim);
    if (fabric_type == "CN") {
      fabric_CN_ = temp_fabric;
    } else if (fabric_type == "PO") {
      fabric_PO_ = temp_fabric;
    } else {
      throw std::runtime_error("Fabric type is not properly assigned!\n");
    }
  }
}

// Get fabric
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, Tdim> mpm::Particle<Tdim>::fabric(
    std::string fabric_type) const {
  Eigen::Matrix<double, Tdim, Tdim> fabric;
  fabric.setZero();
  if (material_id_ != 999 && fabric_type == "CN") fabric = fabric_CN_;
  if (material_id_ != 999 && fabric_type == "PO") fabric = fabric_PO_;
  return fabric;
}

//==============================================================================
//  APPENDIX 5: AFFINE material point method
//==============================================================================
// Compute velocity gradient
template <unsigned Tdim>
inline Eigen::Matrix<double, Tdim, Tdim> mpm::Particle<Tdim>::compute_velocity_gradient(
    const Eigen::MatrixXd& dn_dx, unsigned phase) noexcept {
  // Define strain rate
  Eigen::Matrix<double, Tdim, Tdim> velocity_gradient;
  velocity_gradient.setZero();
  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    Eigen::Matrix<double, Tdim, 1> vel = nodes_[i]->velocity(phase);
    // vel_gradient = partial x / partial X
    velocity_gradient += (dn_dx.row(i).transpose()) * (vel.transpose());
  }
  for (unsigned i = 0; i < Tdim; ++i) {
    for (unsigned j = 0; j < Tdim; ++j) {
      if (std::fabs(velocity_gradient(i, j)) < 1.E-15)
        velocity_gradient(i, j) = 0.;
    }
  }
  return velocity_gradient;
}

// Compute affine matrix 
template <unsigned Tdim>
inline Eigen::Matrix<double, Tdim, Tdim> mpm::Particle<Tdim>::compute_affine_matrix(
    const Eigen::VectorXd& shapefn, unsigned phase) noexcept {
  // Compute B matrix
  Eigen::Matrix<double, Tdim, Tdim> B_matrix, D_matrix, C_matrix;
  B_matrix.setZero();
  D_matrix.setZero();
  C_matrix.setZero(); 
  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    Eigen::Matrix<double, Tdim, 1> vel = nodes_[i]->velocity(phase);
    Eigen::Matrix<double, Tdim, 1> rel_coord = nodes_[i]->coordinates() - this->coordinates_;
    B_matrix += shapefn[i] * vel * rel_coord.transpose(); 
    D_matrix += shapefn[i] * rel_coord * rel_coord.transpose(); 
  }
  C_matrix = B_matrix * D_matrix.inverse();
  return C_matrix;
}

//==============================================================================
//  APPENDIX 5: UNUSED FUNCTIONS
//==============================================================================

// Assign material id of this particle to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::append_material_id_to_nodes() const {
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->append_material_id(material_id_);
}

// Add a neighbour particle and return the status of addition
template <unsigned Tdim>
bool mpm::Particle<Tdim>::add_neighbour(mpm::Index neighbour_id) {
  bool insertion_status = false;
  try {
    // If particle id is not the same as the current particle
    if (neighbour_id != this->id())
      insertion_status = (neighbours_.insert(neighbour_id)).second;
    else
      throw std::runtime_error("Invalid local id of a neighbour particle");

  } catch (std::exception& exception) {
    console_->error("{} {}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return insertion_status;
}

// Assign neighbour particles
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_neighbours(
    const std::vector<mpm::Index>& neighbours_set) {
  bool status = true;
  try {
    // Remove the existing neighbours' id
    neighbours_.clear();

    // Loop over neighbours and check if the id is the same the current particle
    // id
    for (const auto neighbour_id : neighbours_set) {
       if (neighbour_id != this->id()) {
        neighbours_.insert(neighbour_id);
       }    
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}
