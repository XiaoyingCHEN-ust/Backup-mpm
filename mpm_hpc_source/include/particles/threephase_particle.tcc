// Construct a two phase particle with id and coordinates
template <unsigned Tdim>
mpm::ThreePhaseParticle<Tdim>::ThreePhaseParticle(Index id, const VectorDim& coord)
    : mpm::Particle<Tdim>(id, coord) {

    this->initialise_liquid_phase();

    // Set material pointer to null
    liquid_material_ = nullptr;
    // Logger
    std::string logger = "ThreePhaseParticle" + std::to_string(Tdim) + "d::" + std::to_string(id);
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//==============================================================================
// ASSIGN INITIAL CONDITIONS
//==============================================================================

// Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::ThreePhaseParticle<Tdim>::initialise_particle(
    const HDF5Particle& particle) {
    // Derive from particle
    mpm::Particle<Tdim>::initialise_particle(particle);

    // MIXTURE data
    this->liquid_material_id_ = particle.liquid_material_id;
    this->pore_pressure_ = particle.pore_pressure;
    // this->suction_pressure_ = particle.suction_pressure;
    // this->suction_pressure_ = particle.suction_pressure;

    // LIQUID PHASE data
    this->liquid_velocity_[0] = particle.liquid_velocity_x;
    this->liquid_velocity_[1] = particle.liquid_velocity_y;
    this->liquid_velocity_[2] = particle.liquid_velocity_z;
    this->liquid_pressure_ = particle.liquid_pressure;
    this->PIC_liquid_pressure_ = this->liquid_pressure_;
    this->liquid_saturation_ = particle.liquid_saturation;
    this->liquid_fraction_ = particle.liquid_fraction;
    this->liquid_density_ = particle.liquid_density;
    // this->liquid_volume_ = particle.liquid_volume;
    this->liquid_mass_ = particle.liquid_mass;
    this->liquid_permeability_ = particle.liquid_permeability;

    // GAS PHASE data
    this->gas_velocity_[0] = particle.gas_velocity_x;
    this->gas_velocity_[1] = particle.gas_velocity_y;
    this->gas_velocity_[2] = particle.gas_velocity_z;
    this->gas_pressure_ = particle.gas_pressure;
    this->PIC_gas_pressure_ = this->gas_pressure_;
    this->gas_saturation_ = particle.gas_saturation;
    this->gas_fraction_ = particle.gas_fraction;
    this->gas_density_ = particle.gas_density;
    // this->gas_volume_ = particle.gas_volume;
    this->gas_mass_ = particle.gas_mass;
    this->gas_permeability_ = particle.gas_permeability;

    return true;
}

// Initialise particle HDF5 data and material
template <unsigned Tdim>
bool mpm::ThreePhaseParticle<Tdim>::initialise_particle(
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
void mpm::ThreePhaseParticle<Tdim>::initialise_liquid_phase() {
    // MIXTURE
    set_mixture_traction_ = false;
    pore_pressure_ = 0.;
    suction_pressure_ = 0.;
    // The active linear retention law uses dSw/dpw = 1E-6 1/Pa.
    // Initialise it before the first pore-pressure update.
    dSw_dpw_ = 1E-6;
    set_pressure_constraint_ = false;

    // Liquid properties
    liquid_velocity_.setZero();
    liquid_flux_.setZero();
    liquid_strain_.setZero();
    liquid_saturation_ = 0.;
    liquid_fraction_ = 0.;
    liquid_density_ = 0.;
    liquid_volume_ = 0.;
    liquid_mass_ = 0.;
    liquid_mass_density_ = 0.;
    liquid_pressure_ = 0.;
    liquid_pressure_acceleration_ = 0.;
    liquid_volumetric_strain_ = 0.;
    liquid_permeability_ = 1.;
    liquid_source_ = 0.;

    // Gas properties
    gas_velocity_.setZero();
    gas_flux_.setZero();
    gas_strain_.setZero();
    gas_saturation_ = 0.;
    gas_fraction_ = 0.;
    gas_density_ = 0.;
    gas_volume_ = 0.;
    gas_mass_ = 0.;
    gas_mass_density_ = 0.;
    gas_pressure_ = 0.;
    gas_pressure_acceleration_ = 0.;
    gas_volumetric_strain_ = 0.;
    gas_permeability_ = 1.;
    gas_source_ = 0.;

    // Link data with NAME
    this->scalar_property_.insert({
      {"pore_pressures",         [&]() {return this->pore_pressure_;}},
      {"suction_pressures",      [&]() {return this->suction_pressure_;}},
      {"liquid_pressures",       [&]() {return this->liquid_pressure_;}},
      {"PIC_liquid_pressures",   [&]() {return this->PIC_liquid_pressure_;}},
      {"liquid_saturations",     [&]() {return this->liquid_saturation_;}},
      {"liquid_fractions",       [&]() {return this->liquid_fraction_;}},
      {"liquid_chis",            [&]() {return this->liquid_chi_;}},
      {"liquid_densities",       [&]() {return this->liquid_density_;}},
      {"liquid_sources",         [&]() {return this->liquid_source_;}},
      {"liquid_permeabilities",   [&]() {return this->liquid_permeability_;}},
      {"liquid_volumes",         [&]() {return this->liquid_volume_;}},
      {"liquid_masses",          [&]() {return this->liquid_mass_;}},
      {"liquid_volumetric_strains", [&]() {return this->liquid_volumetric_strain_;}},
      {"liquid_vol_strains",     [&]() {return this->liquid_volumetric_strain_;}},
      {"liquid_viscosities",     [&]() {return this->liquid_viscosity_;}},
      {"liquid_critical_times",  [&]() {return this->liquid_critical_time_;}},
      {"gas_pressures",          [&]() {return this->gas_pressure_;}},
      {"PIC_gas_pressures",      [&]() {return this->PIC_gas_pressure_;}},
      {"gas_saturations",        [&]() {return this->gas_saturation_;}},
      {"gas_fractions",          [&]() {return this->gas_fraction_;}},
      {"gas_densities",          [&]() {return this->gas_density_;}}, 
      {"gas_sources",            [&]() {return this->gas_source_;}},
      {"gas_permeabilities",      [&]() {return this->gas_permeability_;}},
      {"gas_volumes",            [&]() {return this->gas_volume_;}},
      {"gas_masses",             [&]() {return this->gas_mass_;}},
      {"gas_mass_densities",     [&]() {return this->gas_mass_density_;}}, 
      {"gas_vol_strains",        [&]() {return this->gas_volumetric_strain_;}},
      {"gas_viscosities",        [&]() {return this->gas_viscosity_;}},
      {"gas_critical_times",     [&]() {return this->gas_critical_time_;}}
    });

    this->vector_property_.insert({
      {"liquid_velocities",      [&]() {return this->liquid_velocity_;}},
      {"liquid_accelerations",   [&]() {return this->liquid_acceleration_;}},
      {"liquid_strains",         [&]() {return this->liquid_strain_rate_;}},
      {"liquid_fluxes",          [&]() {return this->liquid_flux_;}},
      {"liquid_pressure_gradients",[&]() {return this->liquid_pressure_gradient_;}},
      {"gas_velocities",         [&]() {return this->gas_velocity_;}},
      {"gas_accelerations",      [&]() {return this->gas_acceleration_;}},
      {"gas_strains",            [&]() {return this->gas_strain_rate_;}},
      {"gas_fluxes",             [&]() {return this->gas_flux_;}},
      {"gas_pressure_gradients",[&]() {return this->gas_pressure_gradient_;}},
      {"K_matrix",               [&]() {return this->K_matrix_;}}
    }); 
}

// Return particle data in HDF5 format
template <unsigned Tdim>
mpm::HDF5Particle mpm::ThreePhaseParticle<Tdim>::hdf5() {
    // Derive from particle
    auto particle_data = mpm::Particle<Tdim>::hdf5();

    // MIXTURE data
    particle_data.liquid_material_id = this->liquid_material_id_;
    particle_data.pore_pressure = this->pore_pressure_;
    // particle_data.suction_pressure = this->suction_pressure_;

    // LIQUID PHASE data
    particle_data.liquid_velocity_x = this->liquid_velocity_[0];
    particle_data.liquid_velocity_y = this->liquid_velocity_[1];
    particle_data.liquid_velocity_z = this->liquid_velocity_[2];
    particle_data.liquid_pressure = this->liquid_pressure_;
    particle_data.liquid_saturation = this->liquid_saturation_;
    particle_data.liquid_fraction = this->liquid_fraction_;
    particle_data.liquid_density = this->liquid_density_;
    // particle_data.liquid_volume = this->liquid_volume_;
    particle_data.liquid_mass = this->liquid_mass_;
    particle_data.liquid_permeability = this->liquid_permeability_;

    // GAS PHASE data
    particle_data.gas_velocity_x = this->gas_velocity_[0];
    particle_data.gas_velocity_y = this->gas_velocity_[1];
    particle_data.gas_velocity_z = this->gas_velocity_[2];
    particle_data.gas_pressure = this->gas_pressure_;
    particle_data.gas_saturation = this->gas_saturation_;
    particle_data.gas_fraction = this->gas_fraction_;
    particle_data.gas_density = this->gas_density_;
    // particle_data.gas_volume = this->gas_volume_;
    particle_data.gas_mass = this->gas_mass_;
    particle_data.gas_permeability = this->gas_permeability_;

    return particle_data;
}

// Assign a liquid material to particle
template <unsigned Tdim>
bool mpm::ThreePhaseParticle<Tdim>::assign_liquid_material(
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
    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
    return status;
}

// Compute mass of particle (both solid and fluid)
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::compute_mass() {
  mpm::Particle<Tdim>::compute_mass();
  this->assign_initial_properties();
}

// Assign initial properties to particle
template <unsigned Tdim>
bool mpm::ThreePhaseParticle<Tdim>::assign_initial_properties() {
  bool status = true;
  try {
    this->ini_porosity_ = material_->template property<double>(std::string("porosity"));
    this->porosity_ = ini_porosity_; 
    this->intrinsic_permeability_ = 
          material_->template property<double>(std::string("intrinsic_permeability"));
    this->liquid_permeability_ = this->intrinsic_permeability_;
    this->gas_permeability_ = this->intrinsic_permeability_;

    // SOLID PAHSE
    // Constant propeties
    this->solid_thermal_conductivity_ = material_->template
                  property<double>(std::string("thermal_conductivity"));
    this->solid_specific_heat_ = material_->template 
                  property<double>(std::string("specific_heat"));
                  
    this->solid_expansivity_ = material_->template
                  property<double>(std::string("thermal_expansivity"));

    // LIQUID PHASE
    // Constant propeties
    this->liquid_thermal_conductivity_ = liquid_material_->template
                  property<double>(std::string("liquid_thermal_conductivity"));

    this->liquid_specific_heat_ = liquid_material_->template 
                  property<double>(std::string("liquid_specific_heat"));
    this->liquid_expansivity_ = liquid_material_->template 
                  property<double>(std::string("liquid_expansivity"));
    this->liquid_compressibility_ = liquid_material_->template 
                  property<double>(std::string("liquid_compressibility"));
    this->liquid_molar_mass_ = liquid_material_->template 
                  property<double>(std::string("liquid_molar_mass"));
    this->liquid_saturation_res_ = liquid_material_->template 
                  property<double>(std::string("liquid_saturation_res"));
    // Time-dependent propeties
    ini_liquid_density_ = liquid_material_->template
                  property<double>(std::string("density"));
    ini_liquid_saturation_ = liquid_material_->template
                  property<double>(std::string("liquid_saturation"));
    ini_liquid_viscosity_ = liquid_material_->template
                  property<double>(std::string("liquid_viscosity"));
    this->liquid_density_ = ini_liquid_density_;
    this->liquid_saturation_ = ini_liquid_saturation_;
    this->liquid_viscosity_ = ini_liquid_viscosity_;

    // GAS PHASE
    // Constant propeties
    this->gas_thermal_conductivity_ = liquid_material_->template
                  property<double>(std::string("gas_thermal_conductivity"));
    this->gas_specific_heat_ = liquid_material_->template
                  property<double>(std::string("gas_specific_heat"));
    this->gas_constant_ = liquid_material_->template
                  property<double>(std::string("gas_constant"));
    this->gas_molar_mass_ = liquid_material_->template
                  property<double>(std::string("gas_molar_mass"));
    this->gas_saturation_res_ = liquid_material_->template
                  property<double>(std::string("gas_saturation_res"));

    // Time-dependent propeties
    ini_gas_saturation_ = liquid_material_->template
                  property<double>(std::string("gas_saturation"));
    ini_gas_viscosity_ = liquid_material_->template
                  property<double>(std::string("gas_viscosity"));
    double p_ref = material_->template property<double>(std::string("p_ref"));
    this->gas_density_ = gas_molar_mass_ * (gas_pressure_ + p_ref) /
                          gas_constant_ / (PIC_temperature_ + 273.15);
    this->gas_saturation_ = ini_gas_saturation_;
    this->gas_viscosity_ = ini_gas_viscosity_;

    // SOLID PHASE
    this->solid_heat_capacity_ = mass_ * solid_specific_heat_;

    // LIQUID PHASE
    this->ini_liquid_fraction_ = porosity_ * liquid_saturation_;
    this->liquid_fraction_ = this->ini_liquid_fraction_;
    this->liquid_volume_ = liquid_fraction_ * volume_;
    this->liquid_mass_ = liquid_volume_ * liquid_density_;
    this->liquid_mass_density_ = liquid_fraction_ * liquid_density_;
    this->liquid_chi_ = liquid_saturation_ / (liquid_saturation_ + gas_saturation_);
    this->liquid_heat_capacity_ = liquid_mass_ * liquid_specific_heat_;

    // GAS PHASE
    this->ini_gas_fraction_ = porosity_ * gas_saturation_;
    this->gas_fraction_ = this->ini_gas_fraction_;
    this->gas_volume_ = gas_fraction_ * volume_;
    this->gas_mass_ = gas_volume_ * gas_density_;
    this->gas_mass_density_ = gas_fraction_ * gas_density_;
    this->gas_chi_ = 1 - liquid_chi_;
    this->gas_heat_capacity_ = gas_mass_ * gas_specific_heat_;

    // MIXTURE
    this->mixture_mass_ = mass_ + liquid_mass_ + gas_mass_;

    // Initial pore presssure
    const double para_p0 = liquid_material_->template 
                  property<double>(std::string("para_p0"));
    const double para_m = liquid_material_->template 
                  property<double>(std::string("para_m"));
    this->effective_saturation_ = (liquid_saturation_ - this->liquid_saturation_res_) /
                                  (1 - this->liquid_saturation_res_);
    // this->suction_pressure_ = para_p0 * std::pow(effective_saturation_, -1. / para_m);
    // this->suction_pressure_ = para_p0 * std::pow(std::pow(effective_saturation_, -1. / para_m) - 1., 1. - para_m);
    this->suction_pressure_ = (1 - ini_liquid_saturation_) * 1E6;
    this->gas_pressure_ = p_ref;
    this->liquid_pressure_ = this->gas_pressure_ - this->suction_pressure_;
    this->PIC_gas_pressure_ = this->gas_pressure_;
    this->PIC_liquid_pressure_ = this->liquid_pressure_;
    this->pore_pressure_ = liquid_saturation_ * PIC_liquid_pressure_ +
                          gas_saturation_ * PIC_gas_pressure_;

    this->ini_gas_pressure_ = this->gas_pressure_;
    this->ini_liquid_pressure_ = this->liquid_pressure_;
    this->ini_pore_pressure_ = this->pore_pressure_;

  } catch (std::exception& exception) {
    console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                    exception.what());
    status = false;
  }
  return status;
}

//==============================================================================
// MAP PARTICLE INFORMATION TO NODES
//==============================================================================

// Map particle mass and momentum to nodes
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::map_mass_momentum_to_nodes() noexcept {
  if (this->material_id_ != 999) {
    // SOLID PHASE
    mpm::Particle<Tdim>::map_mass_momentum_to_nodes();

    for (unsigned i = 0; i < nodes_.size(); ++i) {
      // Reduce the cost of accessing the shapefn_ array repeatedly
      double shapefn_i = shapefn_[i];

      // LIQUID PHASE
      double liquid_mass = liquid_mass_ * shapefn_i;
      nodes_[i]->update_mass_momentum(true, mpm::ParticlePhase::Liquid,
                                            liquid_mass, 
                                            liquid_mass * liquid_velocity_);
      // GAS PHASE
      double gas_mass = gas_mass_ * shapefn_i;
      nodes_[i]->update_mass_momentum(true, mpm::ParticlePhase::Gas,
                                            gas_mass,
                                            gas_mass * gas_velocity_);
    }
  }
}

// Map particle external force = body force + traction force
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::map_external_force(const VectorDim& pgravity) {
  if (this->material_id_ != 999) {
    try {
      this->pgravity_ = pgravity;
      Eigen::Matrix<double, Tdim, 1> mixture_force, liquid_force, gas_force;
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        
        // External force = body force + boundary traction
        // MIXTURE
        mixture_force.setZero();
        mixture_force = pgravity * mixture_mass_ * shapefn_[i] +
                        this->mixture_traction_ * shapefn_[i];
        nodes_[i]->update_external_force(true, mpm::ParticlePhase::Mixture,
                                                mixture_force);

        // LIQUID PHASE
        liquid_force.setZero(); 
        liquid_force = pgravity * liquid_mass_ * shapefn_[i] +
                      liquid_traction_ * shapefn_[i];
        nodes_[i]->update_external_force(true, mpm::ParticlePhase::Liquid,
                                          liquid_force);

        // GAS PHASE
        gas_force.setZero();
        gas_force = pgravity * gas_mass_ * shapefn_[i] +
                    gas_traction_ * shapefn_[i];
        nodes_[i]->update_external_force(true, mpm::ParticlePhase::Gas, gas_force);
      }
    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

// Map particle internal force 2D
template <>
void mpm::ThreePhaseParticle<2>::map_internal_force() {
  if (this->material_id_ != 999) {
    try {
      Eigen::Matrix<double, 2, 1> mixture_force, liquid_force, gas_force;

      this->total_stress_ = this->stress_;

      total_stress_[0] -= this->pore_pressure_ - this->ini_pore_pressure_;
      total_stress_[1] -= this->pore_pressure_ - this->ini_pore_pressure_;

      // LIQUID PHASE
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        liquid_force[0] = dn_dx_(i, 0) * (liquid_pressure_ - ini_liquid_pressure_);
        liquid_force[1] = dn_dx_(i, 1) * (liquid_pressure_ - ini_liquid_pressure_);
        // liquid_force[0] = dn_dx_(i, 0) * (liquid_fraction_ * liquid_pressure_ - ini_liquid_fraction_ * ini_liquid_pressure_);
        // liquid_force[1] = dn_dx_(i, 1) * (liquid_fraction_ * liquid_pressure_ - ini_liquid_fraction_ * ini_liquid_pressure_);

        liquid_force *= this->volume_ * liquid_fraction_;

        nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, liquid_force);
      }

      // GAS PHASE
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        gas_force[0] = dn_dx_(i, 0) * (gas_pressure_ - ini_gas_pressure_);
        gas_force[1] = dn_dx_(i, 1) * (gas_pressure_ - ini_gas_pressure_);
        // gas_force[0] = dn_dx_(i, 0) * (gas_pressure_ *  gas_fraction_ - ini_gas_pressure_ * ini_gas_fraction_);
        // gas_force[1] = dn_dx_(i, 1) * (gas_pressure_ *  gas_fraction_ - ini_gas_pressure_ * ini_gas_fraction_);

        gas_force *= this->volume_ * gas_fraction_;
        nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Gas, gas_force);
      }

      // MIXTURE
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        mixture_force[0] = dn_dx_(i, 0) * total_stress_[0] +
                          dn_dx_(i, 1) * total_stress_[3];
        mixture_force[1] = dn_dx_(i, 1) * total_stress_[1] +
                          dn_dx_(i, 0) * total_stress_[3];

        mixture_force *= -1. * volume_;
        nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Mixture, mixture_force);
      }

    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    } 
  } 
}

// Map particle internal force 3D
template <>
void mpm::ThreePhaseParticle<3>::map_internal_force() {
  if (this->material_id_ != 999) {
    try {
      Eigen::Matrix<double, 3, 1> mixture_force, liquid_force, gas_force;

      this->total_stress_ = this->stress_;
      total_stress_[0] -= this->pore_pressure_;
      total_stress_[1] -= this->pore_pressure_;
      total_stress_[2] -= this->pore_pressure_;

      for (unsigned i = 0; i < nodes_.size(); ++i) {

        // LIQUID PHASE
        liquid_force.setZero();
        liquid_force[0] = dn_dx_(i, 0) * liquid_pressure_;
        liquid_force[1] = dn_dx_(i, 1) * liquid_pressure_;
        liquid_force[2] = dn_dx_(i, 2) * liquid_pressure_;
        liquid_force *= volume_ * liquid_fraction_;
        nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, liquid_force);

        // GAS PHASE
        gas_force.setZero(); 
        gas_force[0] = dn_dx_(i, 0) * gas_pressure_;
        gas_force[1] = dn_dx_(i, 1) * gas_pressure_;
        gas_force[2] = dn_dx_(i, 2) * gas_pressure_;
        gas_force *= volume_ * gas_fraction_; 
        nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Gas, gas_force);

        // MIXTURE
        //           [0 3 5]
        // [0 1 2] * [3 1 4]
        //           [5 4 2]
        mixture_force.setZero();
        mixture_force[0] = dn_dx_(i, 0) * total_stress_[0] +
                          dn_dx_(i, 1) * total_stress_[3] +
                          dn_dx_(i, 2) * total_stress_[5] ;
        mixture_force[1] = dn_dx_(i, 0) * total_stress_[3] +
                          dn_dx_(i, 1) * total_stress_[1] +
                          dn_dx_(i, 2) * total_stress_[4] ;
        mixture_force[2] = dn_dx_(i, 0) * total_stress_[5] +
                          dn_dx_(i, 1) * total_stress_[4] +
                          dn_dx_(i, 2) * total_stress_[2] ;
        mixture_force *= -1. * volume_;
        nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Mixture, mixture_force);
      }
    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  } 
}

// Compute pressure gradient of the particle
template <unsigned Tdim>
inline Eigen::Matrix<double, Tdim, 1> mpm::ThreePhaseParticle<Tdim>::
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
  return pressure_gradient;
}

// Map drag force coefficient - lumped matrix
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::map_drag_force_coefficient() {
  if (this->material_id_ != 999) {  
    try {
      for (unsigned i = 0; i < nodes_.size(); ++i) {

        // LIQUID PHASE
        double liquid_drag_coeff = liquid_fraction_ * liquid_fraction_ *
                                  liquid_viscosity_ / liquid_permeability_;
        liquid_drag_coeff *= volume_ * shapefn_[i];
        nodes_[i]->update_drag_force_coefficient(true, mpm::ParticlePhase::Liquid,
                                                liquid_drag_coeff);

        // GAS PHASE
        double gas_drag_coeff = gas_fraction_ * gas_fraction_ * 
                                  gas_viscosity_ / gas_permeability_;
        gas_drag_coeff *= volume_ * shapefn_[i];
        nodes_[i]->update_drag_force_coefficient(true, mpm::ParticlePhase::Gas,
                                                gas_drag_coeff);
      }

    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

// Map particle heat to nodes
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::map_heat_to_nodes() {
  if (this->material_id_ != 999) {
    try {
      // Calculate mixture heat capacity  
      double mixture_heat_capacity_ = solid_heat_capacity_ +
                          liquid_heat_capacity_ + gas_heat_capacity_;

      // Map mixture heat capacity & heat
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        double mixture_heat_capacity = mixture_heat_capacity_ * shapefn_[i];
        nodes_[i]->update_heat_capacity(true, mpm::ParticlePhase::Mixture,
                                mixture_heat_capacity);
        nodes_[i]->update_heat(true, mpm::ParticlePhase::Mixture, 
                                mixture_heat_capacity * temperature_);
      }
    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

// Map conductive heat transfer
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::map_heat_conduction() {
  if (this->material_id_ != 999) {
    try {
      // Calculate temperature gradient
      mpm::Particle<Tdim>::compute_temperature_gradient(mpm::ParticlePhase::Solid);

      // Calculate thermal conductivity of mixture
      double mixture_cond = solid_fraction_ * solid_thermal_conductivity_ +
                            liquid_fraction_ * liquid_thermal_conductivity_ +
                            gas_fraction_ * gas_thermal_conductivity_;

      // Map heat conduction to nodes
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        double heat_conduction = 0;
        for (unsigned j = 0; j < Tdim; ++j){
          heat_conduction += dn_dx_(i, j) * this->temperature_gradient_[j];
        }
        heat_conduction *= -1 * this->volume_ * mixture_cond;
        nodes_[i]->update_heat_conduction(true, 
                                mpm::ParticlePhase::Mixture, heat_conduction);
      }
    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

// Map convective heat transfer
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::map_heat_convection() {
  if (this->material_id_ != 999) {
    try {
      // Liquid phase & Gas phase
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        double liquid_heat_convection = 0;
        double gas_heat_convection = 0;
        for (unsigned j = 0; j < Tdim; ++j){
          liquid_heat_convection += shapefn_[i] * 
                temperature_gradient_[j] * (liquid_velocity_[j] - velocity_[j]);
          gas_heat_convection += shapefn_[i] *
                temperature_gradient_[j] * (gas_velocity_[j] - velocity_[j]);
        }
        liquid_heat_convection *= -1 * this->volume_ * liquid_heat_capacity_;
        gas_heat_convection *= -1 * this->volume_ * gas_heat_capacity_;
        nodes_[i]->update_heat_convection(true, mpm::ParticlePhase::Liquid,
                                          liquid_heat_convection);
        nodes_[i]->update_heat_convection(true, mpm::ParticlePhase::Gas,
                                          gas_heat_convection);
      } 
    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

//==============================================================================
//  PART 2: UPDATE PARTICLE INFORMATION
//==============================================================================

// Compute updated velocity of the liquid
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::compute_updated_velocity(
                      double dt, double pic, double damping_factor) {
  mpm::Particle<Tdim>::compute_updated_velocity(dt, pic, damping_factor);

  if (this->material_id_ == 999) {
    this->liquid_velocity_ = this->velocity_;
    this->gas_velocity_ = this->velocity_;
  } else {
    try {
      // Get interpolated nodal acceleration
      Eigen::Matrix<double, Tdim, 1> liquid_acceleration;
      Eigen::Matrix<double, Tdim, 1> gas_acceleration;
      liquid_acceleration.setZero();
      gas_acceleration.setZero();
      //  
      for (unsigned i = 0; i < nodes_.size(); ++i) {
      liquid_acceleration +=
                shapefn_[i] * nodes_[i]->acceleration(mpm::ParticlePhase::Liquid);
      gas_acceleration +=
                shapefn_[i] * nodes_[i]->acceleration(mpm::ParticlePhase::Gas);
      }
      // Particle acceleration
      this->liquid_acceleration_ = liquid_acceleration;
      this->gas_acceleration_ = gas_acceleration;  

      // Get PIC velocity
      Eigen::Matrix<double, Tdim, 1> pic_liquid_velocity;
      Eigen::Matrix<double, Tdim, 1> pic_gas_velocity;
      pic_liquid_velocity.setZero();
      pic_gas_velocity.setZero();
      //                                
      for (unsigned i = 0; i < nodes_.size(); ++i){
        pic_liquid_velocity +=
                shapefn_[i] * nodes_[i]->velocity(mpm::ParticlePhase::Liquid);
        pic_gas_velocity +=
                shapefn_[i] * nodes_[i]->velocity(mpm::ParticlePhase::Gas);
      }

      // Applying particle damping
      liquid_acceleration -= damping_factor * this->liquid_velocity_;
      gas_acceleration -= damping_factor * this->gas_velocity_;

      // Get FLIP velocity
      Eigen::Matrix<double, Tdim, 1> flip_liquid_velocity =
                this->liquid_velocity_ + liquid_acceleration * dt;
      Eigen::Matrix<double, Tdim, 1> flip_gas_velocity =
                this->gas_velocity_ + gas_acceleration * dt;

      // Update particle velocity based on PIC value
      this->liquid_velocity_ = pic * pic_liquid_velocity + 
                              (1. - pic) * flip_liquid_velocity;
      this->gas_velocity_ = pic * pic_gas_velocity + 
                              (1. - pic) * flip_gas_velocity;
    } catch (std::exception& exception) {
        console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                        exception.what());
    }
  }
}

// Compute updated pore pressure
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::compute_pore_pressure(double dt){
  if (this->material_id_ != 999) {
    try {
      const double K_ww = porosity_ * dSw_dpw_ + 
                    porosity_ * liquid_saturation_ * liquid_compressibility_;
      const double K_wg = -porosity_ * dSw_dpw_;
      const double K_gg = porosity_ * gas_saturation_ / gas_density_ * 
                    gas_molar_mass_ / gas_constant_ / (PIC_temperature_ + 273.15) +
                    porosity_ * dSw_dpw_;
      const double K_gw = -porosity_ * dSw_dpw_;

      const double K = porosity_ * liquid_saturation_ * liquid_compressibility_ +
                porosity_ * gas_saturation_ / gas_density_ * 
                gas_molar_mass_ / gas_constant_ / (PIC_temperature_ + 273.15);

      this->liquid_strain_rate_= 
              this->compute_strain_rate(dn_dx_, shapefn_,
                                        mpm::ParticlePhase::Liquid);
      this->gas_strain_rate_= 
              this->compute_strain_rate(dn_dx_, shapefn_,
                                        mpm::ParticlePhase::Gas);

      Eigen::Matrix<double, 6, 1> strain_rate;
      strain_rate = 
              this->compute_strain_rate(dn_dx_, shapefn_,
                                        mpm::ParticlePhase::Solid);

      double solid_strain_rate, liquid_strain_rate, gas_strain_rate;
      if (is_axisymmetric_) {
        solid_strain_rate = strain_rate.head(3).sum();
        liquid_strain_rate = liquid_strain_rate_.head(3).sum();
        gas_strain_rate = gas_strain_rate_.head(3).sum();
      } else {
        solid_strain_rate = strain_rate.head(Tdim).sum();
        liquid_strain_rate = liquid_strain_rate_.head(Tdim).sum();
        gas_strain_rate = gas_strain_rate_.head(Tdim).sum();
      }

      const double beta_w = (1 - porosity_) * liquid_saturation_ * solid_expansivity_ * 3.0 +
                            liquid_fraction_ * liquid_expansivity_ * 3.0 ;

      const double beta_g = (1 - porosity_) * gas_saturation_ * solid_expansivity_ * 3.0 +
                            gas_fraction_ / gas_density_ * gas_pressure_ * gas_molar_mass_ /
                            gas_constant_ / std::pow(PIC_temperature_ + 273.15, 2);

      const double beta_m = solid_fraction_ * solid_expansivity_ * 3.0 +
                      liquid_fraction_ * liquid_expansivity_ * 3.0 +
                      gas_fraction_ / gas_density_ * gas_pressure_ * gas_molar_mass_ /
                      gas_constant_ / std::pow(PIC_temperature_ + 273.15, 2);

      const double f_w = -liquid_saturation_ * solid_strain_rate - 
                    liquid_fraction_ * (liquid_strain_rate - solid_strain_rate) +
                    beta_w * temperature_acceleration_; 

      const double f_g = -(1 - liquid_saturation_) * solid_strain_rate -
                    gas_fraction_ * (gas_strain_rate - solid_strain_rate) +
                    beta_g * temperature_acceleration_; ;

      const double f = beta_m * this->temperature_acceleration_ - solid_strain_rate -
                    liquid_fraction_ * (liquid_strain_rate - solid_strain_rate) -
                    gas_fraction_ * (gas_strain_rate - solid_strain_rate);

      Eigen::Matrix<double, 2, 2> K_matrix, K_matrix_inverse;
      K_matrix(0,0) = K_ww;
      K_matrix(0,1) = K_wg;
      K_matrix(1,0) = K_gw;
      K_matrix(1,1) = K_gg;

      // if (K_gg > 1E-16) {
      //   K_matrix_inverse = K_matrix.inverse();

      //   this->liquid_pressure_acceleration_ = K_matrix_inverse(0, 0) * f_w +
      //                                         K_matrix_inverse(0, 1) * f_g;
      //   this->gas_pressure_acceleration_ = K_matrix_inverse(1, 0) * f_w +
      //                                     K_matrix_inverse(1, 1) * f_g;
      // } else {
          this->liquid_pressure_acceleration_ = f_w / K_ww;
          this->gas_pressure_acceleration_ = 0;
      // }

      // double pressure_acceleration = (f_w + f_g) / K;
      // this->liquid_pressure_acceleration_ = pressure_acceleration;
      // this->gas_pressure_acceleration_ = pressure_acceleration;

      this->liquid_pressure_ += this->liquid_pressure_acceleration_ * dt;
      this->gas_pressure_ += this->gas_pressure_acceleration_ * dt;
    if (this->free_surface()) {

      const double para_p0 = liquid_material_->template 
                              property<double>(std::string("para_p0"));
      const double para_m = liquid_material_->template 
                            property<double>(std::string("para_m"));
      const double suction = liquid_material_->template 
                            property<double>(std::string("suction"));
      const double ini_Se = (ini_liquid_saturation_ - this->liquid_saturation_res_) /
                                  (1 - this->liquid_saturation_res_);
      // this->suction_pressure_ = para_p0 * std::pow(std::pow(ini_Se, -1. / para_m) - 1., 1. - para_m);
      this->suction_pressure_ = (1 - ini_liquid_saturation_) * 1E6;
      this->suction_pressure_ += suction;
      this->gas_pressure_ = material_->template property<double>(std::string("p_ref"));
      this->liquid_pressure_ = this->gas_pressure_ - this->suction_pressure_;
    }

      this->PIC_liquid_pressure_ = this->liquid_pressure_;
      this->PIC_gas_pressure_ = this->gas_pressure_;
      this->pore_pressure_ = this->liquid_saturation_ * this->PIC_liquid_pressure_ +
                            this->gas_saturation_ * this->PIC_gas_pressure_; 

  // // Debug
  // if (this->id_ <= 196 && this->id_ >= 191) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   // std::cout << "    K_ww: " << K_ww << std::endl;
  //   // std::cout << "    thermal_term: " << beta_w * this->temperature_acceleration_ * dt << std::endl;
  //   // std::cout << "    strain_term: " << (liquid_saturation_ * solid_strain_rate + 
  //   //                 liquid_fraction_ * (liquid_strain_rate - solid_strain_rate)) * dt << std::endl;
  //   std::cout << "    liquid_pressure_: " << liquid_pressure_ << std::endl;
  //   std::cout << "    this->temperature_acceleration_: " << this->temperature_acceleration_ << std::endl;
  //   std::cout << "    liquid_strain_rate_.head(2).sum(): " << liquid_strain_rate_.head(2).sum() << std::endl;
  //   std::cout << "    strain_rate.head(2).sum(): " << this->compute_strain_rate(dn_dx_, shapefn_,
  //                                       mpm::ParticlePhase::Solid)<< std::endl;
  // }

    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

// Assign a liquid material to particle
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::update_particle_volume() {
  if (this->material_id_ != 999) {
    try {
      // Solid phase  
      mpm::Particle<Tdim>::update_particle_volume();

      // SOLID PHASE
      this->solid_heat_capacity_ = mass_ * solid_specific_heat_;

      // LIQUID PHASE
      this->liquid_fraction_ = porosity_ * liquid_saturation_;
      this->liquid_volume_ = liquid_fraction_ * volume_;
      this->liquid_mass_ = liquid_volume_ * liquid_density_;
      this->liquid_mass_density_ = liquid_fraction_ * liquid_density_;
      this->liquid_chi_ = liquid_saturation_ / (liquid_saturation_ + gas_saturation_);
      this->liquid_heat_capacity_ = liquid_mass_ * liquid_specific_heat_;

      // GAS PHASE
      this->gas_fraction_ = porosity_ * gas_saturation_;
      this->gas_volume_ = gas_fraction_ * volume_;
      this->gas_mass_ = gas_volume_ * gas_density_;
      this->gas_mass_density_ = gas_fraction_ * gas_density_;
      this->gas_chi_ = 1 - liquid_chi_;
      this->gas_heat_capacity_ = gas_mass_ * gas_specific_heat_;

      //Mixture
      this->mixture_mass_ = mass_ + liquid_mass_ + gas_mass_;

    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

// Calculate absoulute/relative permeability
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::update_permeability() {
  if (this->material_id_ != 999) {
    try {   
      double k_phi = 1;
      double k_a = 1;
      double k_r_liquid = 1.;
      double k_r_gas = 1.; 

      const double para_m = liquid_material_->template 
                    property<double>(std::string("para_m"));

      // // Calculate porosity-dependent permeability, k_phi
      // k_phi = std::pow(porosity_ / ini_porosity_, 1.5) * 
      //         std::pow((1 - ini_porosity_) / (1 - porosity_), 3);

      // k_r_liquid = std::pow(effective_saturation_, (2*para_m+3));
      // k_r_gas = std::pow(1-effective_saturation_, 2) * 
      //           (1 - std::pow(effective_saturation_, (2*para_m+1)));

      // Calculate absolute permeability, k_a
      k_a = intrinsic_permeability_ * k_phi;

      // Calculate permeability
      this->liquid_permeability_ = k_a *  k_r_liquid; 
      this->gas_permeability_ = k_a * k_r_gas;
    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

// Compute updated liquid saturation
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::update_liquid_saturation(double dt){
  if (this->material_id_ != 999) {
    try {
      // this->suction_pressure_ = PIC_gas_pressure_ - PIC_liquid_pressure_;

      // const double para_p0 = liquid_material_->template 
      //               property<double>(std::string("para_p0"));
      // const double para_m = liquid_material_->template 
      //               property<double>(std::string("para_m"));
      // double Se = std::pow(suction_pressure_ / para_p0, -para_m);
      // if (Se > 1.0) Se = 1.0;
      // this->effective_saturation_ = Se;

      // // Compute liquid saturation
      // this->liquid_saturation_ = (Se * (1.0 - liquid_saturation_res_) +
      //                           liquid_saturation_res_);
      // // Compute gas saturation
      // this->gas_saturation_ = 1.0 - this->liquid_saturation_;

      // // Compute dSw_dpc
      // const double dSw_dSe = 1 - liquid_saturation_res_;
      // const double dSe_dpc = -para_m / para_p0 * std::pow(suction_pressure_ / para_p0, -para_m - 1.0);
      // this->dSw_dpw_ = -dSw_dSe * dSe_dpc;
      // if (Se > 1.0) this->dSw_dpw_ = 0.0;


      this->suction_pressure_ = PIC_gas_pressure_ - PIC_liquid_pressure_;
      if (this->suction_pressure_ < 0) this->suction_pressure_ = 0;
      this->liquid_saturation_ = 1 - 1E-6 * this->suction_pressure_;
      this->dSw_dpw_ = 1E-6;

      // this->suction_pressure_ = PIC_gas_pressure_ - PIC_liquid_pressure_;
      // double para_p0 = liquid_material_->template 
      //               property<double>(std::string("para_p0"));
      // double para_m = liquid_material_->template 
      //               property<double>(std::string("para_m"));
      // double A = suction_pressure_ / para_p0;
      // double B = 1.0 + std::pow(abs(A), 1.0 / (1.0 - para_m));
      // double Se = std::pow(B, -para_m);

      // this->liquid_saturation_ = (Se * (1.0 - liquid_saturation_res_) +
      //                   liquid_saturation_res_);

      // double dA_dpw = -1.0 / para_p0;
      // double dB_dA = 1.0 / (1.0 - para_m) * std::pow(abs(A), (1.0 / (1.0 - para_m) - 1.0));
      // double dSe_dB = -para_m * std::pow(B, (-para_m - 1.0));
      // double dSw_dSe = (1 - liquid_saturation_res_ - gas_saturation_res_) * Se; 
      // double dSw_dpw = dSw_dSe * dSe_dB * dB_dA * dA_dpw;

      // this->dSw_dpw_ = dSw_dpw;

      this->gas_saturation_ = 1.0 - this->liquid_saturation_;

      this->liquid_chi_ = this->liquid_saturation_;
      this->gas_chi_ = 1 - this->liquid_chi_;
    } catch (std::exception& exception) {
        console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                        exception.what());
    }
  }
}

// Compute updated density
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::update_particle_density(double dt){
  if (this->material_id_ != 999) {  
    try {
      const double p_ref = material_->template property<double>(std::string("p_ref"));
      this->gas_density_ = gas_molar_mass_ * PIC_gas_pressure_/ 
                            gas_constant_ / (PIC_temperature_ + 273.15);

      // Check if NaN
      if (!(this->gas_density_ > 0.)) gas_density_ = 1.0;

      // this->density_ /= 1 + solid_expansivity_ * temperature_increment_;
      this->liquid_density_ /= 1 + liquid_expansivity_ * temperature_increment_ -
                                  dt * liquid_pressure_acceleration_ * liquid_compressibility_;

    } catch (std::exception& exception) {
        console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                        exception.what());
    }
  }
}

//==============================================================================
//  PART 3: ASSIGN AND APPLY BOUNDARY CONDITIONS
//==============================================================================

// Assign traction
template <unsigned Tdim>
bool mpm::ThreePhaseParticle<Tdim>::assign_particle_traction(unsigned direction,
                                                          double traction) {
  bool status = false;
  if (this->material_id_ != 999) {
    try {
      if (direction >= Tdim * 3  ||
          this->volume_ == std::numeric_limits<double>::max()) {
        throw std::runtime_error(
            "Particle mixture traction property: volume / direction is invalid");
      }
      // Assign mixture traction
      if (direction < Tdim) 
        mixture_traction_(direction) =
            traction * this->volume_ / this->size_(direction);
      else if (direction < Tdim * 2) 
        liquid_traction_(direction - Tdim) =
            -traction * this->volume_ / this->size_(direction - Tdim);
      else 
        gas_traction_(direction - Tdim) =
            -traction * this->volume_ / this->size_(direction - Tdim);

      status = true;
      this->set_mixture_traction_ = true;
    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

// Assign particle liquid phase velocity constraint
template <unsigned Tdim>
bool mpm::ThreePhaseParticle<Tdim>::assign_particle_liquid_velocity_constraint(
      unsigned dir, double velocity) {
  bool status = true;
  try {
    // Constrain directions can take values between 0 and Dim
    if (dir < Tdim * 3)
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
void mpm::ThreePhaseParticle<Tdim>::apply_particle_liquid_velocity_constraints() {
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
bool mpm::ThreePhaseParticle<Tdim>::assign_particle_pore_pressure_constraint(
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
void mpm::ThreePhaseParticle<Tdim>::apply_particle_pore_pressure_constraints(
      double pore_pressure) {
  // Set particle temperature constraint
  // this->pore_pressure_ = pore_pressure;
  // this->liquid_pressure_ = pore_pressure;
  // this->gas_pressure_ = pore_pressure;
  // this->PIC_liquid_pressure_ = pore_pressure;
  // this->PIC_gas_pressure_ = pore_pressure;
  this->set_pressure_constraint_ = true;
}

// Overwrite node velocity to get strain correct
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::map_moving_rigid_velocity_to_nodes(
    unsigned dir, double velocity, double dt) noexcept {
  if (this->material_id_ == 999){  
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      nodes_[i]->assign_velocity_from_rigid(dir, velocity, dt);
    }
  }
}

//==============================================================================
//  PART 4: SMOOTHING TECHNIQUE
//==============================================================================

// Map particle pore liquid pressure to nodes
template <unsigned Tdim>
bool mpm::ThreePhaseParticle<Tdim>::map_pore_pressure_to_nodes(
    double current_time) noexcept {
  if (this->material_id_ != 999) {
    // Check if particle mass is set
    assert(liquid_mass_ != std::numeric_limits<double>::max());

    bool status = true;
    // Map particle liquid mass and pore pressure to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i){
      nodes_[i]->update_mass_pressure(mpm::ParticlePhase::Liquid,
                                      shapefn_[i] * liquid_mass_ * pore_pressure_,
                                      current_time);
    }
    return status;
  }
}

// Map particle pore liquid pressure to nodes
template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::map_mass_pressure_to_nodes() {
  if (this->material_id_ != 999) {  
    // Check if particle mass is set
    assert(liquid_mass_ != std::numeric_limits<double>::max());
    if (this->material_id_ != 999){
      // Map particle liquid mass and pore pressure to nodes
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        if (gas_mass_ > 1E-15) {
          nodes_[i]->update_pressure(true, mpm::ParticlePhase::Gas,
                                          shapefn_[i] * gas_mass_ * gas_pressure_);
        }
        nodes_[i]->update_pressure(true, mpm::ParticlePhase::Liquid,
                                  shapefn_[i] * liquid_mass_ * liquid_pressure_);
      }
    }
  }
}

// Compute pore liquid pressure smoothing based on nodal pressure
template <unsigned Tdim>
bool mpm::ThreePhaseParticle<Tdim>::compute_pore_pressure_smoothing() noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  bool status = true;  
  if (this->material_id_ != 999) {
    double liquid_pressure = 0;
    double gas_pressure = 0;    
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      liquid_pressure += shapefn_(i) * nodes_[i]->pressure(mpm::ParticlePhase::Liquid);
      gas_pressure += shapefn_(i) * nodes_[i]->pressure(mpm::ParticlePhase::Gas);
    }

    // Update pore liquid pressure to interpolated nodal pressure
    if (gas_mass_ > 1E-15) {
      this->gas_pressure_ = gas_pressure;
    }
    // this->liquid_pressure_= liquid_pressure;


    if (this->free_surface()) {
      this->gas_pressure_ = material_->template property<double>(std::string("p_ref"));
      const double para_p0 = liquid_material_->template 
                              property<double>(std::string("para_p0"));
      const double para_m = liquid_material_->template 
                            property<double>(std::string("para_m"));
      const double suction = liquid_material_->template 
                            property<double>(std::string("suction"));
      const double ini_Se = (ini_liquid_saturation_ - this->liquid_saturation_res_) /
                                  (1 - this->liquid_saturation_res_);
      // this->suction_pressure_ = para_p0 * std::pow(std::pow(ini_Se, -1. / para_m) - 1., 1. - para_m);
      this->suction_pressure_ = (1 - ini_liquid_saturation_) * 1E6;
      this->suction_pressure_ += suction;
      this->liquid_pressure_ = this->gas_pressure_ - this->suction_pressure_;
    }

    this->PIC_liquid_pressure_ = this->liquid_pressure_;
    this->PIC_gas_pressure_ = this->gas_pressure_;
    this->pore_pressure_ = this->liquid_chi_ * this->PIC_liquid_pressure_ +
                          this->gas_chi_ * this->PIC_gas_pressure_;

  }
  return status;
}

template <unsigned Tdim>
void mpm::ThreePhaseParticle<Tdim>::output_results(double step_) {
    // Output directory (change if needed)
    const std::string dir_path =
        "/home/chen/mpm/examples/1D_consolidation/MATLABinput_benchmark";

    // Make sure the directory exists
    std::error_code ec;
    if (!std::filesystem::exists(dir_path)) {
        std::filesystem::create_directories(dir_path, ec);
        if (ec) {
            std::cerr << "Failed to create directory: " << dir_path
                      << " (" << ec.message() << ")\n";
            return;
        }
    }

    // Loop over all registered scalar properties
    for (const auto& kv : this->scalar_property_) {
        const std::string& var_name = kv.first;  // e.g. "pore_pressures"
        const double value = kv.second();        // current value from lambda

        // File name pattern: <var_name>_<step>.txt
        const std::string file_path =
            dir_path + "/" + var_name + "_" +
            std::to_string(static_cast<int>(step_)) + ".txt";

        // Append one row per particle
        std::ofstream ofs(file_path, std::ios::out | std::ios::app);
        if (!ofs.is_open()) {
            std::cerr << "Error opening file: " << file_path << "\n";
            continue;
        }

        // Row format:
        // 2D : id x y value
        // 3D : id x y z value
        // 1D : id x value (just in case)
        if (Tdim == 2) {
            ofs << this->id_ << ' '
                << this->coordinates_[0] << ' '
                << this->coordinates_[1] << ' '
                << value << '\n';
        } else if (Tdim == 3) {
            ofs << this->id_ << ' '
                << this->coordinates_[0] << ' '
                << this->coordinates_[1] << ' '
                << this->coordinates_[2] << ' '
                << value << '\n';
        } else {
            ofs << this->id_ << ' '
                << this->coordinates_[0] << ' '
                << value << '\n';
        }
    }
}
