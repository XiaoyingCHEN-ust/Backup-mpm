// Construct a two phase particle with id and coordinates
template <unsigned Tdim>
mpm::MHBSParticle<Tdim>::MHBSParticle(Index id, const VectorDim& coord)
    : mpm::Particle<Tdim>(id, coord) {

    this->initialise_liquid_phase();

    // Set material pointer to null
    liquid_material_ = nullptr;
    // Logger
    std::string logger = "MHBSParticle" + std::to_string(Tdim) + "d::" + std::to_string(id);
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//==============================================================================
// ASSIGN INITIAL CONDITIONS
//==============================================================================

// Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::MHBSParticle<Tdim>::initialise_particle(
    const HDF5Particle& particle) {
    // Derive from particle
    mpm::Particle<Tdim>::initialise_particle(particle);

    // MIXTURE data
    this->liquid_material_id_ = particle.liquid_material_id;
    this->pore_pressure_ = particle.pore_pressure;
    // this->suction_pressure_ = particle.suction_pressure;
    // this->suction_pressure_ = particle.suction_pressure;  

    // HYDRATE PHASE data
    this->hydrate_saturation_ = particle.hydrate_saturation;
    this->hydrate_fraction_ = particle.hydrate_fraction;
    this->hydrate_density_ = particle.hydrate_density;
    // this->hydrate_volume_ = particle.hydrate_volume;
    this->hydrate_mass_ = particle.hydrate_mass;

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
bool mpm::MHBSParticle<Tdim>::initialise_particle(
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
void mpm::MHBSParticle<Tdim>::initialise_liquid_phase() {
    // MIXTURE
    set_mixture_traction_ = false;
    pore_pressure_ = 0.;
    PIC_pore_pressure_ = 0.;
    ini_pore_pressure_ = 0.;
    suction_pressure_ = 0.;
    ini_suction_pressure_ = 0.;
    set_pressure_constraint_ = false;
    PIC_porosity_ = 0.;
    PIC_volume_ = 0.;
    PIC_volumetric_strain_ = 0.;
    PIC_deviatoric_strain_ = 0.;
    PIC_mean_stress_ = 0.;
    PIC_deviatoric_stress_ = 0.;

    // Hydrate properties
    hydrate_saturation_ = 0.;
    PIC_hydrate_saturation_ = 0.;
    hydrate_fraction_ = 0;
    hydrate_density_ = 0.;
    hydrate_volume_ = 0.;
    hydrate_mass_ = 0.;
    hydrate_mass_density_ = 0.;
    hydrate_source_ = 0.;

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
    ini_liquid_pressure_ = 0.;
    PIC_liquid_pressure_ = 0.;
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
    ini_gas_pressure_ = 0.;
    PIC_gas_pressure_ = 0.;
    gas_pressure_acceleration_ = 0.;
    gas_volumetric_strain_ = 0.;
    gas_permeability_ = 1.;
    gas_source_ = 0.;

    // Link data with NAME
    this->scalar_property_.insert({
      {"pore_pressures",         [&]() {return this->pore_pressure_;}},
      {"PIC_pore_pressures",     [&]() {return this->PIC_pore_pressure_;}},
      {"suction_pressures",      [&]() {return this->suction_pressure_;}},
      {"hydrate_saturations",    [&]() {return this->hydrate_saturation_;}},
      {"PIC_hydrate_saturations",[&]() {return this->PIC_hydrate_saturation_;}},
      {"PIC_porosities",         [&]() {return this->PIC_porosity_;}},
      {"PIC_volumes",             [&]() {return this->PIC_volume_;}},
      {"hydrate_fractions",      [&]() {return this->hydrate_fraction_;}},
      {"hydrate_densities",      [&]() {return this->hydrate_density_;}},
      {"hydrate_sources",        [&]() {return this->hydrate_source_;}},
      {"hydrate_volumes",        [&]() {return this->hydrate_volume_;}},
      {"hydrate_masses",         [&]() {return this->hydrate_mass_;}},
      {"hydrate_mass_densities", [&]() {return this->hydrate_mass_density_;}},
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
      {"gas_critical_times",     [&]() {return this->gas_critical_time_;}},
      {"PIC_volumetric_strains", [&]() {return this->PIC_volumetric_strain_;}},
      {"PIC_deviatoric_strains", [&]() {return this->PIC_deviatoric_strain_;}},
      {"mean_stresses",          [&]() {return this->PIC_mean_stress_;}},
      {"deviatoric_stresses",    [&]() {return this->PIC_deviatoric_stress_;}},
      {"deviatoric_strains",     [&]() {return this->deviatoric_strain_;}  },
    });

    this->vector_property_.insert({
      {"liquid_velocities",      [&]() {return this->liquid_velocity_;}},
      {"liquid_accelerations",   [&]() {return this->liquid_acceleration_;}},
      {"liquid_strains",         [&]() {return this->liquid_strain_;}},
      {"liquid_fluxes",          [&]() {return this->liquid_flux_;}},
      {"liquid_pressure_gradients",[&]() {return this->liquid_pressure_gradient_;}},
      {"gas_velocities",         [&]() {return this->gas_velocity_;}},
      {"gas_accelerations",      [&]() {return this->gas_acceleration_;}},
      {"gas_strains",            [&]() {return this->gas_strain_;}},
      {"gas_fluxes",             [&]() {return this->gas_flux_;}},
      {"gas_pressure_gradients",[&]() {return this->gas_pressure_gradient_;}},
      {"K_matrix",               [&]() {return this->K_matrix_;}}
    }); 
}

// Return particle data in HDF5 format
template <unsigned Tdim>
mpm::HDF5Particle mpm::MHBSParticle<Tdim>::hdf5() {
    // Derive from particle
    auto particle_data = mpm::Particle<Tdim>::hdf5();

    // MIXTURE data
    particle_data.liquid_material_id = this->liquid_material_id_;
    particle_data.pore_pressure = this->pore_pressure_;
    // particle_data.suction_pressure = this->suction_pressure_;

    // HYDRATE PHASE data
    particle_data.hydrate_saturation = this->hydrate_saturation_;
    particle_data.hydrate_fraction = this->hydrate_fraction_;
    particle_data.hydrate_density = this->hydrate_density_;
    // particle_data.hydrate_volume = this->hydrate_volume_;
    particle_data.hydrate_mass = this->hydrate_mass_;

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
bool mpm::MHBSParticle<Tdim>::assign_liquid_material(
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
void mpm::MHBSParticle<Tdim>::compute_mass() {
  mpm::Particle<Tdim>::compute_mass();
  this->assign_initial_properties();
}

// Assign initial properties to particle
template <unsigned Tdim>
bool mpm::MHBSParticle<Tdim>::assign_initial_properties() {
  bool status = true;
  try {
    this->ini_porosity_ = material_->template property<double>(std::string("porosity"));
    this->porosity_ = ini_porosity_; 
    this->PIC_porosity_ = ini_porosity_;
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

    // HYDRATE PHASE
    // Constant propeties
    this->hydrate_thermal_conductivity_ = material_->template 
                  property<double>(std::string("hydrate_thermal_conductivity"));
    this->hydrate_specific_heat_ = material_->template 
                  property<double>(std::string("hydrate_specific_heat"));
    this->hydrate_expansivity_ = material_->template 
                  property<double>(std::string("hydrate_expansivity"));
    this->hydrate_molar_mass_ = material_->template 
                  property<double>(std::string("hydrate_molar_mass"));
    this->hydrate_latent_ = material_->template 
                  property<double>(std::string("hydration_latent"));
    this->hydration_number_ = material_->template 
                  property<double>(std::string("hydration_number"));
    // Time-dependent propeties    
    ini_hydrate_density_ = material_->template 
                  property<double>(std::string("hydrate_density"));
    ini_hydrate_saturation_ = material_->template 
                  property<double>(std::string("hydrate_saturation"));
    this->hydrate_density_ = ini_hydrate_density_;   
    this->hydrate_saturation_ = ini_hydrate_saturation_;
    this->PIC_hydrate_saturation_ = ini_hydrate_saturation_;

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
    this->gas_density_ = gas_molar_mass_ * p_ref / 
                          gas_constant_ / (PIC_temperature_ + 273.15);
    this->gas_saturation_ = ini_gas_saturation_;
    this->gas_viscosity_ = ini_gas_viscosity_;

    if (!(gas_density_>0)) gas_density_ = 1.29;

    // SOLID PHASE
    this->solid_heat_capacity_ = mass_ * solid_specific_heat_;

    // HYDRATE PHASE
    this->hydrate_fraction_ = porosity_ * hydrate_saturation_;
    this->hydrate_volume_ = hydrate_fraction_ * volume_;
    this->hydrate_mass_ = hydrate_volume_ * hydrate_density_;    
    this->hydrate_mass_density_ = hydrate_fraction_ * hydrate_density_;
    this->hydrate_heat_capacity_ = hydrate_mass_ * hydrate_specific_heat_;

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

    //MIXTURE
    this->mixture_mass_ = mass_ + hydrate_mass_ + liquid_mass_ + gas_mass_;
    this->PIC_volume_ = this->volume_;

    // // Initial pore presssure
    // const double para_p0 = liquid_material_->template 
    //               property<double>(std::string("para_p0"));
    // const double para_m = liquid_material_->template 
    //               property<double>(std::string("para_m"));
    // this->effective_saturation_ = 
    //                     (liquid_saturation_ / (1. - hydrate_saturation_) -
    //                     liquid_saturation_res_) /
    //                     (1. - liquid_saturation_res_ - gas_saturation_res_);
    // this->suction_pressure_ = para_p0 * 
    //               std::pow(std::pow(effective_saturation_, -1. / para_m) - 1. + 1E-15, 1. - para_m);
    // this->gas_pressure_ = p_ref;
    // this->liquid_pressure_ = this->gas_pressure_ - this->suction_pressure_;
    // this->PIC_gas_pressure_ = this->gas_pressure_;
    // this->PIC_liquid_pressure_ = this->liquid_pressure_;
    // this->pore_pressure_ = liquid_chi_ * PIC_liquid_pressure_ +
    //                       gas_chi_ * PIC_gas_pressure_;

    // this->ini_gas_pressure_ = this->gas_pressure_;
    // this->ini_liquid_pressure_ = this->liquid_pressure_;
    // this->ini_pore_pressure_ = this->pore_pressure_;
    // this->ini_suction_pressure_ = this->suction_pressure_;

    // if (id_ == 1000) {
    //   std::cout << this->effective_saturation_ << "\n"
    //             << this->suction_pressure_ << "\n"
    //             << this->gas_pressure_ << "\n"
    //             << this->liquid_pressure_ << "\n"
    //             << this->pore_pressure_ << "\n";
    // }

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
void mpm::MHBSParticle<Tdim>::map_mass_momentum_to_nodes() noexcept {
  if (this->material_id_ != 999) {
    // SOLID PHASE
    mpm::Particle<Tdim>::map_mass_momentum_to_nodes();

    for (unsigned i = 0; i < nodes_.size(); ++i) {
      // Reduce the cost of accessing the shapefn_ array repeatedly 
      double shapefn_i = shapefn_[i];

      // HYDRATE PHASE
      double hydrate_mass = hydrate_mass_ * shapefn_i;
      nodes_[i]->update_mass_momentum(true, mpm::ParticlePhase::Solid, 
                                            hydrate_mass, 
                                            hydrate_mass * velocity_);
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
void mpm::MHBSParticle<Tdim>::map_external_force(const VectorDim& pgravity) {
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
void mpm::MHBSParticle<2>::map_internal_force() {
  if (this->material_id_ != 999) {
    try {
      Eigen::Matrix<double, 2, 1> mixture_force, liquid_force, gas_force;

      this->total_stress_ = this->stress_;

      total_stress_[0] -= this->pore_pressure_ - this->ini_pore_pressure_;
      total_stress_[1] -= this->pore_pressure_ - this->ini_pore_pressure_;

      if (is_axisymmetric_) total_stress_[2] -= this->pore_pressure_;

      // LIQUID PHASE
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        liquid_force[0] = dn_dx_(i, 0) * (liquid_pressure_ - ini_liquid_pressure_);
        liquid_force[1] = dn_dx_(i, 1) * (liquid_pressure_ - ini_liquid_pressure_);

        // In case of 2D axisymmetric
        if (is_axisymmetric_) {
          liquid_force[0] += shapefn_[i] / this->coordinates_(0) * liquid_pressure_;
        }

        liquid_force *= this->volume_ * liquid_fraction_;
        nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, liquid_force);
      }

      // GAS PHASE
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        gas_force[0] = dn_dx_(i, 0) * (gas_pressure_ - ini_gas_pressure_);
        gas_force[1] = dn_dx_(i, 1) * (gas_pressure_ - ini_gas_pressure_);

        // In case of 2D axisymmetric
        if (is_axisymmetric_) {
          gas_force[0] += shapefn_[i] / this->coordinates_(0) * gas_pressure_;
        }

        gas_force *= this->volume_ * gas_fraction_;
        nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Gas, gas_force);
      }

      // MIXTURE
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        mixture_force[0] = dn_dx_(i, 0) * total_stress_[0] + 
                          dn_dx_(i, 1) * total_stress_[3];
        mixture_force[1] = dn_dx_(i, 1) * total_stress_[1] + 
                          dn_dx_(i, 0) * total_stress_[3];
        if (is_axisymmetric_) 
          mixture_force[0] += shapefn_[i] / this->coordinates_(0) * total_stress_[2];

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
void mpm::MHBSParticle<3>::map_internal_force() {
  if (this->material_id_ != 999) {
    try {
      Eigen::Matrix<double, 3, 1> mixture_force, liquid_force, gas_force;

      this->total_stress_ = this->stress_;
      total_stress_[0] -= this->pore_pressure_ - this->ini_pore_pressure_;
      total_stress_[1] -= this->pore_pressure_ - this->ini_pore_pressure_;
      total_stress_[2] -= this->pore_pressure_ - this->ini_pore_pressure_;

      for (unsigned i = 0; i < nodes_.size(); ++i) {

        // LIQUID PHASE
        liquid_force.setZero();
        liquid_force[0] = dn_dx_(i, 0) * (liquid_pressure_ - ini_liquid_pressure_);
        liquid_force[1] = dn_dx_(i, 1) * (liquid_pressure_ - ini_liquid_pressure_);
        liquid_force[2] = dn_dx_(i, 2) * (liquid_pressure_ - ini_liquid_pressure_);
        liquid_force *= volume_ * liquid_fraction_;
        nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, liquid_force);

        // GAS PHASE
        gas_force.setZero(); 
        gas_force[0] = dn_dx_(i, 0) * (gas_pressure_ - ini_gas_pressure_);
        gas_force[1] = dn_dx_(i, 1) * (gas_pressure_ - ini_gas_pressure_);
        gas_force[2] = dn_dx_(i, 1) * (gas_pressure_ - ini_gas_pressure_);
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

// Map drag force coefficient - lumped matrix
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::map_drag_force_coefficient() {
  if (this->material_id_ != 999) {  
    try {
      for (unsigned i = 0; i < nodes_.size(); ++i) {

        // LIQUID PHASE
        double liquid_drag_coeff = liquid_fraction_ * liquid_fraction_ * 
                                  0.981 * liquid_viscosity_ / liquid_permeability_;
        liquid_drag_coeff *= volume_ * shapefn_[i];
        nodes_[i]->update_drag_force_coefficient(true, mpm::ParticlePhase::Liquid, 
                                                liquid_drag_coeff);

        // GAS PHASE
        double gas_drag_coeff = gas_fraction_ * gas_fraction_ * 
                                  0.981 * gas_viscosity_ / gas_permeability_;
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
void mpm::MHBSParticle<Tdim>::map_heat_to_nodes() {
  if (this->material_id_ != 999) {
    try {
      // Calculate mixture heat capacity  
      double mixture_heat_capacity_ = 
                          solid_heat_capacity_ + hydrate_heat_capacity_ + 
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
void mpm::MHBSParticle<Tdim>::map_heat_conduction() {
  if (this->material_id_ != 999) {
    try {
      // Calculate temperature gradient
      mpm::Particle<Tdim>::compute_temperature_gradient(mpm::ParticlePhase::Solid);

      // Calculate thermal conductivity of mixture
      double mixture_cond = solid_fraction_ * solid_thermal_conductivity_ +
                            hydrate_fraction_ * hydrate_thermal_conductivity_ +
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
void mpm::MHBSParticle<Tdim>::map_heat_convection() {
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
        liquid_heat_convection *= -1 * liquid_heat_capacity_;
        gas_heat_convection *= -1 * gas_heat_capacity_;
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

// Map convective heat transfer
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::map_mass_convection() {
  try {
    // Liquid phase & Gas phase
    this->gas_density_gradient_ = 
                    this->compute_density_gradient(mpm::ParticlePhase::Gas);  
    this->liquid_density_gradient_ =
                    this->compute_density_gradient(mpm::ParticlePhase::Liquid);

    double liquid_coeff = liquid_density_ / liquid_compressibility_;
    double gas_coeff = gas_density_ * gas_density_ * gas_constant_ * 
                      (PIC_temperature_ + 273.15) / gas_molar_mass_;

    for (unsigned i = 0; i < nodes_.size(); ++i) {
      double liquid_mass_convection = 0;
      double gas_mass_convection = 0;
      for (unsigned j = 0; j < Tdim; ++j){
        liquid_mass_convection += shapefn_[i] *
               liquid_density_gradient_[j] * (liquid_velocity_[j] - velocity_[j]);
        gas_mass_convection += shapefn_[i] * 
              gas_density_gradient_[j] * (gas_velocity_[j] - velocity_[j]);
      }
      liquid_mass_convection *= -liquid_fraction_ / liquid_density_ * 
                                  this->volume_ * liquid_coeff;
      gas_mass_convection *= -gas_fraction_ / gas_density_ * 
                                  this->volume_ * gas_coeff;
      nodes_[i]->update_mass_convection(true, mpm::ParticlePhase::Liquid, 
                                        liquid_mass_convection);
      nodes_[i]->update_mass_convection(true, mpm::ParticlePhase::Gas, 
                                        gas_mass_convection);
    } 
  } catch (std::exception& exception) {
    console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__, 
                    exception.what());
  }  
}

// Compute density gradient of the particle
template <unsigned Tdim>
inline Eigen::Matrix<double, Tdim, 1> mpm::MHBSParticle<Tdim>::
                      compute_density_gradient(unsigned phase) noexcept {

  Eigen::Matrix<double, Tdim, 1> density_gradient;
  density_gradient.setZero();
  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    double density = nodes_[i]->density(phase);
    for (unsigned j = 0; j < Tdim; ++j) {
      // density_gradient = partial p / partial X = p_{i,j}
      density_gradient[j] += dn_dx_(i, j) * density;
      if (std::fabs(density_gradient[j]) < 1.E-15)
        density_gradient[j] = 0.;
    }
  }
  return density_gradient;
}

// Map particle pore liquid pressure to nodes
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::map_density_to_nodes() {
  // Check if particle mass is set
  assert(liquid_mass_ != std::numeric_limits<double>::max());
  if (this->material_id_ != 999){
    // Map particle liquid mass and pore pressure to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      nodes_[i]->update_density(true, mpm::ParticlePhase::Gas,
                                      shapefn_[i] * gas_density_);
      nodes_[i]->update_density(true, mpm::ParticlePhase::Liquid,
                                      shapefn_[i] * liquid_density_);
    }
  }
}

// Map source term due to phase transition - lumped matrix
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::map_heat_source() {
  if (this->material_id_ != 999) {  
    try {
      double liquid_coeff = liquid_density_ / liquid_compressibility_;
      double gas_coeff = gas_density_ * gas_density_ * gas_constant_ * 
                                (PIC_temperature_ + 273.15) / gas_molar_mass_;

      for (unsigned i = 0; i < nodes_.size(); ++i) {

        // LIQUID PHASE
        double liquid_source = liquid_source_ / liquid_density_ * 
                                          volume_ * shapefn_[i] * liquid_coeff;
        nodes_[i]->update_mass_source(true, mpm::ParticlePhase::Liquid, 
                                                liquid_source);

        // GAS PHASE
        double gas_source = (gas_source_ / gas_density_ + 
                              hydrate_source_ / hydrate_density_) * 
                              volume_ * shapefn_[i] * gas_coeff;
        nodes_[i]->update_mass_source(true, mpm::ParticlePhase::Gas, 
                                                gas_source);

        // ENERGY
        double heat_source = heat_source_ * volume_ * shapefn_[i] ;
        nodes_[i]->update_heat_source(true, mpm::ParticlePhase::Mixture, 
                                                heat_source);
      }
    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

// Map hydraulic conduction
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::map_hydraulic_conduction() {
  if (this->material_id_ != 999) {
    double liquid_coeff = liquid_density_ / liquid_compressibility_;
    double gas_coeff = gas_density_ * gas_density_ * gas_constant_ * 
                      (PIC_temperature_ + 273.15) / gas_molar_mass_;

    double K_wg = porosity_ * dSw_dpw_;
    double K_gw = porosity_ * dSw_dpw_;

    K_wg *= liquid_coeff;
    K_gw *= gas_coeff;

    double liquid_conductivity = liquid_permeability_ / liquid_viscosity_;
    double gas_conductivity = gas_permeability_ / gas_viscosity_;

    liquid_conductivity *= liquid_coeff;
    gas_conductivity *= gas_coeff;

    // Assign liquid gradient
    this->gas_pressure_gradient_ = 
                      this->compute_pressure_gradient(mpm::ParticlePhase::Gas);  
    this->liquid_pressure_gradient_ = 
                      this->compute_pressure_gradient(mpm::ParticlePhase::Liquid);

    this->liquid_flux_ = -liquid_permeability_ / liquid_viscosity_ * 
            (liquid_pressure_gradient_ + liquid_density_ * liquid_acceleration_);
    this->gas_flux_ = -gas_permeability_ / gas_viscosity_ * (gas_pressure_gradient_ 
                          + gas_density_ * gas_acceleration_);
    // this->liquid_velocity_ = liquid_flux_ / liquid_fraction_;
    // this->gas_velocity_ = gas_flux_ / gas_fraction_;

    // Compute nodal gas conduction
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      double gas_conduction = 0;
      for (unsigned j = 0; j < Tdim; ++j){
        gas_conduction += 
            dn_dx_(i, j) * this->gas_pressure_gradient_[j] * gas_conductivity -
            dn_dx_(i, j) * gas_density_ * 
                    (gas_acceleration_[j] - pgravity_[j]) * gas_conductivity -
            dn_dx_(i, j) * velocity_[j] * gas_coeff * (1 - liquid_saturation_);
      }
      gas_conduction *= -1 * this->volume_;
      nodes_[i]->update_hydraulic_conduction(true, mpm::ParticlePhase::Gas, 
                                              gas_conduction);
      nodes_[i]->update_K_coeff(true, mpm::ParticlePhase::Gas, 
                                              this->volume_ * K_gw);
    }

    // Compute nodal liquid conduction
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      double liquid_conduction = 0;  
      for (unsigned j = 0; j < Tdim; ++j){
        liquid_conduction += 
            dn_dx_(i, j) * this->liquid_pressure_gradient_[j] * liquid_conductivity -
            dn_dx_(i, j) * liquid_density_ * 
                  (liquid_acceleration_[j] - pgravity_[j]) * liquid_conductivity -
            dn_dx_(i, j) * velocity_[j] * liquid_coeff * liquid_saturation_;
        // liquid_conduction += -dn_dx_centroid_(i, j) * (liquid_velocity_[j] - velocity_[j]) * liquid_coeff * liquid_fraction_ -
        //                      dn_dx_centroid_(i, j) * velocity_[j] * liquid_coeff * liquid_saturation_;

      }
      liquid_conduction *= -1 * this->volume_; 
      nodes_[i]->update_hydraulic_conduction(true, mpm::ParticlePhase::Liquid, 
                                              liquid_conduction);
      nodes_[i]->update_K_coeff(true, mpm::ParticlePhase::Liquid, 
                                              this->volume_ * K_wg);
    }
  }
}

// Map volumetric strain
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::map_volumetric_strain() {
  double liquid_coeff = liquid_density_ / liquid_compressibility_;
  double gas_coeff = gas_density_ * gas_density_ * gas_constant_ * 
                              (PIC_temperature_ + 273.15) / gas_molar_mass_;
  // Compute nodal liquid conduction
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    double volumetric_strain = 0;  

    volumetric_strain = this->volume_ * shapefn_[i] * dvolumetric_strain_;

    nodes_[i]->update_volumetric_strain(true, mpm::ParticlePhase::Liquid, 
                      volumetric_strain * liquid_saturation_ * liquid_coeff);
    nodes_[i]->update_volumetric_strain(true, mpm::ParticlePhase::Gas, 
                      volumetric_strain * (1 - liquid_saturation_) * gas_coeff);
  }

}

// Compute pressure gradient of the particle
template <unsigned Tdim>
inline Eigen::Matrix<double, Tdim, 1> mpm::MHBSParticle<Tdim>::
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

//==============================================================================
//  PART 2: UPDATE PARTICLE INFORMATION
//==============================================================================

// Compute updated pore_pressure of the particle
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::update_particle_pore_pressure(
    double dt, double pic_t) noexcept {
  if (this->material_id_ != 999) {
    // Check if particle has a valid cell ptr and pic_t value
    assert(cell_ != nullptr);
    assert((-1.E-15) <= pic_t && pic_t <= (1 + 1.E-15));

    // Get PIC pressure
    double PIC_liquid_pressure = 0;
    double PIC_gas_pressure = 0;

    for (unsigned i = 0; i < nodes_.size(); ++i){
      PIC_liquid_pressure +=
          shapefn_[i] * nodes_[i]->pressure(mpm::ParticlePhase::Liquid);
      PIC_gas_pressure +=
          shapefn_[i] * nodes_[i]->pressure(mpm::ParticlePhase::Gas);
    }

    // pressure increment
    this->liquid_pressure_increment_ = 
                              PIC_liquid_pressure - this->PIC_liquid_pressure_; 
    this->gas_pressure_increment_ = PIC_gas_pressure - this->PIC_gas_pressure_;     

    this->PIC_liquid_pressure_ = PIC_liquid_pressure;
    this->PIC_gas_pressure_ = PIC_gas_pressure;

    // Get interpolated nodal pressure acceleration
    double liquid_pressure_acceleration = 0.;
    double gas_pressure_acceleration = 0.;

    for (unsigned i = 0; i < nodes_.size(); ++i){
      liquid_pressure_acceleration += shapefn_[i] * 
                  nodes_[i]->pressure_acceleration(mpm::ParticlePhase::Liquid);
      gas_pressure_acceleration += shapefn_[i] * 
                  nodes_[i]->pressure_acceleration(mpm::ParticlePhase::Gas);
    }

    this->liquid_pressure_acceleration_ = liquid_pressure_acceleration;
    this->gas_pressure_acceleration_ = gas_pressure_acceleration;

    // Get FLIP pressure
    double FLIP_liquid_pressure =
        this->liquid_pressure_ + this->liquid_pressure_acceleration_ * dt;
    double FLIP_gas_pressure =
        this->gas_pressure_ + this->gas_pressure_acceleration_ * dt;

    this->FLIP_liquid_pressure_ = FLIP_liquid_pressure;
    this->FLIP_gas_pressure_ = FLIP_gas_pressure;

    // Update particle pressure based on PIC value
    this->liquid_pressure_ = pic_t * PIC_liquid_pressure + 
                            (1. - pic_t) * FLIP_liquid_pressure;
    this->gas_pressure_ = pic_t * PIC_gas_pressure + 
                            (1. - pic_t) * FLIP_gas_pressure;

    this->pore_pressure_ = this->liquid_pressure_ * liquid_chi_ + 
                            this->gas_pressure_ * gas_chi_;
    // this->suction_pressure_ = PIC_gas_pressure_ - PIC_liquid_pressure_;
  }
}

// Compute updated velocity of the liquid
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::compute_updated_velocity(
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

// Compute updated hydrate saturation
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::compute_pore_pressure(double dt){
  if (this->material_id_ != 999) {
    try {
      const double K_0 = material_->template property<double>(std::string("K_0"));
      const double p_ref = material_->template property<double>(std::string("p_ref"));
  
      // Calculate hydrate stability pressure pe:
      // const double A_1 = 1000;
      // const double A_2 = 38.98;
      // const double A_3 = 8533.8;
      // const double p_e = A_1 * std::exp(A_2 - A_3 / (PIC_temperature_ + 273.15));
      // // Calculate kinetic dissociation rate Kd: 
      // const double A_geo = 1E6;
      // const double K_d = K_0;

      // Calculate hydrate stability pressure pe:  
      const double A_1 = 1.15;
      const double A_2 = 49.3185; 
      const double A_3 = 9459;
      const double p_e = A_1 * std::exp(A_2 - A_3 / (PIC_temperature_ + 273.15));
      // Calculate kinetic dissociation rate Kd:
      const double Delta_E = 78151;
      const double A_geo = 7.5E5;
      const double K_d = K_0 * std::exp(-Delta_E / gas_constant_ / (PIC_temperature_ + 273.15));

      double Rr_coeff_gas = 0;
      double Rr_coeff_liquid = 0;
      if ((PIC_pore_pressure_ + p_ref) < p_e) {
        // Calculate specific surface area As:
        const double A_s = porosity_ * A_geo * hydrate_saturation_;
        // Calculate reaction rate Rr
        Rr_coeff_gas = K_d * A_s * (
                    gas_molar_mass_ /gas_density_ - hydrate_molar_mass_ / hydrate_density_);
        Rr_coeff_liquid = K_d * A_s * (hydration_number_ * liquid_molar_mass_ / liquid_density_);
      }

      // if ((pore_pressure_ +  p_ref) > p_e) {
      //   // Calculate specific surface area As:
      //   const double A_s = liquid_saturation_ * gas_saturation_ * A_geo * (1 - hydrate_saturation_);
      //   // Calculate reaction rate Rr
      //   Rr_coeff_gas = K_d * A_s * (
      //               gas_molar_mass_ /gas_density_ - hydrate_molar_mass_ / hydrate_density_);
      //   Rr_coeff_liquid = K_d * A_s * (hydration_number_ * liquid_molar_mass_ / liquid_density_);   
      // }

      double K_ww = porosity_ * dSw_dpw_ + 
                    porosity_ * liquid_saturation_ * liquid_compressibility_ + 
                    dt * (liquid_chi_ * Rr_coeff_liquid);
      double K_wg = -porosity_ * dSw_dpw_ + 
                    dt * (gas_chi_ * Rr_coeff_liquid);
      double K_gg = porosity_ * gas_saturation_ / gas_density_ *
                    gas_molar_mass_ / gas_constant_ / (PIC_temperature_ + 273.15) +
                    porosity_ * dSw_dpw_ +
                    dt * (gas_chi_ * Rr_coeff_gas);
      double K_gw = -porosity_ * dSw_dpw_ +
                    dt * (liquid_chi_ * Rr_coeff_gas);

      double K = porosity_ * liquid_saturation_ * liquid_compressibility_ + 
                porosity_ * gas_saturation_ / gas_density_ * 
                gas_molar_mass_ / gas_constant_ / (PIC_temperature_ + 273.15) + 
                dt * (Rr_coeff_liquid + Rr_coeff_gas);

      this->liquid_strain_rate_= 
              this->compute_strain_rate(dn_dx_centroid_, shapefn_centroid_, 
                                        mpm::ParticlePhase::Liquid);
      this->gas_strain_rate_= 
              this->compute_strain_rate(dn_dx_centroid_, shapefn_centroid_, 
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

      double beta_m = solid_fraction_ * solid_expansivity_ +
                      liquid_fraction_ * liquid_expansivity_ +
                      hydrate_fraction_ * hydrate_expansivity_ +
                      gas_fraction_ / gas_density_ * PIC_gas_pressure_ * gas_molar_mass_ /
                      gas_constant_/ std::pow(PIC_temperature_ + 273.15, 2);

      double f_w = liquid_source_ / liquid_density_ -
                    liquid_saturation_ * solid_strain_rate - 
                    liquid_fraction_ * (liquid_strain_rate - solid_strain_rate); 

      double f_g = gas_source_ / gas_density_ + 
                    hydrate_source_ / hydrate_density_ -
                    (1 - liquid_saturation_) * solid_strain_rate -
                    gas_fraction_ * (gas_strain_rate - solid_strain_rate);

      double f = liquid_source_ / liquid_density_ + gas_source_ / gas_density_ + 
                    hydrate_source_ / hydrate_density_ +
                    beta_m * this->temperature_acceleration_ -
                    solid_strain_rate - 
                    liquid_fraction_ * (liquid_strain_rate - solid_strain_rate) -
                    gas_fraction_ * (gas_strain_rate - solid_strain_rate);


      Eigen::Matrix<double, 2, 2> K_matrix, K_matrix_inverse;
      K_matrix(0,0) = K_ww;
      K_matrix(0,1) = K_wg;
      K_matrix(1,0) = K_gw;
      K_matrix(1,1) = K_gg;

    if (K_gg > 1E-16) {
      K_matrix_inverse = K_matrix.inverse();

      this->liquid_pressure_acceleration_ = K_matrix_inverse(0, 0) * f_w + 
                                            K_matrix_inverse(0, 1) * f_g;
      this->gas_pressure_acceleration_ = K_matrix_inverse(1, 0) * f_w + 
                                        K_matrix_inverse(1, 1) * f_g;
    } else {
        this->liquid_pressure_acceleration_ = f_w / K_ww;
        this->gas_pressure_acceleration_ = 0;
    }

      // double pressure_acceleration = f/K;
      // this->liquid_pressure_acceleration_ = pressure_acceleration;
      // this->gas_pressure_acceleration_ = pressure_acceleration;

      this->liquid_pressure_ += this->liquid_pressure_acceleration_ * dt;
      this->gas_pressure_ += this->gas_pressure_acceleration_ * dt;

      // if (this->free_surface()) {
      //   const double suction = liquid_material_->template 
      //                       property<double>(std::string("suction"));
      //   this->suction_pressure_ = this->ini_suction_pressure_ + suction;
      //   this->gas_pressure_ = material_->template property<double>(std::string("p_ref"));
      //   this->liquid_pressure_ = this->gas_pressure_ - this->suction_pressure_;
      // }

      this->pore_pressure_ = this->liquid_chi_ * this->liquid_pressure_ +
                            this->gas_chi_ * this->gas_pressure_;
      this->PIC_liquid_pressure_ = this->liquid_pressure_;
      this->PIC_gas_pressure_ = this->gas_pressure_;
      this->PIC_pore_pressure_ = this->pore_pressure_; 

    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__, 
                      exception.what());
    }
  }
}

// Assign a liquid material to particle
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::update_particle_volume() {
  if (this->material_id_ != 999) {
    try {
      // Solid phase  
      mpm::Particle<Tdim>::update_particle_volume();

      // SOLID PHASE
      this->solid_heat_capacity_ = mass_ * solid_specific_heat_;

      // HYDRATE PHASE
      this->hydrate_fraction_ = PIC_porosity_ * PIC_hydrate_saturation_;
      this->hydrate_volume_ = hydrate_fraction_ * PIC_volume_;
      this->hydrate_mass_ = hydrate_volume_ * hydrate_density_;
      this->hydrate_mass_density_ = hydrate_fraction_ * hydrate_density_;
      this->hydrate_heat_capacity_ = hydrate_mass_ * hydrate_specific_heat_;

      // LIQUID PHASE
      this->liquid_fraction_ = PIC_porosity_ * liquid_saturation_;
      this->liquid_volume_ = liquid_fraction_ * PIC_volume_;
      this->liquid_mass_ = liquid_volume_ * liquid_density_;
      this->liquid_mass_density_ = liquid_fraction_ * liquid_density_;
      this->liquid_chi_ = liquid_saturation_ / (liquid_saturation_ + gas_saturation_);
      this->liquid_heat_capacity_ = liquid_mass_ * liquid_specific_heat_;

      // GAS PHASE
      this->gas_fraction_ = PIC_porosity_ * gas_saturation_;
      this->gas_volume_ = gas_fraction_ * PIC_volume_;
      this->gas_mass_ = gas_volume_ * gas_density_;
      this->gas_mass_density_ = gas_fraction_ * gas_density_;
      this->gas_chi_ = 1 - liquid_chi_;
      this->gas_heat_capacity_ = gas_mass_ * gas_specific_heat_;

      //Mixture
      this->mixture_mass_ = mass_ + hydrate_mass_ + liquid_mass_ + gas_mass_;

    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__, 
                      exception.what());
    }
  }
}

// Calculate absoulute/relative permeability
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::update_permeability() {
  if (this->material_id_ != 999) {
    try {   
      double k_phi = 1;
      double k_sh = 1;
      double k_a = 1;
      double k_r_liquid = 1.;
      double k_r_gas = 1.; 

      // Calculate porosity-dependent permeability, k_phi
      k_phi = std::pow(porosity_ / ini_porosity_, 1.5) * 
                      std::pow((1 - ini_porosity_) / (1 - PIC_porosity_), 3);

      // Calculate hydrate saturation-dependent permeability, k_sh
      // k_sh = std::pow(1 - hydrate_saturation_, 19./6.);
      k_sh = 0.01;// 0.0075;
      double Shc = 1E-4;
      if (hydrate_saturation_ < Shc) {
        double k_sh = 0.01 + (1. - 0.01) * (Shc - PIC_hydrate_saturation_) / Shc;
      }

      // Calculate absolute permeability, k_a
      k_a = intrinsic_permeability_ * k_phi * k_sh;

      // Calculate relative permeability
      // k_r_liquid = std::pow((liquid_saturation_ - 0.1) / (1 - 0.1), 0.82); 
      // k_r_gas = std::pow(gas_saturation_ / (1 - 0.1), 2.11);

      k_r_liquid = std::pow((liquid_chi_ - liquid_saturation_res_) / 
                      (1 - liquid_saturation_res_ - gas_saturation_res_), 4);
      k_r_gas = std::pow((gas_chi_ ) / 
                      (1 - liquid_saturation_res_ - gas_saturation_res_), 2);

      // Calculate permeability
      this->liquid_permeability_ = k_a *  k_r_liquid; 
      this->gas_permeability_ = k_a * k_r_gas;

      if (nodes_.size() >=2 ) {
        double cell_size = abs(nodes_[0]->coordinates()[0] - 
                                nodes_[1]->coordinates()[0]);
        this->liquid_critical_time_ = cell_size * cell_size * liquid_fraction_ * 
            liquid_compressibility_ / liquid_permeability_ * liquid_viscosity_;
        this->gas_critical_time_ = cell_size * cell_size * gas_fraction_ * 
            gas_molar_mass_ / gas_density_ / gas_constant_ / 
            (PIC_temperature_ + 273.15) / liquid_permeability_ * liquid_viscosity_;
      }
    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__, 
                      exception.what());
    }
  }
}

// Calculate source term in balance equations
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::update_source_term() {
  if (this->material_id_ != 999) {
    try {
      // Calculate hydrate stability pressure pe:
      const double A_1 = 1.15;
      const double A_2 = 49.3185; 
      const double A_3 = 9459;
      // const double A_1 = 1000;
      // const double A_2 = 38.98; 
      // const double A_3 = 8533.8;

      const double p_e = A_1 * std::exp(A_2 - A_3 / (PIC_temperature_ + 273.15));
      double Reaction_rate = 0;
      const double p_ref = material_->template property<double>(std::string("p_ref"));

      const double K_0 = material_->template property<double>(std::string("K_0"));
      const double Delta_E = 78151;
      const double A_geo = 7.5E5;
      // Calculate kinetic dissociation rate Kd:     
      const double K_d = K_0 * std::exp(-Delta_E / gas_constant_ / (PIC_temperature_ + 273.15));

      // const double A_geo = 1E6;
      // // Calculate kinetic dissociation rate Kd:     
      // const double K_d = K_0;

      if ((PIC_pore_pressure_ + p_ref) < p_e) {
        // Calculate specific surface area As:
        const double A_s = PIC_porosity_ * A_geo * PIC_hydrate_saturation_;
        // Calculate reaction rate Rr
        Reaction_rate = K_d * A_s * (p_e - (PIC_pore_pressure_ +  p_ref));
      }

      // if ((pore_pressure_ +  p_ref) > p_e) {
      //   // Calculate specific surface area As:
      //   const double A_s = liquid_saturation_ * gas_saturation_ * A_geo * (1 - hydrate_saturation_);
      //   // Calculate reaction rate Rr
      //   Reaction_rate = K_d * A_s * (p_e - (pore_pressure_ +  p_ref));
      // }

      // Calculate mass source term
      hydrate_source_ = -hydrate_molar_mass_ * Reaction_rate;
      liquid_source_ = hydration_number_ * liquid_molar_mass_ * Reaction_rate;
      gas_source_ = gas_molar_mass_ * Reaction_rate;  
      // Calculate energy source term
      const double C0 = 446120;
      const double C1 = 132.638;
      // const double C0 = 473631.8;
      // const double C1 = 140.117;
      heat_source_ = (C0 - C1 * (PIC_temperature_ + 273.15)) * hydrate_source_;


    } catch (std::exception& exception) {
        console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__, 
                        exception.what());
    }
  }
}

// Compute updated hydrate saturation
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::update_hydrate_saturation(double dt){
  if (this->material_id_ != 999) {
    try {
      // // Forward euler
      // double Sh = (dt * hydrate_source_ / porosity_ / hydrate_density_ + 
      //             hydrate_saturation_) /
      //             (1 - dt * ((1. / porosity_ - 1. ) * solid_expansivity_ + 
      //             hydrate_expansivity_) * temperature_acceleration_- 
      //             1. / porosity_ * strain_rate_.head(Tdim).sum());

      // Backward euler
      double Sh_rate = (hydrate_source_ / PIC_porosity_ / hydrate_density_ - 
                      PIC_hydrate_saturation_ / PIC_porosity_ * strain_rate_.head(Tdim).sum()) /
                      (1 + PIC_porosity_ * strain_rate_.head(Tdim).sum() * dt);
      this->hydrate_saturation_ += Sh_rate * dt;
      this->PIC_hydrate_saturation_ = this->hydrate_saturation_;
    } catch (std::exception& exception) {
        console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__, 
                        exception.what());
    }
  }
}

// Compute updated liquid saturation
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::update_liquid_saturation(double dt){
  if (this->material_id_ != 999) {
    try {
      this->suction_pressure_ = PIC_gas_pressure_ - PIC_liquid_pressure_;
      // if (this->suction_pressure_ < 0) this->suction_pressure_ = 0;

      double para_p0 = liquid_material_->template 
                    property<double>(std::string("para_p0"));
      double para_m = liquid_material_->template 
                    property<double>(std::string("para_m"));
      double A = suction_pressure_ / para_p0;
      double B = 1.0 + std::pow(abs(A), 1.0 / (1.0 - para_m));
      double Se = std::pow(B, -para_m);

      this->liquid_saturation_ = (Se * (1.0 - liquid_saturation_res_ - gas_saturation_res_) +
                        liquid_saturation_res_) * (1 - PIC_hydrate_saturation_);

      double dA_dpw = -1.0 / para_p0;
      double dB_dA = 1.0 / (1.0 - para_m) * std::pow(abs(A), (1.0 / (1.0 - para_m) - 1.0));
      double dSe_dB = -para_m * std::pow(B, (-para_m - 1.0));
      double dSw_dSe = (1 - PIC_hydrate_saturation_) * (1 - liquid_saturation_res_ - 
                                                    gas_saturation_res_) * Se; 
      double dSw_dpw = dSw_dSe * dSe_dB * dB_dA * dA_dpw;

      this->dSw_dpw_ = dSw_dpw;

      this->gas_saturation_ = 
                  1.0 - this->liquid_saturation_ - this->PIC_hydrate_saturation_;

      this->liquid_chi_ = this->liquid_saturation_ / ( 1 -this->PIC_hydrate_saturation_);
      this->gas_chi_ = 1 - this->liquid_chi_;

    } catch (std::exception& exception) {
        console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__, 
                        exception.what());
    }
  }
}

// Compute updated density
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::update_particle_density(double dt){
  if (this->material_id_ != 999) {  
    try {
      double p_ref = material_->template property<double>(std::string("p_ref"));
      this->gas_density_ = gas_molar_mass_ * (PIC_gas_pressure_ + p_ref) / 
                            gas_constant_ / (PIC_temperature_ + 273.15);

      // Check if NaN
      if (!(this->gas_density_ > 0.)) gas_density_ = 1.0;

      // this->density_ /= 1 + dt * solid_expansivity_ * temperature_acceleration_; 
      // this->hydrate_density_ /= 1 + dt * hydrate_expansivity_ * temperature_acceleration_; 
      this->liquid_density_ /= 1 + liquid_expansivity_ * temperature_increment_ - 
                                  liquid_pressure_increment_ * liquid_compressibility_;

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
bool mpm::MHBSParticle<Tdim>::assign_particle_traction(unsigned direction,
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
bool mpm::MHBSParticle<Tdim>::assign_particle_liquid_velocity_constraint(
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
void mpm::MHBSParticle<Tdim>::apply_particle_liquid_velocity_constraints() {
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
bool mpm::MHBSParticle<Tdim>::assign_particle_pore_pressure_constraint(
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
void mpm::MHBSParticle<Tdim>::apply_particle_pore_pressure_constraints(
      double pore_pressure) {
  // Set particle temperature constraint
  this->pore_pressure_ = pore_pressure;
  this->liquid_pressure_ = pore_pressure;
  this->gas_pressure_ = pore_pressure;
  this->PIC_liquid_pressure_ = pore_pressure;
  this->PIC_gas_pressure_ = pore_pressure;  
  this->set_pressure_constraint_ = true;
}

// Overwrite node velocity to get strain correct
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::map_moving_rigid_velocity_to_nodes(
    unsigned dir, double velocity, double dt) noexcept {   
  if (this->material_id_ == 999){  
    for (unsigned i = 0; i < 4; ++i) {
      nodes_[i]->assign_velocity_from_rigid(dir, velocity, dt);
    }       
  }
}


//==============================================================================
//  PART 4: SMOOTHING TECHNIQUE
//==============================================================================

// Map particle pore liquid pressure to nodes
template <unsigned Tdim>
bool mpm::MHBSParticle<Tdim>::map_pore_pressure_to_nodes(
    double current_time) noexcept {
  if (this->material_id_ != 999) {
    // Check if particle mass is set
    assert(liquid_mass_ != std::numeric_limits<double>::max());

    bool status = true;
    // Map particle liquid mass and pore pressure to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i){
      nodes_[i]->update_mass_pressure(mpm::ParticlePhase::Liquid,
                                      shapefn_[i] * liquid_mass_ * liquid_pressure_,
                                      current_time);
      nodes_[i]->update_mass_pressure(mpm::ParticlePhase::Liquid,
                                      shapefn_[i] * liquid_mass_ * gas_pressure_,
                                      current_time);
    }
    return status;
  }
}

// Map particle pore liquid pressure to nodes
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::map_mass_pressure_to_nodes() {
  if (this->material_id_ != 999) {  
    // Check if particle mass is set
    assert(liquid_mass_ != std::numeric_limits<double>::max());
    if (this->material_id_ != 999){
      // Map particle liquid mass and pore pressure to nodes
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        nodes_[i]->update_pressure(true, mpm::ParticlePhase::Gas,
                                        shapefn_[i] * 1.0 * gas_pressure_);
        nodes_[i]->update_pressure(true, mpm::ParticlePhase::Liquid,
                                  shapefn_[i] *  1.0 * liquid_pressure_);
      }
    }
  }
}

// Compute pore liquid pressure smoothing based on nodal pressure
template <unsigned Tdim>
bool mpm::MHBSParticle<Tdim>::compute_pore_pressure_smoothing() noexcept {
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
    this->PIC_liquid_pressure_= liquid_pressure; 
    this->PIC_gas_pressure_ = gas_pressure;

    // if (this->free_surface()) {
    //   this->liquid_pressure_ = 0.0;
    //   this->gas_pressure_ = 0.0;
    // }

    // this->PIC_liquid_pressure_ = this->liquid_pressure_;
    // this->PIC_gas_pressure_ = this->gas_pressure_;
    this->PIC_pore_pressure_ = this->liquid_chi_ * this->PIC_liquid_pressure_ +
                          this->gas_chi_ * this->PIC_gas_pressure_;

  }
  return status;
}

// Map particle pore liquid pressure to nodes
template <unsigned Tdim>
void mpm::MHBSParticle<Tdim>::map_scalers_to_nodes() {
  if (this->material_id_ != 999) {  
    // Check if particle mass is set
    assert(liquid_mass_ != std::numeric_limits<double>::max());
    if (this->material_id_ != 999){
      if (Tdim == 2) {
        // Compute deviatoric strain
        const double volumetric_strain = this->strain_.head(2).sum();
        Eigen::Matrix<double, 6, 1> strain_dev = this->strain_;
        strain_dev.head(2).noalias() += 
                          -Eigen::Vector2d::Constant(volumetric_strain / 2);
        const double deviatoric_strain = std::pow(2. * 
                            (strain_dev.dot(strain_dev) + 
                            strain_dev.tail(3).dot(strain_dev.tail(3))), 0.5);
        // Compute deviatoric stress
        const double mean_stress = this->stress_.head(2).sum() / 2.0;
        Eigen::Matrix<double, 6, 1> stress_dev = this->stress_;
        stress_dev.head(2).noalias() += 
                              -Eigen::Vector2d::Constant(mean_stress);
        const double deviatoric_stress = std::pow(0.5 * 
                            (stress_dev.dot(stress_dev) + 
                            stress_dev.tail(3).dot(stress_dev.tail(3))), 0.5);
        this->deviatoric_strain_ = deviatoric_strain;
        // Map particle liquid mass and pore pressure to nodes
        for (unsigned i = 0; i < nodes_.size(); ++i) {
          nodes_[i]->update_scalers(true, 0, shapefn_[i] * 1.0);
          nodes_[i]->update_scalers(true, 1, shapefn_[i] * 1.0 * porosity_);
          nodes_[i]->update_scalers(true, 2, shapefn_[i] * 1.0 * volume_);
          nodes_[i]->update_scalers(true, 3, shapefn_[i] * 1.0 * hydrate_saturation_);
          nodes_[i]->update_scalers(true, 4, shapefn_[i] * 1.0 * volumetric_strain);
          nodes_[i]->update_scalers(true, 5, shapefn_[i] * 1.0 * deviatoric_strain);
          nodes_[i]->update_scalers(true, 6, shapefn_[i] * 1.0 * mean_stress);
          nodes_[i]->update_scalers(true, 7, shapefn_[i] * 1.0 * deviatoric_stress);
        }
      }
    }
  }
}

// Compute pore liquid pressure smoothing based on nodal pressure
template <unsigned Tdim>
bool mpm::MHBSParticle<Tdim>::compute_scalers_smoothing() noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  bool status = true;  
  if (this->material_id_ != 999) {
    double porosity = 0;
    double volume = 0;
    double hydrate_saturation = 0;
    double volumetric_strain = 0;
    double deviatoric_strain = 0;
    double mean_stress = 0;
    double deviatoric_stress = 0;
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      porosity += shapefn_(i) * nodes_[i]->smoothed_scalers(1);
      volume += shapefn_(i) * nodes_[i]->smoothed_scalers(2);
      hydrate_saturation += shapefn_(i) * nodes_[i]->smoothed_scalers(3);
      volumetric_strain += shapefn_(i) * nodes_[i]->smoothed_scalers(4);
      deviatoric_strain += shapefn_(i) * nodes_[i]->smoothed_scalers(5);
      mean_stress += shapefn_(i) * nodes_[i]->smoothed_scalers(6);
      deviatoric_stress += shapefn_(i) * nodes_[i]->smoothed_scalers(7);
    }

    // Update pore liquid pressure to interpolated nodal pressure
    this->PIC_porosity_ = porosity; 
    this->PIC_volume_ = volume; 
    this->PIC_hydrate_saturation_ = hydrate_saturation;
    this->PIC_volumetric_strain_ = volumetric_strain; 
    this->PIC_deviatoric_strain_ = deviatoric_strain; 
    this->PIC_mean_stress_ = mean_stress;
    this->PIC_deviatoric_stress_ = deviatoric_stress;
  }
  return status;
}