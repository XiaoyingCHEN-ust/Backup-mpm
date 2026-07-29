// Construct a three phase particle with id and coordinates
template <unsigned Tdim>
mpm::ThreePhaseParticleNew<Tdim>::ThreePhaseParticleNew(Index id,
                                                  const VectorDim& coord)
  : mpm::Particle<Tdim>(id, coord) {
  this->initialise_liquid_gas_phases();
  // Set material pointer to null
  liquid_material_ = nullptr;
  // Logger
  std::string logger =
      "ThreePhaseParticleNew" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//==============================================================================
// ASSIGN INITIAL CONDITIONS
//==============================================================================

// Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::ThreePhaseParticleNew<Tdim>::initialise_particle(const HDF5Particle& particle) {
    // Derive from particle
    mpm::Particle<Tdim>::initialise_particle(particle);

    // MIXTURE data
    this->liquid_material_id_ = particle.liquid_material_id;
    this->pore_pressure_ = particle.pore_pressure;
    this->suction_pressure_ = particle.suction_pressure;
    this->mixture_mass_ = particle.mixture_mass;
    this->dSw_dpw_ = particle.dSw_dpw;
    this->total_stress_[0] = particle.total_stress_xx;
    this->total_stress_[1] = particle.total_stress_yy;
    this->total_stress_[2] = particle.total_stress_zz;
    this->total_stress_[3] = particle.total_stress_tau_xy;
    this->total_stress_[4] = particle.total_stress_tau_yz;
    this->total_stress_[5] = particle.total_stress_tau_xz;

    // Vector properties
    Eigen::Vector3d liquid_velocity, liquid_acceleration;
    liquid_velocity << particle.liquid_velocity_x, particle.liquid_velocity_y, particle.liquid_velocity_z;
    liquid_acceleration  << particle.liquid_acceleration_x, particle.liquid_acceleration_y, particle.liquid_acceleration_z;
    for (unsigned i = 0; i < Tdim; ++i) {
      this->liquid_velocity_[i] = liquid_velocity[i];
      this->liquid_acceleration_[i] = liquid_acceleration[i];
    } 

    // LIQUID PHASE data
    this->liquid_strain_[0] = particle.liquid_strain_xx;
    this->liquid_strain_[1] = particle.liquid_strain_yy;
    this->liquid_strain_[2] = particle.liquid_strain_zz;
    this->liquid_strain_[3] = particle.liquid_strain_gamma_xy;
    this->liquid_strain_[4] = particle.liquid_strain_gamma_yz;
    this->liquid_strain_[5] = particle.liquid_strain_gamma_xz;
    this->liquid_strain_rate_[0] = particle.liquid_strain_rate_xx;
    this->liquid_strain_rate_[1] = particle.liquid_strain_rate_yy;
    this->liquid_strain_rate_[2] = particle.liquid_strain_rate_zz;
    this->liquid_strain_rate_[3] = particle.liquid_strain_rate_xy;
    this->liquid_strain_rate_[4] = particle.liquid_strain_rate_yz;
    this->liquid_strain_rate_[5] = particle.liquid_strain_rate_xz;
    this->liquid_saturation_        = particle.liquid_saturation;
    this->effective_saturation_     = particle.effective_saturation;
    this->liquid_chi_               = particle.liquid_chi;
    this->liquid_fraction_          = particle.liquid_fraction;
    this->liquid_density_           = particle.liquid_density;
    this->liquid_volume_            = particle.liquid_volume;
    this->liquid_mass_              = particle.liquid_mass;
    this->liquid_mass_density_      = particle.liquid_mass_density;
    this->liquid_pressure_          = particle.liquid_pressure;
    this->liquid_pressure_acceleration_ = particle.liquid_pressure_acceleration;
    this->liquid_volumetric_strain_ = particle.liquid_volumetric_strain;
    this->liquid_permeability_      = particle.liquid_permeability;
    this->PIC_liquid_pressure_      = particle.PIC_liquid_pressure;
    this->FLIP_liquid_pressure_     = particle.FLIP_liquid_pressure;

    // Vector properties
    Eigen::Vector3d gas_velocity, gas_acceleration, pgravity;
    gas_velocity << particle.gas_velocity_x, particle.gas_velocity_y, particle.gas_velocity_z;
    gas_acceleration  << particle.gas_acceleration_x, particle.gas_acceleration_y, particle.gas_acceleration_z;
    pgravity  << particle.pgravity_x, particle.pgravity_y, particle.pgravity_z;
    for (unsigned i = 0; i < Tdim; ++i) {
      this->gas_velocity_[i] = gas_velocity[i];
      this->gas_acceleration_[i] = gas_acceleration[i];
      this->pgravity_[i] = pgravity[i];
    } 

    // GAS PHASE data
    this->gas_strain_[0] = particle.gas_strain_xx;
    this->gas_strain_[1] = particle.gas_strain_yy;
    this->gas_strain_[2] = particle.gas_strain_zz;
    this->gas_strain_[3] = particle.gas_strain_gamma_xy;
    this->gas_strain_[4] = particle.gas_strain_gamma_yz;
    this->gas_strain_[5] = particle.gas_strain_gamma_xz;
    this->gas_strain_rate_[0] = particle.gas_strain_rate_xx;
    this->gas_strain_rate_[1] = particle.gas_strain_rate_yy;
    this->gas_strain_rate_[2] = particle.gas_strain_rate_zz;
    this->gas_strain_rate_[3] = particle.gas_strain_rate_xy;
    this->gas_strain_rate_[4] = particle.gas_strain_rate_yz;
    this->gas_strain_rate_[5] = particle.gas_strain_rate_xz;
    this->gas_saturation_            = particle.gas_saturation;
    this->gas_fraction_              = particle.gas_fraction;
    this->gas_chi_                   = particle.gas_chi;
    this->gas_density_               = particle.gas_density;
    this->gas_volume_                = particle.gas_volume;
    this->gas_mass_                  = particle.gas_mass;
    this->gas_mass_density_          = particle.gas_mass_density;
    this->gas_pressure_              = particle.gas_pressure;
    this->gas_pressure_acceleration_ = particle.gas_pressure_acceleration;
    this->PIC_gas_pressure_          = particle.PIC_gas_pressure;
    this->FLIP_gas_pressure_         = particle.FLIP_gas_pressure;
    this->gas_pressure_increment_    = particle.gas_pressure_increment;
    this->gas_volumetric_strain_     = particle.gas_volumetric_strain;
    this->gas_permeability_          = particle.gas_permeability;
        
    // WAVE data
    this->wave_pressure_ = particle.wave_pressure;
    this->gamma_b_ = particle.gamma_b;

    return true;
}

// Initialise particle HDF5 data and material
template <unsigned Tdim>
bool mpm::ThreePhaseParticleNew<Tdim>::initialise_particle(
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

// Initialise liquid and gas phase particle properties
template <unsigned Tdim>
void mpm::ThreePhaseParticleNew<Tdim>::initialise_liquid_gas_phases() {
  // Mixture
    set_mixture_traction_ = false;
    set_pressure_constraint_ = false;
    pore_pressure_ = 0.;
    suction_pressure_ = 0.;

  // Liquid
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
    liquid_critical_time_ = 0.;

  // Gas
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
    gas_critical_time_ = 0.;
    dSw_dpw_ = 0.;

  // Wave
    seabed_surface_x_.assign(this->Nx_, 0.0);
    water_depth_.assign(this->Nx_, 0.0);
    beta_k_.assign(this->Nx_, 0.0);
    wave_height_.assign(this->Nx_, 0.0);
    wave_P0_.assign(this->Nx_, 0.0);
    wave_phi_.assign(this->Nx_, 0.0);
    wave_pf_.assign(this->Nx_, 0.0);
    gamma_b_ = 0.78;

    
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
      {"liquid_volumetric_strains",[&]() {return this->liquid_volumetric_strain_;}},
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
mpm::HDF5Particle mpm::ThreePhaseParticleNew<Tdim>::hdf5() {
    // Derive from particle
    auto particle_data = mpm::Particle<Tdim>::hdf5();

    // MIXTURE data
    particle_data.liquid_material_id = this->liquid_material_id_;
    particle_data.pore_pressure = this->pore_pressure_;
    particle_data.suction_pressure = this->suction_pressure_;
    particle_data.mixture_mass = this->mixture_mass_;
    particle_data.dSw_dpw = this->dSw_dpw_;
    particle_data.total_stress_xx = this->total_stress_[0];
    particle_data.total_stress_yy = this->total_stress_[1];
    particle_data.total_stress_zz = this->total_stress_[2];
    particle_data.total_stress_tau_xy = this->total_stress_[3];
    particle_data.total_stress_tau_yz = this->total_stress_[4];
    particle_data.total_stress_tau_xz = this->total_stress_[5];

    // Vector properties
    Eigen::Vector3d liquid_velocity, liquid_acceleration;
    liquid_acceleration.setZero();
    liquid_velocity.setZero();
    for (unsigned i = 0; i < Tdim; ++i) {
      liquid_velocity[i] = this->liquid_velocity_[i];
      liquid_acceleration[i] = this->liquid_acceleration_[i];
    } 

    // LIQUID PHASE data
    particle_data.liquid_velocity_x = liquid_velocity[0];
    particle_data.liquid_velocity_y = liquid_velocity[1];
    particle_data.liquid_velocity_z = liquid_velocity[2];
    particle_data.liquid_acceleration_x = liquid_acceleration[0];
    particle_data.liquid_acceleration_y = liquid_acceleration[1];
    particle_data.liquid_acceleration_z = liquid_acceleration[2];
    particle_data.liquid_strain_xx = this->liquid_strain_[0];
    particle_data.liquid_strain_yy = this->liquid_strain_[1];
    particle_data.liquid_strain_zz = this->liquid_strain_[2];
    particle_data.liquid_strain_gamma_xy = this->liquid_strain_[3];
    particle_data.liquid_strain_gamma_yz = this->liquid_strain_[4];
    particle_data.liquid_strain_gamma_xz = this->liquid_strain_[5];
    particle_data.liquid_strain_rate_xx = this->liquid_strain_rate_[0];
    particle_data.liquid_strain_rate_yy = this->liquid_strain_rate_[1];
    particle_data.liquid_strain_rate_zz = this->liquid_strain_rate_[2];
    particle_data.liquid_strain_rate_xy = this->liquid_strain_rate_[3];
    particle_data.liquid_strain_rate_yz = this->liquid_strain_rate_[4];
    particle_data.liquid_strain_rate_xz = this->liquid_strain_rate_[5];
    particle_data.liquid_saturation = this->liquid_saturation_;
    particle_data.effective_saturation = this->effective_saturation_;
    particle_data.liquid_chi = this->liquid_chi_;
    particle_data.liquid_fraction = this->liquid_fraction_;
    particle_data.liquid_density = this->liquid_density_;
    particle_data.liquid_volume = this->liquid_volume_;
    particle_data.liquid_mass = this->liquid_mass_;
    particle_data.liquid_mass_density = this->liquid_mass_density_;
    particle_data.liquid_pressure = this->liquid_pressure_;
    particle_data.liquid_pressure_acceleration = this->liquid_pressure_acceleration_;
    particle_data.liquid_volumetric_strain = this->liquid_volumetric_strain_;
    particle_data.liquid_permeability = this->liquid_permeability_;
    particle_data.PIC_liquid_pressure = this->PIC_liquid_pressure_;
    particle_data.FLIP_liquid_pressure = this->FLIP_liquid_pressure_;

    // Vector properties
    Eigen::Vector3d gas_velocity, gas_acceleration, pgravity;
    gas_acceleration.setZero();
    gas_velocity.setZero();
    pgravity.setZero();
    for (unsigned i = 0; i < Tdim; ++i) {
      gas_velocity[i] = this->gas_velocity_[i];
      gas_acceleration[i] = this->gas_acceleration_[i];
      pgravity[i] = this->pgravity_[i];
    } 

    // GAS PHASE data
    particle_data.gas_velocity_x = gas_velocity[0];
    particle_data.gas_velocity_y = gas_velocity[1];
    particle_data.gas_velocity_z = gas_velocity[2];
    particle_data.gas_acceleration_x = gas_acceleration[0];
    particle_data.gas_acceleration_y = gas_acceleration[1];
    particle_data.gas_acceleration_z = gas_acceleration[2];
    particle_data.pgravity_x = pgravity[0];
    particle_data.pgravity_y = pgravity[1];
    particle_data.pgravity_z = pgravity[2];
    particle_data.gas_strain_xx = this->gas_strain_[0];
    particle_data.gas_strain_yy = this->gas_strain_[1];
    particle_data.gas_strain_zz = this->gas_strain_[2];
    particle_data.gas_strain_gamma_xy = this->gas_strain_[3];
    particle_data.gas_strain_gamma_yz = this->gas_strain_[4];
    particle_data.gas_strain_gamma_xz = this->gas_strain_[5];
    particle_data.gas_strain_rate_xx = this->gas_strain_rate_[0];
    particle_data.gas_strain_rate_yy = this->gas_strain_rate_[1];
    particle_data.gas_strain_rate_zz = this->gas_strain_rate_[2];
    particle_data.gas_strain_rate_xy = this->gas_strain_rate_[3];
    particle_data.gas_strain_rate_yz = this->gas_strain_rate_[4];
    particle_data.gas_strain_rate_xz = this->gas_strain_rate_[5];
    particle_data.gas_saturation = this->gas_saturation_;
    particle_data.gas_fraction = this->gas_fraction_;
    particle_data.gas_chi = this->gas_chi_;
    particle_data.gas_density = this->gas_density_;
    particle_data.gas_volume = this->gas_volume_;
    particle_data.gas_mass = this->gas_mass_;
    particle_data.gas_mass_density = this->gas_mass_density_;
    particle_data.gas_pressure = this->gas_pressure_;
    particle_data.gas_pressure_acceleration = this->gas_pressure_acceleration_;
    particle_data.PIC_gas_pressure = this->PIC_gas_pressure_;
    particle_data.FLIP_gas_pressure = this->FLIP_gas_pressure_;
    particle_data.gas_pressure_increment = this->gas_pressure_increment_;
    particle_data.gas_volumetric_strain = this->gas_volumetric_strain_;
    particle_data.gas_permeability = this->gas_permeability_;
        
    // WAVE data
    particle_data.wave_pressure = this->wave_pressure_;
    particle_data.gamma_b = this->gamma_b_;
    

    return particle_data;
}

// Assign a liquid material to particle
template <unsigned Tdim>
bool mpm::ThreePhaseParticleNew<Tdim>::assign_liquid_material(
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

// Compute mass of particle (solid, fluid and gas)
template <unsigned Tdim>
void mpm::ThreePhaseParticleNew<Tdim>::compute_mass() {
  mpm::Particle<Tdim>::compute_mass();
  this->assign_initial_properties();
}

// Assign initial properties to particle
template <unsigned Tdim>
bool mpm::ThreePhaseParticleNew<Tdim>::assign_initial_properties() {
  bool status = true;
  try {
    this->ini_porosity_ = 
          material_->template property<double>(std::string("porosity"));
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

    // Wave properties
    this->wave_pressure_ = liquid_material_->template
                  property<bool>(std::string("wave_pressure"));
    this->domain_length_x_ = liquid_material_->template
                  property<double>(std::string("domain_length_x"));
    this->Nx_ = liquid_material_->template
                  property<double>(std::string("Nx"));
    this->sea_level_ = liquid_material_->template
                  property<double>(std::string("sea_level"));   
    this->depth_left_ = liquid_material_->template
                  property<double>(std::string("depth_left"));
    this->depth_right_ = liquid_material_->template
                  property<double>(std::string("depth_right"));
    this->wave_height_ini_ = liquid_material_->template
                  property<double>(std::string("wave_height_ini_"));
    this->wave_period_ = liquid_material_->template
                  property<double>(std::string("wave_period"));       

    // Time-dependent propeties
    ini_gas_saturation_ = liquid_material_->template
                  property<double>(std::string("gas_saturation"));
    ini_gas_viscosity_ = liquid_material_->template
                  property<double>(std::string("gas_viscosity"));
    double p_ref = material_->template property<double>(std::string("p_ref"));
    this->gas_saturation_ = ini_gas_saturation_;
    this->gas_pressure_ = 0.0;
    this->gas_density_ = gas_molar_mass_ * (gas_pressure_ + p_ref) /
                          gas_constant_ / (PIC_temperature_ + 273.15);
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
    const double initial_mobile_saturation =
        1.0 - liquid_saturation_res_ - gas_saturation_res_;
    if (initial_mobile_saturation <= 0.0)
      throw std::runtime_error(
          "Initial state requires positive mobile saturation range");
    this->effective_saturation_ = std::min(
        1.0, std::max(0.0,
            (this->liquid_saturation_ - liquid_saturation_res_) /
            initial_mobile_saturation));
    // Preserve the historical 20 psi initial suction unless the input file
    // explicitly supplies a value or requests consistency with the SWRC.
    this->suction_pressure_ = liquid_material_->template property_or<double>(
        std::string("initial_suction"), 137340.0);

    const bool initial_suction_from_swrc =
        liquid_material_->template property_or<bool>(
            std::string("initial_suction_from_swrc"), false);
    if (initial_suction_from_swrc) {
      if (para_p0 <= 0.0 || para_m <= 0.0 || para_m >= 1.0)
        throw std::runtime_error(
            "Initial SWRC suction requires para_p0 > 0 and 0 < para_m < 1");

      const double saturation_range =
          1.0 - liquid_saturation_res_ - gas_saturation_res_;
      if (saturation_range <= 0.0)
        throw std::runtime_error(
            "Initial SWRC suction requires positive mobile saturation range");

      double initial_effective_saturation =
          (liquid_saturation_ - liquid_saturation_res_) / saturation_range;
      initial_effective_saturation = std::min(
          1.0 - 1.0e-12,
          std::max(1.0e-12, initial_effective_saturation));

      this->suction_pressure_ = para_p0 * std::pow(
          std::pow(initial_effective_saturation, -1.0 / para_m) - 1.0,
          1.0 - para_m);
    }
    this->suction_pressure_ = std::abs(this->suction_pressure_);

    // Initialise the retention-curve tangent before the first pressure
    // update. Leaving dSw_dpw_ undefined makes the first coupling matrix
    // depend on uninitialised memory.
    this->dSw_dpw_ = 0.0;
    if (para_p0 > 0.0 && para_m > 0.0 && para_m < 1.0) {
      const double A = this->suction_pressure_ / para_p0;
      const double exponent = 1.0 / (1.0 - para_m);
      const double abs_A = std::max(std::abs(A), 1.0e-12);
      const double B = 1.0 + std::pow(abs_A, exponent);
      const double dB_dA = exponent * std::pow(abs_A, exponent - 1.0);
      const double dSe_dB = -para_m * std::pow(B, -para_m - 1.0);
      const double mobile_saturation =
          1.0 - liquid_saturation_res_ - gas_saturation_res_;
      this->dSw_dpw_ = mobile_saturation * dSe_dB * dB_dA *
                       (-1.0 / para_p0);
    }
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

// pre-compute pf for wave
template <unsigned Tdim>
bool mpm::ThreePhaseParticleNew<Tdim>::build_wave_pf_context(){
  bool status = false;
  try {
      const double pgravity = -this->pgravity_[1];

      // 1) k(x) by dispersion
      auto solve_beta_k = [&] (double omega, double depth)->double{
        if (depth <= 0.0) return 0.0;
        double beta_k = std::max((omega * omega) / pgravity, 1e-8);
        if(beta_k * depth < 0.1) beta_k = omega /std::sqrt(std::max(pgravity * depth, 1e-12));
        for(int it = 0; it < 60; ++it){
          double kd = beta_k * depth, th = std::tanh(kd), ch = std::cosh(kd);
          double f = pgravity * beta_k * th - (omega * omega);
          if(std::abs(f) < 1e-12) break;
          double df = pgravity * (th + kd / (ch * ch));
          double kn = beta_k - f / df;
          if(!( kn > 0.0) || !std::isfinite(kn)) kn = std::max(0.5 * beta_k, 1e-12);
          beta_k = kn;
        }
      return std::max(beta_k, 1e-12);
      };

      const double omega = 2.0 * M_PI / this->wave_period_;
      this->seabed_surface_x_.resize(this->Nx_);
      this->water_depth_.resize(this->Nx_);
      this->beta_k_.assign(this->Nx_,0.0);
      for(int i = 0; i < this->Nx_; ++i){
          this->seabed_surface_x_[i] = this->domain_length_x_ * i / (this->Nx_ - 1);
          this->water_depth_[i] = this->depth_left_ + (this->depth_right_ - this->depth_left_) 
                            * (this->seabed_surface_x_[i] / this->domain_length_x_);
          this->beta_k_[i] = solve_beta_k(omega, this->water_depth_[i]);
      }

      // 2) H(x) by shoaling + breaker cap
      this->wave_height_.assign(this->Nx_,0.0);
      std::vector<double> cg(this->Nx_,0.0);
      for(int i = 0; i < this->Nx_; ++i) if(this->beta_k_[i] > 0){
          double c = omega / this->beta_k_[i], kd = this->beta_k_[i] * this->water_depth_[i];
          double s2 = std::sinh(2.0 * kd);
          double term = 1.0 + 2.0 * kd / (std::abs(s2) > 1e-12 ? s2 : (kd > 0 ? 1e-12 : -1e-12));
          cg[i]= 0.5 * c * term;
      }
      
      int iref = 0; 
      for(int i = 1; i < this->Nx_; ++i) if(this->water_depth_[i] > this->water_depth_[iref]) iref=i;
      double cg0 = std::max(cg[iref], 1e-12);
      for(int i = 0; i < this->Nx_; ++i) if(this->water_depth_[i] > 0 && this->beta_k_[i]>0){
          double Hx = this->wave_height_ini_ * std::sqrt(cg0/std::max(cg[i],1e-12));
          Hx = std::min(Hx, this->gamma_b_ * this->water_depth_[i]);
          this->wave_height_[i] = Hx;
      }

      // 3) P0(x) = rho g H / (2 cosh(k d))
      this->wave_P0_.assign(this->Nx_,0.0);
      for(int i = 0; i < this->Nx_; ++i) if(this->beta_k_[i]>0 && this->water_depth_[i]>0){
          double kd = this->beta_k_[i] * this->water_depth_[i];
          this->wave_P0_[i]=(this->liquid_density_ * pgravity * this->wave_height_[i])
                            /(2.0*std::max(std::cosh(kd),1e-12));
      }

      // 4) phi(x) = ∫ k dx
      this->wave_phi_.assign(this->Nx_, 0.0);
      for(int i = 1; i < this->Nx_; ++i){
          const double dx = seabed_surface_x_[i] - seabed_surface_x_[i-1];
          this->wave_phi_[i] = this->wave_phi_[i-1] + 0.5 * (this->beta_k_[i] + this->beta_k_[i-1]) * dx;
      }

      status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                    exception.what());
  }
  return status;
}

//==============================================================================
// PART 1: MAP PARTICLE INFORMATION TO NODES
//==============================================================================

//! Map particle mass and momentum to nodes
template <unsigned Tdim>
void mpm::ThreePhaseParticleNew<Tdim>::map_mass_momentum_to_nodes() noexcept {
  if (this->material_id_ != 999) {
    // SOLID PHASE
    mpm::Particle<Tdim>::map_mass_momentum_to_nodes();

    for (unsigned i = 0; i < nodes_.size(); ++i) {
      // Reduce the cost of accessing the shapefn_ array repeatedly
      double shapefn_i = shapefn_[i];

      // LIQUID PHASE
      double liquid_mass = liquid_mass_ * shapefn_i;
      nodes_[i]->update_mass_momentum(true, mpm::ParticlePhase::Liquid, liquid_mass,
                                                    liquid_mass * liquid_velocity_);

      // GAS PHASE
      double gas_mass = gas_mass_ * shapefn_i;
      nodes_[i]->update_mass_momentum(true, mpm::ParticlePhase::Gas, gas_mass,
                                                    gas_mass * gas_velocity_); 
      }
  }
}

// Map particle external force = body force + traction force
template <unsigned Tdim>
void mpm::ThreePhaseParticleNew<Tdim>::map_external_force(const VectorDim& pgravity) {
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

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   std::cout << "    mixture_force: " << mixture_force[0] << std::endl;
  //   std::cout << "    liquid_force: " << liquid_force[0] << std::endl;
  //   std::cout << "    gas_force: " << gas_force[0] << std::endl;
  // }

    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

// Map particle internal force 2D
template <>
void mpm::ThreePhaseParticleNew<2>::map_internal_force() {
  if (this->material_id_ != 999) {
    try {
      Eigen::Matrix<double, 2, 1> mixture_force, liquid_force, gas_force;

      this->total_stress_ = this->stress_;
      total_stress_[0] -= this->pore_pressure_;
      total_stress_[1] -= this->pore_pressure_;      

      // LIQUID PHASE
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        liquid_force[0] = dn_dx_(i, 0) * this->PIC_liquid_pressure_;
        liquid_force[1] = dn_dx_(i, 1) * this->PIC_liquid_pressure_;

        liquid_force *= this->volume_ * this->liquid_fraction_;

        nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, liquid_force);
      }

      // GAS PHASE
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        gas_force[0] = dn_dx_(i, 0) * (this->PIC_gas_pressure_);
        gas_force[1] = dn_dx_(i, 1) * (this->PIC_gas_pressure_);

        gas_force *= this->volume_ * this->gas_fraction_;
        nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Gas, gas_force);
      }

      // MIXTURE
      // [0 1] * [0 3] = 1*2 Matrix [mixture_force[0],mixture_force[1]]
      //         [3 1]
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        mixture_force[0] = dn_dx_(i, 0) * total_stress_[0] +
                           dn_dx_(i, 1) * total_stress_[3];
        mixture_force[1] = dn_dx_(i, 1) * total_stress_[1] +
                           dn_dx_(i, 0) * total_stress_[3];

        mixture_force *= -1. * volume_;
        nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Mixture, mixture_force);
      }

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   // std::cout << "    pore_pressure_: " << this->pore_pressure_ << std::endl;
  //   // std::cout << "    ini_pore_pressure_: " << this->ini_pore_pressure_ << std::endl;
  //   std::cout << "    liquid_pressure_: " << liquid_pressure_ << std::endl;
  //       //std::cout << "    liquid_force/v: " << (liquid_pressure_ - ini_liquid_pressure_)* liquid_fraction_ << std::endl;
  //   //std::cout << "    dn_dx_(0,0): " << dn_dx_(0,0) << std::endl;
  //   //std::cout << "    dn_dx_(0,1): " << dn_dx_(0,1) << std::endl;
  //   //std::cout << "    dn_dx_(1,0): " << dn_dx_(1,0) << std::endl;
  //   //std::cout << "    dn_dx_(1,1): " << dn_dx_(1,1) << std::endl;
  // }

    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    } 
  } 
}

// Map particle internal force 3D
template <>
void mpm::ThreePhaseParticleNew<3>::map_internal_force() {
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
        //           [0 3 5]   At nodali = 1
        // [0 1 2] * [3 1 4] = 1*3 Matrix [mixture_force[0],mixture_force[1],mixture_force[2]]
        //           [5 4 2]
        // For nodali = 4
        // [0 1 2]   [0 3 5]              [mixture_force[0],mixture_force[1],mixture_force[2]]
        // [1 1 2] * [3 1 4] = 4*3 Matrix [mixture_force[0],mixture_force[1],mixture_force[2]]
        // [2 1 2]   [5 4 2]              [mixture_force[0],mixture_force[1],mixture_force[2]]
        // [3 1 2]                        [mixture_force[0],mixture_force[1],mixture_force[2]]
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
inline Eigen::Matrix<double, Tdim, 1> mpm::ThreePhaseParticleNew<Tdim>::
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
void mpm::ThreePhaseParticleNew<Tdim>::map_drag_force_coefficient() {
  if (this->material_id_ != 999) {  
    try {
      for (unsigned i = 0; i < nodes_.size(); ++i) {

        // // T-viscosity
        // double muG = (0.0075 / (0.999 * this->PIC_temperature_ + 392.877)) *
        //                std::pow(((273.15 + this->PIC_temperature_) / 291.15) , 1.5);
        // this->gas_viscosity_ = muG;
        
        // double muL = 1.984e-6 * std::exp(1825.85 / (273.15 + this->PIC_temperature_));
        // this->liquid_viscosity_ = muL;

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
void mpm::ThreePhaseParticleNew<Tdim>::map_heat_to_nodes() {
  if (this->material_id_ != 999) {
    try {
      // Calculate mixture heat capacity  
      double mixture_heat_capacity_ = this->solid_heat_capacity_ +
                          this->liquid_heat_capacity_ + this->gas_heat_capacity_;

      // Map mixture heat capacity & heat
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        double mixture_heat_capacity = mixture_heat_capacity_ * shapefn_[i];
        nodes_[i]->update_heat_capacity(true, mpm::ParticlePhase::Mixture,
                                mixture_heat_capacity);
        nodes_[i]->update_heat(true, mpm::ParticlePhase::Mixture, 
                                mixture_heat_capacity * this->temperature_);
      }
    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

// Map conductive heat transfer
template <unsigned Tdim>
void mpm::ThreePhaseParticleNew<Tdim>::map_heat_conduction() {
  if (this->material_id_ != 999) {
    try {
      // Calculate temperature gradient
      mpm::Particle<Tdim>::compute_temperature_gradient(mpm::ParticlePhase::Solid);

      // Calculate thermal conductivity of mixture
      double mixture_cond = this->solid_fraction_ * this->solid_thermal_conductivity_ +
                            this->liquid_fraction_ * this->liquid_thermal_conductivity_ +
                            this->gas_fraction_ * this->gas_thermal_conductivity_;

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

  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   //std::cout << "    Mheat_conduction: " << -1 * this->volume_ * mixture_cond << std::endl;
  //   //std::cout << "    volume_: " <<this->volume_;
  //   //std::cout << "    liquid_fraction_: " << this->liquid_fraction_ << std::endl;
  //   std::cout << "    temperature_gradient_y: " << this->temperature_gradient_[1] << std::endl;
  // }

    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

// Map convective heat transfer
template <unsigned Tdim>
void mpm::ThreePhaseParticleNew<Tdim>::map_heat_convection() {
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
        liquid_heat_convection *= -1 * liquid_heat_capacity_* this->volume_;
        gas_heat_convection *= -1 * gas_heat_capacity_ * this->volume_;
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
void mpm::ThreePhaseParticleNew<Tdim>::compute_updated_velocity(
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
  // // Debug
  // if (this->id_ <= 200 && this->id_ >= 189) {
  //   std::cout << "particle id: " << this->id_ << std::endl;
  //   //std::cout << "    K_ww: " << 1 / (porosity_ * liquid_compressibility_) << std::endl;
  //   //std::cout << "    thermal_term: " << thermal_term << std::endl;
  //   //std::cout << "    strain_term: " << strain_term << std::endl;
  //   std::cout << "    gas_velocity_: " << this->gas_velocity_ << std::endl;
  //   std::cout << "    liquid_velocity_: " << this->liquid_velocity_ * dt<< std::endl;
  // }

    } catch (std::exception& exception) {
        console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                        exception.what());
    }
  }
}

// Return frozen coefficients for a global semi-implicit pressure solve
template <unsigned Tdim>
mpm::ThreePhasePressureState<Tdim>
mpm::ThreePhaseParticleNew<Tdim>::semi_implicit_pressure_state() {
  ThreePhasePressureState<Tdim> state;
  if (this->material_id_ == 999) return state;

  state.active_node_ids.reserve(nodes_.size());
  for (const auto& node : nodes_) state.active_node_ids.emplace_back(node->active_id());
  state.shapefn = shapefn_;
  state.dn_dx = dn_dx_;
  state.volume = volume_;
  state.porosity = porosity_;
  state.liquid_saturation = liquid_saturation_;
  state.gas_saturation = gas_saturation_;
  state.liquid_pressure = PIC_liquid_pressure_;
  state.gas_pressure = PIC_gas_pressure_;
  state.saturation_pressure_tangent = std::max(dSw_dpw_, 0.0);
  state.liquid_compressibility = liquid_compressibility_;
  state.gas_compressibility =
      gas_saturation_ / std::max(gas_density_, 1.0e-12) * gas_molar_mass_ /
      gas_constant_ / std::max(PIC_temperature_ + 273.15, 1.0);
  // The three-phase momentum equation uses viscosity / permeability in its
  // drag coefficient.  Use the corresponding Darcy mobility here so the
  // semi-implicit and explicit paths have the same steady response.
  state.liquid_conductivity =
      liquid_permeability_ / std::max(liquid_viscosity_, 1.0e-30);
  state.gas_conductivity =
      gas_permeability_ / std::max(gas_viscosity_, 1.0e-30);
  state.liquid_density = liquid_density_;
  state.gas_density = gas_density_;
  state.fixed_gas_pressure =
      liquid_material_->template property_or<bool>(
          std::string("fixed_gas_pressure"), false);

  const auto solid_strain_rate = this->compute_strain_rate(
      dn_dx_, shapefn_, mpm::ParticlePhase::Solid);
  const double solid_volumetric_strain_rate =
      is_axisymmetric_ ? solid_strain_rate.head(3).sum()
                       : solid_strain_rate.head(Tdim).sum();
  const double beta_w =
      3.0 * (1.0 - porosity_) * liquid_saturation_ * solid_expansivity_ +
      3.0 * liquid_fraction_ * liquid_expansivity_;
  const double beta_g =
      3.0 * (1.0 - porosity_) * gas_saturation_ * solid_expansivity_ +
      gas_fraction_ / std::max(gas_density_, 1.0e-12) * gas_pressure_ *
          gas_molar_mass_ / gas_constant_ /
          std::pow(std::max(PIC_temperature_ + 273.15, 1.0), 2);
  state.liquid_mass_source =
      beta_w * temperature_acceleration_ -
      liquid_saturation_ * solid_volumetric_strain_rate;
  state.gas_mass_source =
      beta_g * temperature_acceleration_ -
      gas_saturation_ * solid_volumetric_strain_rate;

  state.pressure_boundary = set_pressure_constraint_ || this->free_surface();
  if (set_pressure_constraint_) {
    state.boundary_liquid_pressure = 0.0;
    state.boundary_gas_pressure = 0.0;
  } else if (this->free_surface()) {
    if ((coordinates_[1] <= sea_level_) && (coordinates_[0] >= 25.0) &&
        wave_pressure_) {
      const double pgravity = -pgravity_[1];
      state.boundary_liquid_pressure =
          liquid_density_ * pgravity * (sea_level_ - coordinates_[1]) +
          this->pf_seabed_surface_particle();
      state.boundary_gas_pressure = state.boundary_liquid_pressure;
    } else {
      state.boundary_liquid_pressure =
          liquid_material_->template property_or<double>(
              std::string("surface_liquid_pressure"), 25000.0);
      state.boundary_gas_pressure =
          state.fixed_gas_pressure
              ? gas_pressure_
              : liquid_material_->template property_or<double>(
                    std::string("surface_gas_pressure_ratio"), 1.0) *
                    state.boundary_liquid_pressure;
    }
  }
  return state;
}

// Increment particle pressure from a global semi-implicit nodal solution.
// Incremental transfer avoids repeated particle-node-particle smoothing of
// the absolute pressure when the physical pressure increment is zero.
template <unsigned Tdim>
int mpm::ThreePhaseParticleNew<Tdim>::update_semi_implicit_pressure(
    const Eigen::VectorXd& nodal_liquid_pressure_increment,
    const Eigen::VectorXd& nodal_gas_pressure_increment,
    const Eigen::VectorXd& nodal_liquid_pressure,
    const Eigen::VectorXd& nodal_gas_pressure,
    bool reconstruct_pressure_gradient, bool bounded_transfer, double dt) {
  if (this->material_id_ == 999) return 0;
  try {
    const double old_liquid_pressure = liquid_pressure_;
    const double old_gas_pressure = gas_pressure_;
    double minimum_liquid_pressure = std::numeric_limits<double>::max();
    double maximum_liquid_pressure = std::numeric_limits<double>::lowest();
    double minimum_gas_pressure = std::numeric_limits<double>::max();
    double maximum_gas_pressure = std::numeric_limits<double>::lowest();
    double minimum_capillary_pressure = std::numeric_limits<double>::max();
    double maximum_capillary_pressure = std::numeric_limits<double>::lowest();
    const bool reconstruct_gradient =
        reconstruct_pressure_gradient || bounded_transfer;
    if (reconstruct_gradient) {
      liquid_pressure_gradient_.setZero();
      gas_pressure_gradient_.setZero();
    }
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      const auto active_id = nodes_[i]->active_id();
      if (active_id >=
              static_cast<Index>(nodal_liquid_pressure_increment.size()) ||
          active_id >= static_cast<Index>(nodal_gas_pressure_increment.size()) ||
          active_id >= static_cast<Index>(nodal_liquid_pressure.size()) ||
          active_id >= static_cast<Index>(nodal_gas_pressure.size()))
        throw std::runtime_error("Invalid active node id in pressure update");
      liquid_pressure_ +=
          shapefn_[i] * nodal_liquid_pressure_increment(active_id);
      gas_pressure_ +=
          shapefn_[i] * nodal_gas_pressure_increment(active_id);
      if (reconstruct_gradient) {
        liquid_pressure_gradient_ +=
            dn_dx_.row(i).transpose() * nodal_liquid_pressure(active_id);
        gas_pressure_gradient_ +=
            dn_dx_.row(i).transpose() * nodal_gas_pressure(active_id);
      } else {
        liquid_pressure_gradient_ += dn_dx_.row(i).transpose() *
                                     nodal_liquid_pressure_increment(active_id);
        gas_pressure_gradient_ += dn_dx_.row(i).transpose() *
                                  nodal_gas_pressure_increment(active_id);
      }
      if (bounded_transfer) {
        minimum_liquid_pressure = std::min(
            minimum_liquid_pressure, nodal_liquid_pressure(active_id));
        maximum_liquid_pressure = std::max(
            maximum_liquid_pressure, nodal_liquid_pressure(active_id));
        minimum_gas_pressure = std::min(
            minimum_gas_pressure, nodal_gas_pressure(active_id));
        maximum_gas_pressure = std::max(
            maximum_gas_pressure, nodal_gas_pressure(active_id));
        const double capillary_pressure =
            nodal_gas_pressure(active_id) - nodal_liquid_pressure(active_id);
        minimum_capillary_pressure =
            std::min(minimum_capillary_pressure, capillary_pressure);
        maximum_capillary_pressure =
            std::max(maximum_capillary_pressure, capillary_pressure);
      }
    }

    const bool fixed_gas_pressure =
        liquid_material_->template property_or<bool>(
            std::string("fixed_gas_pressure"), false);
    if (fixed_gas_pressure) gas_pressure_ = old_gas_pressure;
    unsigned limiter_flags = 0;
    const bool pressure_boundary = set_pressure_constraint_ || this->free_surface();
    if (bounded_transfer && !pressure_boundary) {
      // Bound the FLIP pressure transfer by the local resolved nodal field.
      // Pure incremental transfer preserves particle modes that are invisible
      // to the pressure grid.  Gravity amplified one such two-particle mode
      // into alternating wet and dry layers.  Clipping only new extrema keeps
      // resolved FLIP increments unchanged and is distinct from mechanical
      // PIC damping.
      const double candidate_liquid_pressure = liquid_pressure_;
      const double candidate_gas_pressure = gas_pressure_;
      if (!fixed_gas_pressure) {
        gas_pressure_ = std::max(
            minimum_gas_pressure,
            std::min(maximum_gas_pressure, gas_pressure_));
      }
      // Preserve the pore-air pressure where possible and project liquid
      // pressure onto the intersection of its phase and capillary bounds.
      const double feasible_liquid_lower = std::max(
          minimum_liquid_pressure,
          gas_pressure_ - maximum_capillary_pressure);
      const double feasible_liquid_upper = std::min(
          maximum_liquid_pressure,
          gas_pressure_ - minimum_capillary_pressure);
      if (feasible_liquid_lower <= feasible_liquid_upper) {
        liquid_pressure_ = std::max(
            feasible_liquid_lower,
            std::min(feasible_liquid_upper, liquid_pressure_));
      } else {
        // The independent gas clamp can very rarely miss the coupled feasible
        // set.  A local nodal pressure pair is feasible by construction.
        double closest_distance = std::numeric_limits<double>::max();
        for (unsigned i = 0; i < nodes_.size(); ++i) {
          const auto active_id = nodes_[i]->active_id();
          const double nodal_gas = fixed_gas_pressure
              ? old_gas_pressure
              : nodal_gas_pressure(active_id);
          const double distance =
              std::pow(nodal_liquid_pressure(active_id) -
                           candidate_liquid_pressure, 2) +
              std::pow(nodal_gas - candidate_gas_pressure, 2);
          if (distance < closest_distance) {
            closest_distance = distance;
            liquid_pressure_ = nodal_liquid_pressure(active_id);
            gas_pressure_ = nodal_gas;
          }
        }
      }
      if (liquid_pressure_ != candidate_liquid_pressure ||
          gas_pressure_ != candidate_gas_pressure)
        limiter_flags |= 1U;
      if (gas_pressure_ - liquid_pressure_ !=
          candidate_gas_pressure - candidate_liquid_pressure)
        limiter_flags |= 2U;
    }
    if (set_pressure_constraint_) {
      liquid_pressure_ = 0.0;
      gas_pressure_ = 0.0;
    } else if (this->free_surface()) {
      if ((coordinates_[1] <= sea_level_) && (coordinates_[0] >= 25.0) &&
          wave_pressure_) {
        const double pgravity = -pgravity_[1];
        liquid_pressure_ =
            liquid_density_ * pgravity * (sea_level_ - coordinates_[1]) +
            this->pf_seabed_surface_particle();
        gas_pressure_ = liquid_pressure_;
      } else {
        liquid_pressure_ = liquid_material_->template property_or<double>(
            std::string("surface_liquid_pressure"), 25000.0);
        if (!fixed_gas_pressure)
          gas_pressure_ =
              liquid_material_->template property_or<double>(
                  std::string("surface_gas_pressure_ratio"), 1.0) *
              liquid_pressure_;
      }
    }

    liquid_pressure_acceleration_ =
        (liquid_pressure_ - old_liquid_pressure) / dt;
    gas_pressure_acceleration_ =
        fixed_gas_pressure ? 0.0 : (gas_pressure_ - old_gas_pressure) / dt;
    PIC_liquid_pressure_ = liquid_pressure_;
    PIC_gas_pressure_ = gas_pressure_;
    suction_pressure_ = std::max(gas_pressure_ - liquid_pressure_, 0.0);
    pore_pressure_ = liquid_saturation_ * liquid_pressure_ +
                     gas_saturation_ * gas_pressure_;
    if (!std::isfinite(liquid_pressure_) || !std::isfinite(gas_pressure_))
      return -1;
    return static_cast<int>(limiter_flags);
  } catch (std::exception& exception) {
    console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__,
                    __func__, exception.what());
    return -1;
  }
}

// Compute updated pore pressure
template <unsigned Tdim>
void mpm::ThreePhaseParticleNew<Tdim>::compute_pore_pressure(double dt){
    if (this->material_id_ != 999) {
    try {

// Update particle porosity
  bool two_phase = liquid_material_->template 
                    property<bool>(std::string("two_phase"));
  if (two_phase) {

    // Compute at centroid
    // get liquid phase strain rate at cell centre
    auto liquid_strain_rate_centroid =
        this->compute_strain_rate(dn_dx_, shapefn_, mpm::ParticlePhase::Liquid);

    auto strain_rate_centroid =
        this->compute_strain_rate(dn_dx_centroid_, shapefn_centroid_, mpm::ParticlePhase::Solid);

    // Compute mass-convection-induced pressure
    double convective_term = 0;

    // Compute thermal-expansion-induced pressure
    double thermal_term = (3 * porosity_ * liquid_expansivity_ + 
                            3 * (1 - porosity_) * solid_expansivity_)
                          * this->temperature_increment_cent_;

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
        1 / (porosity_ * liquid_compressibility_) * (thermal_term - strain_term);
      this->liquid_pressure_acceleration_ = pore_pressure_increments_ / dt;
      this->liquid_pressure_ += this->liquid_pressure_acceleration_ * dt;

      this->PIC_gas_pressure_ = this->gas_pressure_;
      this->PIC_liquid_pressure_ = this->liquid_pressure_;
      this->pore_pressure_ = this->liquid_saturation_ * this->PIC_liquid_pressure_ +
                            this->gas_saturation_ * this->PIC_gas_pressure_; 
    }    

  } else {

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

      const double beta_w = 3.0 * (1 - porosity_) * liquid_saturation_ * solid_expansivity_ +
                            3.0 * liquid_fraction_ * liquid_expansivity_;

      const double beta_g = 3.0 * (1 - porosity_) * gas_saturation_ * solid_expansivity_ +
                            gas_fraction_ / gas_density_ * gas_pressure_ * gas_molar_mass_ /
                            gas_constant_ / std::pow(PIC_temperature_ + 273.15, 2);

      const double beta_m = 3 * solid_fraction_ * solid_expansivity_ +
                      3 * liquid_fraction_ * liquid_expansivity_ +
                      gas_fraction_ / gas_density_ * gas_pressure_ * gas_molar_mass_ /
                      gas_constant_ / std::pow(PIC_temperature_ + 273.15, 2);

      const double f_w = beta_w * this->temperature_acceleration_-liquid_saturation_ * solid_strain_rate - 
                    liquid_fraction_ * (liquid_strain_rate - solid_strain_rate); 

      const double f_g = beta_g * this->temperature_acceleration_-(1 - liquid_saturation_) * solid_strain_rate -
                    gas_fraction_ * (gas_strain_rate - solid_strain_rate);

      const double f = beta_m * this->temperature_acceleration_ - solid_strain_rate -
                    liquid_fraction_ * (liquid_strain_rate - solid_strain_rate) -
                    gas_fraction_ * (gas_strain_rate - solid_strain_rate);

      Eigen::Matrix<double, 2, 2> K_matrix, K_matrix_inverse;
      K_matrix(0,0) = K_ww;
      K_matrix(0,1) = K_wg;
      K_matrix(1,0) = K_gw;
      K_matrix(1,1) = K_gg;

      const bool fixed_gas_pressure =
          liquid_material_->template property_or<bool>(
              std::string("fixed_gas_pressure"), false);

      // Fixed-gas formulation: solve only the liquid mass balance.
      if (fixed_gas_pressure) {
          this->liquid_pressure_acceleration_ = f_w / K_ww;
          this->gas_pressure_acceleration_ = 0;
      } else if (K_gg > 1E-16) {
        K_matrix_inverse = K_matrix.inverse();

        this->liquid_pressure_acceleration_ = K_matrix_inverse(0, 0) * f_w +
                                              K_matrix_inverse(0, 1) * f_g;
        this->gas_pressure_acceleration_ = K_matrix_inverse(1, 0) * f_w +
                                          K_matrix_inverse(1, 1) * f_g;
      } else {
          this->liquid_pressure_acceleration_ = f_w / K_ww;
          this->gas_pressure_acceleration_ = 0;
      }

      this->liquid_pressure_ += this->liquid_pressure_acceleration_ * dt;
      if (!fixed_gas_pressure) {
        this->gas_pressure_ += this->gas_pressure_acceleration_ * dt;
      } else {
        this->gas_pressure_acceleration_ = 0.0;
      }

    if (this->set_pressure_constraint_){
      // this->gas_pressure_ = 0.0;
      // this->suction_pressure_ = this->liquid_density_ * 9.81 * this->coordinates_[1];
      // this->liquid_pressure_ = this->gas_pressure_ - this->suction_pressure_;
      this->gas_pressure_ = 0.0;
      this->suction_pressure_ = 0.0;
      this->liquid_pressure_ = 0.0;
      // std::cout << "    id_: " << id_ << std::endl;
    }
    if (this->free_surface()) {
      if ((this->coordinates_[1]<=this->sea_level_) && (this->coordinates_[0]>=25.0) && this->wave_pressure_){
        const double pgravity = -this->pgravity_[1];

        double wave_pressure_surface = this->liquid_density_ * pgravity 
                         *(this->sea_level_ - this->coordinates_[1]);
        wave_pressure_surface += this->pf_seabed_surface_particle();

        this->liquid_pressure_ = wave_pressure_surface;
        this->gas_pressure_ = this->liquid_pressure_;
        this->suction_pressure_ = this->gas_pressure_ - this->liquid_pressure_;
        if (this->suction_pressure_ < 0) this->suction_pressure_ = 0.0;
      } else {
        // Slope surface: hydraulic head and gas venting are controlled by JSON.
        const double surface_liquid_pressure =
            liquid_material_->template property_or<double>(
                std::string("surface_liquid_pressure"), 25000.0);
        const double surface_gas_pressure_ratio =
            liquid_material_->template property_or<double>(
                std::string("surface_gas_pressure_ratio"), 1.0);
        this->liquid_pressure_ = surface_liquid_pressure;
        if (!fixed_gas_pressure) {
          this->gas_pressure_ =
              surface_gas_pressure_ratio * surface_liquid_pressure;
        }
        this->suction_pressure_ =
            std::max(this->gas_pressure_ - this->liquid_pressure_, 0.0);
      }
    }
      this->PIC_gas_pressure_ = this->gas_pressure_;
      this->PIC_liquid_pressure_ = this->liquid_pressure_;

      this->pore_pressure_ = this->liquid_saturation_ * this->PIC_liquid_pressure_ +
                            this->gas_saturation_ * this->PIC_gas_pressure_;
      
    }
    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

// Compute updated update porosity of the particle
template <unsigned Tdim>
bool mpm::ThreePhaseParticleNew<Tdim>::update_particle_porosity(double dt) {
  bool status = true;
  if (this->material_id_ != 999) {
    try {
      // Update particle porosity
      bool two_phase = liquid_material_->template 
                    property<bool>(std::string("two_phase"));
      if (!two_phase) {
      this->porosity_ = 
          1 - (1 - this->porosity_) / (1 + this->dvolumetric_strain_);
      }

      // Clamp porosity to valid range [0, 1]
      constexpr double porosity_min = 1.0e-6;
      constexpr double porosity_max = 1.0 - 1.0e-6;

      if (this->porosity_ <= 0) {
        console_->warn("Particle #{}: porosity {} clamped to {}", 
                       this->id_, this->porosity_, porosity_min);
        this->porosity_ = porosity_min;
      } else if (this->porosity_ >= 1) {
        console_->warn("Particle #{}: porosity {} clamped to {}", 
                       this->id_, this->porosity_, porosity_max);
        this->porosity_ = porosity_max;
      }

      this->solid_fraction_ = 1.0 - this->porosity_;

    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

// Update particle volume
template <unsigned Tdim>
void mpm::ThreePhaseParticleNew<Tdim>::update_particle_volume() {
  if (this->material_id_ != 999) {
    try {
      // SOLID PHASE
      // Strain rate for reduced integration
      bool two_phase = liquid_material_->template 
                    property<bool>(std::string("two_phase"));
      if (!two_phase) {
      this->volume_ *= (1. + dvolumetric_strain_);
      this->mass_density_ = this->mass_density_ / (1. + dvolumetric_strain_);
      }

      // mpm::Particle<Tdim>::update_particle_volume();

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
void mpm::ThreePhaseParticleNew<Tdim>::update_permeability() {
  if (this->material_id_ != 999) {
    try {
      // SWRC is Linearized

      double k_phi = 1;
      double k_a = 1;
      double k_r_liquid = 1.;
      double k_r_gas = 1.; 

      const double para_m = liquid_material_->template 
                    property<double>(std::string("para_m"));
      const double relative_permeability_floor =
          liquid_material_->template property_or<double>(
              std::string("relative_permeability_floor"), 1.0e-12);
      if (relative_permeability_floor <= 0.0 ||
          relative_permeability_floor >= 1.0)
        throw std::runtime_error(
            "relative_permeability_floor must be between zero and one");

      // Calculate porosity-dependent permeability, k_phi
      k_phi = std::pow(porosity_ / ini_porosity_, 1.5) * 
              std::pow((1 - ini_porosity_) / (1 - porosity_), 3);
      // this->effective_saturation_ = (this->liquid_saturation_ - this->liquid_saturation_res_) /
      //                               (1 - this->liquid_saturation_res_);
      k_r_liquid = std::pow(this->effective_saturation_, (2*para_m+3));
      k_r_gas = std::pow(1-this->effective_saturation_, 2) * 
                (1 - std::pow(this->effective_saturation_, (2*para_m+1)));
      k_r_liquid = std::max(k_r_liquid, relative_permeability_floor);
      k_r_gas = std::max(k_r_gas, relative_permeability_floor);

      // Calculate absolute permeability, k_a
      k_a = this->intrinsic_permeability_ * k_phi;

      // Calculate permeability
      this->liquid_permeability_ = k_a *  k_r_liquid; 
      this->gas_permeability_ = k_a * k_r_gas;


      // SWRC of VGM
      // // Porosity parameter, using updated porosity_
      // double porosity_0 = material_->template property<double>("porosity");
      // // Porosity parameter
      // double k_phi = std::pow(this->porosity_, 3) / std::pow((1. - this->porosity_), 2) /
      //           (std::pow(porosity_0, 3) / std::pow((1. - porosity_0), 2));
      // double k_a = 1;
      // double k_r_liquid = 1.;
      // double k_r_gas = 1.; 

      // const double constant_mw = liquid_material_->template 
      //               property<double>(std::string("constant_mw"));
      // const double constant_mg = liquid_material_->template 
      //               property<double>(std::string("constant_mg"));

      // // Relative permeability
      // k_r_liquid = compute_relative_permeability(
      //     this->liquid_saturation_, this->liquid_saturation_res_,
      //     this->liquid_saturation_res_,this->gas_saturation_res_, constant_mw);
      // k_r_gas = compute_relative_permeability(
      //     this->gas_saturation_, this->gas_saturation_res_, this->liquid_saturation_res_,
      //     this->gas_saturation_res_, constant_mg);

      // // Calculate absolute permeability, k_a
      // k_a = intrinsic_permeability_ * k_phi;

      // // Calculate permeability
      // this->liquid_permeability_ = k_a *  k_r_liquid; 
      // this->gas_permeability_ = k_a * k_r_gas;
    } catch (std::exception& exception) {
      console_->error("{} #{}: Function: {}, {}\n", __FILE__, __LINE__, __func__,
                      exception.what());
    }
  }
}

// Calculate relative permeability
template <unsigned Tdim>
double mpm::ThreePhaseParticleNew<Tdim>::compute_relative_permeability(
    double Spi, double Srpi, double Srw, double Srg, double m_pi) {
  double denominator = 1. - Srw - Srg;
  double numerator = Spi - Srpi;
  if (denominator <= 0.) return 0.;
  double value = numerator / denominator;
  value = std::max(value, 0.);
  return std::pow(value, m_pi);
}

// Compute updated liquid and gas saturation
template <unsigned Tdim>
void mpm::ThreePhaseParticleNew<Tdim>::update_liquid_gas_saturation(double dt){
  if (this->material_id_ != 999) {
    try {
      this->suction_pressure_ = PIC_gas_pressure_ - PIC_liquid_pressure_;
      if (this->suction_pressure_ < 0) this->suction_pressure_ = 0.0;
      // this->liquid_saturation_ = 1 - 1E-6 * this->suction_pressure_;
      // this->dSw_dpw_ = 1E-6;

      // // slope-T
      // const double ini_temperature = 
      //     material_->template property<double>(std::string("initial_temperature"));
      // const double a0 = 15000.0;
      // const double beta = -72.93;
      // const double SLr0 = 0.37;
      // const double lam_r = -5.93e-4;
      // const double nv0 = 3.48;
      // const double kn = -0.052;
      
      // double a = a0 * (beta + this->PIC_temperature_) / (beta + ini_temperature);
      // double SLr = SLr0 + lam_r * (this->PIC_temperature_ - ini_temperature);
      // double nv = nv0 + kn * (this->PIC_temperature_ - ini_temperature);

      // double ratio = this->suction_pressure_ / a;
      // double ratio_nv = std::pow(ratio, nv);
      // double g = 1.0 + ratio_nv;
      // double expo = 1.0 / nv - 1.0;

      // this->liquid_saturation_ = SLr + (1.0 - SLr) * std::pow(g, expo);
      // if (this->liquid_saturation_<0.0) this->liquid_saturation_ = 0.0;
      // if (this->liquid_saturation_>1.0) this->liquid_saturation_ = 1.0;

      // double dSw_dPc = (1.0 - SLr) * (1.0 - nv) / a *
      //                 std::pow(ratio, nv - 1.0) * std::pow(g, 1.0 / nv - 2.0);
      // this->dSw_dpw_ = -dSw_dPc;

      // this->suction_pressure_ = PIC_gas_pressure_ - PIC_liquid_pressure_;
      double para_p0 = liquid_material_->template 
                    property<double>(std::string("para_p0"));
      double para_m = liquid_material_->template 
                    property<double>(std::string("para_m"));
      bool two_phase = liquid_material_->template 
                    property<bool>(std::string("two_phase"));
      double A = suction_pressure_ / para_p0;
      double B = 1.0 + std::pow(std::abs(A), 1.0 / (1.0 - para_m));
      double Se = std::pow(B, -para_m);
      Se = std::min(1.0, std::max(0.0, Se));
      
      // if (two_phase) {
      //   // For two phase cases
      //   this->liquid_saturation_ = 1.0;
      // }
      // else {
      // this->liquid_saturation_ = (Se * (1.0 - liquid_saturation_res_ -
      //               gas_saturation_res_) + liquid_saturation_res_);
      // }

      this->liquid_saturation_ = (Se * (1.0 - liquid_saturation_res_ -
              gas_saturation_res_) + liquid_saturation_res_);
      const double mobile_saturation =
          1.0 - liquid_saturation_res_ - gas_saturation_res_;
      if (mobile_saturation > 0.0) {
        this->effective_saturation_ = std::min(
            1.0, std::max(0.0,
                (this->liquid_saturation_ - liquid_saturation_res_) /
                mobile_saturation));
      } else {
        this->effective_saturation_ = 0.0;
      }

      // Some terms For dB_dA
      double sgnA = (A >= 0.0 ? 1.0 : -1.0);
      const double eps = 1e-12;
      double absA_eps = std::max(std::abs(A), eps);
      
      // For dSw_dpw
      double dA_dpw = -1.0 / para_p0;
      double dB_dA = 1.0 / (1.0 - para_m) * std::pow(absA_eps, (1.0 / (1.0 - para_m) - 1.0)) * sgnA;
      double dSe_dB = -para_m * std::pow(B, (-para_m - 1.0));
      double dSw_dSe = (1 - this->liquid_saturation_res_ - this->gas_saturation_res_); 
      double dSw_dpw = dSw_dSe * dSe_dB * dB_dA * dA_dpw;

      this->dSw_dpw_ = dSw_dpw;

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
void mpm::ThreePhaseParticleNew<Tdim>::update_particle_density(double dt){
  if (this->material_id_ != 999) {  
    try {
        // Update particle porosity
      bool two_phase = liquid_material_->template 
                    property<bool>(std::string("two_phase"));
      if (!two_phase) {
      double p_ref = material_->template property<double>(std::string("p_ref"));
      this->gas_density_ = gas_molar_mass_ * (PIC_gas_pressure_ + p_ref)/ 
                            gas_constant_ / (this->PIC_temperature_ + 273.15);

      // Check if NaN
      if (!(this->gas_density_ > 0.)) gas_density_ = 1.2;

      // this->density_ /= 1 + solid_expansivity_ * temperature_increment_;
      this->liquid_density_ /= 1 + liquid_expansivity_ * temperature_increment_ -
                                  dt * liquid_pressure_acceleration_ * liquid_compressibility_;
       }


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
bool mpm::ThreePhaseParticleNew<Tdim>::assign_particle_traction(unsigned direction,
                                                          double traction) {
  bool status = false;
  if (this->material_id_ != 999) {
    try {
      if (direction >= Tdim * 3 ||
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
        gas_traction_(direction - 2 * Tdim) =
            -traction * this->volume_ / this->size_(direction - 2 * Tdim);

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
bool mpm::ThreePhaseParticleNew<Tdim>::assign_particle_liquid_velocity_constraint(
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
void mpm::ThreePhaseParticleNew<Tdim>::apply_particle_liquid_velocity_constraints() {
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
bool mpm::ThreePhaseParticleNew<Tdim>::assign_particle_pore_pressure_constraint(
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
void mpm::ThreePhaseParticleNew<Tdim>::apply_particle_pore_pressure_constraints(
      double pore_pressure) {

  // Set particle temperature constraint
  // this->pore_pressure_ = pore_pressure;
  // this->liquid_pressure_ = this->liquid_density_ * 9.81
  //           *(0.0 - this->coordinates_[1]);
  // this->gas_pressure_ = pore_pressure;
  // this->PIC_liquid_pressure_ = pore_pressure;
  // this->PIC_gas_pressure_ = pore_pressure;
  this->set_pressure_constraint_ = true;
}

// Overwrite node velocity to get strain correct
template <unsigned Tdim>
void mpm::ThreePhaseParticleNew<Tdim>::map_moving_rigid_velocity_to_nodes(
    unsigned dir, double velocity, double dt) noexcept {
  if (this->material_id_ == 999){  
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      nodes_[i]->assign_velocity_from_rigid(dir, velocity, dt);
    }
  }
}

//! Initial pore pressure
template <unsigned Tdim>
bool mpm::ThreePhaseParticleNew<Tdim>::initialise_pore_pressure_watertable(
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
        this->liquid_pressure_ =
            ((h0_right - h0_left) / (right_boundary - left_boundary) *
                 (position - left_boundary) +
             h0_left - this->coordinates_(dir_v)) *
            1000 * 9.81;
      } else
        // Particle with only left boundary
        this->liquid_pressure_ =
            (h0_left - this->coordinates_(dir_v)) * 1000 * 9.81;
    }
    // Particle with only right boundary
    else if (right_boundary != std::numeric_limits<double>::max())
      this->liquid_pressure_ =
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
//  PART 4: SMOOTHING TECHNIQUE
//==============================================================================

// Map particle pore liquid pressure to nodes
template <unsigned Tdim>
bool mpm::ThreePhaseParticleNew<Tdim>::map_pore_pressure_to_nodes(
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
      nodes_[i]->update_mass_pressure(mpm::ParticlePhase::Gas,
                                      shapefn_[i] * gas_mass_ * gas_pressure_,
                                      current_time);                           
    }
    return status;
  }
}

// Map particle pore liquid pressure to nodes
template <unsigned Tdim>
void mpm::ThreePhaseParticleNew<Tdim>::map_mass_pressure_to_nodes() {
  if (this->material_id_ != 999) {  
    // Check if particle mass is set
    assert(liquid_mass_ != std::numeric_limits<double>::max());
    if (this->material_id_ != 999){
      // Map particle liquid mass and pore pressure to nodes
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        nodes_[i]->update_scalers(true, 0, shapefn_[i] * 1.0);
        nodes_[i]->update_pressure(true, mpm::ParticlePhase::Gas,
                                        shapefn_[i] * 1.0 * this->gas_pressure_);
        nodes_[i]->update_pressure(true, mpm::ParticlePhase::Liquid,
                                  shapefn_[i] * 1.0 * this->liquid_pressure_);
      }
    }
  }
}


// Compute pore liquid pressure smoothing based on nodal pressure
template <unsigned Tdim>
bool mpm::ThreePhaseParticleNew<Tdim>::compute_pore_pressure_smoothing() noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  bool status = true;  
  if (this->material_id_ != 999) {
    const bool fixed_gas_pressure =
        liquid_material_->template property_or<bool>(
            std::string("fixed_gas_pressure"), false);
    double liquid_pressure = 0;
    double gas_pressure = 0;    
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      liquid_pressure += shapefn_(i) * nodes_[i]->pressure(mpm::ParticlePhase::Liquid);
      gas_pressure += shapefn_(i) * nodes_[i]->pressure(mpm::ParticlePhase::Gas);
    }
    if (fixed_gas_pressure) gas_pressure = this->gas_pressure_;

    if (this->set_pressure_constraint_){
      // gas_pressure = 0.0;
      // liquid_pressure = gas_pressure - this->liquid_density_ * 9.81 * this->coordinates_[1];
      gas_pressure = 0.0;
      liquid_pressure = 0.0;
    }
    if (this->free_surface()) {
      if ((this->coordinates_[1]<=this->sea_level_) && (this->coordinates_[0]>=25.0) && this->wave_pressure_){
        const double pgravity = -this->pgravity_[1];

        double wave_pressure_surface = this->liquid_density_ * pgravity 
                         *(this->sea_level_ - this->coordinates_[1]);
        wave_pressure_surface += this->pf_seabed_surface_particle();

        liquid_pressure = wave_pressure_surface;
        gas_pressure = liquid_pressure;

      } else {
        // Slope surface: hydraulic head and gas venting are controlled by JSON.
        const double surface_liquid_pressure =
            liquid_material_->template property_or<double>(
                std::string("surface_liquid_pressure"), 25000.0);
        const double surface_gas_pressure_ratio =
            liquid_material_->template property_or<double>(
                std::string("surface_gas_pressure_ratio"), 1.0);
        liquid_pressure = surface_liquid_pressure;
        if (!fixed_gas_pressure) {
          gas_pressure =
              surface_gas_pressure_ratio * surface_liquid_pressure;
        }
      }
    }
    
    // Update pore liquid pressure to interpolated nodal pressure
    this->PIC_liquid_pressure_= liquid_pressure; 
    this->PIC_gas_pressure_ = gas_pressure;
    this->liquid_pressure_ = this->PIC_liquid_pressure_;
    this->gas_pressure_ = this->PIC_gas_pressure_;
    this->pore_pressure_ = this->liquid_saturation_ * this->PIC_liquid_pressure_ +
                            this->gas_saturation_ * this->PIC_gas_pressure_;

  }
  return status;
}

// Compute wave pressure on the seabed surface
template <unsigned Tdim>
double mpm::ThreePhaseParticleNew<Tdim>::pf_seabed_surface_particle(){
    const double omega = 2.0 * M_PI / this->wave_period_;
    this->wave_pf_.assign(this->Nx_,0.0);
    for (int i = 0; i < this->Nx_; ++i) {
        this->wave_pf_[i] = this->wave_P0_[i] * std::cos(this->wave_phi_[i] - omega * this->current_time_);
    }

    // interp1(xg, yg, xq-need)
    const std::size_t n = seabed_surface_x_.size();
    if (n == 0) return 0.0;
    if (n == 1) return this->wave_pf_[0];

    if (this->coordinates_[0] <= this->seabed_surface_x_.front()) return this->wave_pf_.front();
    if (this->coordinates_[0] >= this->seabed_surface_x_.back())  return this->wave_pf_.back();

    auto it = std::upper_bound(seabed_surface_x_.begin(), seabed_surface_x_.end(), this->coordinates_[0]);
    std::size_t j = std::size_t(it - seabed_surface_x_.begin());
    std::size_t i = j - 1;

    const double x1 = seabed_surface_x_[i], y1 = this->wave_pf_[i];
    const double x2 = seabed_surface_x_[j], y2 = this->wave_pf_[j];
    const double t  = (this->coordinates_[0] - x1) / (x2 - x1);
    return y1 + this->current_time_ * (y2 - y1);
}


template <unsigned Tdim>
void mpm::ThreePhaseParticleNew<Tdim>::output_results(double step_) {
    // Output directory (change if needed)
    const std::string dir_path =
        "/home/chen/mpm/examples/1D_consolidation/MATLABinput";

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
