#ifndef MPM_MHBS_PARTICLE_WITH_PHASE_CHANGE_H_
#define MPM_MHBS_PARTICLE_WITH_PHASE_CHANGE_H_

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "logger.h"
#include "particle.h"

namespace mpm {

// MHBSParticle class
template <unsigned Tdim>
class MHBSParticle : public mpm::Particle<Tdim> {

public:
  // Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //============================================================================
  // CONSTRUCT AND DESTRUCT A PARTICLE

  // Construct a particle with id and coordinates
  MHBSParticle(Index id, const VectorDim& coord);

  // Destructor
  ~MHBSParticle() override{};

  // Delete copy constructor
  MHBSParticle(const MHBSParticle<Tdim>&) = delete;

  // Delete assignment operator
  MHBSParticle& operator = (const MHBSParticle<Tdim>&) = delete;

  //============================================================================
  // ASSIGN INITIAL CONDITIONS

  // Initialise particle from HDF5 data
  bool initialise_particle(const HDF5Particle& particle) override;

  // Initialise particle HDF5 data and material
  bool initialise_particle(
      const HDF5Particle& particle,
      const std::shared_ptr<Material<Tdim>>& material) override;

  // Initialise liquid phase particle properties
  void initialise_liquid_phase() override;

  // Retrun particle data as HDF5
  HDF5Particle hdf5() override;

  // Assign material
  bool assign_liquid_material(
      const std::shared_ptr<Material<Tdim>>& material) override;

  // Compute both solid and liquid mass
  void compute_mass() override;

  // Assign initial properties to particle
  bool assign_initial_properties();

  // Initial pore pressure
  void initial_pore_pressure(double pore_pressure) override {
    this->pore_pressure_ = pore_pressure;
    this->liquid_pressure_ = pore_pressure;
    this->gas_pressure_ = pore_pressure;
    this->PIC_liquid_pressure_ = pore_pressure;
    this->PIC_gas_pressure_ = pore_pressure;

  }

  //============================================================================
  // APPLY BOUNDARY CONDITIONS

  // Assign traction to the particle
  bool assign_particle_traction(unsigned direction, double traction) override;

  // // Assign contact to the particle
  // bool assign_particle_contact(unsigned dir, double normal) override;    

  // Assign particle liquid phase velocity constraints
  bool assign_particle_liquid_velocity_constraint(unsigned dir,
                                                  double velocity) override;

  // Apply particle liquid phase velocity constraints
  void apply_particle_liquid_velocity_constraints() override;

  // Assign particle pressure constraints
  bool assign_particle_pore_pressure_constraint(double pressure) override;

  // Apply particle pore pressure constraints
  void apply_particle_pore_pressure_constraints(double pore_pressure) override;

  // Overwrite node velocity to get strain correct
  void map_moving_rigid_velocity_to_nodes(unsigned dir, 
                        double velocity, double dt) noexcept override;

  //============================================================================
  // MAP PARTICLE INFORMATION TO NODES

  // Map particle mass and momentum to nodes (both solid and liquid)
  void map_mass_momentum_to_nodes() noexcept override;

  // Map particle mass and pressure to nodes
  void map_mass_pressure_to_nodes() override;

  // Map particle scaler properties to nodes
  void map_scalers_to_nodes() override;

  // Assign density to nodes
  void map_density_to_nodes() override;    

  // Map body force
  void map_external_force(const VectorDim& pgravity) override;

  // Map internal force
  void map_internal_force () override;

  // Map drag force coefficient
  void map_drag_force_coefficient() override;

  // Map particle heat capacity and heat to nodes
  void map_heat_to_nodes() override;

  // Map particle heat conduction to node
  void map_heat_conduction() override;

  // Map particle hydraulic conduction to node
  void map_hydraulic_conduction() override;

  // Map heat convection of mixture
  void map_heat_convection() override;

  // Map heat convection of mixture
  void map_mass_convection() override;

  // Map heat source
  void map_heat_source() override;

  // Map volumetric strain
  void map_volumetric_strain() override; 

  // Assign pore pressure to nodes
  bool map_pore_pressure_to_nodes(double current_time = 0.) noexcept override;

  //==========================================================================
  // UPDATE PARTICLE INFORMATION

  // compute updated pore pressure
  void update_particle_pore_pressure(double dt, double pic_t) noexcept override;

  // Compute updated velocity and position of the particle
  void compute_updated_velocity(double dt, double pic = 0,
                                double damping_factor = 0) override;

  // Map nodal pore pressure to particles
  void compute_pore_pressure(double dt) override;

  // Update mass of the particle
  void update_particle_volume() override;

  // Update particle permeability
  void update_permeability() override;

  // Update source term
  void update_source_term() override;

  // Update hydrate saturation
  void update_hydrate_saturation(double dt) override;

  // Update liquid water saturation
  void update_liquid_saturation(double dt) override;

  // Update density of the particle
  void update_particle_density(double dt) override;

  // Compute pore pressure somoothening by interpolating nodal pressure
  bool compute_pore_pressure_smoothing() noexcept override;

  // Compute smoothed scaler values by interpolating nodal values
  bool compute_scalers_smoothing() noexcept override;

//============================================================================
// RETURN PARTICLE DATA

  // Return liquid pore pressure
  double pore_pressure() const override { return pore_pressure_; }

  // Return liquid_saturation
  double liquid_saturation() const override { return liquid_saturation_; }

  // Return liquid_saturation
  double hydrate_saturation() const override { return hydrate_saturation_; }  

  // Return liquid_saturation
  double ini_hydrate_saturation() const override { return ini_hydrate_saturation_; }    

protected:

  // Map particle pressure gradient to node
  inline Eigen::Matrix<double, Tdim, 1> compute_pressure_gradient(
    unsigned phase) noexcept;

  // Map particle density gradient to node
  inline Eigen::Matrix<double, Tdim, 1> compute_density_gradient(
    unsigned phase) noexcept;    


protected:

  // Inherit properties from ParticleBase class 
  using ParticleBase<Tdim>::id_;
  using ParticleBase<Tdim>::coordinates_;
  using ParticleBase<Tdim>::cell_;
  using ParticleBase<Tdim>::nodes_;

  // Inherit properties from Particle class 
  using Particle<Tdim>::material_;
  using Particle<Tdim>::shapefn_;
  using Particle<Tdim>::shapefn_centroid_;
  using Particle<Tdim>::dn_dx_;  
  using Particle<Tdim>::dn_dx_centroid_; 
  using Particle<Tdim>::is_axisymmetric_; 
  using Particle<Tdim>::volume_;
  using Particle<Tdim>::density_;
  using Particle<Tdim>::mass_;  
  using Particle<Tdim>::mass_density_;
  using Particle<Tdim>::porosity_;
  using Particle<Tdim>::solid_fraction_;  
  using Particle<Tdim>::temperature_;
  using Particle<Tdim>::PIC_temperature_; 
  using Particle<Tdim>::temperature_increment_;
  using Particle<Tdim>::heat_source_;  
  using Particle<Tdim>::temperature_gradient_;
  using Particle<Tdim>::temperature_acceleration_;
  using Particle<Tdim>::dvolumetric_strain_;  
  using Particle<Tdim>::dthermal_strain_;
  using Particle<Tdim>::dthermal_volumetric_strain_;
  using Particle<Tdim>::thermal_strain_;
  using Particle<Tdim>::thermal_volumetric_strain_;
  using Particle<Tdim>::stress_;
  using Particle<Tdim>::strain_rate_;
  using Particle<Tdim>::velocity_;
  using Particle<Tdim>::acceleration_;
  // To be deleted
  using Particle<Tdim>::contact_normal_;
  using Particle<Tdim>::contact_tangential_;
  using Particle<Tdim>::set_contact_;

  // Liquid Material
  std::shared_ptr<Material<Tdim>> liquid_material_;
  unsigned liquid_material_id_{std::numeric_limits<unsigned>::max()};
  double pore_pressure_;
  double PIC_pore_pressure_;
  double ini_pore_pressure_;
  double suction_pressure_;
  double ini_suction_pressure_;
  double mixture_mass_;
  double ini_porosity_;
  double intrinsic_permeability_;
  double dSw_dpw_;
  Eigen::Matrix<double, 6, 1> total_stress_;
  double PIC_porosity_;
  double PIC_volume_;
  double PIC_volumetric_strain_;
  double PIC_deviatoric_strain_;
  double PIC_mean_stress_;
  double PIC_deviatoric_stress_;

  // Solid properties
  double solid_thermal_conductivity_;
  double solid_specific_heat_;
  double solid_heat_capacity_;  
  double solid_expansivity_;

  // Hydrate properties
  double hydrate_saturation_;
  double PIC_hydrate_saturation_;
  double hydrate_fraction_;
  double hydrate_density_;
  double hydrate_volume_;
  double hydrate_mass_;
  double hydrate_mass_density_;
  double hydrate_source_;

  double hydrate_thermal_conductivity_;
  double hydrate_specific_heat_;
  double hydrate_heat_capacity_;  
  double hydrate_expansivity_;
  double hydrate_molar_mass_;
  double hydrate_latent_;
  double hydration_number_;
  double ini_hydrate_density_;
  double ini_hydrate_saturation_;

  // Liquid property
  Eigen::Matrix<double, Tdim, 1> liquid_velocity_;
  Eigen::Matrix<double, Tdim, 1> liquid_acceleration_;
  Eigen::Matrix<double, Tdim, 1> liquid_flux_;
  Eigen::Matrix<double, Tdim, 1> liquid_traction_;
  Eigen::Matrix<double, Tdim, 1> liquid_pressure_gradient_;
  Eigen::Matrix<double, Tdim, 1> liquid_density_gradient_;
  Eigen::Matrix<double, 6, 1> liquid_strain_;
  Eigen::Matrix<double, 6, 1> liquid_strain_rate_;
  double liquid_saturation_;
  double liquid_chi_;
  double liquid_fraction_;
  double liquid_density_;
  double liquid_volume_;
  double liquid_mass_;
  double liquid_mass_density_;
  double liquid_pressure_;
  double ini_liquid_pressure_;
  double liquid_pressure_acceleration_;
  double liquid_volumetric_strain_;
  double liquid_permeability_;
  double liquid_source_;
  double PIC_liquid_pressure_;
  double FLIP_liquid_pressure_;
  double liquid_pressure_increment_; 
  double liquid_critical_time_; 

  double liquid_thermal_conductivity_;
  double liquid_specific_heat_;
  double liquid_heat_capacity_;
  double liquid_expansivity_;
  double liquid_compressibility_;
  double liquid_viscosity_;
  double liquid_molar_mass_;
  double liquid_saturation_res_;
  double ini_liquid_density_;
  double ini_liquid_saturation_;
  double ini_liquid_viscosity_;
  double effective_saturation_;

  // Gas properties
  Eigen::Matrix<double, Tdim, 1> gas_velocity_;
  Eigen::Matrix<double, Tdim, 1> gas_acceleration_;
  Eigen::Matrix<double, Tdim, 1> gas_flux_;
  Eigen::Matrix<double, Tdim, 1> gas_traction_; 
  Eigen::Matrix<double, Tdim, 1> gas_pressure_gradient_;
  Eigen::Matrix<double, Tdim, 1> gas_density_gradient_;
  Eigen::Matrix<double, Tdim, 1> pgravity_;
  Eigen::Matrix<double, 6, 1> gas_strain_;
  Eigen::Matrix<double, 6, 1> gas_strain_rate_;
  double gas_saturation_;
  double gas_fraction_; 
  double gas_chi_;
  double gas_density_;
  double gas_volume_;
  double gas_mass_;
  double gas_mass_density_;

  double gas_pressure_;
  double ini_gas_pressure_;
  double gas_pressure_acceleration_;
  double PIC_gas_pressure_;
  double FLIP_gas_pressure_;
  double gas_pressure_increment_;
  double gas_volumetric_strain_;
  double gas_permeability_;
  double gas_source_;
  double gas_molar_mass_;
  double gas_thermal_conductivity_;
  double gas_specific_heat_;
  double gas_heat_capacity_;
  double gas_viscosity_;
  double gas_constant_;
  double gas_saturation_res_;
  double ini_gas_saturation_;
  double ini_gas_viscosity_;
  double gas_critical_time_; 

  // Boundary conditions
  bool set_mixture_traction_;
  bool set_pressure_constraint_;
  Eigen::Matrix<double, Tdim, 1> mixture_traction_;
  std::map<unsigned, double> liquid_velocity_constraints_;
  double pore_pressure_constraint_{std::numeric_limits<unsigned>::max()};

  // Logger
  std::unique_ptr<spdlog::logger> console_;
  bool debug_{false};

  Eigen::Matrix<double, 6, 1> K_matrix_;

};  // MHBSParticle class
}  // namespace mpm

#include "MHBS_particle.tcc"

#endif  // MPM_MHBS_PARTICLE_WITH_PHASE_CHANGE_H__
