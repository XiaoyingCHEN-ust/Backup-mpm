#ifndef MPM_FROZEN_PARTICLE_UNSATURATED_H_
#define MPM_FROZEN_PARTICLE_UNSATURATED_H_

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include "logger.h"
#include "particle.h"

namespace mpm {

// UnsatFrozenParticle class
template <unsigned Tdim>
class UnsatFrozenParticle : public mpm::Particle<Tdim> {

public:
  // Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //============================================================================
  // CONSTRUCT AND DESTRUCT A PARTICLE

  // Construct a particle with id and coordinates
  UnsatFrozenParticle(Index id, const VectorDim& coord);

  // Destructor
  ~UnsatFrozenParticle() override{};

  // Delete copy constructor
  UnsatFrozenParticle(const UnsatFrozenParticle<Tdim>&) = delete;

  // Delete assignment operator
  UnsatFrozenParticle& operator = (const UnsatFrozenParticle<Tdim>&) = delete;

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

  // Assign degree of liquid water saturation
  bool assign_liquid_saturation_degree();

  // Compute both solid and liquid mass
  void compute_mass() override;

  // Initial pore pressure
  void initial_pore_pressure(double pore_pressure) override {
    this->pore_pressure_ = pore_pressure;
  }

  // Assign particles initial pore pressure by watertable
  bool initialise_pore_pressure_watertable(
      const unsigned dir_v, const unsigned dir_h,
      std::map<double, double>& refernece_points);

  //============================================================================
  // APPLY BOUNDARY CONDITIONS

  // Assign traction to the particle
  bool assign_particle_traction(unsigned direction, double traction) override;

  // Assign contact to the particle
  bool assign_particle_contact(unsigned dir, double normal) override; 

  // Assign particle liquid phase velocity constraints
  bool assign_particle_liquid_velocity_constraint(unsigned dir,
                                                  double velocity) override;

  // Apply particle liquid phase velocity constraints
  void apply_particle_liquid_velocity_constraints() override;

  // Assign particle pressure constraints
  bool assign_particle_pore_pressure_constraint(double pressure) override;

  // Apply particle pore pressure constraints
  void apply_particle_pore_pressure_constraints(double pore_pressure) override;

  //============================================================================
  // MAP PARTICLE INFORMATION TO NODES

  // Map particle mass and momentum to nodes (both solid and liquid)
  void map_mass_momentum_to_nodes() noexcept override;

  // Map body force
  void map_external_force(const VectorDim& pgravity) override;

  // Map traction force
  void map_traction_force() noexcept override;

  // Map internal force
  void map_internal_force_semi (double beta = 1) noexcept override;

  // Assign pore pressure to nodes
  bool map_pore_pressure_to_nodes(double current_time = 0.) noexcept override;

  // Assign gas saturation to nodes
  bool map_gas_saturation_to_nodes() noexcept override; 

  // Map drag force coefficient
  void map_drag_force_coefficient() override;

  // Map particle heat capacity and heat to nodes
  void map_heat_to_nodes() override;

  // Map particle heat conduction to node
  void map_heat_conduction() override;

  // Map heat convection of mixture
  void map_heat_convection() override;

  // Map latent heat of ice phase
  void map_latent_heat(double dt) noexcept override; 

  // Map heat convection of mixture
  void map_covective_heat_flux(double current_time) noexcept override;

  // Map reaction direction to nodes
  void map_moving_rigid_velocity_to_nodes(unsigned dir, double velocity, 
                                          double dt) noexcept override;  

  // map moving rigid velocity and surface normal to nodes
  void map_rigid_mass_momentum_to_nodes() noexcept override;  

  //------------------------------------------------------------
  // Semi-implict mpm

  // Map K inter element matrix to cell
  bool map_K_inter_to_cell() override;

  // Map laplacian element matrix to cell
  bool map_L_to_cell(double dt, double alpha) override;

  // Map F element matrix to cell (used in poisson equation RHS)
  bool map_F_to_cell() override;

  // Map T element matrix to cell (used in poisson equation RHS)
  bool map_T_to_cell() override;

  // // Map T element matrix to cell (used in poisson equation RHS)
  // bool map_P_to_cell(double beta) override;  

  // Map K_cor element matrix to cell
  bool map_K_cor_to_cell(double dt, double alpha) override;

  //==========================================================================
  // UPDATE PARTICLE INFORMATION

  // Compute updated velocity and position of the particle
  void compute_updated_velocity(double dt, double pic = 0,
                                double damping_factor = 0) override;

  // Map nodal pore pressure to particles
  bool compute_updated_pore_pressure(double beta) override;

  // Compute pore pressure somoothening by interpolating nodal pressure
  bool compute_pore_pressure_smoothing() noexcept override;

  // Compute ga saturation somoothening
  bool compute_gas_saturation_smoothing() noexcept override;  

  // //! Compute pore liquid pressure and pore ice pressure at particles
  // bool compute_liquid_ice_pore_pressure(double beta) override; 

  //! Compute thermal strain
  void compute_ice_thermal_strain() noexcept override;

  // //! Compute frost heave
  // void compute_frost_heave_strain() noexcept override;

  // Update gas saturation
  void update_gas_saturation(double dt) override;  

  // Update liquid water saturation
  void update_liquid_saturation(double dt) override;

  // Update viscosity
  void update_viscosity() noexcept override;

  // Update particle permeability
  void update_permeability() override;

  // Update density of the particle
  void update_particle_density(double dt) override;

  // Update mass of the particle
  void update_particle_volume()  override;    

//============================================================================
// RETURN PARTICLE DATA

  // Return liquid pore pressure
  double pore_pressure() const override { return pore_pressure_; }

  // Return liquid_saturation
  double liquid_saturation() const override { return liquid_saturation_; }  

  // Return ice_saturation
  double ice_saturation() const override {
    double sign = 1; 
    if (ice_saturation_rate_ < 0) sign = -1; 
    return ice_saturation_ * sign;
    }    

  //! Return liquid pore pressure
  double pore_ice_pressure() const override { return pore_ice_pressure_; }        

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
  using Particle<Tdim>::temperature_gradient_;
  using Particle<Tdim>::temperature_acceleration_;
  using Particle<Tdim>::dthermal_strain_;
  using Particle<Tdim>::dthermal_volumetric_strain_;
  using Particle<Tdim>::thermal_strain_;
  using Particle<Tdim>::thermal_volumetric_strain_;
  using Particle<Tdim>::stress_;
  using Particle<Tdim>::strain_rate_;
  using Particle<Tdim>::velocity_;
  // To be deleted
  using Particle<Tdim>::contact_normal_;
  using Particle<Tdim>::contact_tangential_;
  using Particle<Tdim>::set_contact_; 

  // Liquid Material
  std::shared_ptr<Material<Tdim>> liquid_material_;
  unsigned liquid_material_id_{std::numeric_limits<unsigned>::max()};

  // Liquid property
  double liquid_density_;
  double liquid_fraction_;
  double liquid_mass_density_;
  double liquid_mass_;
  double liquid_saturation_;
  double smoothed_liquid_saturation_;  
  double liquid_saturation_min_;  
  double gas_density_;
  double gas_fraction_;
  double gas_mass_density_;
  double gas_mass_;
  double gas_saturation_;
  double smoothed_gas_saturation_; 
  double pore_pressure_;
  double pore_pressure_increment_;  
  double permeability_;
  double viscosity_;
  double dSl_dT_;
  double ice_density_;
  double ice_fraction_;
  double smoothed_ice_saturation_;
  double ice_saturation_;    
  double ice_mass_density_;
  double ice_mass_;
  double ice_saturation_rate_;

  // Vector property
  Eigen::Matrix<double, Tdim, 1> liquid_velocity_;
  Eigen::Matrix<double, Tdim, 1> liquid_flux_;
  Eigen::Matrix<double, 6, 1> liquid_strain_rate_;
  Eigen::Matrix<double, 6, 1> liquid_strain_;

  // Boundary conditions
  bool set_mixture_traction_;  
  Eigen::Matrix<double, Tdim, 1> mixture_traction_;
  std::map<unsigned, double> liquid_velocity_constraints_;
  double pore_pressure_constraint_{std::numeric_limits<unsigned>::max()};

  // Logger
  std::unique_ptr<spdlog::logger> console_;
  
  // To be deleted
  double pore_liquid_pressure_;
  double pore_liquid_pressure_increment_;
  double pore_ice_pressure_;
  double pore_ice_pressure_increment_;
  double ice_fraction_increment_;
  double dheave_volumetric_strain_;
  Eigen::Matrix<double, 6, 1> dheave_strain_;

private:

  // Assign particle permeability
  virtual bool assign_permeability();

  // Compute liquid mass
  virtual void compute_liquid_mass() noexcept;    

  // Assign mixture traction
  virtual bool assign_mixture_traction(unsigned direction, double traction);

  // Assign liquid mass and momentum to nodes
  virtual void map_liquid_mass_momentum_to_nodes() noexcept;

  // Assign ice mass and momentum to nodes
  virtual void map_ice_mass_momentum_to_nodes() noexcept; 

  // Map two phase mixture body force
  virtual void map_mixture_body_force(unsigned mixture,
                                      const VectorDim& pgravity) noexcept;

  // Map liquid body force
  virtual void map_liquid_body_force(const VectorDim& pgravity) noexcept;

  // Map two phase mixture traction force
  virtual void map_mixture_traction_force(unsigned mixture) noexcept;

  // Map liquid internal force
  virtual void map_liquid_internal_force(double beta = 1);

  // Map two phase mixture internal force
  virtual void map_mixture_internal_force(unsigned mixture, double beta = 1);

  // Map liquid heat capacity and heat to nodes
  virtual void map_liquid_heat_to_nodes() noexcept;

  // Map liquid phase heat conduction 
  virtual void map_liquid_heat_conduction() noexcept;

  // Compute updated liquid velocity
  virtual void compute_updated_liquid_velocity(
      double dt, double pic = 0., double damping_factor = 0.);

  // update liquid density of the particle
  virtual void update_liquid_density() noexcept;

  // update liquid mass of particle
  virtual void update_liquid_mass() noexcept;

};  // UnsatFrozenParticle class
}  // namespace mpm

#include "Frozen_particle_unsaturated.tcc"

#endif  // MPM_TWOPHASE_PARTICLE_UNSATURATED_H__
