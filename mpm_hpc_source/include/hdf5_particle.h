#ifndef MPM_HDF5_H_
#define MPM_HDF5_H_

#include <cstddef>

// HDF5
#include "hdf5.h"
#include "hdf5_hl.h"

#include "data_types.h"

namespace mpm {
// Define a struct of particle
typedef struct HDF5Particle {
  // Index
  mpm::Index id;
  // Index 
  mpm::Index cell_id;
  // Material id
  unsigned material_id;
  // Liquid material id
  unsigned liquid_material_id;
  // Current time
  double current_time;
  // Mass
  double mass;
  // Liquid mass
  double liquid_mass;
  // Ice mass
  double ice_mass; 
  // Hydrate mass
  double hydrate_mass;
  // Gas mass
  double gas_mass;   
  // Volume
  double volume;
  // Density
  double density;
  // Liquid density
  double liquid_density;  
  // Ice density
  double ice_density; 
  // Hydrate density
  double hydrate_density;  
  // Gas density
  double gas_density;   
  // Porosity
  double porosity;
  // Liquid water saturation
  double liquid_saturation;
  // Ice saturation
  double ice_saturation;  
  // Hydrate saturation
  double hydrate_saturation;
  // Gas saturation
  double gas_saturation;  
  // Liquid fraction
  double liquid_fraction;
  // Ice fraction
  double ice_fraction;
  // Hydrate fraction
  double hydrate_fraction;
  // Gas fraction
  double gas_fraction;
  // viscosity
  double viscosity;
  // Permeability
  double permeability;
  // Permeability
  double liquid_permeability;
   // Permeability
  double gas_permeability;     
  // Pressure
  double pressure;
  // Pore pressure
  double pore_pressure;
  // Pore pressure
  double pore_liquid_pressure;
  // Pore pressure
  double pore_ice_pressure;
  // Pore pressure
  double liquid_pressure;
  // Pore pressure
  double gas_pressure;
  // Temperature
  double temperature;
  // Temperature
  double PIC_temperature;  
  // Coordinates
  double coord_x, coord_y, coord_z;
  // Displacement
  double displacement_x, displacement_y, displacement_z;
  // Natural particle size
  double nsize_x, nsize_y, nsize_z;
  // Velocity
  double velocity_x, velocity_y, velocity_z;
  // Liquid velocity
  double liquid_velocity_x, liquid_velocity_y, liquid_velocity_z;
  // Gas velocity
  double gas_velocity_x, gas_velocity_y, gas_velocity_z;  
  // Stresses
  double stress_xx, stress_yy, stress_zz;
  double tau_xy, tau_yz, tau_xz;
  // Strains
  double strain_xx, strain_yy, strain_zz;
  double gamma_xy, gamma_yz, gamma_xz;
  // Thermal strains
  double thermal_strain_xx, thermal_strain_yy, thermal_strain_zz;
  double thermal_gamma_xy, thermal_gamma_yz, thermal_gamma_xz;
  // Liquid strain
  double liquid_strain_xx, liquid_strain_yy, liquid_strain_zz;
  double liquid_gamma_xy, liquid_gamma_yz, liquid_gamma_xz;
  // Volumetric strain centroid
  double epsilon_v;
  // Thermal volumetric strain centroid
  double thermal_epsilon_v;
  // Deformation gradient
  double fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz;
  // Fabric
  double fabric;
  // Rotation
  double rotation_xy, rotation_yz, rotation_xz;
  // State variables (init to zero)
  double svars[20] = {0};
  // Status
  bool status;
  // Number of state variables
  unsigned nstate_vars;
  
  // ===== Mixture (new) =====
  double suction_pressure;                 // matric suction
  double mixture_mass;                     // total mixture mass
  double dSw_dpw;                          // d(S_l)/d(p_l)
  double total_stress_xx, total_stress_yy, total_stress_zz;
  double total_stress_tau_xy, total_stress_tau_yz, total_stress_tau_xz;

  // ===== Liquid: kinematics & rates (new) =====
  double liquid_acceleration_x, liquid_acceleration_y, liquid_acceleration_z;
  double liquid_strain_gamma_xy, liquid_strain_gamma_yz, liquid_strain_gamma_xz;
  double liquid_strain_rate_xx, liquid_strain_rate_yy, liquid_strain_rate_zz;
  double liquid_strain_rate_xy, liquid_strain_rate_yz, liquid_strain_rate_xz;

  // ===== Liquid: scalars (new) =====
  double effective_saturation;             // S_e
  double liquid_chi;                       // Bishop's χ for liquid
  double liquid_volume;
  double liquid_mass_density;              // ρ_l (per mixture vol)
  double liquid_pressure_acceleration;     // ṗ_l contribution
  double liquid_volumetric_strain;         // ε_v^l
  double PIC_liquid_pressure;              // PIC update of p_l
  double FLIP_liquid_pressure;             // FLIP update of p_l

  // ===== Gas: kinematics & rates (new) =====
  double gas_acceleration_x, gas_acceleration_y, gas_acceleration_z;
  double pgravity_x, pgravity_y, pgravity_z;            // particle gravity vector
  double gas_strain_xx, gas_strain_yy, gas_strain_zz;
  double gas_strain_gamma_xy, gas_strain_gamma_yz, gas_strain_gamma_xz;
  double gas_strain_rate_xx, gas_strain_rate_yy, gas_strain_rate_zz;
  double gas_strain_rate_xy, gas_strain_rate_yz, gas_strain_rate_xz;

  // ===== Gas: scalars (new) =====
  double gas_chi;                         // Bishop's χ for gas
  double gas_volume;                      // gas_volume
  double gas_mass_density;                // ρ_g (per mixture vol)
  double gas_pressure_acceleration;       // ṗ_g contribution
  double PIC_gas_pressure;                // PIC update of p_g
  double FLIP_gas_pressure;               // FLIP update of p_g
  double gas_pressure_increment;          // Δp_g
  double gas_volumetric_strain;           // ε_v^g

  // ===== Wave (new) =====
  bool wave_pressure;                   // wave-induced pore pressure (value)
  double gamma_b;                         // breaker index or related coefficient
} HDF5Particle;

namespace hdf5::particle {
// Keep the compound-table schema append-only so older checkpoints retain the
// same field positions. Fourteen SANISAND state-variable fields are appended
// after the original 159 fields in hdf5_particle.cc.
const hsize_t LEGACY_NFIELDS = 159;
const hsize_t NFIELDS = 173;
const unsigned LEGACY_NSTATE_VARS = 6;

const size_t dst_size = sizeof(HDF5Particle);

// Destination offset
extern const size_t dst_offset[NFIELDS];

// Destination size
extern const size_t dst_sizes[NFIELDS];

// Define particle field information
extern const char* field_names[NFIELDS];

// Initialise field types
extern const hid_t field_type[NFIELDS];

//! Validate checkpoint table dimensions before reading into HDF5Particle
void validate_table_metadata(hsize_t field_count, hsize_t record_count,
                             hsize_t expected_records);

//! Validate restart history and zero state variables absent from legacy tables
void prepare_state_variables_for_restart(hsize_t field_count,
                                         HDF5Particle* particles,
                                         std::size_t particle_count);

}  // namespace hdf5::particle

}  // namespace mpm

#endif  // MPM_HDF5_H_
