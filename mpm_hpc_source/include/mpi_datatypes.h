#ifndef MPM_MPI_HDF5_PARTICLE_H_
#define MPM_MPI_HDF5_PARTICLE_H_

// Offset macro
#include <cstddef>

// MPI
#ifdef USE_MPI
#include "mpi.h"

#include "hdf5_particle.h"

namespace mpm {
//! Initialise MPI particle data types
inline MPI_Datatype register_mpi_particle_type(const HDF5Particle& particle) {
  // Number of blocks to create
  const unsigned nblocks = 64;
  // Array containing the length of each block
  int lengths[nblocks] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 20};
  // Array containing the displacement for each block, expressed in bytes. The
  // displacement is the distance between the start of the MPI datatype created
  // and the start of the block
  const std::size_t size_ull = sizeof(unsigned long long);

  const MPI_Aint displacements[nblocks] = {
      offsetof(HDF5Particle, id),                    // id
      offsetof(HDF5Particle, mass),                  // mass
      offsetof(HDF5Particle, liquid_mass),           // liquid_mass
      offsetof(HDF5Particle, volume),                // volume
      offsetof(HDF5Particle, pressure),              // pressure
      offsetof(HDF5Particle, pore_pressure),         // pore_pressure
      offsetof(HDF5Particle, liquid_pressure),       // liquid pressure
      offsetof(HDF5Particle, gas_pressure),          // gas pressure      
      offsetof(HDF5Particle, coord_x),               // coord_x
      offsetof(HDF5Particle, coord_y),               // coord_y
      offsetof(HDF5Particle, coord_z),               // coord_z
      offsetof(HDF5Particle, displacement_x),        // disp_x
      offsetof(HDF5Particle, displacement_y),        // disp_y
      offsetof(HDF5Particle, displacement_z),        // disp_z
      offsetof(HDF5Particle, nsize_x),               // nsize_x
      offsetof(HDF5Particle, nsize_y),               // nsize_y
      offsetof(HDF5Particle, nsize_z),               // nsize_z
      offsetof(HDF5Particle, velocity_x),            // vel_x
      offsetof(HDF5Particle, velocity_y),            // vel_y
      offsetof(HDF5Particle, velocity_z),            // vel_z
      offsetof(HDF5Particle, liquid_velocity_x),     // liquid_vel_x
      offsetof(HDF5Particle, liquid_velocity_y),     // liquid_vel_y
      offsetof(HDF5Particle, liquid_velocity_z),     // liquid_vel_z
      offsetof(HDF5Particle, gas_velocity_y),        // gas_vel_y
      offsetof(HDF5Particle, gas_velocity_x),        // gas_vel_x      
      offsetof(HDF5Particle, gas_velocity_z),        // gas_vel_z      
      offsetof(HDF5Particle, stress_xx),             // stress_xx
      offsetof(HDF5Particle, stress_yy),             // stress_yy
      offsetof(HDF5Particle, stress_zz),             // stress_zz
      offsetof(HDF5Particle, tau_xy),                // tau_xy
      offsetof(HDF5Particle, tau_yz),                // tau_yz
      offsetof(HDF5Particle, tau_xz),                // tau_xz
      offsetof(HDF5Particle, strain_xx),             // strain_xx
      offsetof(HDF5Particle, strain_yy),             // strain_yy
      offsetof(HDF5Particle, strain_zz),             // strain_zz
      offsetof(HDF5Particle, gamma_xy),              // gamma_xy
      offsetof(HDF5Particle, gamma_yz),              // gamma_yz
      offsetof(HDF5Particle, gamma_xz),              // gamma_xz
      offsetof(HDF5Particle, thermal_strain_xx),     // thermal_strain_xx
      offsetof(HDF5Particle, thermal_strain_yy),     // thermal_strain_yy
      offsetof(HDF5Particle, thermal_strain_zz),     // thermal_strain_zz
      offsetof(HDF5Particle, thermal_gamma_xy),      // thermal_gamma_xy
      offsetof(HDF5Particle, thermal_gamma_yz),      // thermal_gamma_yz
      offsetof(HDF5Particle, thermal_gamma_xz),      // thermal_gamma_xz
      offsetof(HDF5Particle, liquid_strain_xx),      // liquid_strain_xx
      offsetof(HDF5Particle, liquid_strain_yy),      // liquid_strain_yy
      offsetof(HDF5Particle, liquid_strain_zz),      // liquid_strain_zz
      offsetof(HDF5Particle, liquid_gamma_xy),       // liquid_gamma_xy
      offsetof(HDF5Particle, liquid_gamma_yz),       // liquid_gamma_yz
      offsetof(HDF5Particle, liquid_gamma_xz),       // liquid_gamma_xz
      offsetof(HDF5Particle, epsilon_v),             // epsv
      offsetof(HDF5Particle, cell_id),               // cellid
      offsetof(HDF5Particle, thermal_epsilon_v),     // thermal_epsv
      offsetof(HDF5Particle, status),                // status
      offsetof(HDF5Particle, material_id),           // material_id
      offsetof(HDF5Particle, liquid_material_id),    // liquid_material_id
      offsetof(HDF5Particle, temperature),           // temperature
      offsetof(HDF5Particle, temperature_increment), // temperature_increment
      offsetof(HDF5Particle, gradT_x),               // gradT_x
      offsetof(HDF5Particle, gradT_y),               // gradT_y
      offsetof(HDF5Particle, gradT_z),               // gradT_z
      offsetof(HDF5Particle, density),                // density
      offsetof(HDF5Particle, nstate_vars),           // nstate_vars
      offsetof(HDF5Particle, svars)                  // state_vars
      // ===== MIXTURE: new fields (9) =====
      offsetof(HDF5Particle, suction_pressure),     // suction_pressure
      offsetof(HDF5Particle, mixture_mass),         // mixture_mass
      offsetof(HDF5Particle, dSw_dpw),              // dSw_dpw
      offsetof(HDF5Particle, total_stress_xx),      // total_stress_xx
      offsetof(HDF5Particle, total_stress_yy),      // total_stress_yy
      offsetof(HDF5Particle, total_stress_zz),      // total_stress_zz
      offsetof(HDF5Particle, total_stress_tau_xy),  // total_stress_tau_xy
      offsetof(HDF5Particle, total_stress_tau_yz),  // total_stress_tau_yz
      offsetof(HDF5Particle, total_stress_tau_xz),  // total_stress_tau_xz

      // ===== LIQUID: new kinematics & rates (12) =====
      offsetof(HDF5Particle, liquid_acceleration_x), // liquid_acceleration_x
      offsetof(HDF5Particle, liquid_acceleration_y), // liquid_acceleration_y
      offsetof(HDF5Particle, liquid_acceleration_z), // liquid_acceleration_z
      offsetof(HDF5Particle, liquid_strain_gamma_xy),// liquid_strain_gamma_xy
      offsetof(HDF5Particle, liquid_strain_gamma_yz),// liquid_strain_gamma_yz
      offsetof(HDF5Particle, liquid_strain_gamma_xz),// liquid_strain_gamma_xz
      offsetof(HDF5Particle, liquid_strain_rate_xx), // liquid_strain_rate_xx
      offsetof(HDF5Particle, liquid_strain_rate_yy), // liquid_strain_rate_yy
      offsetof(HDF5Particle, liquid_strain_rate_zz), // liquid_strain_rate_zz
      offsetof(HDF5Particle, liquid_strain_rate_xy), // liquid_strain_rate_xy
      offsetof(HDF5Particle, liquid_strain_rate_yz), // liquid_strain_rate_yz
      offsetof(HDF5Particle, liquid_strain_rate_xz), // liquid_strain_rate_xz

      // ===== LIQUID: new scalars (8) =====
      offsetof(HDF5Particle, effective_saturation),      // effective_saturation
      offsetof(HDF5Particle, liquid_chi),                // liquid_chi
      offsetof(HDF5Particle, liquid_volume),
      offsetof(HDF5Particle, liquid_mass_density),       // liquid_mass_density
      offsetof(HDF5Particle, liquid_pressure_acceleration), // liquid_pressure_acceleration
      offsetof(HDF5Particle, liquid_volumetric_strain),  // liquid_volumetric_strain
      offsetof(HDF5Particle, PIC_liquid_pressure),       // PIC_liquid_pressure
      offsetof(HDF5Particle, FLIP_liquid_pressure),      // FLIP_liquid_pressure

      // ===== GAS: new kinematics & rates (18) =====
      offsetof(HDF5Particle, gas_acceleration_x),    // gas_acceleration_x
      offsetof(HDF5Particle, gas_acceleration_y),    // gas_acceleration_y
      offsetof(HDF5Particle, gas_acceleration_z),    // gas_acceleration_z
      offsetof(HDF5Particle, pgravity_x),            // pgravity_x
      offsetof(HDF5Particle, pgravity_y),            // pgravity_y
      offsetof(HDF5Particle, pgravity_z),            // pgravity_z
      offsetof(HDF5Particle, gas_strain_xx),         // gas_strain_xx
      offsetof(HDF5Particle, gas_strain_yy),         // gas_strain_yy
      offsetof(HDF5Particle, gas_strain_zz),         // gas_strain_zz
      offsetof(HDF5Particle, gas_strain_gamma_xy),   // gas_strain_gamma_xy
      offsetof(HDF5Particle, gas_strain_gamma_yz),   // gas_strain_gamma_yz
      offsetof(HDF5Particle, gas_strain_gamma_xz),   // gas_strain_gamma_xz
      offsetof(HDF5Particle, gas_strain_rate_xx),    // gas_strain_rate_xx
      offsetof(HDF5Particle, gas_strain_rate_yy),    // gas_strain_rate_yy
      offsetof(HDF5Particle, gas_strain_rate_zz),    // gas_strain_rate_zz
      offsetof(HDF5Particle, gas_strain_rate_xy),    // gas_strain_rate_xy
      offsetof(HDF5Particle, gas_strain_rate_yz),    // gas_strain_rate_yz
      offsetof(HDF5Particle, gas_strain_rate_xz),    // gas_strain_rate_xz

      // ===== GAS: new scalars (8) =====
      offsetof(HDF5Particle, gas_chi),                   // gas_chi
      offsetof(HDF5Particle, gas_volume),                // gas_volume
      offsetof(HDF5Particle, gas_mass_density),          // gas_mass_density
      offsetof(HDF5Particle, gas_pressure_acceleration), // gas_pressure_acceleration
      offsetof(HDF5Particle, PIC_gas_pressure),          // PIC_gas_pressure
      offsetof(HDF5Particle, FLIP_gas_pressure),         // FLIP_gas_pressure
      offsetof(HDF5Particle, gas_pressure_increment),    // gas_pressure_increment
      offsetof(HDF5Particle, gas_volumetric_strain),     // gas_volumetric_strain

      // ===== WAVE: new fields (2) =====
      offsetof(HDF5Particle, wave_pressure), // wave_pressure
      offsetof(HDF5Particle, gamma_b)        // gamma_b
      
  };

  // Array containing the MPI datatypes to replicate to make each block.
  MPI_Datatype types[nblocks] = {
      MPI_UNSIGNED_LONG_LONG,  // id
      MPI_DOUBLE,              // mass
      MPI_DOUBLE,              // liquid_mass
      MPI_DOUBLE,              // volume
      MPI_DOUBLE,              // pressure
      MPI_DOUBLE,              // pore_pressure
      MPI_DOUBLE,              // lqiuid pressure
      MPI_DOUBLE,              // gas pressure      
      MPI_DOUBLE,              // coord_x
      MPI_DOUBLE,              // coord_y
      MPI_DOUBLE,              // coord_z
      MPI_DOUBLE,              // disp_x
      MPI_DOUBLE,              // disp_y
      MPI_DOUBLE,              // disp_z
      MPI_DOUBLE,              // nsize_x
      MPI_DOUBLE,              // nsize_y
      MPI_DOUBLE,              // nsize_z
      MPI_DOUBLE,              // vel_x
      MPI_DOUBLE,              // vel_y
      MPI_DOUBLE,              // vel_z
      MPI_DOUBLE,              // liquid_vel_x
      MPI_DOUBLE,              // liquid_vel_y
      MPI_DOUBLE,              // liquid_vel_z
      MPI_DOUBLE,              // gas_vel_x
      MPI_DOUBLE,              // gas_vel_y
      MPI_DOUBLE,              // gas_vel_z      
      MPI_DOUBLE,              // stress_xx
      MPI_DOUBLE,              // stress_yy
      MPI_DOUBLE,              // stress_zz
      MPI_DOUBLE,              // tau_xx
      MPI_DOUBLE,              // tau_yy
      MPI_DOUBLE,              // tau_zz
      MPI_DOUBLE,              // strain_xx
      MPI_DOUBLE,              // strain_yy
      MPI_DOUBLE,              // strain_zz
      MPI_DOUBLE,              // gamma_xy
      MPI_DOUBLE,              // gamma_yz
      MPI_DOUBLE,              // gamma_zx
      MPI_DOUBLE,              // thermal_strain_xx
      MPI_DOUBLE,              // thermal_strain_yy
      MPI_DOUBLE,              // thermal_strain_zz
      MPI_DOUBLE,              // thermal_gamma_xy
      MPI_DOUBLE,              // thermal_gamma_yz
      MPI_DOUBLE,              // thermal_gamma_zx
      MPI_DOUBLE,              // liquid_strain_xx
      MPI_DOUBLE,              // liquid_strain_yy
      MPI_DOUBLE,              // liquid_strain_zz
      MPI_DOUBLE,              // liquid_gamma_xy
      MPI_DOUBLE,              // liquid_gamma_yz
      MPI_DOUBLE,              // liquid_gamma_zx
      MPI_DOUBLE,              // epsv
      MPI_DOUBLE,              // thermal_epsv
      MPI_UNSIGNED_LONG_LONG,  // cell_id
      MPI_C_BOOL,              // status
      MPI_UNSIGNED,            // material_id
      MPI_UNSIGNED,            // liquid_material_id
      MPI_DOUBLE,              // temperature
      MPI_DOUBLE,              // temperature_increment
      MPI_DOUBLE,              // gradT_x
      MPI_DOUBLE,              // gradT_y
      MPI_DOUBLE,              // gradT_z
      MPI_DOUBLE,              // density
      MPI_UNSIGNED,            // nstate_vars
      MPI_DOUBLE,              // state variables
      // ===== MIXTURE: new fields (9) =====
      MPI_DOUBLE, // suction_pressure
      MPI_DOUBLE, // mixture_mass
      MPI_DOUBLE, // dSw_dpw
      MPI_DOUBLE, // total_stress_xx
      MPI_DOUBLE, // total_stress_yy
      MPI_DOUBLE, // total_stress_zz
      MPI_DOUBLE, // total_stress_tau_xy
      MPI_DOUBLE, // total_stress_tau_yz
      MPI_DOUBLE, // total_stress_tau_xz

      // ===== LIQUID: new kinematics & rates (12) =====
      MPI_DOUBLE, // liquid_acceleration_x
      MPI_DOUBLE, // liquid_acceleration_y
      MPI_DOUBLE, // liquid_acceleration_z
      MPI_DOUBLE, // liquid_strain_gamma_xy
      MPI_DOUBLE, // liquid_strain_gamma_yz
      MPI_DOUBLE, // liquid_strain_gamma_xz
      MPI_DOUBLE, // liquid_strain_rate_xx
      MPI_DOUBLE, // liquid_strain_rate_yy
      MPI_DOUBLE, // liquid_strain_rate_zz
      MPI_DOUBLE, // liquid_strain_rate_xy
      MPI_DOUBLE, // liquid_strain_rate_yz
      MPI_DOUBLE, // liquid_strain_rate_xz

      // ===== LIQUID: new scalars (8) =====
      MPI_DOUBLE, // effective_saturation
      MPI_DOUBLE, // liquid_chi
      MPI_DOUBLE, // liquid_volume
      MPI_DOUBLE, // liquid_mass_density
      MPI_DOUBLE, // liquid_pressure_acceleration
      MPI_DOUBLE, // liquid_volumetric_strain
      MPI_DOUBLE, // PIC_liquid_pressure
      MPI_DOUBLE, // FLIP_liquid_pressure

      // ===== GAS: new kinematics & rates (18) =====
      MPI_DOUBLE, // gas_acceleration_x
      MPI_DOUBLE, // gas_acceleration_y
      MPI_DOUBLE, // gas_acceleration_z
      MPI_DOUBLE, // pgravity_x
      MPI_DOUBLE, // pgravity_y
      MPI_DOUBLE, // pgravity_z
      MPI_DOUBLE, // gas_strain_xx
      MPI_DOUBLE, // gas_strain_yy
      MPI_DOUBLE, // gas_strain_zz
      MPI_DOUBLE, // gas_strain_gamma_xy
      MPI_DOUBLE, // gas_strain_gamma_yz
      MPI_DOUBLE, // gas_strain_gamma_xz
      MPI_DOUBLE, // gas_strain_rate_xx
      MPI_DOUBLE, // gas_strain_rate_yy
      MPI_DOUBLE, // gas_strain_rate_zz
      MPI_DOUBLE, // gas_strain_rate_xy
      MPI_DOUBLE, // gas_strain_rate_yz
      MPI_DOUBLE, // gas_strain_rate_xz

      // ===== GAS: new scalars (8) =====
      MPI_DOUBLE, // gas_chi
      MPI_DOUBLE, // gas_volume
      MPI_DOUBLE, // gas_mass_density
      MPI_DOUBLE, // gas_pressure_acceleration
      MPI_DOUBLE, // PIC_gas_pressure
      MPI_DOUBLE, // FLIP_gas_pressure
      MPI_DOUBLE, // gas_pressure_increment
      MPI_DOUBLE, // gas_volumetric_strain

      // ===== WAVE: new fields (2) =====
      MPI_HBOOL, // wave_pressure
      MPI_DOUBLE, // gamma_b
  };

  // Create the Particle datatype
  MPI_Datatype MPIParticle;

  // Create particle data types
  MPI_Type_create_struct(nblocks, lengths, displacements, types, &MPIParticle);
  MPI_Type_commit(&MPIParticle);
  return MPIParticle;
}

//! Deregister MPI particle data type
inline void deregister_mpi_particle_type(MPI_Datatype& particle_type) {
  MPI_Type_free(&particle_type);
}
}  // namespace mpm

#endif  // MPI

#endif  // MPM_MPI_HDF5_PARTICLE_H_
