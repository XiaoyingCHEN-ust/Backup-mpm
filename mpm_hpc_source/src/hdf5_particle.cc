#include "hdf5_particle.h"

#include <stdexcept>
#include <string>

namespace mpm {
namespace hdf5::particle {
const size_t dst_offset[NFIELDS] = {
    HOFFSET(HDF5Particle, id),
    HOFFSET(HDF5Particle, cell_id),
    HOFFSET(HDF5Particle, material_id),
    HOFFSET(HDF5Particle, liquid_material_id),
    HOFFSET(HDF5Particle, current_time),   
    HOFFSET(HDF5Particle, mass),
    HOFFSET(HDF5Particle, liquid_mass),
    HOFFSET(HDF5Particle, ice_mass),
    HOFFSET(HDF5Particle, hydrate_mass),
    HOFFSET(HDF5Particle, gas_mass),                    
    HOFFSET(HDF5Particle, volume),
    HOFFSET(HDF5Particle, density),
    HOFFSET(HDF5Particle, liquid_density),
    HOFFSET(HDF5Particle, ice_density),
    HOFFSET(HDF5Particle, hydrate_density),
    HOFFSET(HDF5Particle, gas_density),                
    HOFFSET(HDF5Particle, porosity),
    HOFFSET(HDF5Particle, liquid_saturation),
    HOFFSET(HDF5Particle, ice_saturation),
    HOFFSET(HDF5Particle, hydrate_saturation),
    HOFFSET(HDF5Particle, gas_saturation),              
    HOFFSET(HDF5Particle, liquid_fraction),
    HOFFSET(HDF5Particle, ice_fraction),
    HOFFSET(HDF5Particle, hydrate_fraction),
    HOFFSET(HDF5Particle, gas_fraction),        
    HOFFSET(HDF5Particle, viscosity), 
    HOFFSET(HDF5Particle, permeability),
    HOFFSET(HDF5Particle, liquid_permeability),
    HOFFSET(HDF5Particle, gas_permeability),                               
    HOFFSET(HDF5Particle, pressure),
    HOFFSET(HDF5Particle, pore_pressure),
    HOFFSET(HDF5Particle, pore_liquid_pressure),
    HOFFSET(HDF5Particle, pore_ice_pressure),         
    HOFFSET(HDF5Particle, liquid_pressure),
    HOFFSET(HDF5Particle, gas_pressure),
    HOFFSET(HDF5Particle, temperature),
    HOFFSET(HDF5Particle, PIC_temperature),
    HOFFSET(HDF5Particle, coord_x),
    HOFFSET(HDF5Particle, coord_y),
    HOFFSET(HDF5Particle, coord_z),
    HOFFSET(HDF5Particle, displacement_x),
    HOFFSET(HDF5Particle, displacement_y),
    HOFFSET(HDF5Particle, displacement_z),
    HOFFSET(HDF5Particle, nsize_x),
    HOFFSET(HDF5Particle, nsize_y),
    HOFFSET(HDF5Particle, nsize_z),
    HOFFSET(HDF5Particle, velocity_x),
    HOFFSET(HDF5Particle, velocity_y),
    HOFFSET(HDF5Particle, velocity_z),
    HOFFSET(HDF5Particle, liquid_velocity_x),
    HOFFSET(HDF5Particle, liquid_velocity_y),
    HOFFSET(HDF5Particle, liquid_velocity_z),
    HOFFSET(HDF5Particle, gas_velocity_x),
    HOFFSET(HDF5Particle, gas_velocity_y),
    HOFFSET(HDF5Particle, gas_velocity_z),    
    HOFFSET(HDF5Particle, stress_xx),
    HOFFSET(HDF5Particle, stress_yy),
    HOFFSET(HDF5Particle, stress_zz),
    HOFFSET(HDF5Particle, tau_xy),
    HOFFSET(HDF5Particle, tau_yz),
    HOFFSET(HDF5Particle, tau_xz),
    HOFFSET(HDF5Particle, strain_xx),
    HOFFSET(HDF5Particle, strain_yy),
    HOFFSET(HDF5Particle, strain_zz),
    HOFFSET(HDF5Particle, gamma_xy),
    HOFFSET(HDF5Particle, gamma_yz),
    HOFFSET(HDF5Particle, gamma_xz),
    HOFFSET(HDF5Particle, thermal_strain_xx),
    HOFFSET(HDF5Particle, thermal_strain_yy),
    HOFFSET(HDF5Particle, thermal_strain_zz),
    HOFFSET(HDF5Particle, thermal_gamma_xy),
    HOFFSET(HDF5Particle, thermal_gamma_yz),
    HOFFSET(HDF5Particle, thermal_gamma_xz),
    HOFFSET(HDF5Particle, liquid_strain_xx),
    HOFFSET(HDF5Particle, liquid_strain_yy),
    HOFFSET(HDF5Particle, liquid_strain_zz),
    HOFFSET(HDF5Particle, liquid_gamma_xy),
    HOFFSET(HDF5Particle, liquid_gamma_yz),
    HOFFSET(HDF5Particle, liquid_gamma_xz),
    HOFFSET(HDF5Particle, epsilon_v),
    HOFFSET(HDF5Particle, thermal_epsilon_v),      
    HOFFSET(HDF5Particle, fxx),
    HOFFSET(HDF5Particle, fxy),
    HOFFSET(HDF5Particle, fxz),
    HOFFSET(HDF5Particle, fyx),
    HOFFSET(HDF5Particle, fyy),
    HOFFSET(HDF5Particle, fyz),
    HOFFSET(HDF5Particle, fzx),
    HOFFSET(HDF5Particle, fzy),
    HOFFSET(HDF5Particle, fzz),
    HOFFSET(HDF5Particle, fabric),
    HOFFSET(HDF5Particle, rotation_xy),
    HOFFSET(HDF5Particle, rotation_yz),
    HOFFSET(HDF5Particle, rotation_xz),
    HOFFSET(HDF5Particle, svars[0]),
    HOFFSET(HDF5Particle, svars[1]),
    HOFFSET(HDF5Particle, svars[2]),
    HOFFSET(HDF5Particle, svars[3]),
    HOFFSET(HDF5Particle, svars[4]),
    HOFFSET(HDF5Particle, svars[5]),
    HOFFSET(HDF5Particle, status),
    HOFFSET(HDF5Particle, nstate_vars),

    // ===== MIXTURE: new fields =====
    HOFFSET(HDF5Particle, suction_pressure),
    HOFFSET(HDF5Particle, mixture_mass),
    HOFFSET(HDF5Particle, dSw_dpw),
    HOFFSET(HDF5Particle, total_stress_xx),
    HOFFSET(HDF5Particle, total_stress_yy),
    HOFFSET(HDF5Particle, total_stress_zz),
    HOFFSET(HDF5Particle, total_stress_tau_xy),
    HOFFSET(HDF5Particle, total_stress_tau_yz),
    HOFFSET(HDF5Particle, total_stress_tau_xz),

    // ===== LIQUID: new kinematics & rates =====
    HOFFSET(HDF5Particle, liquid_acceleration_x),
    HOFFSET(HDF5Particle, liquid_acceleration_y),
    HOFFSET(HDF5Particle, liquid_acceleration_z),
    HOFFSET(HDF5Particle, liquid_strain_gamma_xy),
    HOFFSET(HDF5Particle, liquid_strain_gamma_yz),
    HOFFSET(HDF5Particle, liquid_strain_gamma_xz),
    HOFFSET(HDF5Particle, liquid_strain_rate_xx),
    HOFFSET(HDF5Particle, liquid_strain_rate_yy),
    HOFFSET(HDF5Particle, liquid_strain_rate_zz),
    HOFFSET(HDF5Particle, liquid_strain_rate_xy),
    HOFFSET(HDF5Particle, liquid_strain_rate_yz),
    HOFFSET(HDF5Particle, liquid_strain_rate_xz),

    // ===== LIQUID: new scalars =====
    HOFFSET(HDF5Particle, effective_saturation),
    HOFFSET(HDF5Particle, liquid_chi),
    HOFFSET(HDF5Particle, liquid_volume),
    HOFFSET(HDF5Particle, liquid_mass_density),
    HOFFSET(HDF5Particle, liquid_pressure_acceleration),
    HOFFSET(HDF5Particle, liquid_volumetric_strain),
    HOFFSET(HDF5Particle, PIC_liquid_pressure),
    HOFFSET(HDF5Particle, FLIP_liquid_pressure),

    // ===== GAS: new kinematics & rates =====
    HOFFSET(HDF5Particle, gas_acceleration_x),
    HOFFSET(HDF5Particle, gas_acceleration_y),
    HOFFSET(HDF5Particle, gas_acceleration_z),
    HOFFSET(HDF5Particle, pgravity_x),
    HOFFSET(HDF5Particle, pgravity_y),
    HOFFSET(HDF5Particle, pgravity_z),
    HOFFSET(HDF5Particle, gas_strain_xx),
    HOFFSET(HDF5Particle, gas_strain_yy),
    HOFFSET(HDF5Particle, gas_strain_zz),
    HOFFSET(HDF5Particle, gas_strain_gamma_xy),
    HOFFSET(HDF5Particle, gas_strain_gamma_yz),
    HOFFSET(HDF5Particle, gas_strain_gamma_xz),
    HOFFSET(HDF5Particle, gas_strain_rate_xx),
    HOFFSET(HDF5Particle, gas_strain_rate_yy),
    HOFFSET(HDF5Particle, gas_strain_rate_zz),
    HOFFSET(HDF5Particle, gas_strain_rate_xy),
    HOFFSET(HDF5Particle, gas_strain_rate_yz),
    HOFFSET(HDF5Particle, gas_strain_rate_xz),

    // ===== GAS: new scalars =====
    HOFFSET(HDF5Particle, gas_chi),
    HOFFSET(HDF5Particle, gas_volume),
    HOFFSET(HDF5Particle, gas_mass_density),
    HOFFSET(HDF5Particle, gas_pressure_acceleration),
    HOFFSET(HDF5Particle, PIC_gas_pressure),
    HOFFSET(HDF5Particle, FLIP_gas_pressure),
    HOFFSET(HDF5Particle, gas_pressure_increment),
    HOFFSET(HDF5Particle, gas_volumetric_strain),

    // ===== WAVE: new fields =====
    HOFFSET(HDF5Particle, wave_pressure),
    HOFFSET(HDF5Particle, gamma_b),

    // ===== STATE VARIABLES: append-only checkpoint extension =====
    HOFFSET(HDF5Particle, svars[6]),
    HOFFSET(HDF5Particle, svars[7]),
    HOFFSET(HDF5Particle, svars[8]),
    HOFFSET(HDF5Particle, svars[9]),
    HOFFSET(HDF5Particle, svars[10]),
    HOFFSET(HDF5Particle, svars[11]),
    HOFFSET(HDF5Particle, svars[12]),
    HOFFSET(HDF5Particle, svars[13]),
    HOFFSET(HDF5Particle, svars[14]),
    HOFFSET(HDF5Particle, svars[15]),
    HOFFSET(HDF5Particle, svars[16]),
    HOFFSET(HDF5Particle, svars[17]),
    HOFFSET(HDF5Particle, svars[18]),
    HOFFSET(HDF5Particle, svars[19])
};

// Get size of particle
HDF5Particle particle;
const size_t dst_sizes[NFIELDS] = {
    sizeof(particle.id),
    sizeof(particle.cell_id),
    sizeof(particle.material_id),
    sizeof(particle.liquid_material_id),    
    sizeof(particle.current_time),    
    sizeof(particle.mass),
    sizeof(particle.liquid_mass),  
    sizeof(particle.ice_mass),
    sizeof(particle.hydrate_mass),  
    sizeof(particle.gas_mass),          
    sizeof(particle.volume),
    sizeof(particle.density),  
    sizeof(particle.liquid_density),  
    sizeof(particle.ice_density),
    sizeof(particle.hydrate_density),  
    sizeof(particle.gas_density),              
    sizeof(particle.porosity),  
    sizeof(particle.liquid_saturation),
    sizeof(particle.ice_saturation),
    sizeof(particle.hydrate_saturation),
    sizeof(particle.gas_saturation),   
    sizeof(particle.liquid_fraction),    
    sizeof(particle.ice_fraction),
    sizeof(particle.hydrate_fraction),
    sizeof(particle.gas_fraction),    
    sizeof(particle.viscosity), 
    sizeof(particle.permeability),
    sizeof(particle.liquid_permeability),
    sizeof(particle.gas_permeability),               
    sizeof(particle.pressure),
    sizeof(particle.pore_pressure),
    sizeof(particle.pore_liquid_pressure),
    sizeof(particle.pore_ice_pressure),
    sizeof(particle.liquid_pressure),
    sizeof(particle.gas_pressure),            
    sizeof(particle.temperature),
    sizeof(particle.PIC_temperature),      
    sizeof(particle.coord_x),
    sizeof(particle.coord_y),
    sizeof(particle.coord_z),
    sizeof(particle.displacement_x),
    sizeof(particle.displacement_y),
    sizeof(particle.displacement_z),
    sizeof(particle.nsize_x),
    sizeof(particle.nsize_y),
    sizeof(particle.nsize_z),
    sizeof(particle.velocity_x),
    sizeof(particle.velocity_y),
    sizeof(particle.velocity_z),
    sizeof(particle.liquid_velocity_x),
    sizeof(particle.liquid_velocity_y),
    sizeof(particle.liquid_velocity_z),
    sizeof(particle.gas_velocity_x),
    sizeof(particle.gas_velocity_y),
    sizeof(particle.gas_velocity_z),    
    sizeof(particle.stress_xx),
    sizeof(particle.stress_yy),
    sizeof(particle.stress_zz),
    sizeof(particle.tau_xy),
    sizeof(particle.tau_yz),
    sizeof(particle.tau_xz),
    sizeof(particle.strain_xx),
    sizeof(particle.strain_yy),
    sizeof(particle.strain_zz),
    sizeof(particle.gamma_xy),
    sizeof(particle.gamma_yz),
    sizeof(particle.gamma_xz),
    sizeof(particle.thermal_strain_xx),
    sizeof(particle.thermal_strain_yy),
    sizeof(particle.thermal_strain_zz),
    sizeof(particle.thermal_gamma_xy),
    sizeof(particle.thermal_gamma_yz),
    sizeof(particle.thermal_gamma_xz),
    sizeof(particle.liquid_strain_xx),
    sizeof(particle.liquid_strain_yy),
    sizeof(particle.liquid_strain_zz),
    sizeof(particle.liquid_gamma_xy),
    sizeof(particle.liquid_gamma_yz),
    sizeof(particle.liquid_gamma_xz),
    sizeof(particle.epsilon_v),
    sizeof(particle.thermal_epsilon_v),
    sizeof(particle.fxx),
    sizeof(particle.fxy),
    sizeof(particle.fxz),
    sizeof(particle.fyx),
    sizeof(particle.fyy),
    sizeof(particle.fyz),
    sizeof(particle.fzx),
    sizeof(particle.fzy),
    sizeof(particle.fzz),
    sizeof(particle.fabric),  
    sizeof(particle.rotation_xy),
    sizeof(particle.rotation_yz),
    sizeof(particle.rotation_xz),
    sizeof(particle.svars[0]),
    sizeof(particle.svars[1]),
    sizeof(particle.svars[2]),
    sizeof(particle.svars[3]),
    sizeof(particle.svars[4]),
    sizeof(particle.svars[5]),
    sizeof(particle.status),
    sizeof(particle.nstate_vars),
    // ===== MIXTURE: new fields =====
    sizeof(particle.suction_pressure),
    sizeof(particle.mixture_mass),
    sizeof(particle.dSw_dpw),
    sizeof(particle.total_stress_xx),
    sizeof(particle.total_stress_yy),
    sizeof(particle.total_stress_zz),
    sizeof(particle.total_stress_tau_xy),
    sizeof(particle.total_stress_tau_yz),
    sizeof(particle.total_stress_tau_xz),
    // ===== LIQUID: new kinematics & rates =====
    sizeof(particle.liquid_acceleration_x),
    sizeof(particle.liquid_acceleration_y),
    sizeof(particle.liquid_acceleration_z),
    sizeof(particle.liquid_strain_gamma_xy),
    sizeof(particle.liquid_strain_gamma_yz),
    sizeof(particle.liquid_strain_gamma_xz),
    sizeof(particle.liquid_strain_rate_xx),
    sizeof(particle.liquid_strain_rate_yy),
    sizeof(particle.liquid_strain_rate_zz),
    sizeof(particle.liquid_strain_rate_xy),
    sizeof(particle.liquid_strain_rate_yz),
    sizeof(particle.liquid_strain_rate_xz),
    // ===== LIQUID: new scalars =====
    sizeof(particle.effective_saturation),
    sizeof(particle.liquid_chi),
    sizeof(particle.liquid_volume),
    sizeof(particle.liquid_mass_density),
    sizeof(particle.liquid_pressure_acceleration),
    sizeof(particle.liquid_volumetric_strain),
    sizeof(particle.PIC_liquid_pressure),
    sizeof(particle.FLIP_liquid_pressure),
    // ===== GAS: new kinematics & rates =====
    sizeof(particle.gas_acceleration_x),
    sizeof(particle.gas_acceleration_y),
    sizeof(particle.gas_acceleration_z),
    sizeof(particle.pgravity_x),
    sizeof(particle.pgravity_y),
    sizeof(particle.pgravity_z),
    sizeof(particle.gas_strain_xx),
    sizeof(particle.gas_strain_yy),
    sizeof(particle.gas_strain_zz),
    sizeof(particle.gas_strain_gamma_xy),
    sizeof(particle.gas_strain_gamma_yz),
    sizeof(particle.gas_strain_gamma_xz),
    sizeof(particle.gas_strain_rate_xx),
    sizeof(particle.gas_strain_rate_yy),
    sizeof(particle.gas_strain_rate_zz),
    sizeof(particle.gas_strain_rate_xy),
    sizeof(particle.gas_strain_rate_yz),
    sizeof(particle.gas_strain_rate_xz),
    // ===== GAS: new scalars =====
    sizeof(particle.gas_chi),
    sizeof(particle.gas_volume),
    sizeof(particle.gas_mass_density),
    sizeof(particle.gas_pressure_acceleration),
    sizeof(particle.PIC_gas_pressure),
    sizeof(particle.FLIP_gas_pressure),
    sizeof(particle.gas_pressure_increment),
    sizeof(particle.gas_volumetric_strain),
    // ===== WAVE: new fields =====
    sizeof(particle.wave_pressure),
    sizeof(particle.gamma_b),

    // ===== STATE VARIABLES: append-only checkpoint extension =====
    sizeof(particle.svars[6]),
    sizeof(particle.svars[7]),
    sizeof(particle.svars[8]),
    sizeof(particle.svars[9]),
    sizeof(particle.svars[10]),
    sizeof(particle.svars[11]),
    sizeof(particle.svars[12]),
    sizeof(particle.svars[13]),
    sizeof(particle.svars[14]),
    sizeof(particle.svars[15]),
    sizeof(particle.svars[16]),
    sizeof(particle.svars[17]),
    sizeof(particle.svars[18]),
    sizeof(particle.svars[19])
};

// Define particle field information
const char* field_names[NFIELDS] = {
    "id",
    "cell_id",
    "material_id",
    "liquid_material_id",    
    "current_time",
    "mass",
    "liquid_mass",
    "ice_mass",
    "hydrate_mass",
    "gas_mass",        
    "volume",
    "density",
    "liquid_density",
    "ice_density",
    "hydrate_density",
    "gas_density",    
    "porosity",
    "liquid_saturation",
    "ice_saturation",
    "hydrate_saturation",
    "gas_saturation",            
    "liquid_fraction",   
    "ice_fraction",
    "hydrate_fraction",   
    "gas_fraction",     
    "viscosity",   
    "permeability",
    "liquid_permeability",
    "gas_permeability",                       
    "pressure",
    "pore_pressure",
    "pore_liquid_pressure",   
    "pore_ice_pressure",
    "liquid_pressure",
    "gas_pressure",                          
    "temperature",
    "PIC_temperature",    
    "coord_x",
    "coord_y",
    "coord_z",
    "displacement_x",
    "displacement_y",
    "displacement_z",
    "nsize_x",
    "nsize_y",
    "nsize_z",
    "velocity_x",
    "velocity_y",
    "velocity_z",
    "liquid_velocity_x",
    "liquid_velocity_y",
    "liquid_velocity_z",
    "gas_velocity_x",
    "gas_velocity_y",
    "gas_velocity_z",    
    "stress_xx",
    "stress_yy",
    "stress_zz",
    "tau_xy",
    "tau_yz",
    "tau_xz",
    "strain_xx",
    "strain_yy",
    "strain_zz",
    "gamma_xy",
    "gamma_yz",
    "gamma_xz",
    "thermal_strain_xx",
    "thermal_strain_yy",
    "thermal_strain_zz",
    "thermal_gamma_xy",
    "thermal_gamma_yz",
    "thermal_gamma_xz",
    "liquid_strain_xx",
    "liquid_strain_yy",
    "liquid_strain_zz",
    "liquid_gamma_xy",
    "liquid_gamma_yz",
    "liquid_gamma_xz",
    "epsilon_v",
    "thermal_epsilon_v",
    "fxx",
    "fxy",
    "fxz",
    "fyx",
    "fyy",
    "fyz",
    "fzx",
    "fzy",
    "fzz",
    "fabric",    
    "rotation_xy",
    "rotation_yz",
    "rotation_xz",
    "svars_0",
    "svars_1",
    "svars_2",
    "svars_3",
    "svars_4",
    "svars_5",
    "status",
    "nstate_vars",
    // ===== MIXTURE: new fields =====
    "suction_pressure",
    "mixture_mass",
    "dSw_dpw",
    "total_stress_xx",
    "total_stress_yy",
    "total_stress_zz",
    "total_stress_tau_xy",
    "total_stress_tau_yz",
    "total_stress_tau_xz",
    // ===== LIQUID: new kinematics & rates =====
    "liquid_acceleration_x",
    "liquid_acceleration_y",
    "liquid_acceleration_z",
    "liquid_strain_gamma_xy",
    "liquid_strain_gamma_yz",
    "liquid_strain_gamma_xz",
    "liquid_strain_rate_xx",
    "liquid_strain_rate_yy",
    "liquid_strain_rate_zz",
    "liquid_strain_rate_xy",
    "liquid_strain_rate_yz",
    "liquid_strain_rate_xz",
    // ===== LIQUID: new scalars =====
    "effective_saturation",
    "liquid_chi",
    "liquid_volume",
    "liquid_mass_density",
    "liquid_pressure_acceleration",
    "liquid_volumetric_strain",
    "PIC_liquid_pressure",
    "FLIP_liquid_pressure",
    // ===== GAS: new kinematics & rates =====
    "gas_acceleration_x",
    "gas_acceleration_y",
    "gas_acceleration_z",
    "pgravity_x",
    "pgravity_y",
    "pgravity_z",
    "gas_strain_xx",
    "gas_strain_yy",
    "gas_strain_zz",
    "gas_strain_gamma_xy",
    "gas_strain_gamma_yz",
    "gas_strain_gamma_xz",
    "gas_strain_rate_xx",
    "gas_strain_rate_yy",
    "gas_strain_rate_zz",
    "gas_strain_rate_xy",
    "gas_strain_rate_yz",
    "gas_strain_rate_xz",
    // ===== GAS: new scalars =====
    "gas_chi",
    "gas_volume",
    "gas_mass_density",
    "gas_pressure_acceleration",
    "PIC_gas_pressure",
    "FLIP_gas_pressure",
    "gas_pressure_increment",
    "gas_volumetric_strain",
    // ===== WAVE: new fields =====
    "wave_pressure",
    "gamma_b",

    // ===== STATE VARIABLES: append-only checkpoint extension =====
    "svars_6",
    "svars_7",
    "svars_8",
    "svars_9",
    "svars_10",
    "svars_11",
    "svars_12",
    "svars_13",
    "svars_14",
    "svars_15",
    "svars_16",
    "svars_17",
    "svars_18",
    "svars_19"
    };

// Initialise field types
const hid_t field_type[NFIELDS] = {
    H5T_NATIVE_LLONG, 
    H5T_NATIVE_LLONG,           
    H5T_NATIVE_UINT,            
    H5T_NATIVE_UINT,     
    H5T_NATIVE_DOUBLE,              
    H5T_NATIVE_DOUBLE, 
    H5T_NATIVE_DOUBLE,             
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,              
    H5T_NATIVE_DOUBLE, 
    H5T_NATIVE_DOUBLE,             
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,    
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,                     
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,         
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE, 
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,          
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE,         
    H5T_NATIVE_DOUBLE,     
    H5T_NATIVE_HBOOL,                      
    H5T_NATIVE_UINT,
    // ===== MIXTURE: new fields (9) =====
    H5T_NATIVE_DOUBLE, // suction_pressure
    H5T_NATIVE_DOUBLE, // mixture_mass
    H5T_NATIVE_DOUBLE, // dSw_dpw
    H5T_NATIVE_DOUBLE, // total_stress_xx
    H5T_NATIVE_DOUBLE, // total_stress_yy
    H5T_NATIVE_DOUBLE, // total_stress_zz
    H5T_NATIVE_DOUBLE, // total_stress_tau_xy
    H5T_NATIVE_DOUBLE, // total_stress_tau_yz
    H5T_NATIVE_DOUBLE, // total_stress_tau_xz
    // ===== LIQUID: new kinematics & rates (12) =====
    H5T_NATIVE_DOUBLE, // liquid_acceleration_x
    H5T_NATIVE_DOUBLE, // liquid_acceleration_y
    H5T_NATIVE_DOUBLE, // liquid_acceleration_z
    H5T_NATIVE_DOUBLE, // liquid_strain_gamma_xy
    H5T_NATIVE_DOUBLE, // liquid_strain_gamma_yz
    H5T_NATIVE_DOUBLE, // liquid_strain_gamma_xz
    H5T_NATIVE_DOUBLE, // liquid_strain_rate_xx
    H5T_NATIVE_DOUBLE, // liquid_strain_rate_yy
    H5T_NATIVE_DOUBLE, // liquid_strain_rate_zz
    H5T_NATIVE_DOUBLE, // liquid_strain_rate_xy
    H5T_NATIVE_DOUBLE, // liquid_strain_rate_yz
    H5T_NATIVE_DOUBLE, // liquid_strain_rate_xz
    // ===== LIQUID: new scalars (7) =====
    H5T_NATIVE_DOUBLE, // effective_saturation
    H5T_NATIVE_DOUBLE, // liquid_chi
    H5T_NATIVE_DOUBLE, // liquid_volume
    H5T_NATIVE_DOUBLE, // liquid_mass_density
    H5T_NATIVE_DOUBLE, // liquid_pressure_acceleration
    H5T_NATIVE_DOUBLE, // liquid_volumetric_strain
    H5T_NATIVE_DOUBLE, // PIC_liquid_pressure
    H5T_NATIVE_DOUBLE, // FLIP_liquid_pressure
    // ===== GAS: new kinematics & rates (18) =====
    H5T_NATIVE_DOUBLE, // gas_acceleration_x
    H5T_NATIVE_DOUBLE, // gas_acceleration_y
    H5T_NATIVE_DOUBLE, // gas_acceleration_z
    H5T_NATIVE_DOUBLE, // pgravity_x
    H5T_NATIVE_DOUBLE, // pgravity_y
    H5T_NATIVE_DOUBLE, // pgravity_z
    H5T_NATIVE_DOUBLE, // gas_strain_xx
    H5T_NATIVE_DOUBLE, // gas_strain_yy
    H5T_NATIVE_DOUBLE, // gas_strain_zz
    H5T_NATIVE_DOUBLE, // gas_strain_gamma_xy
    H5T_NATIVE_DOUBLE, // gas_strain_gamma_yz
    H5T_NATIVE_DOUBLE, // gas_strain_gamma_xz
    H5T_NATIVE_DOUBLE, // gas_strain_rate_xx
    H5T_NATIVE_DOUBLE, // gas_strain_rate_yy
    H5T_NATIVE_DOUBLE, // gas_strain_rate_zz
    H5T_NATIVE_DOUBLE, // gas_strain_rate_xy
    H5T_NATIVE_DOUBLE, // gas_strain_rate_yz
    H5T_NATIVE_DOUBLE, // gas_strain_rate_xz
    // ===== GAS: new scalars (7) =====
    H5T_NATIVE_DOUBLE, // gas_chi
    H5T_NATIVE_DOUBLE, // gas_volume
    H5T_NATIVE_DOUBLE, // gas_mass_density
    H5T_NATIVE_DOUBLE, // gas_pressure_acceleration
    H5T_NATIVE_DOUBLE, // PIC_gas_pressure
    H5T_NATIVE_DOUBLE, // FLIP_gas_pressure
    H5T_NATIVE_DOUBLE, // gas_pressure_increment
    H5T_NATIVE_DOUBLE, // gas_volumetric_strain
    // ===== WAVE: new fields (2) =====
    H5T_NATIVE_HBOOL, // wave_pressure
    H5T_NATIVE_DOUBLE, // gamma_b

    // ===== STATE VARIABLES: append-only checkpoint extension (14) =====
    H5T_NATIVE_DOUBLE, // svars_6
    H5T_NATIVE_DOUBLE, // svars_7
    H5T_NATIVE_DOUBLE, // svars_8
    H5T_NATIVE_DOUBLE, // svars_9
    H5T_NATIVE_DOUBLE, // svars_10
    H5T_NATIVE_DOUBLE, // svars_11
    H5T_NATIVE_DOUBLE, // svars_12
    H5T_NATIVE_DOUBLE, // svars_13
    H5T_NATIVE_DOUBLE, // svars_14
    H5T_NATIVE_DOUBLE, // svars_15
    H5T_NATIVE_DOUBLE, // svars_16
    H5T_NATIVE_DOUBLE, // svars_17
    H5T_NATIVE_DOUBLE, // svars_18
    H5T_NATIVE_DOUBLE, // svars_19
    };         

void validate_table_metadata(hsize_t field_count, hsize_t record_count,
                             hsize_t expected_records) {
  if (field_count != LEGACY_NFIELDS && field_count != NFIELDS)
    throw std::runtime_error(
        "Unsupported HDF5 particle schema: field_count=" +
        std::to_string(field_count) + ", expected " +
        std::to_string(LEGACY_NFIELDS) + " (legacy) or " +
        std::to_string(NFIELDS));
  if (record_count != expected_records)
    throw std::runtime_error(
        "HDF5 particle record count mismatch: file=" +
        std::to_string(record_count) +
        ", mesh=" + std::to_string(expected_records));
}

void validate_legacy_state_variables(hsize_t field_count,
                                     const HDF5Particle* particles,
                                     std::size_t particle_count) {
  if (field_count == NFIELDS) return;
  if (field_count != LEGACY_NFIELDS)
    throw std::runtime_error(
        "Cannot validate state variables for an unsupported HDF5 particle "
        "schema");
  if (particles == nullptr && particle_count != 0)
    throw std::invalid_argument("HDF5 particle buffer is null");

  for (std::size_t i = 0; i < particle_count; ++i) {
    if (particles[i].nstate_vars > LEGACY_NSTATE_VARS)
      throw std::runtime_error(
          "Legacy HDF5 particle schema stores only svars_0..svars_5, but "
          "particle " +
          std::to_string(particles[i].id) + " declares " +
          std::to_string(particles[i].nstate_vars) +
          " state variables; refusing a lossy restart");
  }
}
}  // namespace hdf5::particle
}  // namespace mpm





  
