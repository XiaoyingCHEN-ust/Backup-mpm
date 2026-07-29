#ifndef MPM_PARTICLEBASE_H_
#define MPM_PARTICLEBASE_H_

#include <array>
#include <limits>
#include <memory>
#include <vector>

#include "cell.h"
#include "data_types.h"
#include "function_base.h"
#include "hdf5_particle.h"
#include "material.h"

namespace mpm {

// Forward declaration of Material
template <unsigned Tdim>
class Material;

// Particle phases
enum ParticlePhase : unsigned int {
  Mixture = 0,
  Solid = 0,
  Hydrate = 0,
  Ice = 0,
  Liquid = 1,
  Gas = 2
};

//! Force and diagnostic data returned by a circular rigid-contact operation.
template <unsigned Tdim>
struct RigidCircleContactResult {
  Eigen::Matrix<double, Tdim, 1> force =
      Eigen::Matrix<double, Tdim, 1>::Zero();
  double moment{0.};
  double maximum_penetration{0.};
  unsigned contacts{0};
};

// ParticleBase class
template <unsigned Tdim>
class ParticleBase {

public:
  // Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //============================================================================
  // CONSTRUCT AND DESTRUCT A PARTICLE

  // Constructor with id and coordinates
  ParticleBase(Index id, const VectorDim& coord);

  // Constructor with id, coordinates and status
  ParticleBase(Index id, const VectorDim& coord, bool status);

  // Destructor
  virtual ~ParticleBase(){};

  // Delete copy constructor
  ParticleBase(const ParticleBase<Tdim>&) = delete;

  // Delete assignement operator
  ParticleBase& operator=(const ParticleBase<Tdim>&) = delete;

  //============================================================================
  // ASSIGN INITIAL CONDITIONS

  // Initialise particle HDF5 data
  virtual bool initialise_particle(const HDF5Particle& particle) = 0;

  // Initialise particle HDF5 data and material
  virtual bool initialise_particle(const HDF5Particle& particle,
      const std::shared_ptr<Material<Tdim>>& material) = 0;

  // Initialise properties
  virtual void initialise() = 0;

  // Retrun particle data as HDF5
  virtual HDF5Particle hdf5() = 0;

  // Assign material
  virtual bool assign_material(
      const std::shared_ptr<Material<Tdim>>& material) = 0;

  // Assign affine particle-grid velocity transfer
  virtual void assign_affine_mpm(bool affine_mpm) { (void) affine_mpm; };

  // Assign volume
  virtual bool assign_initial_volume(double volume) = 0;

  // Compute volume of particle
  virtual void compute_volume(bool is_axisymmetric) noexcept = 0;

  // Compute mass of particle
  virtual void compute_mass() = 0;

  // Initial stress
  virtual void assign_initial_stress(const Eigen::Matrix<double, 6, 1>& stress) = 0;

  // Initial temperature
  virtual void assign_initial_temperature(double temperature) = 0;

  // Initialise liquid phase
  virtual void initialise_liquid_phase() {
    throw std::runtime_error(
        "Calling the base class function (initialise_liquid_phase) in "
        "ParticleBase:: illegal operation!");
  };

  // Initialise liquid phase
  virtual void initialise_liquid_gas_phases() {
    throw std::runtime_error(
        "Calling the base class function (initialise_liquid_gas_phases) in "
        "ParticleBase:: illegal operation!");
  };

  // Assign material
  virtual bool assign_liquid_material(
      const std::shared_ptr<Material<Tdim>>& material) {
    throw std::runtime_error(
        "Calling the base class function (assign_liquid_material) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Initial pore pressure
  virtual void initial_pore_pressure(double pore_pressure) {
    throw std::runtime_error(
        "Calling the base class function "
        "(initial_pore_pressure) in "
        "ParticleBase:: illegal operation!");
  };

  // Initialise particle pore pressure by watertable
  virtual bool initialise_pore_pressure_watertable(
                    const unsigned dir_v, const unsigned dir_h,
                    std::map<double, double>& refernece_points) {
    throw std::runtime_error(
        "Calling the base class function "
        "(initial_pore_pressure_watertable) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  // Assign saturation degree
  virtual bool assign_liquid_saturation_degree() {
    throw std::runtime_error(
        "Calling the base class function ( assign_liquid_saturation_degree) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Output particle results
  virtual void output_results(double step_) {
    throw std::runtime_error(
        "Calling the base class function ( assign_liquid_saturation_degree) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Assign externally prescribed liquid and gas pressures
  virtual bool assign_prescribed_phase_pressures(double liquid_pressure,
                                                 double gas_pressure) {
    return false;
  };

  //============================================================================
  // APPLY BOUNDARY CONDITIONS

  // Apply particle velocity constraints
  virtual void apply_particle_velocity_constraints(unsigned dir,
                                                  double velocity) = 0;

  // Apply particle temperature constraints
  virtual void apply_particle_temperature_constraints(double temperature) = 0;

  // Assign initial free surface particle
  virtual void assign_particle_free_surface(bool initial_free_surface) = 0;

  // Assign initial nonfree surface particle
  virtual void assign_particle_nonfree_surface(bool initial_nonfree_surface) = 0;

  // Assign traction
  virtual bool assign_particle_traction(unsigned direction, double traction) = 0;

  // Assign heat source
  virtual bool assign_particle_heat_source(double heat_source, double dt) = 0;

  // Compute free surface
  virtual bool compute_particle_free_surface() = 0;

  // Assign traction
  virtual bool assign_particle_contact(unsigned dir, double normal) = 0; 

  // Assign free surface particle manually
  virtual bool assign_particle_free_surfaces() = 0;

  // Assign particle liquid phase velocity constraints
  virtual bool assign_particle_liquid_velocity_constraint(unsigned dir,
                                                          double velocity) {
    throw std::runtime_error(
        "Calling the base class function "
        "(assign_particle_liquid_velocity_constraint) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Apply particle liquid phase velocity constraints
  virtual void apply_particle_liquid_velocity_constraints() {
    throw std::runtime_error(
        "Calling the base class function "
        "(apply_particle_liquid_velocity_constraints) in "
        "ParticleBase:: illegal operation!");
  };

  // Assign particle pressure constraints
  virtual bool assign_particle_pore_pressure_constraint(double pressure) {
    throw std::runtime_error(
        "Calling the base class function "
        "(assign_particle_pore_pressure_constraint) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Apply particle pore pressure constraints
  virtual void apply_particle_pore_pressure_constraints(double pore_pressure)  {
    throw std::runtime_error(
        "Calling the base class function (apply_particle_pore_pressure_constraints) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  //==========================================================================
  // LOCATE PARTICLES & COMPUTE SHAPE FUNCTIONS

  // Assign coordinates
  void assign_coordinates(const VectorDim& coord) { coordinates_ = coord; }

  // Assign cell
  virtual bool assign_cell(const std::shared_ptr<Cell<Tdim>>& cellptr) = 0;

  // Assign cell and xi
  virtual bool assign_cell_xi(const std::shared_ptr<Cell<Tdim>>& cellptr,
                              const Eigen::Matrix<double, Tdim, 1>& xi) = 0;

  // Assign cell id
  virtual bool assign_cell_id(Index id) = 0;

  // Remove cell
  virtual void remove_cell() = 0;

  // Compute reference coordinates in a cell
  virtual bool compute_reference_location() = 0;

  // Compute shape functions
  virtual void compute_shapefn() noexcept = 0;

  //==========================================================================
  // MAP PARTICLE INFORMATION TO NODES

  // Map particle mass and momentum to nodes
  virtual void map_mass_momentum_to_nodes() noexcept = 0;

  // Map particle mass and pressure to nodes
  virtual void map_mass_pressure_to_nodes() {
    throw std::runtime_error(
        "Calling the base class function (map_mass_pressure_to_nodes) in "
        "ParticleBase:: illegal operation!");
  };

  // Map particle PIC pressure to nodes
  virtual void map_pic_mass_pressure_to_nodes() {
    throw std::runtime_error(
        "Calling the base class function (map_pic_mass_pressure_to_nodes) in "
        "ParticleBase:: illegal operation!");
  };

  // Compute particle pressure gradient from the current nodal pressure field
  // and cache it on the particle for output.
  virtual void compute_pressure_gradient_to_cache(unsigned phase) {
    throw std::runtime_error(
        "Calling the base class function "
        "(compute_pressure_gradient_to_cache) in "
        "ParticleBase:: illegal operation!");
  };

  // Map cached particle pressure gradients to nodes for output smoothing.
  virtual void map_pressure_gradient_to_nodes() {
    throw std::runtime_error(
        "Calling the base class function (map_pressure_gradient_to_nodes) in "
        "ParticleBase:: illegal operation!");
  };

  // Map particle scaler propertiesto nodes
  virtual void map_scalers_to_nodes() {
    throw std::runtime_error(
        "Calling the base class function (map_scalers_to_nodes) in "
        "ParticleBase:: illegal operation!");
  };

  // Map particle density to nodes
  virtual void map_density_to_nodes() {
    throw std::runtime_error(
        "Calling the base class function (map_density_to_nodes) in "
        "ParticleBase:: illegal operation!");
  };  

  // Map body force
  virtual void map_external_force(const VectorDim& pgravity) = 0;

  //! Map soil contact force from a circular rigid pipe to mixture nodes.
  //! Non-supporting particle types return an inactive contact by default.
  virtual RigidCircleContactResult<Tdim> map_rigid_circle_contact_force(
      const VectorDim&, const VectorDim&, double, double, double, double,
      double, double, double) {
    return RigidCircleContactResult<Tdim>();
  }

  // Map internal force
  virtual void map_internal_force() = 0;

  // Map traction force
  virtual void map_traction_force() noexcept = 0;

  // Map particle pressure to nodes
  virtual bool map_pressure_to_nodes(double current_time = 0.) noexcept = 0;

  // Map particle heat capacity and heat to nodes
  virtual void map_heat_to_nodes() = 0;

  // Map heat conduction
  virtual void map_heat_conduction() = 0;

   // Map heat conduction
  virtual void map_virtual_heat_flux(bool convective, 
                                    const double para_1,
                                    const double para_2) = 0; 

  // Map hydraulic convection of mixture
  virtual void map_hydraulic_conduction() {
    throw std::runtime_error(
        "Calling the base class function (map_hydraulic_conduction) in "
        "ParticleBase:: illegal operation!");
  };

  // Map heat source force
  virtual void map_heat_source() = 0;

  // Map plastic work
  virtual void map_plastic_work(double dt) noexcept = 0;

  // Map rigid particle velocity to nodes
  virtual void map_moving_rigid_velocity_to_nodes(unsigned dir, 
                                  double velocity, double dt) noexcept = 0;

  // map moving rigid velocity and surface normal to nodes
  virtual void map_rigid_mass_momentum_to_nodes() noexcept = 0; 

  // Map internal force
  virtual void map_internal_force_semi(double beta){
    throw std::runtime_error(
        "Calling the base class function (map_internal_force_semi) in "
        "ParticleBase:: illegal operation!");};

  // Map drag force coefficient
  virtual void map_drag_force_coefficient() {
    throw std::runtime_error(
        "Calling the base class function (map_drag_force_coefficient) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Map heat convection of mixture
  virtual void map_heat_convection() {
    throw std::runtime_error(
        "Calling the base class function (map_heat_convection) in "
        "ParticleBase:: illegal operation!");
  };

  // Map mass convection of mixture
  virtual void map_mass_convection() {
    throw std::runtime_error(
        "Calling the base class function (map_mass_convection) in "
        "ParticleBase:: illegal operation!");
  };

  // Map latent heat of ice phase
  virtual void map_latent_heat(double dt)  {
    throw std::runtime_error(
        "Calling the base class function (map_latent_heat) in "
        "ParticleBase:: illegal operation!");
  };

  // Map heat convection of mixture
  virtual void map_covective_heat_flux(double current_time) {
    throw std::runtime_error(
        "Calling the base class function (map_covective_heat_flux) in "
        "ParticleBase:: illegal operation!");
  };

  // Assign pore pressure to nodes
  virtual bool map_pore_pressure_to_nodes(double current_time = 0.) {
    throw std::runtime_error(
        "Calling the base class function (map_pore_pressure_to_nodes) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Assign pore pressure to nodes
  virtual bool map_gas_saturation_to_nodes() {
    throw std::runtime_error(
        "Calling the base class function (map_gas_saturation_to_nodes) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Map volumetric strain
  virtual void map_volumetric_strain() {
    throw std::runtime_error(
        "Calling the base class function (map_volumetric_strain) in "
        "ParticleBase:: illegal operation!");
  };

  //------------------------------------------------------------
  // Implict mpm
  // Map heat laplacian materix KTT to cell
  virtual bool map_KTT_to_cell() {
    throw std::runtime_error(
        "Calling the base class function (map_KTT_to_cell) in ParticleBase:: "
        "illegal operation!");
    return 0;
  };

  // Map heat laplacian materix KTT to cell
  virtual bool map_MTT_to_cell() {
    throw std::runtime_error(
        "Calling the base class function (map_MTT_to_cell) in ParticleBase:: "
        "illegal operation!");
    return 0;
  };

  //------------------------------------------------------------
  // Semi-implict mpm
  // Map K inter element matrix to cell
  virtual bool map_K_inter_to_cell() {
    throw std::runtime_error(
        "Calling the base class function (map_K_inter_to_cell) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Map laplacian element matrix to cell
  virtual bool map_L_to_cell(double dt, double alpha) {
    throw std::runtime_error(
        "Calling the base class function (map_L_to_cell) in ParticleBase:: "
        "illegal operation!");
    return 0;
  };

  // Map F element matrix to cell (used in poisson equation RHS)
  virtual bool map_F_to_cell() {
    throw std::runtime_error(
        "Calling the base class function (map_F_to_cell) in ParticleBase:: "
        "illegal operation!");
    return 0;
  };  

  // Map T element matrix to cell (used in poisson equation RHS)
  virtual bool map_T_to_cell() {
    throw std::runtime_error(
        "Calling the base class function (map_T_to_cell) in ParticleBase:: "
        "illegal operation!");
    return 0;
  }; 

  // Map T element matrix to cell (used in poisson equation RHS)
  virtual bool map_P_to_cell(double beta) {
    throw std::runtime_error(
        "Calling the base class function (map_P_to_cell) in ParticleBase:: "
        "illegal operation!");
    return 0;
  }; 


  // Map K_cor element matrix to cell
  virtual bool map_K_cor_to_cell(double dt, double alpha) {
    throw std::runtime_error(
        "Calling the base class function (map_K_cor_to_cell) in ParticleBase:: "
        "illegal operation!");
    return 0;
  };

  //==========================================================================
  // UPDATE PARTICLE INFORMATION

  // Compute updated position
  virtual void compute_updated_velocity(double dt, 
                                        double pic = 0,
                                        double damping_factor = 0) = 0;

  // Compute updated temperature and temperature changing rate 
  virtual void update_particle_temperature(double dt, double pic_t) noexcept = 0;

  // compute updated pore pressure
  virtual void update_particle_pore_pressure(double dt, double pic_t) noexcept {
    throw std::runtime_error(
        "Calling the base class function (update_particle_pore_pressure) in "
        "ParticleBase:: illegal operation!");
  };

  // Compute strain
  virtual void update_particle_strain(double dt) noexcept = 0;

  // Compute thermal strain of the particle
  virtual void update_particle_thermal_strain() noexcept  = 0;

  // Compute stress
  virtual void update_particle_stress() noexcept = 0;

  // Update volume based on centre volumetric strain rate
  virtual void update_particle_volume()  = 0;

  // Update porosity
  virtual bool update_particle_porosity(double dt) = 0;

  // update material density of particle (considering thermal expansion)
  virtual void update_particle_density(double dt) = 0;

  // Compute pressure smoothing of the particle based on nodal pressure
  virtual bool compute_pressure_smoothing() noexcept = 0;

  // Compute pore pressure
  virtual void compute_pore_pressure(double dt) {
    throw std::runtime_error(
        "Calling the base class function (compute_pore_pressure) in "
        "ParticleBase:: illegal operation!");
  };

  // Compute pore pressure
  virtual double compute_relative_permeability(
    double Spi, double Srpi, double Srw, double Srg, double m_pi) {
    throw std::runtime_error(
        "Calling the base class function (compute_relative_permeability) in "
        "ParticleBase:: illegal operation!");
  };

  // Compute pre-compute pf for wave
  virtual bool build_wave_pf_context() {
    throw std::runtime_error(
        "Calling the base class function (compute_pre-compute pf for wave) in "
        "ParticleBase:: illegal operation!");
  };

  // Compute pore pressure (thermal)
  virtual void compute_thermal_pore_pressure (double dt) {
    throw std::runtime_error(
        "Calling the base class function (compute_thermal_pore_pressure) in "
        "ParticleBase:: illegal operation!");
  };

  // Compute updated pore pressure
  virtual bool compute_updated_pore_pressure(double beta) {
    throw std::runtime_error(
        "Calling the base class function (compute_updated_pore_pressure) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Compute pore pressure somoothening by interpolating nodal pressure
  virtual bool compute_pore_pressure_smoothing() {
    throw std::runtime_error(
        "Calling the base class function (compute_pore_pressure_smoothing) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Compute pore pressure smoothing for PIC output without state write-back
  virtual bool compute_pore_pressure_smoothing_for_pic() {
    throw std::runtime_error(
        "Calling the base class function "
        "(compute_pore_pressure_smoothing_for_pic) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Compute scalers somoothening by interpolating nodal values
  virtual bool compute_scalers_smoothing() {
    throw std::runtime_error(
        "Calling the base class function (compute_scalers_smoothing) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Smooth cached pressure gradients by interpolating nodal averaged values.
  virtual bool compute_pressure_gradient_smoothing() {
    throw std::runtime_error(
        "Calling the base class function "
        "(compute_pressure_gradient_smoothing) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Compute pore pressure somoothening by interpolating nodal pressure
  virtual bool compute_gas_saturation_smoothing() {
    throw std::runtime_error(
        "Calling the base class function (compute_gas_saturation_smoothing) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Return gas saturation
  virtual void update_gas_saturation(double dt) {
    throw std::runtime_error(
        "Calling the base class function (gas_saturation) in "
        "ParticleBase:: illegal operation!");
  };

  // Return liquid water saturation
  virtual void update_liquid_saturation(double dt) {
    throw std::runtime_error(
        "Calling the base class function (liquid_saturation) in "
        "ParticleBase:: illegal operation!");
  };

  // Return liquid water saturation
  virtual void update_liquid_gas_saturation(double dt) {
    throw std::runtime_error(
        "Calling the base class function (liquid_saturation) in "
        "ParticleBase:: illegal operation!");
  };

  // Update particle permeability
  virtual void update_permeability() {
    throw std::runtime_error(
        "Calling the base class function (update_permeability) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Compute viscosity
  virtual void update_viscosity() noexcept {
    throw std::runtime_error(
        "Calling the base class function (update_viscosity) in "
        "ParticleBase:: illegal operation!");
  };

  // Compute pore liquid pressure and pore ice pressure at particles
  virtual bool compute_liquid_ice_pore_pressure(double beta) {
    throw std::runtime_error(
        "Calling the base class function (compute_liquid_ice_pore_pressure) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

    // Compute thermal strain
  virtual void compute_ice_thermal_strain() noexcept {
    throw std::runtime_error(
        "Calling the base class function (compute_ice_thermal_strain) in "
        "ParticleBase:: illegal operation!");
  };

    // Compute frost heave
  virtual void compute_frost_heave_strain() noexcept {
    throw std::runtime_error(
        "Calling the base class function (compute_frost_heave_strain) in "
        "ParticleBase:: illegal operation!");
  };

    // Compute source term
  virtual void update_source_term() {
    throw std::runtime_error(
        "Calling the base class function (update_source_term) in "
        "ParticleBase:: illegal operation!");
  };

    // Compute hydrate saturation
  virtual void update_hydrate_saturation(double dt) {
    throw std::runtime_error(
        "Calling the base class function (update_hydrate_saturation) in "
        "ParticleBase:: illegal operation!");
  };  


  //==========================================================================
  // MULTISCALE FUNCTIONS

  // Compute displacement gradient
  virtual void compute_displacement_gradient(double dt, bool thermal) noexcept = 0;

  // Set stress
  virtual void set_stress(const Eigen::MatrixXd& stresses,
                          bool increment) noexcept = 0;

  // Set porosity
  virtual void set_porosity(const Eigen::MatrixXd& porosities) noexcept = 0;

  // Set fabric
  virtual void set_fabric(std::string fabric_type,
                          const Eigen::MatrixXd& fabrics) = 0;

  // Set rotation
  virtual void set_rotation(const Eigen::MatrixXd& rotations) noexcept = 0;

  // Return fabric
  virtual Eigen::Matrix<double, Tdim, Tdim> fabric(
      std::string fabric_type) const = 0;

  //============================================================================
  // RETURN PARTICLE DATA

  // Return id of the particleBase
  Index id() const { return id_; }

  // Return material id
  unsigned material_id() const { return material_id_; }

  // Status
  bool status() const { return status_; }

  // Return coordinates
  VectorDim coordinates() const { return coordinates_; }

  // Return cell id
  virtual Index cell_id() const = 0;

  // Return cell ptr status
  virtual bool cell_ptr() const = 0;

  // Return vector data of particles
  virtual Eigen::VectorXd vector_data(const std::string& property) = 0;

  // Return scalar data of particles
  virtual double scalar_data(const std::string& property) = 0;

  // Return a state variable
  virtual double state_variable(const std::string& var) const = 0;

  // Return pressure
  virtual double pressure() const = 0;

  // Return volume
  virtual double volume() const = 0;

  // dvolumetric strain
  virtual double dvolumetric_strain() const = 0;

  // Return pic temperature of the particle
  virtual double temperature() const = 0;

  // Return temperature increment of the particle
  virtual double temperature_increment() const = 0;

  // Return displacement of the particle
  virtual VectorDim displacement() const = 0;

  // Strain rate
  virtual Eigen::Matrix<double, 6, 1> strain_rate() const = 0;

  // Return displacement gradient of particle
  virtual Eigen::Matrix<double, Tdim, Tdim> displacement_gradient() const = 0;

  // Return initital free surface particle
  virtual bool initial_free_surface() = 0;

  // Return free surface particle
  virtual bool free_surface() = 0;

  // Return liquid pore pressure
  virtual double pore_pressure() const {
    throw std::runtime_error(
        "Calling the base class function (pore_pressure) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Return liquid ice pressure
  virtual double pore_ice_pressure() const {
    throw std::runtime_error(
        "Calling the base class function (pore_ice_pressure) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // return pore pressure increment
  virtual double pore_pressure_increments(){
    throw std::runtime_error(
        "Calling the base class function (pore_pressure_increments) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Return liquid saturation
  virtual double liquid_saturation() const {
    throw std::runtime_error(
        "Calling the base class function (liquid_saturation) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Return liquid saturation
  virtual double hydrate_saturation() const {
    throw std::runtime_error(
        "Calling the base class function (hydrate_saturation) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Return liquid saturation
  virtual double ini_hydrate_saturation() const {
    throw std::runtime_error(
        "Calling the base class function (ini_hydrate_saturation) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  // Return liquid saturation
  virtual double ice_saturation() const {
    throw std::runtime_error(
        "Calling the base class function (ice_saturation) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };  

  // Record current time
  virtual void record_time(double current_time) noexcept = 0;

  //============================================================================
  // APPENDIX 4: UNUSED FUNCTIONS

  // Assign material id of this particle to nodes
  virtual void append_material_id_to_nodes() const = 0;

  // Return the number of neighbour particles
  virtual unsigned nneighbours() const = 0;

  // Add a neighbour particle
  virtual bool add_neighbour(mpm::Index neighbour_id) = 0;

  // Assign neighbour particles
  virtual bool assign_neighbours(const std::vector<mpm::Index>& neighbours) = 0;

  // Return neighbour ids
  virtual std::set<mpm::Index> neighbours() const = 0;

protected:

  // particleBase id
  Index id_{std::numeric_limits<Index>::max()};
  // coordinates
  VectorDim coordinates_;
  // Cell id
  Index cell_id_{std::numeric_limits<Index>::max()};
  // Status
  bool status_{true};
  // Reference coordinates (in a cell)
  Eigen::Matrix<double, Tdim, 1> xi_;
  // Cell
  std::shared_ptr<Cell<Tdim>> cell_;
  // Container of nodeal pointers
  std::vector<std::shared_ptr<NodeBase<Tdim>>> nodes_;
  // Material
  std::shared_ptr<Material<Tdim>> material_;
  // Unsigned material id
  unsigned material_id_{std::numeric_limits<unsigned>::max()};
  // Material state history variables
  mpm::dense_map state_variables_;

  };  // ParticleBase class
}  // namespace mpm

#include "particle_base.tcc"

#endif  // MPM_PARTICLEBASE_H__
