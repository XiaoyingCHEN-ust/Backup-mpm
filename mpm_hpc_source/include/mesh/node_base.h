#ifndef MPM_NODE_BASE_H_
#define MPM_NODE_BASE_H_

#include <array>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include "data_types.h"
#include "function_base.h"

namespace mpm {

// NodeBase base class for nodes
//! \brief Base class that stores the information about nodes
//! \details Node class: id_ and coordinates.
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim>
class NodeBase {
  
public:
    // Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //============================================================================
  // CONSTRUCT AND DESTRUCT A NODE

  // Constructor with id and coordinates
  NodeBase(mpm::Index id, const VectorDim& coords){};

  // Destructor
  virtual ~NodeBase(){};

  // Delete copy constructor
  NodeBase(const NodeBase<Tdim>&) = delete;

  // Delete assignement operator
  NodeBase& operator=(const NodeBase<Tdim>&) = delete;

  //============================================================================
  // INITIALISE NODE

  // Initialise properties
  virtual void initialise() noexcept = 0;

  // Assign status
  virtual void assign_status(bool status) = 0;

  // Assign coordinates
  virtual void assign_coordinates(const VectorDim& coord) = 0;

  // Assign active id
  virtual void assign_active_id(Index id) = 0;

  // Set ghost id
  virtual void ghost_id(Index gid) = 0; 
  
  // append material id
  virtual void append_material_id(unsigned id) = 0;

  // Assign reference step
  virtual void assign_reference_step(Index ref_step) = 0;

  // Assign rotation matrix
  virtual void assign_rotation_matrix(
      const Eigen::Matrix<double, Tdim, Tdim>& rotation_matrix) = 0;

  // Assign MPI rank to node
  virtual bool mpi_rank(unsigned rank) = 0;

  // Assign MPI rank to node
  virtual std::set<unsigned> mpi_ranks() const = 0;

  // Clear MPI ranks on node
  virtual void clear_mpi_ranks() = 0;

  //============================================================================
  // UPDATE NODAL INFORMATION FROM PARTICLES

  // Update mass at the nodes from particle
  virtual void update_mass(bool update, unsigned phase,
                          double mass) noexcept = 0;

  // Update nodal momentum
  virtual void update_momentum(bool update, unsigned phase,
                              const VectorDim& momentum) noexcept = 0;

  // Update mass at the nodes from particle
  virtual void update_mass_momentum(bool update, unsigned phase, double mass, 
                            const VectorDim& momentum) noexcept = 0;

  // Update external force (body force / traction force)
  virtual void update_external_force(bool update, unsigned phase,
                                    const VectorDim& force) noexcept = 0;

  // Update internal force 
  virtual void update_internal_force(bool update, unsigned phase,
                                    const VectorDim& force) noexcept = 0;

  // update heat capacity at nodes from particle
  virtual void update_heat_capacity(bool update, unsigned phase, 
      const double heat_capacity) noexcept = 0;
  
  // Update heat at the nodes
  virtual void update_heat(bool update, unsigned phase, 
      const double heat) noexcept = 0;

  // Update heat capacity and heat at the nodes
  virtual void update_energy(bool update, unsigned phase, 
      const double heat_capacity, const double heat) noexcept = 0;

  // Update heat conduction at the nodes
  virtual void update_heat_conduction(bool update, unsigned phase, 
      const double heat_conduction) noexcept = 0;

  // Update hydraulic conduction at the nodes
  virtual void update_hydraulic_conduction(bool update, unsigned phase, 
      const double hydraulic_conduction) noexcept = 0;      

  // Update heat convection at the nodes
  virtual void update_heat_convection(bool update, unsigned phase, 
      const double heat_convection) noexcept = 0;

  // Update mass convection at the nodes
  virtual void update_mass_convection(bool update, unsigned phase, 
      const double mass_convection) noexcept = 0;

  // Update covective heat flux_ at the nodes
  virtual void update_covective_heat_flux(bool update, unsigned phase, 
      const double covective_heat_flux) noexcept = 0;

  // Update latent heat at the nodes
  virtual void update_latent_heat(bool update, unsigned phase, 
      const double latent_heat) noexcept = 0;

  // Update heat source at the nodes
  virtual void update_heat_source(bool update, unsigned phase, 
      const double heat_source) noexcept = 0;

  // Update mass source at the nodes
  virtual void update_mass_source(bool update, unsigned phase, 
      const double mass_source) noexcept = 0;

  // Update volumetric strain at the nodes
  virtual void update_volumetric_strain(bool update, unsigned phase, 
      const double volumetric_strain) noexcept = 0;      

  // Update plastic work at the nodes
  virtual void update_plastic_work(bool update, unsigned phase, 
      const double plastic_work) noexcept = 0;

  // Update drag force
  virtual void update_drag_force(const VectorDim& drag_force) = 0;

  // Update internal force (body force / traction force)
  virtual bool update_drag_force_coefficient(bool update, unsigned phase,
                                            const double drag_force) = 0;

  // Update pressure at the nodes from particle
  virtual void update_mass_pressure(unsigned phase, double mass_pressure,
                                    double current_time = 0.) noexcept = 0;

  // Update gas saturation at the nodes from particle
  virtual void update_gas_saturation(unsigned phase, double gas_saturation) noexcept = 0;                                

  // Update pressure at the nodes from particle
  virtual void update_pressure(bool update, unsigned phase, 
                                      double mass_pressure) noexcept = 0;

  // Update density at the nodes from particle
  virtual void update_density(bool update, unsigned phase, 
                                      double density) noexcept = 0;                                       

  // Update K_coeff at the nodes from particle
  virtual void update_K_coeff(bool update, unsigned phase, 
                                      double K_coeff) noexcept = 0; 

  // Update volume at the nodes from particle
  virtual void update_volume(bool update, unsigned phase,
                            double volume, double mean_length) noexcept = 0;

  // Update scalers at the nodes from particle
  virtual void update_scalers(bool update, unsigned scaler_id, 
                                      double scaler) noexcept = 0;

  //============================================================================
  // COMPUTE NODAL VARIABLES

  // Compute velocity from the momentum
  virtual void compute_velocity(double dt) = 0;

  // Compute nodal temperature from heat
  virtual void compute_temperature(unsigned phase) = 0;

  // Compute nodal pressure
  virtual void compute_pressure(unsigned phase) = 0;

  // Compute nodal scaler
  virtual void compute_scalers() = 0;

  // Compute acceleration
  virtual bool compute_acceleration_velocity(unsigned phase,
                                            double dt) noexcept = 0;

  // Compute temperature acceleration and temperature
  virtual bool compute_acceleration_temperature(unsigned phase, 
                                              double dt) noexcept = 0;

  // Compute pressure acceleration and pressure
  virtual bool compute_acceleration_pressure(unsigned phase, 
                                              double dt) noexcept = 0;

  // Compute pressure acceleration and pressure
  virtual bool compute_acceleration_pressure_threephase(unsigned pore_fluid, 
                                              unsigned pore_gas, 
                                              double dt) noexcept = 0;                                              

  // Compute acceleration and velocity for two phase
  virtual bool compute_acc_vel_twophase_explicit(unsigned soil_skeleton, 
      unsigned pore_fluid, unsigned mixture, double dt) = 0;

  // Compute acceleration and velocity for two phase
  virtual bool compute_acc_vel_threephase_explicit(unsigned soil_skeleton, 
      unsigned pore_fluid, unsigned pore_gas, unsigned mixture, 
      double dt) = 0;      

  // Compute acceleration and velocity for two phase
  virtual bool compute_inter_acc_vel_twophase_semi(unsigned soil_skeleton,
      unsigned pore_fluid, unsigned mixture, double dt) = 0;

  // Compute semi-implicit acceleration and velocity
  virtual bool compute_acc_vel_twophase_semi(unsigned phase, double dt) = 0;

  // Update intermediate velocity at the node
  virtual void update_intermediate_velocity(
      const unsigned phase, const Eigen::MatrixXd& velocity_inter,
      double dt) = 0;

  // Update update correction force
  virtual void update_correction_force(const unsigned phase, const unsigned dim,
                                      const double force_cor) = 0;

  // Compute intermediate force
  virtual bool compute_intermediate_force(const double dt) = 0;

  // Compute nodal corrected force
  virtual bool compute_nodal_corrected_force(
      VectorDim& force_cor_part_solid, VectorDim& force_cor_part_water) = 0;

  // Update pore pressure increment at the node
  virtual void update_pore_pressure_increment(
      const Eigen::VectorXd& pore_pressure_increment,
      double current_time = 0.) = 0;

  // Update temperature rate at the node
  virtual void update_temperature_rate(
      const Eigen::VectorXd& temperature_rate, double dt,
      double current_time = 0.) = 0;

  // Compute nodal density
  virtual void compute_density() = 0;

  //============================================================================
  // ASSIGN AND APPLY NODAL BOUDARY CONDITIONS

  // Assign velocity constraint
  virtual bool assign_velocity_constraint(unsigned dir, double velocity) = 0;

  // Apply velocity constraints
  virtual void apply_velocity_constraints() = 0;

  // Assign temperature constraint
  virtual bool assign_temperature_constraint(
      const unsigned phase, const double temperature,
      const std::shared_ptr<FunctionBase>& function) = 0;

  // Apply temperature constraints
  virtual bool apply_temperature_constraints(
    const unsigned phase, const double current_time) = 0;

  // Assign convective_heat constraint
  virtual bool assign_convective_heat_constraint(
      const unsigned phase, const double temperature,
      const std::shared_ptr<FunctionBase>& function, 
      const double coeff, const int set_id) = 0;

  // Apply convective_heat constraints
  virtual bool apply_convective_heat_constraints(
    const unsigned phase, const double current_time) = 0;

  // Assign concentrated force to the node
  virtual bool assign_concentrated_force(
      unsigned phase, unsigned direction, double traction,
      const std::shared_ptr<FunctionBase>& function) = 0;

  // Apply concentrated force to external force
  virtual void apply_concentrated_force(unsigned phase,
                                        double current_time) = 0;

  // Assign friction constraint
  virtual bool assign_friction_constraint(unsigned dir, int sign,
                                          double friction) = 0;

  // Apply friction constraints
  virtual void apply_friction_constraints(double dt) = 0;

  // Assign heat source
  virtual bool assign_heat_source(
      const unsigned phase, const double heat_source,
      const std::shared_ptr<FunctionBase>& function) = 0;

  // Apply heat source
  virtual bool apply_heat_source(
    const unsigned phase, const double current_time) = 0;

  // Assign pressure constraint
  virtual bool assign_pressure_constraint(
      const unsigned phase, double pressure,
      const std::shared_ptr<FunctionBase>& function) = 0;

  // Apply pressure constraints
  virtual bool apply_pressure_constraints(const unsigned phase,
                                          const double current_time) = 0;

  // Assign free surface
  virtual void assign_free_surface(bool free_surface) = 0;

  // Assign water table
  virtual bool assign_water_table(const std::shared_ptr<FunctionBase>& function,
                                  const unsigned dir, const double h0) = 0;

  // Apply water table
  virtual bool apply_water_table(const double current_time) = 0;

  // Assign pressure at the nodes from particle
  virtual void assign_pressure(unsigned phase, double mass_pressure) = 0;

  // Assign nodal pressure increment constraints
  virtual bool assign_nodal_pressure_increment_constraints(const Index step,
                                                          const double dt) = 0;

  //============================================================================
  // RETURN NODE DATA

  // Return id of the nodebase
  virtual Index id() const = 0;

  // Return active id
  virtual Index active_id() = 0;

  // Return ghost id
  virtual Index ghost_id() const = 0;

  // Return status
  virtual bool status() const = 0;

  // Return degrees of freedom
  virtual unsigned dof() const = 0;

  // Return mass at a given node for a given phase
  virtual double mass(unsigned phase) const = 0;

  // Return volume at a given node for a given phase
  virtual double volume(unsigned phase) const = 0;

  // Return real density at a given node for a given phase
  virtual double density(unsigned phase) = 0;

  // Return pressure at a given node for a given phase
  virtual double pressure(unsigned phase) const = 0;

  // Return pressure at a given node for a given phase
  virtual double pressure_acceleration(unsigned phase) const = 0;  

  // Return pore pressure increment
  virtual double pore_pressure_increment() const = 0;

   // Return pore pressure increment
  virtual double gas_saturation() const = 0; 

  // Return temperature
  virtual double temperature(unsigned phase) const = 0;

  // Return heat at a given node for a given phase
  virtual double heat(unsigned phase) const = 0;

  // Return heat conduction at a given node for a given phase
  virtual double heat_conduction(unsigned phase) const = 0;

  // Return covective heat flux_ at a given node for a given phase
  virtual double covective_heat_flux(unsigned phase) const = 0;

  // Return latent heat at a given node for a given phase
  virtual double latent_heat(unsigned phase) const = 0;  

  // Return heat capacity at a given node for a given phase
  virtual double heat_capacity(unsigned phase) const = 0;

  // Return heat source at a given node for a given phase
  virtual double heat_source(unsigned phase) const = 0;

  // Return heat convection at a given node for a given phase
  virtual double heat_convection(unsigned phase) const = 0;

  // Return temperature acceleration
  virtual double temperature_acceleration(unsigned phase) const = 0;

  // Return temperature constraint
  virtual double temperature_constraint(const unsigned phase,
                                const double current_time) = 0; 

  // Return heat source constraint
  virtual double heat_source_constraint(const unsigned phase,
                                const double current_time) = 0; 

  // Return pressure constraint
  virtual double pressure_constraint(const unsigned phase,
                                    const double current_time) = 0;

  // Return velocity constraint
  virtual std::map<unsigned, double>& velocity_constraints() = 0;

  // Return coordinates
  virtual VectorDim coordinates() const = 0;

  // Return velocity
  virtual VectorDim velocity(unsigned phase) const = 0;

  // Return acceleration
  virtual VectorDim acceleration(unsigned phase) const = 0;

  // Return momentum
  virtual VectorDim momentum(unsigned phase) const = 0;

  // Return external force
  virtual VectorDim external_force(unsigned phase) const = 0;

  // Return internal force
  virtual VectorDim internal_force(unsigned phase) const = 0;

  // Return intermediate velocity at the node
  virtual VectorDim intermediate_velocity(const unsigned phase) = 0;

  // Return intermediate acceleration at the node
  virtual VectorDim intermediate_acceleration(const unsigned phase) = 0;

  // Return total intermediate force
  virtual VectorDim force_total_inter() = 0;

  // Return fluid intermediate force
  virtual VectorDim force_fluid_inter() = 0;

  // Return drag force at a given node
  virtual double drag_force_coefficient(unsigned phase) const = 0;

  //! Return external force
  virtual VectorDim reaction_force() const = 0;

  // Return free surface
  virtual bool free_surface() = 0;

  // Return free surface
  virtual bool set_velocity_constraints() = 0;  

  // Return material ids in node
  virtual std::set<unsigned> material_ids() const = 0;

  // Return vector data of particles
  virtual Eigen::VectorXd nodal_vector_data(const std::string& property) = 0;

  // Return scalar data of particles
  virtual double nodal_scalar_data(const std::string& property) = 0;

  // Return scalars of particles
  virtual double smoothed_scalers(unsigned scaler_id) = 0;

  //============================================================================
  // CONTACT： TO BE OPTIMIZED

  // overwrite nodal velocity by rigid particle velocity
  virtual void assign_velocity_from_rigid(unsigned dir,
                                          const double velocity, double dt) = 0;

  // overwrite nodal intermediate velocity by rigid particle velocity
  virtual void assign_intermediate_velocity_from_rigid(double dt) = 0;

  // overwrite nodal corrected  velocity by rigid particle velocity
  virtual void assign_corrected_velocity_from_rigid(double dt) = 0;

  // assign contact normal for rigid particles
  virtual void update_contact_normal(double mass,  
              const Eigen::Matrix<double, Tdim, 1>& mass_contact_normal) = 0;

  // assign contact normal for rigid particles
  virtual void update_rigid_mass_momentum(double mass, 
              const Eigen::MatrixXd& mass_contact_velocity) = 0;

  // Compute contact normal 
  virtual void compute_rigid_velocity(double dt) = 0; 

};  // NodeBase class
}  // namespace mpm

#endif  // MPM_NODE_BASE_H_
