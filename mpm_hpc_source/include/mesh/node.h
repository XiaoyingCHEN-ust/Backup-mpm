#ifndef MPM_NODE_H_
#define MPM_NODE_H_

#include "logger.h"
#include "node_base.h"
#include <iostream>

namespace mpm {

// Node class
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
class Node : public NodeBase<Tdim> {

public:
  // Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //============================================================================
  // CONSTRUCT AND DESTRUCT A NODE

  // Constructor with id, coordinates and dof
  Node(Index id, const VectorDim& coord);

  // Virtual destructor
  ~Node() override{};

  // Delete copy constructor
  Node(const Node<Tdim, Tdof, Tnphases>&) = delete;

  // Delete assignement operator
  Node& operator=(const Node<Tdim, Tdof, Tnphases>&) = delete;

  //============================================================================
  // INITIALISE NODE

  // Initialise nodal properties
  void initialise() noexcept override;

  // Assign status
  void assign_status(bool status) override { status_ = status; }

  // Assign coordinates
  void assign_coordinates(const VectorDim& coord) override {
    coordinates_ = coord;
  }

  // Assign active id
  void assign_active_id(Index id) override { active_id_ = id; }

  // Set ghost id
  void ghost_id(Index gid) override { ghost_id_ = gid; }

  // Add material id from material points to list of materials in materials_
  void append_material_id(unsigned id) override;

  // Assign reference step
  void assign_reference_step(Index ref_step) override { ref_step_ = ref_step; }

  // Assign rotation matrix
  void assign_rotation_matrix(
      const Eigen::Matrix<double, Tdim, Tdim>& rotation_matrix) override {
    rotation_matrix_ = rotation_matrix;
    generic_boundary_constraints_ = true;
  }

  // Assign MPI rank to node
  bool mpi_rank(unsigned rank) override;

  // Assign MPI rank to node
  std::set<unsigned> mpi_ranks() const override { return mpi_ranks_; }

  // Clear MPI rank
  void clear_mpi_ranks() override { mpi_ranks_.clear(); }

  //============================================================================
  // UPDATE NODAL INFORMATION FROM PARTICLES

  // Update mass at the nodes from particle
  void update_mass(bool update, unsigned phase, double mass) noexcept override;

  // Update momentum at the nodes
  void update_momentum(bool update, unsigned phase,
                      const VectorDim& momentum) noexcept override;

  // Update mass and momentum at the nodes from particle
  void update_mass_momentum(bool update, unsigned phase, double mass, 
                            const VectorDim& momentum) noexcept override;

  // Update external force (body force / traction force)
  void update_external_force(bool update, unsigned phase,
                              const VectorDim& force) noexcept override;

  // Update internal force 
  void update_internal_force(bool update, unsigned phase,
                              const VectorDim& force) noexcept override;

  // update heat capacity at nodes from particle
  void update_heat_capacity(bool update, unsigned phase, 
      const double heat_capacity) noexcept override;
  
  // Update heat at the nodes
  void update_heat(bool update, unsigned phase, 
      const double heat) noexcept override;

  // Update heat capacity and heat at the nodes
  void update_energy(bool update, unsigned phase, 
      const double heat_capacity, const double heat) noexcept override;

  // Update heat conduction at the nodes
  void update_heat_conduction(bool update, unsigned phase, 
      const double heat_conduction) noexcept override;

  // Update hydraulic conduction at the nodes
  void update_hydraulic_conduction(bool update, unsigned phase, 
      const double hydraulic_conduction) noexcept override;      

  // Update heat convection at the nodes
  void update_heat_convection(bool update, unsigned phase, 
      const double heat_convection) noexcept override;

  // Update mass convection at the nodes
  void update_mass_convection(bool update, unsigned phase, 
      const double mass_convection) noexcept override;      

  // Update covective heat flux at the nodes
  void update_covective_heat_flux(bool update, unsigned phase, 
      const double covective_heat_flux) noexcept override;

  // Update latent heat at the nodes
  void update_latent_heat(bool update, unsigned phase, 
      const double latent_heat) noexcept override;

  // Update heat source at the nodes
  void update_heat_source(bool update, unsigned phase, 
      const double heat_source) noexcept override;

  // Update volumetric strain at the nodes
  void update_volumetric_strain(bool update, unsigned phase, 
      const double volumetric_strain) noexcept override;  

   // Update heat source at the nodes
  void update_mass_source(bool update, unsigned phase, 
      const double mass_source) noexcept override;  
               

  // Update plastic work at the nodes
  void update_plastic_work(bool update, unsigned phase, 
      const double plastic_work) noexcept override;      

  // Update drag force
  void update_drag_force(const VectorDim& drag_force) {
    drag_force_ += drag_force;
  }

  // Update drag force coefficient
  bool update_drag_force_coefficient(bool update, unsigned phase,
                                      const double drag_force) override;

  // Update pressure at the nodes from particle
  void update_mass_pressure(unsigned phase, double mass_pressure,
                            double current_time = 0.) noexcept override;

  // Update gas saturation at the nodes from particle
  void update_gas_saturation(unsigned phase, double gas_saturation) noexcept override;                           

  // Update pressure at the nodes from particle
  void update_pressure(bool update, unsigned phase, 
                                double mass_pressure) noexcept override;

  // Update density at the nodes from particle
  void update_density(bool update, unsigned phase, 
                                double density) noexcept override;                                                            

  // Update K_coeff at the nodes from particle
  void update_K_coeff(bool update, unsigned phase, 
                                double K_coeff) noexcept override;    

  // Update volume at the nodes from particle
  void update_volume(bool update, unsigned phase,
                      double volume, double mean_length) noexcept override;

  // Update scalers at the nodes from particle
  void update_scalers(bool update, unsigned scaler_id,
                      double scaler) noexcept override;

  //============================================================================
  // COMPUTE NODAL VARIABLES

  // Compute velocity from the momentum
  void compute_velocity(double dt) override;

  // Compute nodal temperature from heat
  void compute_temperature(unsigned phase) override;

  // Compute nodal pressure
  void compute_pressure(unsigned phase) override;

  // Compute nodal scaler
  void compute_scalers() override;
  
  // Compute acceleration and velocity
  bool compute_acceleration_velocity(unsigned phase,
                                      double dt) noexcept override;

  // Compute temperature acceleration and temperature
  bool compute_acceleration_temperature(unsigned phase, 
      double dt) noexcept override;

  // Compute pressure acceleration and pressure
  bool compute_acceleration_pressure_threephase(unsigned pore_fluid, unsigned pore_gas, 
      double dt) noexcept override;

  // Compute pressure acceleration and pressure
  bool compute_acceleration_pressure(unsigned phase, 
      double dt) noexcept override;      

  // Compute acceleration and velocity for two phase
  bool compute_acc_vel_twophase_explicit(unsigned soil_skeleton,
        unsigned pore_fluid, unsigned mixture, double dt) override;

  // Compute acceleration and velocity for two phase
  bool compute_acc_vel_threephase_explicit(unsigned soil_skeleton,
        unsigned pore_fluid, unsigned pore_gas, unsigned mixture, 
        double dt) override;

  // Compute acceleration and velocity for two phase
  bool compute_inter_acc_vel_twophase_semi(unsigned soil_skeleton,
        unsigned pore_fluid, unsigned mixture, double dt) override;

  // Compute semi-implicit acceleration and velocity
  bool compute_acc_vel_twophase_semi(unsigned phase, double dt);

  // Update intermediate elocity at the node
  void update_intermediate_velocity(const unsigned phase,
                                    const Eigen::MatrixXd& acceleration_inter,
                                    double dt) {
    // index
    acceleration_inter_.col(phase) =
        acceleration_inter.row(active_id_).transpose();
    velocity_inter_.col(phase) =
        velocity_.col(phase) +
        dt * acceleration_inter.row(active_id_).transpose();         
  }

  // Update update correction force
  void update_correction_force(const unsigned phase, const unsigned dim,
                                const double force_cor) {
    force_cor_.col(phase)(dim) = force_cor;
  }

  // Compute intermediate force
  bool compute_intermediate_force(const double dt) override;

  // Compute nodal corrected force
  bool compute_nodal_corrected_force(VectorDim& force_cor_part_solid,
                                    VectorDim& force_cor_part_water) override;

  // Update pore pressure increment at the node
  void update_pore_pressure_increment(
      const Eigen::VectorXd& pore_pressure_increment,
      double current_time = 0.) override;

  // Update temperature rate at the node
  void update_temperature_rate(
      const Eigen::VectorXd& temperature_rate, double dt,
      double current_time = 0.) override;

  // Compute nodal density
  void compute_density() override;

  //============================================================================
  // ASSIGN AND APPLY NODAL BOUDARY CONDITIONS

  // Assign velocity constraint
  bool assign_velocity_constraint(unsigned dir, double velocity) override;

  // Apply velocity constraints
  void apply_velocity_constraints() override;

  // Assign temperature constraint
  bool assign_temperature_constraint(
      const unsigned phase, const double temperature,
      const std::shared_ptr<FunctionBase>& function) override;   
  
  // Apply temperature constraints
  bool apply_temperature_constraints(
    const unsigned phase, const double current_time) override;

  // Assign convective_heat constraint
  bool assign_convective_heat_constraint(
      const unsigned phase, const double temperature,
      const std::shared_ptr<FunctionBase>& function, 
      const double coeff, const int set_id) override;   

  // Apply convective_heat constraints
  bool apply_convective_heat_constraints(
    const unsigned phase, const double current_time) override;

  // Assign concentrated force to the node
  bool assign_concentrated_force(
      unsigned phase, unsigned direction, double force,
      const std::shared_ptr<FunctionBase>& function) override;

  // Apply concentrated force to external force
  void apply_concentrated_force(unsigned phase, double current_time) override;

  // Assign friction constraint
  bool assign_friction_constraint(unsigned dir, int sign,
                                  double friction) override;

  // Apply friction constraints
  void apply_friction_constraints(double dt) override;

  // Assign heat source
  bool assign_heat_source(
      const unsigned phase, const double heat_source,
      const std::shared_ptr<FunctionBase>& function) override;   
  
  // Apply heat source
  bool apply_heat_source(
    const unsigned phase, const double current_time) override;

  // Assign pressure constraint
  bool assign_pressure_constraint(
      const unsigned phase, const double pressure,
      const std::shared_ptr<FunctionBase>& function) override;

  // Apply pressure constraints
  bool apply_pressure_constraints(const unsigned phase,
                                  const double current_time) override;

  // Assign free surface
  void assign_free_surface(bool free_surface) override {
    free_surface_ = free_surface;
  }

  // Assign water table
  bool assign_water_table(const std::shared_ptr<FunctionBase>& function,
                          const unsigned dir, const double h0) override;

  // Apply water table
  bool apply_water_table(const double current_time) override;

  // Assign pressure at the nodes from particle
  void assign_pressure(unsigned phase, double mass_pressure) override;

  // Assign nodal pressure increment constraints
  bool assign_nodal_pressure_increment_constraints(const Index step,
                                                    const double dt) override;

  //============================================================================
  // RETURN NODE DATA

  // Return id of the nodebase
  Index id() const override { return id_; }

  // Return active id
  Index active_id() override { return active_id_; }

  // Return ghost id
  Index ghost_id() const override { return ghost_id_; }

  // Return status
  bool status() const override { return status_; }

  // Return degrees of freedom
  unsigned dof() const override { return dof_; }

  // Return mass at a given node for a given phase
  double mass(unsigned phase) const override { return mass_(phase); }

  // Return volume at a given node for a given phase
  double volume(unsigned phase) const override { return volume_(phase); }

  // Return real density at a given node for a given phase
  double density(unsigned phase) override { return density_(phase); }

  // Return pressure at a given node for a given phase
  double pressure(unsigned phase) const override { return pressure_(phase); }

  // Return pressure at a given node for a given phase
  double gas_saturation() const override { return gas_saturation_; }

  // Return pore pressure increment
  double pore_pressure_increment() const override {
    return pore_pressure_increment_;
  }

  // Return pore pressure acceleration
  double pressure_acceleration(unsigned phase) const override {
    return  pressure_acceleration_(phase);
  }  

  // Return temperature at a given node for a given phase
  double temperature(unsigned phase) const override {
    return temperature_(phase);
  }

  // Return momentum at a given node for a given phase
  double heat(unsigned phase) const override { return heat_(phase); }

  // Return heat conduction at a given node for a given phase
  double heat_conduction(unsigned phase) const override {
    return heat_conduction_(phase);
  }

  // Return covective_heat_flux at a given node for a given phase
  double covective_heat_flux(unsigned phase) const override {
    return convective_heat_flux_(phase);
  }

  // Return latent heat at a given node for a given phase
  double latent_heat(unsigned phase) const override {
    return latent_heat_(phase);
  }  

  // Return heat source at a given node for a given phase
  double heat_source(unsigned phase) const override {
    return heat_source_(phase);
  }

  // Return heat capacity at a given node for a given phase
  double heat_capacity(unsigned phase) const override { 
    return heat_capacity_(phase); 
  }
  
  // Return heat convection at a given node for a given phase
  double heat_convection(unsigned phase) const override {
    return heat_convection_(phase);
  }

  // Return temperature acceleration at a given node for a given phase
  double temperature_acceleration(unsigned phase) const override {
    return temperature_acceleration_(phase);
  }

  // Return temperature constraint
  double temperature_constraint(const unsigned phase,
                              const double current_time) override {
    if (temperature_constraints_.find(phase) != temperature_constraints_.end()) {
      const double scalar =
          (temperature_function_.find(phase) != temperature_function_.end())
              ? temperature_function_[phase]->value(current_time)
              : 1.0;

      return scalar * temperature_constraints_[phase];
    } else
      return std::numeric_limits<double>::max();
  }

  // Return heat source constraint
  double heat_source_constraint(const unsigned phase,
                              const double current_time) override {
    if (heat_source_constraints_.find(phase) != heat_source_constraints_.end()) {
      const double scalar =
          (heat_source_function_.find(phase) != heat_source_function_.end())
              ? heat_source_function_[phase]->value(current_time)
              : 1.0;

      return scalar * heat_source_constraints_[phase];
    } else
      return std::numeric_limits<double>::max();
  }

  // Return pressure constraint
  double pressure_constraint(const unsigned phase,
                              const double current_time) override {
    if (pressure_constraints_.find(phase) != pressure_constraints_.end()) {
      const double scalar =
          (pressure_function_.find(phase) != pressure_function_.end())
              ? pressure_function_[phase]->value(current_time)
              : 1.0;

      return scalar * pressure_constraints_[phase];
    } else
      return std::numeric_limits<double>::max();
  }

  // Return velocity constraint
  std::map<unsigned, double>& velocity_constraints() override {
    return velocity_constraints_;
  }

  // Return coordinates
  VectorDim coordinates() const override { return coordinates_; }

  // Return velocity at a given node for a given phase
  VectorDim velocity(unsigned phase) const override {
    return velocity_.col(phase);
  }

  // Return acceleration at a given node for a given phase
  VectorDim acceleration(unsigned phase) const override {
    return acceleration_.col(phase);
  }

  // Return momentum at a given node for a given phase
  VectorDim momentum(unsigned phase) const override {
    return momentum_.col(phase);
  }

  // Return external force at a given node for a given phase
  VectorDim external_force(unsigned phase) const override {
    return external_force_.col(phase);
  }

  // Return internal force at a given node for a given phase
  VectorDim internal_force(unsigned phase) const override {
    return internal_force_.col(phase);
  }

  // Update intermediate elocity at the node
  VectorDim intermediate_velocity(const unsigned phase) {  
    return velocity_inter_.col(phase);
  }

  VectorDim intermediate_acceleration(const unsigned phase) {
    return acceleration_inter_.col(phase);
  }

  // Return total intermediate force
  VectorDim force_total_inter() override { return force_total_inter_; }

  // Return fluid intermediate force
  VectorDim force_fluid_inter() override { return force_fluid_inter_; }

  // Return drag force at a given node
  double drag_force_coefficient(unsigned phase) const override {
    return drag_force_coefficient_(phase);
  }

  // Return external force
  VectorDim reaction_force() const override { return reaction_force_; }

  // Return free surface bool
  bool free_surface() override { return free_surface_; }

  // Return material ids in node
  std::set<unsigned> material_ids() const override { return material_ids_; }

  // Return particle vector data
  Eigen::VectorXd nodal_vector_data(const std::string& property) override {
    return this->nodal_vector_property_.at(property)();      
  }

  // Return particle scalar data
  double nodal_scalar_data(const std::string& property) override {
    return this->nodal_scalar_property_.at(property)();
  }

  bool set_velocity_constraints() override {return set_velocity_constraints_; }; 

  double smoothed_scalers(unsigned scaler_id) override { 
    return scaler_(scaler_id);
  }

  //============================================================================
  // CONTACT： TO BE OPTIMIZED

  // overwrite nodal velocity by rigid particle velocity
  void assign_velocity_from_rigid(unsigned dir,
                                  const double velocity, double dt) override;

  // overwrite nodal intermediate velocity by rigid particle velocity
  void assign_intermediate_velocity_from_rigid(double dt) override;

  // overwrite nodal corrected velocity by rigid particle velocity
  void assign_corrected_velocity_from_rigid(double dt) override;  

  // assign contact normal for rigid particles  
  void update_contact_normal(double mass, 
            const Eigen::Matrix<double, Tdim, 1>& mass_contact_normal) override;

  // assign contact normal for rigid particles  
  void update_rigid_mass_momentum(double mass,
            const Eigen::MatrixXd& mass_contact_velocity) override;                      

  // Compute contact normal                                                 
  void compute_rigid_velocity(double dt) override; 

private:

  // Mutex
  std::mutex node_mutex_;
  // nodebase id
  Index id_{std::numeric_limits<Index>::max()};
  // Global index for active node
  Index active_id_{std::numeric_limits<Index>::max()};  
  // shared ghost id
  Index ghost_id_{std::numeric_limits<Index>::max()};
  // nodal coordinates
  VectorDim coordinates_;
  // Degrees of freedom
  unsigned dof_{std::numeric_limits<unsigned>::max()};
  // Status if the node is actived
  bool status_{false};
  // current time
  double current_time_{0};
  // MPI ranks
  std::set<unsigned> mpi_ranks_;

  //  Reaction force
  Eigen::Matrix<double, 1, Tnphases> mass_;
  Eigen::Matrix<double, 1, Tnphases> volume_;
  Eigen::Matrix<double, 1, Tnphases> volume_fraction_;   
  Eigen::Matrix<double, 1, Tnphases> density_;
  Eigen::Matrix<double, 1, Tnphases> mass_density_; 
  Eigen::Matrix<double, 1, Tnphases> pressure_;
  Eigen::Matrix<double, 1, Tnphases> pressure_acceleration_;  
  Eigen::Matrix<double, 1, Tnphases> mass_pressure_;  
  Eigen::Matrix<double, 1, Tnphases> temperature_;
  Eigen::Matrix<double, 1, Tnphases> temperature_acceleration_;  
  Eigen::Matrix<double, 1, Tnphases> heat_;
  Eigen::Matrix<double, 1, Tnphases> heat_capacity_;
  Eigen::Matrix<double, 1, Tnphases> heat_conduction_;
  Eigen::Matrix<double, 1, Tnphases> hydraulic_conduction_;  
  Eigen::Matrix<double, 1, Tnphases> heat_source_;
  Eigen::Matrix<double, 1, Tnphases> mass_source_;  
  Eigen::Matrix<double, 1, Tnphases> plastic_work_;
  Eigen::Matrix<double, 1, Tnphases> latent_heat_;
  Eigen::Matrix<double, 1, Tnphases> heat_convection_;
  Eigen::Matrix<double, 1, Tnphases> mass_convection_;  
  Eigen::Matrix<double, 1, Tnphases> convective_heat_flux_;
  Eigen::Matrix<double, 1, Tnphases> drag_force_coefficient_;
  Eigen::Matrix<double, 1, Tnphases> K_coeff_; 
  Eigen::Matrix<double, 1, Tnphases> volumetric_strain_;
  Eigen::Matrix<double, 1, 8> scaler_;
  // p^(t+1) - beta * p^(t)
  double pore_pressure_increment_;
  double mean_length_;
  double gas_saturation_;

  // Vector properties
  Eigen::Matrix<double, Tdim, Tnphases> velocity_;
  Eigen::Matrix<double, Tdim, Tnphases> velocity_inter_;  
  Eigen::Matrix<double, Tdim, Tnphases> acceleration_;
  Eigen::Matrix<double, Tdim, Tnphases> acceleration_inter_;
  Eigen::Matrix<double, Tdim, Tnphases> momentum_;
  Eigen::Matrix<double, Tdim, Tnphases> external_force_;
  Eigen::Matrix<double, Tdim, Tnphases> internal_force_;
  Eigen::Matrix<double, Tdim, Tnphases> concentrated_force_;
  Eigen::Matrix<double, Tdim, Tnphases> force_cor_;


  VectorDim force_fluid_inter_;
  VectorDim force_total_inter_;
  VectorDim drag_force_;
  VectorDim drag_force_liquid_;
  VectorDim drag_force_gas_;

  Eigen::Matrix<double, Tdim, 1> reaction_force_;
  Eigen::Matrix<unsigned, Tdim, 1> reaction_dir_;
  // for general velocity constraints
  Eigen::Matrix<double, Tdim, Tdim> rotation_matrix_;

  // Boundary condition constraints
  std::set<unsigned> material_ids_;
  std::map<unsigned, double> velocity_constraints_;
  std::map<unsigned, double> pressure_constraints_;
  std::tuple<unsigned, int, double> friction_constraint_;
  std::map<unsigned, double> temperature_constraints_;
  std::map<unsigned, double> heat_source_constraints_;
  std::map<unsigned, double> convective_heat_constraints_;  
  // Mathematical function for force
  std::map<unsigned, std::shared_ptr<FunctionBase>> pressure_function_;
  std::map<unsigned, std::shared_ptr<FunctionBase>> temperature_function_;
  std::map<unsigned, std::shared_ptr<FunctionBase>> heat_source_function_;
  std::shared_ptr<FunctionBase> force_function_{nullptr};
  std::shared_ptr<FunctionBase> water_table_function_{nullptr};

  // A general velocity (non-Cartesian/inclined) constraint is specified at the node
  bool generic_boundary_constraints_{false};
  bool set_velocity_constraints_{false};  
  // Frictional constraints
  bool friction_{false};
  // Contact
  bool contact_{false}; 
  // Free surface
  bool free_surface_{false};
  // Water table
  bool water_table_{false};
  bool convective_heat_boundary_{false};
  // Water table
  std::pair<unsigned, double> h0_;
  // Logger
  std::unique_ptr<spdlog::logger> console_;  
  // Refence step for pressure increment constraint
  Index ref_step_{std::numeric_limits<Index>::max()};

  std::map<std::string, std::function<double()>> nodal_scalar_property_;  
  std::map<std::string, std::function<Eigen::MatrixXd()>> nodal_vector_property_;

  // Contact - TO BE OPTIMIZED 
  double rigid_mass_, rigid_mass2_;
  Eigen::Matrix<double, Tdim, Tnphases> rigid_velocity_;
  Eigen::Matrix<double, Tdim, Tnphases> rigid_momentum_;
  Eigen::Matrix<double, Tdim, Tnphases> mass_rigid_velocity_;
  Eigen::Matrix<double, Tdim, Tnphases> rigid_acceleration_;
  Eigen::Matrix<double, Tdim, 1> mass_contact_normal_;
  Eigen::Matrix<double, Tdim, 1> contact_normal_;
  Eigen::Matrix<double, Tdim, 1> contact_tangential_;

};  // Node class
}  // namespace mpm

#include "node.tcc"
#include "node_twophase.tcc"
#include "node_threephase.tcc"
#endif  // MPM_NODE_H_

