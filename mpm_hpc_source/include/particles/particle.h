#ifndef MPM_PARTICLE_H_
#define MPM_PARTICLE_H_

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <cmath>

#include "cell.h"
#include "logger.h"
#include "particle_base.h"

namespace mpm {

// Particle class
// Base class that stores the information about particles
template <unsigned Tdim>
class Particle : public ParticleBase<Tdim> {

  public:

    // Define a vector of size dimension
    using VectorDim = Eigen::Matrix<double, Tdim, 1>;

    // Define DOFs ？？？
    static const unsigned Tdof = (Tdim == 1) ? 1 : 3 * (Tdim - 1);

    //============================================================================
    // CONSTRUCT AND DESTRUCT A PARTICLE

    // Constructor with id and coordinates
    Particle(Index id, const VectorDim& coord);

    // Constructor with id, coordinates and status    
    Particle(Index id, const VectorDim& coord, bool status);

    // Destructor    
    ~Particle() override{};

    // Delete copy constructor
    Particle(const Particle<Tdim>&) = delete;

    // Delete assignement operator
    Particle& operator=(const Particle<Tdim>&) = delete; 

    //============================================================================
    // ASSIGN INITIAL CONDITIONS

    // Initialise particle HDF5 data
    bool initialise_particle(const HDF5Particle& particle) override;

    // Initialise particle HDF5 data and material
    bool initialise_particle(const HDF5Particle& particle,
          const std::shared_ptr<Material<Tdim>>& material) override;

    // Initialise properties
    void initialise() override;

    // Retrun particle data as HDF5
    HDF5Particle hdf5() override;

    // Assign a material to particle
    bool assign_material(
            const std::shared_ptr<Material<Tdim>>& material) override;

    // Assign affine particle-grid velocity transfer
    void assign_affine_mpm(bool affine_mpm) override {
        this->affine_mpm_ = affine_mpm;
    }

    // Assign initial volume, porosity, density, and mass to the particle
    bool assign_initial_volume(double volume) override;
    void compute_volume(bool is_axisymmetric) noexcept override;
    void compute_mass() override;


    // Assign initial stress
    void assign_initial_stress(
          const Eigen::Matrix<double, 6, 1>& stress) override {
        this->stress_ = stress; 
    }

    // Assign initial temperature
    void assign_initial_temperature(double temperature) override {
        this->temperature_ = temperature;
        this->PIC_temperature_ = temperature;
    }

    //============================================================================
    // APPLY BOUNDARY CONDITIONS

    // Apply particle velocity constraints
    void apply_particle_velocity_constraints(unsigned dir,
                                            double velocity) override {
      this->velocity_(dir) = velocity;
    };

    // Apply particle temperature constraints
    void apply_particle_temperature_constraints(double temperature) override {
      this->temperature_ = temperature;
      this->temperature_acceleration_ = 0;
      this->temperature_increment_ = 0;
    };
    
    // Assign initial free surface
    void assign_particle_free_surface(bool initial_free_surface) override {
      this->initial_free_surface_ = initial_free_surface;
    };

    // Assign initial nonfree surface
    void assign_particle_nonfree_surface(bool initial_nonfree_surface) override {
      this->initial_nonfree_surface_ = initial_nonfree_surface;
    };

    // Assign traction to the particle
    bool assign_particle_traction(unsigned direction, double traction) override;

    // Assign heat source to the particle
    bool assign_particle_heat_source(double heat_source, double dt) override;

    // Assign free surface particle manually
    bool assign_particle_free_surfaces() override;

    // Compute free surface particle based on density
    bool compute_particle_free_surface() override;

    // Assign contact to the particle // TODO: TO BE OPTIMIZED
    bool assign_particle_contact(unsigned dir, double normal) override;
    void map_moving_rigid_velocity_to_nodes(unsigned dir, 
                                            double velocity, 
                                            double dt) noexcept override;
    void map_rigid_mass_momentum_to_nodes() noexcept override; 

    //==========================================================================
    // LOCATE PARTICLES & COMPUTE SHAPE FUNCTIONS

    // Assign a cell to particle
    bool assign_cell(const std::shared_ptr<Cell<Tdim>>& cellptr) override;

    // Assign a cell to particle  
    bool assign_cell_xi(const std::shared_ptr<Cell<Tdim>>& cellptr,
                        const Eigen::Matrix<double, Tdim, 1>& xi) override;

    // Assign a cell id to particle
    bool assign_cell_id(Index id) override;

    // Remove cell for the particle
    void remove_cell() override;

    // Compute reference location cell to particle
    bool compute_reference_location() noexcept override;

    // Compute shape functions and gradients
    void compute_shapefn() noexcept override;

    //==========================================================================
    // MAP PARTICLE INFORMATION TO NODES

    // Map particle mass and momentum to nodes
    void map_mass_momentum_to_nodes() noexcept override;

    // Map body force to nodes
    void map_external_force(const VectorDim& pgravity) override;

    // Map internal force 1D to nodes
    void map_internal_force() override;

    // Map traction force to nodes
    void map_traction_force() noexcept override;

    // Map pressure to nodes
    bool map_pressure_to_nodes(double current_time = 0.) noexcept override;

    // Map particle heat capacity and heat to nodes
    void map_heat_to_nodes() override;

    // Map heat conduction to nodes
    void map_heat_conduction() override;

    // Map heat conduction to nodes
    void map_virtual_heat_flux(bool convective, 
                              const double para_1,
                              const double para_2) override;    

    // Map heat source to nodes
    void map_heat_source() override;

    // Map plastic work to nodes
    void map_plastic_work(double dt) noexcept override;

    //==========================================================================
    // UPDATE PARTICLE INFORMATION

    // Compute updated velocity and position of the particle
    void compute_updated_velocity(double dt, 
                                  double pic = 0,
                                  double damping_factor = 0) override;

    // Compute updated temperature of the particle
    void update_particle_temperature(double dt, double pic_t) noexcept override;

    // Update particle strain
    void update_particle_strain(double dt) noexcept override;

    // Compute thermal strain of the particle    
    void update_particle_thermal_strain() noexcept override;

    // Update particle stress
    void update_particle_stress() noexcept override;

    // Compute updated update porosity of the particle
    void update_particle_volume()  override;

    // Update particle volume of the particle
    void update_particle_density(double dt) override;

    // Compute updated update porosity of the particle
    bool update_particle_porosity(double dt) override;

    // Compute pressure smoothing
    bool compute_pressure_smoothing() noexcept override;    

    //==========================================================================
    // MULTISCALE FUNCTIONS

    // Set stress
    void set_stress(const Eigen::MatrixXd& stresses,
                    bool increment) noexcept override;

    // Set porosity
    void set_porosity(const Eigen::MatrixXd& porosities) noexcept override {
        if (material_id_ != 999) {
            porosity_ = porosities(id_);
        }
    }

    // Set fabric
    void set_fabric(std::string fabric_type,
                    const Eigen::MatrixXd& fabrics) override;

    // Set rotation
    void set_rotation(const Eigen::MatrixXd& rotations) noexcept override {
        if (material_id_ != 999) {
            this->rotation_ = rotations.col(id_);
        }
    }

    // Set fabric
    Eigen::Matrix<double, Tdim, Tdim> fabric(
                  std::string fabric_type) const override;

    // Compute displacement gradient
    void compute_displacement_gradient(double dt, bool thermal) noexcept override;

    //============================================================================
    // RETURN PARTICLE DATA

    // Return particle vector data
    Eigen::VectorXd vector_data(const std::string& property) override {
      return this->vector_property_.at(property)();      
    }

    // Return particle scalar data
    double scalar_data(const std::string& property) override {
      return this->scalar_property_.at(property)();
    } 

    // Return a state variable
    double state_variable(const std::string& var) const override {
      return (state_variables_.find(var) != state_variables_.end())
                ? state_variables_.at(var)
                : std::numeric_limits<double>::quiet_NaN();
    }

    // Return pressure of the particles
    double pressure() const override {
      return (state_variables_.find("pressure") != state_variables_.end())
                ? state_variables_.at("pressure")
                : std::numeric_limits<double>::quiet_NaN();
    }    

    // Return volume
    double volume() const override { 
      return volume_; 
    }

    // Return incremental volumetric strain
    double dvolumetric_strain() const override { 
      return dvolumetric_strain_; 
    }

    // Return PIC temperature
    double temperature() const override { 
      return PIC_temperature_; 
    }

    // Return temperature increment
    double temperature_increment() const override { 
      return temperature_increment_; 
    } 

    // Return displacement
    VectorDim displacement() const override { 
      return displacement_; 
    }

    // Return strain rate
    Eigen::Matrix<double, 6, 1> strain_rate() const override {
      return strain_rate_; 
    }

    // Return displacement gradient
    Eigen::Matrix<double, Tdim, Tdim> displacement_gradient() const override { 
      return displacement_gradient_; 
    }

    // Return cell ID
    Index cell_id() const override { 
      return cell_id_; 
    }

    // Check if this is an active cell
    bool cell_ptr() const override { 
      return cell_ != nullptr; 
    } 

    //! Check if is an initial free surface particle
    bool initial_free_surface() override { 
      return initial_free_surface_; 
    }

    // Check if this is a free surface particle
    bool free_surface() override {
      return free_surface_; 
    };

    void record_time(double current_time) noexcept override {
      current_time_ = current_time;
    }

  //============================================================================
  // APPENDIX 4: UNUSED FUNCTIONS

    // Assign material id of this particle to nodes
    void append_material_id_to_nodes() const override;

    // Add a neighbour particle
    bool add_neighbour(mpm::Index neighbour_id) override;

    // Assign neighbour particles
    bool assign_neighbours(
        const std::vector<mpm::Index>& neighbours_set) override;

    // Return the number of neighbour particles
    unsigned nneighbours() const override { 
      return neighbours_.size(); 
    };

    // Return neighbour ids
    std::set<mpm::Index> neighbours() const override { 
      return neighbours_; 
    };

  protected:

    // Compute strain rate of the particle
    inline Eigen::Matrix<double, 6, 1> compute_strain_rate(
        const Eigen::MatrixXd& grad_shapefn, Eigen::VectorXd& shapefn, 
        unsigned phase) noexcept;

    // Compute temperature gradient of the particle
    inline Eigen::Matrix<double, Tdim, 1> compute_temperature_gradient(
      unsigned phase) noexcept; 

    // Compute mass gradient of the particle
    inline Eigen::Matrix<double, Tdim, 1> compute_mass_gradient(
      unsigned phase) noexcept;       

    // Compute velocity gradient of the particle
    inline Eigen::Matrix<double, Tdim, Tdim> compute_velocity_gradient(
        const Eigen::MatrixXd& grad_shapefn, unsigned phase) noexcept;

    // Compute velocity gradient of the particle
    inline Eigen::Matrix<double, Tdim, Tdim> compute_affine_matrix(
        const Eigen::VectorXd& shapefn, unsigned phase) noexcept;

    // Compute jaumann stress of the particle
    inline Eigen::Matrix<double, 6, 1> compute_jaumann_stress() noexcept; 

    // Reshape a 3 * 3 tensor to a vector
    inline Eigen::Matrix<double, 9, 1> reshape_tensor(
                        Eigen::Matrix<double, Tdim, Tdim> tensor) {
        Eigen::Matrix<double, 9, 1> reshaped_tensor = Eigen::VectorXd::Zero(9);
        for (unsigned i = 0; i < Tdim; ++i)
            for (unsigned j = 0; j < Tdim; ++j) {
                reshaped_tensor(3 * i + j) = tensor(i, j);
            }
        return reshaped_tensor;
    };

  protected:
  
    // Inherit properties from class ParticleBase
    using ParticleBase<Tdim>::id_;
    using ParticleBase<Tdim>::coordinates_;
    using ParticleBase<Tdim>::xi_; // Reference coordinates (in a cell)
    using ParticleBase<Tdim>::cell_;
    using ParticleBase<Tdim>::cell_id_;
    using ParticleBase<Tdim>::nodes_;
    using ParticleBase<Tdim>::status_;
    using ParticleBase<Tdim>::material_;
    using ParticleBase<Tdim>::material_id_;
    using ParticleBase<Tdim>::state_variables_;

    // Shape functions
    Eigen::VectorXd shapefn_;
    Eigen::VectorXd shapefn_centroid_;
    Eigen::MatrixXd dn_dx_;
    Eigen::MatrixXd dn_dx_centroid_;
    Eigen::VectorXd gradient_shapefn_;;

    // Scalar properties
    double current_time_{0}; 
    double volume_{std::numeric_limits<double>::max()};
    double density_{1.0};
    double porosity_;
    double mass_{1.0};
    double solid_fraction_{1.0}; 
    double mass_density_{1.0};
    double dvolumetric_strain_;  
    double volumetric_strain_;
    double ddeviatoric_strain_;  
    double deviatoric_strain_;
    double temperature_;
    double PIC_temperature_;
    double FLIP_temperature_;
    double temperature_acceleration_;
    double temperature_increment_;
    double temperature_increment_cent_;    
    double heat_capacity_{1.0};  
    double heat_source_;
    double dthermal_volumetric_strain_; 
    double thermal_volumetric_strain_;
    // To be optimized
    double plastic_work_;

    // Vector properties 
    Eigen::Matrix<double, Tdim, 1> natural_size_;
    Eigen::Matrix<double, Tdim, 1> size_;
    Eigen::Matrix<double, Tdim, 1> velocity_;
    Eigen::Matrix<double, Tdim, 1> acceleration_;
    Eigen::Matrix<double, Tdim, 1> displacement_;
    Eigen::Matrix<double, Tdim, 1> mass_gradient_; 
    Eigen::Matrix<double, Tdim, 1> outward_normal_;        
    Eigen::Matrix<double, Tdim, 1> temperature_gradient_;    
    Eigen::Matrix<double, Tdim, 1> heat_flux_;
    Eigen::Matrix<double, Tdim, 1> traction_;
    Eigen::Matrix<double, Tdim, 1> contact_normal_;
    Eigen::Matrix<double, Tdim, 1> contact_tangential_;
    Eigen::Matrix<double, 3, 1> rotation_;

    // Tensor properties
    Eigen::Matrix<double, 6, 1> stress_;
    Eigen::Matrix<double, 6, 1> jaumann_stress_;
    Eigen::Matrix<double, 6, 1> strain_rate_;
    Eigen::Matrix<double, 6, 1> dstrain_;
    Eigen::Matrix<double, 6, 1> strain_;
    Eigen::Matrix<double, 6, 1> dthermal_strain_;
    Eigen::Matrix<double, 6, 1> thermal_strain_;
    Eigen::Matrix<double, Tdim, Tdim> displacement_gradient_;
    Eigen::Matrix<double, Tdim, Tdim> deformation_gradient_;
    Eigen::Matrix<double, Tdim, Tdim> fabric_CN_; // Contact normal based
    Eigen::Matrix<double, Tdim, Tdim> fabric_PO_; // Particle orientation based
    Eigen::Matrix<double, Tdim, Tdim> C_matrix_;
    Eigen::Matrix<double, Tdim, Tdim> velocity_gradient_;        

    // Bool properties
    bool set_traction_{false};
    bool set_contact_{false};  
    bool set_heat_source_{false};
    bool initial_free_surface_{false};
    bool initial_nonfree_surface_{false};
    bool free_surface_{false};
    bool is_axisymmetric_{false};
    bool affine_mpm_{false};

    // Map
    std::map<unsigned, double> particle_velocity_constraints_;
    std::map<unsigned, double> particle_temperature_constraints_;
    std::map<std::string, std::function<double()>> scalar_property_;  
    std::map<std::string, std::function<Eigen::MatrixXd()>> vector_property_;
    std::set<mpm::Index> neighbours_;
    std::unique_ptr<spdlog::logger> console_;

  };  // Particle class
}  // namespace mpm

#include "particle.tcc"

#endif  // MPM_PARTICLE_H__
