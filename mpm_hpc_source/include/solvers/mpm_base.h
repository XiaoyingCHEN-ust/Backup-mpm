#ifndef MPM_MPM_BASE_H_
#define MPM_MPM_BASE_H_

#include <numeric>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <Eigen/Sparse>
#include "cg_eigen.h"

// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif
#include "tbb/task_group.h"

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "container.h"
#include "loads_bcs/friction_constraint.h"
#include "loads_bcs/velocity_constraint.h"
#include "loads_bcs/temperature_constraint.h"
#include "loads_bcs/pore_pressure_constraint.h"
#include "particle.h"
#include "solvers/mpm.h"
#include "velocity_constraint.h"
#include "temperature_constraint.h"
#include "pore_pressure_constraint.h"

namespace mpm {

//! Stress update method
//! USF: Update Stress First
//! USL: Update Stress Last
//! MUSL: Modified Stress Last
enum class StressUpdate { USF, USL, MUSL };
extern std::map<std::string, StressUpdate> stress_update;

//! Damping type
//! None: No damping is specified
//! Cundall: Cundall damping
//! Viscous: Viscous damping
enum class Damping { None, Cundall, Viscous };

//! MPMBase class
//! \brief A class that implements the fully base one phase mpm
//! \details A Base MPM class
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMBase : public MPM {
 public:
  //! Default constructor
  MPMBase(const std::shared_ptr<IO>& io);

  //! Initialise mesh
  bool initialise_mesh() override;

  //! Initialise particles
  bool initialise_particles() override;

  //! Initialise materials
  bool initialise_materials() override;

  //! Initialise loading
  bool initialise_loads() override;

  //! Initialise math functions
  bool initialise_math_functions(const Json&) override;

  //! Initialise vtk output for single phase
  bool initialise_vtk() override;

  //! Initialise vtk output for two phase
  bool initialise_vtk_twophase() override;

  //! Solve
  bool solve() override { return true; }

  //! Preprocessing
  bool pre_process() override { return true; }

  //! Get displcament gradient
  bool get_deformation_task() override { return true; }

  //! Get analysis information
  void get_info(unsigned& dim, bool& resume,
                unsigned& checkpoint_step) override{};

  //! get step, time
  void get_status(double& dt, unsigned& step, unsigned& nsteps,
                  unsigned& output_steps) override {}

  //! send deformation task
  bool send_deformation_task(
      std::vector<unsigned>& id,
      std::vector<Eigen::MatrixXd>& displacement_gradients) override {
    return true;
  }

  //! send temperature task
  bool send_temperature_task(
      std::vector<unsigned>& id,
      std::vector<double>& particle_temperature) override {
    return true;
  }

  //! send temperature task
  bool send_saturation_task(
      std::vector<unsigned>& id,
      std::vector<double>& hydrate_saturation) override {
    return true;
  }  

  //! Set particle stress
  bool set_stress_task(const Eigen::MatrixXd& stresses,
                       bool increment) override {
    return true;
  }

  //! Set particle porosity
  bool set_porosity_task(const Eigen::MatrixXd& porosities) override {
    return true;
  }

  // Set particle fabric
  bool set_fabric_task(std::string fabric_type,
                       const Eigen::MatrixXd& fabrics) override {
    return true;
  }

  // Set particle rotation
  bool set_rotation_task(const Eigen::MatrixXd& rotations) override {
    return true;
  }

  bool update_state_task() override { return true; }

  //! Checkpoint resume
  bool checkpoint_resume() override;

#ifdef USE_VTK
  //! Write VTK files
  void write_vtk(mpm::Index step, mpm::Index max_steps) override;
#endif

#ifdef USE_PARTIO
  //! Write PARTIO files
  void write_partio(mpm::Index step, mpm::Index max_steps) override;
#endif

  //! Write HDF5 files
  void write_hdf5(mpm::Index step, mpm::Index max_steps) override;

  //! Write reaction force
  void write_reaction_force(bool overwrite, mpm::Index step,
                            mpm::Index max_steps) override;

 private:
  //! Return if a mesh will be isoparametric or not
  //! \retval isoparametric Status of mesh type
  bool is_isoparametric();

  //! Node entity sets
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] check Check duplicates
  void node_entity_sets(const Json& mesh_prop, bool check);

  //! Node Euler angles
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] mesh_io Mesh IO handle
  void node_euler_angles(const Json& mesh_prop,
                         const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io);

  //! Nodal velocity constraints
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] mesh_io Mesh IO handle
  void nodal_velocity_constraints(
      const Json& mesh_prop, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io);

  //! Nodal frictional constraints
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] mesh_io Mesh IO handle
  void nodal_frictional_constraints(
      const Json& mesh_prop, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io);

  //! Cell entity sets
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] check Check duplicates
  void cell_entity_sets(const Json& mesh_prop, bool check);

  //! Particles cells
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] particle_io Particle IO handle
  void particles_cells(const Json& mesh_prop,
                       const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io);

  //! Particles volumes
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] particle_io Particle IO handle
  void particles_volumes(const Json& mesh_prop,
                         const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io);

  //! Particle velocity constraints
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] particle_io Particle IO handle
  void particle_velocity_constraints(
      const Json& mesh_prop,
      const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io);

  //! Particles stresses
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] particle_io Particle IO handle
  void particles_stresses(
      const Json& mesh_prop,
      const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io);

  //! Nodal pore pressure constraints
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] mesh_io Mesh IO handle
  void nodal_pore_pressure_constraints(
      const Json& mesh_prop, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io);

  // Particles pore pressures
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] particle_io Particle IO handle
  void particles_pore_pressures(
      const Json& mesh_prop,
      const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io);

  // Particles pore pressure constraints
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] particle_io Particle IO handle
  void particle_pore_pressure_constraints(
      const Json& mesh_prop,
      const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io);  

  //! Particle entity sets
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] check Check duplicates
  void particle_entity_sets(const Json& mesh_prop, bool check);

  //! Particle material assignment by particle sets
  //! \param[in] mesh_prop Mesh properties
  void particle_materials(const Json& mesh_prop);

  //! Initialise damping
  //! \param[in] damping_props Damping properties
  bool initialise_damping(const Json& damping_props);

  //! Initialise nodal water table
  //! \param[in] water_table Water table properties
  bool initialise_nodal_water_table(const Json& water_table);

  //! Pressure smoothing
  //! \param[in] phase Phase to smooth pressure
  virtual void pressure_smoothing(unsigned phase) {
    throw std::runtime_error(
        "Calling the base class function (pressure_smoothing) in MPMBase:: "
        "illegal operation!");
  };
/////////////////////////////////////////////////////////////////////
//                     THERMAL PART                      ////////////
/////////////////////////////////////////////////////////////////////
  //! Nodal temperature constraints
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] mesh_io Mesh IO handle
  void nodal_temperature_constraints(
      const Json& mesh_prop, 
      const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io);

  //! Nodal_convective_heat_constraints
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] mesh_io Mesh IO handle
  void nodal_convective_heat_constraints(
      const Json& mesh_prop, 
      const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io);

  //! Nodal heat source
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] mesh_io Mesh IO handle
  void nodal_heat_source(
      const Json& mesh_prop, 
      const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io);

  // Particles temperature constraints
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] particle_io Particle IO handle
  void particle_temperature_constraints(
      const Json& mesh_prop,
      const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io);    

  // Particles temperatures
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] particle_io Particle IO handle
  void particles_temperatures(
      const Json& mesh_prop,
      const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io);     

 protected:
  // Generate a unique id for the analysis
  using mpm::MPM::uuid_;
  //! Time step size
  using mpm::MPM::dt_;
  //! Current step
  using mpm::MPM::step_;
  //! Number of steps
  using mpm::MPM::nsteps_;
  //! Output steps
  using mpm::MPM::output_steps_;
  //! Current_time
  using mpm::MPM::current_time_;
  //! A unique ptr to IO object
  using mpm::MPM::io_; 
  //! JSON analysis object
  using mpm::MPM::analysis_;
  //! JSON post-process object
  using mpm::MPM::post_process_;
  //! Logger
  using mpm::MPM::console_;

  //! Stress update method (default USF = 0, USL = 1, MUSL = 2)
  mpm::StressUpdate stress_update_{mpm::StressUpdate::USF};
  //! PIC value: PIC update(pic = 0.), FLIP(pic = 1.)
  double pic_{0.};
  //! PIC_T value: PIC_T update(PIC_T = 0.), FLIP(PIC_T = 1.)
  double pic_t_{0.};
  //! Use affine particle-grid velocity transfer for the solid phase
  bool affine_mpm_{false};
  //! Gravity
  Eigen::Matrix<double, Tdim, 1> gravity_;
  //! Mesh object
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
  //! Materials
  std::map<unsigned, std::shared_ptr<mpm::Material<Tdim>>> materials_;
  //! Mathematical functions
  std::map<unsigned, std::shared_ptr<mpm::FunctionBase>> math_functions_;
  //! VTK attributes for solid phase
  std::vector<std::string> vtk_attributes_;
  //! VTK state variables
  std::vector<std::string> vtk_statevars_;
  //! Set node concentrated force
  bool set_node_concentrated_force_{false};
  //! Damping type
  mpm::Damping damping_type_{mpm::Damping::None};
  //! Damping factor
  double damping_factor_{0.};
  //! Output hdf5
  bool write_hdf5_{false};
  //! Output vtk
  bool write_vtk_{false};

  // Time step matrix
  Eigen::VectorXd dt_matrix_;
  // Steps
  Eigen::VectorXd nsteps_matrix_;
  // Time step matrix size
  int dt_matrix_size;
  //is_axisymmetric
  bool is_axisymmetric_{false};

#ifdef USE_GRAPH_PARTITIONING
  // graph pass the address of the container of cell
  std::shared_ptr<Graph<Tdim>> graph_{nullptr};
#endif
};  // MPMBase class
}  // namespace mpm

#include "mpm_base.tcc"

#endif  // MPM_MPM_BASE_H_
