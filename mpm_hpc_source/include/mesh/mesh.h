#ifndef MPM_MESH_H_
#define MPM_MESH_H_

#include <iostream>
#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

#include <Eigen/Sparse>

// Eigen
#include "Eigen/Dense"
// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif
// TBB
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>

#include <tsl/robin_map.h>
// JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "cell.h"
#include "container.h"
#include "factory.h"
#include "friction_constraint.h"
#include "function_base.h"
#include "geometry.h"
#include "hdf5_particle.h"
#include "io.h"
#include "io_mesh.h"
#include "loads_bcs/friction_constraint.h"
#include "loads_bcs/traction.h"
#include "loads_bcs/contact.h"
#include "loads_bcs/heat_source.h"
#include "loads_bcs/velocity_constraint.h"
#include "loads_bcs/temperature_constraint.h"
#include "loads_bcs/pore_pressure_constraint.h"
#include "logger.h"
#include "material/material.h"
#include "mpi_datatypes.h"
#include "node.h"
#include "particle.h"
#include "particle_base.h"

namespace mpm {

//! Mesh class
//! \brief Base class that stores the information about meshes
//! \details Mesh class which stores the particles, nodes, cells and neighbours
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Mesh {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  // Construct a mesh with a global unique id
  //! \param[in] id Global mesh id
  //! \param[in] isoparametric Mesh is isoparametric
  Mesh(unsigned id, bool isoparametric = true);

  //! Default destructor
  ~Mesh() = default;

  //! Delete copy constructor
  Mesh(const Mesh<Tdim>&) = delete;

  //! Delete assignement operator
  Mesh& operator=(const Mesh<Tdim>&) = delete;

  //! Return id of the mesh
  unsigned id() const { return id_; }

  //! Return if a mesh is isoparametric
  bool is_isoparametric() const { return isoparametric_; }

  //! Create nodes from coordinates
  //! \param[in] gnid Global node id
  //! \param[in] node_type Node type
  //! \param[in] coordinates Nodal coordinates
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval status Create node status
  bool create_nodes(mpm::Index gnid, const std::string& node_type,
                    const std::vector<VectorDim>& coordinates,
                    bool check_duplicates = true);

  //! Add a node to the mesh
  //! \param[in] node A shared pointer to node
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval insertion_status Return the successful addition of a node
  bool add_node(const std::shared_ptr<mpm::NodeBase<Tdim>>& node,
                bool check_duplicates = true);

  //! Remove a node from the mesh
  //! \param[in] node A shared pointer to node
  //! \retval insertion_status Return the successful addition of a node
  bool remove_node(const std::shared_ptr<mpm::NodeBase<Tdim>>& node);

  //! Return the number of nodes
  mpm::Index nnodes() const { return nodes_.size(); }

  //! Iterate over nodes
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_nodes(Toper oper);

  //! Iterate over nodes with predicate
  //! \tparam Toper Callable object typically a baseclass functor
  //! \tparam Tpred Predicate
  template <typename Toper, typename Tpred>
  void iterate_over_nodes_predicate(Toper oper, Tpred pred);

  //! Create a list of active nodes in mesh
  void find_active_nodes();

  //! Iterate over active nodes
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_active_nodes(Toper oper);

#ifdef USE_MPI
  //! All reduce over nodal property
  //! \tparam Ttype Type of property to accumulate
  //! \tparam Tnparam Size of individual property
  //! \tparam Tgetfunctor Functor for getter
  //! \tparam Tsetfunctor Functor for setter
  //! \param[in] getter Getter function
  template <typename Ttype, unsigned Tnparam, typename Tgetfunctor,
            typename Tsetfunctor>
  void nodal_halo_exchange(Tgetfunctor getter, Tsetfunctor setter);
#endif

  //! Create cells from list of nodes
  //! \param[in] gcid Global cell id
  //! \param[in] element Element type
  //! \param[in] cells Node ids of cells
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval status Create cells status
  bool create_cells(mpm::Index gnid,
                    const std::shared_ptr<mpm::Element<Tdim>>& element,
                    const std::vector<std::vector<mpm::Index>>& cells,
                    bool check_duplicates = true);

  //! Add a cell from the mesh
  //! \param[in] cell A shared pointer to cell
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval insertion_status Return the successful addition of a cell
  bool add_cell(const std::shared_ptr<mpm::Cell<Tdim>>& cell,
                bool check_duplicates = true);

  //! Remove a cell from the mesh
  //! \param[in] cell A shared pointer to cell
  //! \retval insertion_status Return the successful addition of a cell
  bool remove_cell(const std::shared_ptr<mpm::Cell<Tdim>>& cell);

  //! Number of cells in the mesh
  mpm::Index ncells() const { return cells_.size(); }

  //! Iterate over cells
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_cells(Toper oper);

  //! Create particles from coordinates
  //! \param[in] particle_type Particle type
  //! \param[in] coordinates Nodal coordinates
  //! \param[in] material_id ID of the material
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \param[in] liquid_material_id ID of the liquid material
  //! \retval status Create particle status
  bool create_particles(
      const std::string& particle_type,
      const std::vector<VectorDim>& coordinates, unsigned material_id,
      bool check_duplicates = true,
      unsigned liquid_material_id = std::numeric_limits<unsigned>::max());

  //! Add a particle to the mesh
  //! \param[in] particle A shared pointer to particle
  //! \param[in] checks Parameter to check duplicates and addition
  //! \retval insertion_status Return the successful addition of a particle
  bool add_particle(const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle,
                    bool checks = true);

  //! Remove a particle from the mesh
  //! \param[in] particle A shared pointer to particle
  //! \retval insertion_status Return the successful addition of a particle
  bool remove_particle(
      const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle);

  //! Remove a particle by id
  bool remove_particle_by_id(mpm::Index id);

  //! Remove a particle from the mesh
  //! \param[in] pids Vector of particle ids
  void remove_particles(const std::vector<mpm::Index>& pids);

  //! Remove all particles in a cell in nonlocal rank
  void remove_all_nonrank_particles();

  //! Transfer particles to different ranks in nonlocal rank cells
  void transfer_nonrank_particles();

  //! Find shared nodes across MPI domains in the mesh
  void find_domain_shared_nodes();

  //! Number of particles in the mesh
  mpm::Index nparticles() const { return particles_.size(); }

  //! Locate particles in a cell
  //! Iterate over all cells in a mesh to find the cell in which particles
  //! are located.
  //! \retval particles Particles which cannot be located in the mesh
  std::vector<std::shared_ptr<mpm::ParticleBase<Tdim>>> locate_particles_mesh();

  //! Iterate over particles
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_particles(Toper oper);

  //! Iterate over particle set
  //! \tparam Toper Callable object typically a baseclass functor
  //! \param[in] set_id particle set id
  template <typename Toper>
  void iterate_over_particle_set(int set_id, Toper oper);

  //! Iterate over node set
  //! \tparam Toper Callable object typically a baseclass functor
  //! \param[in] set_id node set id
  template <typename Toper>
  void iterate_over_node_set(int set_id, Toper oper);

  //! Return coordinates of particles
  std::vector<Eigen::Matrix<double, 3, 1>> particle_coordinates();

  //! Return particles tensor data
  //! \param[in] attribute Name of the tensor data attribute
  template <unsigned Tsize>
  std::vector<Eigen::Matrix<double, Tsize, 1>> particles_vector_data(
      const std::string& attribute);

  //! Return particles tensor data
  //! \param[in] attribute Name of the tensor data attribute
  std::vector<double> particles_scalar_data(const std::string& attribute);

  //! Return particles tensor data
  //! \param[in] attribute Name of the tensor data attribute
  template <unsigned Tsize>
  std::vector<Eigen::Matrix<double, Tsize, 1>> nodal_vector_data(
      const std::string& attribute);

  //! Return particles tensor data
  //! \param[in] attribute Name of the tensor data attribute
  std::vector<double> nodal_scalar_data(const std::string& attribute);

  //! Return particles scalar data
  //! \param[in] attribute Name of the scalar data attribute
  std::vector<double> particles_statevars_data(const std::string& attribute);

  //! Compute and assign rotation matrix to nodes
  //! \param[in] euler_angles Map of node number and respective euler_angles
  bool compute_nodal_rotation_matrices(
      const std::map<mpm::Index, Eigen::Matrix<double, Tdim, 1>>& euler_angles);

  //! Assign particles volumes
  //! \param[in] particle_volumes Volume at dir on particle
//   bool assign_particles_volumes(
//       const std::vector<std::tuple<mpm::Index, double>>& particle_volumes);
  bool assign_particles_volumes(
      const std::vector<double>& particle_volumes);


  //! Create particles tractions
  //! \param[in] mfunction Math function if defined
  //! \param[in] setid Particle set id
  //! \param[in] dir Direction of traction load
  //! \param[in] traction Particle traction
  bool create_particles_tractions(
      const std::shared_ptr<FunctionBase>& mfunction, int set_id,
      unsigned facet, unsigned dir, double traction);

  //! Apply traction to particles
  //! \param[in] current_time Current time
  void apply_traction_on_particles(double current_time);

  //! Create particles contacts
  bool create_particles_contacts(int set_id, unsigned dir, double normal);

  //! Apply contact to particles
  void apply_contact_on_particles();

  //! Create particle velocity constraints 
  //! \param[in] setid Node set id
  bool create_particle_velocity_constraint(
      int set_id, const std::shared_ptr<mpm::VelocityConstraint>& constraint);

  //! Create particle velocity constraints
  //! \param[in] setid Node set id
  bool create_nodal_velocity_constraint(
      int set_id, const std::shared_ptr<mpm::VelocityConstraint>& constraint);

  //! Apply particles and nodal velocity constraints
  void apply_velocity_constraints(double current_time);

  //! Assign rigid particle velocity to nodes
  void apply_moving_rigid_boundary(double current_time, double dt);

  // //! Assign nodal pressure reference step
  // //! \param[in] setid Node set id
  // //! \param[in] ref_step Reference step
  // bool assign_nodal_pressure_reference_step(int set_id, const Index ref_step);

  //! Assign nodal velocity constraints
  //! \param[in] mfunction Math function if defined
  //! \param[in] setid Node set id
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] pconstraint Pressure constraint at node
  bool assign_nodal_pressure_constraint(
      const std::shared_ptr<FunctionBase>& mfunction, int set_id,
      const unsigned phase, double pconstraint);

  //! Assign nodal velocity constraints
  //! \param[in] setid Node set id
  //! \param[in] velocity_constraints Velocity constraint at node, dir, velocity
  bool assign_nodal_velocity_constraint(
      int set_id, const std::shared_ptr<mpm::VelocityConstraint>& constraint);

  //! Assign nodal frictional constraints
  //! \param[in] setid Node set id
  //! \param[in] friction_constraints Constraint at node, dir, sign, friction
  bool assign_nodal_frictional_constraint(
      int nset_id,
      const std::shared_ptr<mpm::FrictionConstraint>& fconstraints);

  //! Assign velocity constraints to nodes
  //! \param[in] velocity_constraints Constraint at node, dir, and velocity
  bool assign_nodal_velocity_constraints(
      const std::vector<std::tuple<mpm::Index, unsigned, double>>&
          velocity_constraints);

  //! Assign friction constraints to nodes
  //! \param[in] friction_constraints Constraint at node, dir, sign, and
  //! friction
  bool assign_nodal_friction_constraints(
      const std::vector<std::tuple<mpm::Index, unsigned, int, double>>&
          friction_constraints);

  // //! Assign nodal concentrated force
  // //! \param[in] nodal_forces Force at dir on nodes
  // bool assign_nodal_concentrated_forces(
  //     const std::vector<std::tuple<mpm::Index, unsigned, double>>&
  //         nodal_forces);

  // //! Assign nodal concentrated force
  // //! \param[in] mfunction Math function if defined
  // //! \param[in] setid Node set id
  // //! \param[in] dir Direction of force
  // //! \param[in] node_forces Concentrated force at dir on nodes
  // bool assign_nodal_concentrated_forces(
  //     const std::shared_ptr<FunctionBase>& mfunction, int set_id, unsigned dir,
  //     double force);

  //! Assign particles stresses
  //! \param[in] particle_stresses Initial stresses of particle
  bool assign_particles_stresses(
      const std::vector<Eigen::Matrix<double, 6, 1>>& particle_stresses);

  //! Assign particles cells
  //! \param[in] particles_cells Particles and cells
  bool assign_particles_cells(
      const std::vector<std::array<mpm::Index, 2>>& particles_cells);

  //! Return particles cells
  //! \retval particles_cells Particles and cells
  std::vector<std::array<mpm::Index, 2>> particles_cells() const;

  //! Return status of the mesh. A mesh is active, if at least one particle is
  //! present
  bool status() const { return particles_.size(); }

  //! Generate points
  //! \param[in] nquadratures Number of points per direction in cell
  //! \param[in] particle_type Particle type
  //! \param[in] material_id ID of the material
  //! \param[in] cset_id Set ID of the cell [-1 for all cells]
  //! \param[in] liquid_material_id ID of the liquid material
  //! \retval point Material point coordinates
  bool generate_material_points(
      unsigned nquadratures, const std::string& particle_type,
      unsigned material_id, int cset_id,
      unsigned liquid_material_id = std::numeric_limits<unsigned>::max());

  //! Initialise material models
  //! \param[in] materials Material models
  void initialise_material_models(
      const std::map<unsigned, std::shared_ptr<mpm::Material<Tdim>>>&
          materials) {
    materials_ = materials;
  }

  //! Find cell neighbours
  void compute_cell_neighbours();

  //! Find particle neighbours for all particle
  void compute_particle_neighbours();

  //! Find particle neighbours
  //! \param[in] cell of interest
  void compute_particle_neighbours(const Cell<Tdim>& cell);

  //! Find particle neighbours
  //! \param[in] particle of interest
  void compute_particle_neighbours(const ParticleBase<Tdim>& particle);

  //! Add a neighbour mesh, using the local id for the new mesh and a mesh
  //! pointer
  //! \param[in] local_id local id of the mesh
  //! \param[in] neighbour A shared pointer to the neighbouring mesh
  //! \retval insertion_status Return the successful addition of a node
  bool add_neighbour(unsigned local_id,
                     const std::shared_ptr<Mesh<Tdim>>& neighbour);

  //! Return the number of neighbouring meshes
  unsigned nneighbours() const { return neighbour_meshes_.size(); }

  //! Find ghost boundary cells
  void find_ghost_boundary_cells();

  //! Write HDF5 particles
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] filename Name of HDF5 file to write particles data
  //! \retval status Status of writing HDF5 output
  bool write_particles_hdf5(unsigned phase, const std::string& filename);

  //! Read HDF5 particles
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] filename Name of HDF5 file to write particles data
  //! \retval status Status of reading HDF5 output
  bool read_particles_hdf5(unsigned phase, const std::string& filename);

  //! Return HDF5 particles
  //! \retval particles_hdf5 Vector of HDF5 particles
  std::vector<mpm::HDF5Particle> particles_hdf5();

  //! Return nodal coordinates
  std::vector<Eigen::Matrix<double, 3, 1>> nodal_coordinates() const;

  //! Return node pairs
  std::vector<std::array<mpm::Index, 2>> node_pairs() const;

  //! Create map of container of particles in sets
  //! \param[in] map of particles ids in sets
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval status Status of  create particle sets
  bool create_particle_sets(
      const tsl::robin_map<mpm::Index, std::vector<mpm::Index>>& particle_sets,
      bool check_duplicates);

  //! Assign a material to particles in a particle set
  //! \param[in] set_id Particle set id
  //! \param[in] material_id Material id
  //! \retval status Status of assigning materials
  bool assign_material_to_particle_set(int set_id, unsigned material_id);

  //! Create map of container of nodes in sets
  //! \param[in] map of nodes ids in sets
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval status Status of  create node sets
  bool create_node_sets(
      const tsl::robin_map<mpm::Index, std::vector<mpm::Index>>& node_sets,
      bool check_duplicates);

  //! Create map of container of cells in sets
  //! \param[in] map of cells ids in sets
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval status Status of create cell sets
  bool create_cell_sets(
      const tsl::robin_map<mpm::Index, std::vector<mpm::Index>>& cell_sets,
      bool check_duplicates);

  //! Get the container of cell
  mpm::Container<Cell<Tdim>> cells();

  //! Return particle cell ids
  std::map<mpm::Index, mpm::Index>* particles_cell_ids();

  //! Return nghost cells
  unsigned nghost_cells() const { return ghost_cells_.size(); }

  //! Return nlocal ghost cells
  unsigned nlocal_ghost_cells() const { return local_ghost_cells_.size(); }

  //! Generate particles
  //! \param[in] io IO object handle
  //! \param[in] generator Point generator object
  bool generate_particles(const std::shared_ptr<mpm::IO>& io,
                          const Json& generator);
  //--------------------------------------------------
  //! TWO-PHASE functions

  //! Assign velocity constraints to cells
  //! \param[in] velocity_constraints Constraint at cell id, face id, dir, and
  //! velocity
  bool assign_cell_velocity_constraints(
      const std::vector<std::tuple<mpm::Index, unsigned, unsigned, double>>&
          velocity_constraints);

  //! Assign nodal water table
  //! \param[in] wfunction Math function if defined
  //! \param[in] set_id Node set id
  //! \param[in] dir Water table direction
  //! \param[in] h0 Zero pore pressure height
  bool assign_nodal_water_table(const std::shared_ptr<FunctionBase>& wfunction,
                                const int set_id, const unsigned dir,
                                const double h0);

  //! Assign nodal pressure constraints to nodes
  //! \param[in] pressure_constraints Constraint at node, pressure
  bool assign_nodal_pressure_constraints(
      const unsigned phase,
      const std::vector<std::tuple<mpm::Index, double>>& pressure_constraints);

  //! Assign particles pore pressures
  //! \param[in] particle_pore_pressure Initial pore pressure of particle
  bool assign_particles_pore_pressures(
      const std::vector<double>& particle_pore_pressures);

  //! Create particle pore pressure constraints
  bool create_particle_pore_pressure_constraint(
    int set_id, const std::shared_ptr<mpm::PorepressureConstraint>& pconstraint);

  //! Apply particles pore pressure constraints
  void apply_particle_pore_pressure_constraints(double current_time);

  //! Get global node indices
  //! \retval global_node_indices Global node indices
  std::vector<Eigen::VectorXi> global_node_indices() const;

  bool assign_nodal_K_cor(const unsigned dim, Eigen::VectorXd& force_cor);

  bool compute_nodal_corrected_force(Eigen::SparseMatrix<double>& K_cor_matrix,
                                     Eigen::VectorXd& pore_pressure_increment,
                                     double dt);

  mpm::Container<NodeBase<Tdim>> nodes() { return nodes_; }

  //! Compute free surface
  //! \param[in] map_nodal_volume_for_density Map cell volume when detection
  //! requires nodal density. Set false when the solver already mapped cell
  //! volume during the current nodal initialisation cycle.
  bool compute_free_surface(
      std::string free_surface_particle,
      double tolerance = std::numeric_limits<unsigned>::epsilon(),
      bool map_nodal_volume_for_density = true);

  //! Get free surface node set
  std::set<mpm::Index> free_surface_nodes();

  //! Get free surface cell set
  std::set<mpm::Index> free_surface_cells();

  //! Get free surface particle set
  std::set<mpm::Index> free_surface_particles();

  //! Assign id for active node
  unsigned assign_active_node_id();

  //! Return container of active nodes
  mpm::Container<NodeBase<Tdim>> active_nodes() { return active_nodes_; }

  //! Assign initial free surface particles
  bool assign_free_surface_particles(const std::shared_ptr<mpm::IO>& io);

  //! Assign initial nonfree surface particles
  bool assign_nonfree_surface_particles(const std::shared_ptr<mpm::IO>& io);

  //! Get deformation gradient
  void get_displacement_gradient(
      std::vector<unsigned>& id,
      std::vector<Eigen::MatrixXd>& displacement_gradients);

  //! Get particle temperature
  void get_particle_temperature(
      std::vector<unsigned>& id,
      std::vector<double>& particle_temperature);      

  //! Get particle temperature
  void get_hydrate_saturation(
      std::vector<unsigned>& id,
      std::vector<double>& hydrate_saturation);      


  //! Get total reaction force for rigid material point (material_id =999)
  void get_reaction_force(Eigen::Matrix<double, Tdim, 1>& disp,
                          Eigen::Matrix<double, Tdim, 1>& reaction_force);

 private:
  // Read particles from file
  bool read_particles_file(const std::shared_ptr<mpm::IO>& io,
                           const Json& generator);

  // Locate a particle in mesh cells
  bool locate_particle_cells(
      const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle);

 public:
/////////////////////////////////////////////////////////////////////
//                     THERMAL PART                      ////////////
/////////////////////////////////////////////////////////////////////
//! Assign particles temperatures
  //! \param[in] particle_temperature Initial temperature of particle
  bool assign_particles_temperatures(
      const std::vector<double>& particle_temperature);

  //! Create particle temperature constraints
  bool create_particle_temperature_constraint(
    int set_id, const std::shared_ptr<mpm::TemperatureConstraint>& Tconstraint);


  //! Assign nodal temperature constraint
  //! \param[in] mfunction Math function if defined
  //! \param[in] setid Node set id
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] Tconstraint temperature constraint at node
  bool assign_nodal_temperature_constraint(
      const std::shared_ptr<FunctionBase>& mfunction, int set_id,
      const unsigned phase, const double Tconstraint); 

  //! Assign nodal temperature constraint
  bool assign_nodal_convective_heat_constraint(
      const std::shared_ptr<FunctionBase>& mfunction, int set_id,
      const unsigned phase, const double Tconstraint, const double coeff);

  //! Assign nodal temperature constraints to nodes
  //! \param[in] temperature_constraints Constraint at node, pressure
   bool assign_nodal_temperature_constraints(
      const unsigned phase,
      const std::vector<std::tuple<mpm::Index, double>>& temperature_constraints); 

  //! Apply particles temperature constraints
  void apply_particle_temperature_constraints(double current_time);

  //! Apply nodal temperature constraints
  void apply_nodal_temperature_constraints(unsigned phase, double current_time);

  //! Apply nodal convective heat constraints
  void apply_nodal_convective_heat_constraints(unsigned phase, double current_time);

  //! Create particles heat source
  //! \param[in] mfunction Math function if defined
  //! \param[in] setid Particle set id
  //! \param[in] heat_source Particle heat_source
  bool create_particles_heat_sources(
      const std::shared_ptr<FunctionBase>& mfunction, int set_id, double heat_source);

  //! Apply heat source to particles
  //! \param[in] current_time Current time
  void apply_heat_source_on_particles(double current_time, double dt);
  
  //! Assign nodal heat source
  //! \param[in] mfunction Math function if defined
  //! \param[in] setid Node set id
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] heat_source temperature constraint at node
  bool assign_nodal_heat_source(
      const std::shared_ptr<FunctionBase>& mfunction, int set_id,
      const unsigned phase, const unsigned heat_source); 

  //! Assign nodal temperature constraints to nodes
  //! \param[in] heat_source Constraint at node, pressure
   bool assign_nodal_heat_sources(
      const unsigned phase,
      const std::vector<std::tuple<mpm::Index, double>>& heat_sources); 

  //! Apply heat source on nodes
  void apply_heat_source_on_nodes(unsigned phase, double current_time);


  private:
  //! mesh id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  // TODO: Optimize this and may be remove it later
  //! Isoparametric mesh
  bool isoparametric_{true};
  //! Container of mesh neighbours
  Map<Mesh<Tdim>> neighbour_meshes_;
  //! Container of particles
  Container<ParticleBase<Tdim>> particles_;
  //! Container of particles ids and cell ids
  std::map<mpm::Index, mpm::Index> particles_cell_ids_;
  //! Container of particle sets
  tsl::robin_map<unsigned, tbb::concurrent_vector<mpm::Index>> particle_sets_;
  //! Material id overrides by particle id
  tsl::robin_map<mpm::Index, unsigned> particle_material_ids_;
  //! Map of particles for fast retrieval
  Map<ParticleBase<Tdim>> map_particles_;
  //! Container of nodes
  Container<NodeBase<Tdim>> nodes_;
  //! Container of domain shared nodes
  Container<NodeBase<Tdim>> domain_shared_nodes_;
  //! Boundary nodes
  Container<NodeBase<Tdim>> boundary_nodes_;
  //! Container of node sets
  tsl::robin_map<unsigned, Container<NodeBase<Tdim>>> node_sets_;
  //! Container of active nodes
  Container<NodeBase<Tdim>> active_nodes_;
  //! Map of nodes for fast retrieval
  Map<NodeBase<Tdim>> map_nodes_;
  //! Map of cells for fast retrieval
  Map<Cell<Tdim>> map_cells_;
  //! Container of cells
  Container<Cell<Tdim>> cells_;
  //! Container of ghost cells sharing the current MPI rank
  Container<Cell<Tdim>> ghost_cells_;
  //! Container of local ghost cells
  Container<Cell<Tdim>> local_ghost_cells_;
  //! Container of cell sets
  tsl::robin_map<unsigned, Container<Cell<Tdim>>> cell_sets_;
  //! Map of ghost cells to the neighbours ranks
  std::map<unsigned, std::vector<unsigned>> ghost_cells_neighbour_ranks_;
  //! Faces and cells
  std::multimap<std::vector<mpm::Index>, mpm::Index> faces_cells_;
  //! Materials
  std::map<unsigned, std::shared_ptr<mpm::Material<Tdim>>> materials_;
  //! Loading (Particle tractions)
  std::vector<std::shared_ptr<mpm::Traction>> particle_tractions_;
  //! Particle contacts)
  std::vector<std::shared_ptr<mpm::Contact>> particle_contacts_;
  //! Loading (Particle heat sources)
  std::vector<std::shared_ptr<mpm::Heat_source>> particle_heat_sources_;
  //! Particle velocity constraints
  std::vector<std::shared_ptr<mpm::VelocityConstraint>>
      particle_velocity_constraints_;
  //! Nodal velocity constraints
  std::vector<std::shared_ptr<mpm::VelocityConstraint>>
      nodal_velocity_constraints_;
  //! Particle temperature constraints
  std::vector<std::shared_ptr<mpm::TemperatureConstraint>>
      particle_temperature_constraints_;  
  //! Particle pore pressure constraints 
  std::vector<std::shared_ptr<mpm::PorepressureConstraint>>
      particle_pore_pressure_constraints_;
      
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
  //! TBB grain size
  int tbb_grain_size_{100};
  //! Maximum number of halo nodes
  unsigned nhalo_nodes_{0};
  //! Maximum number of halo nodes
  unsigned ncomms_{0};
};  // Mesh class
}  // namespace mpm

#include "mesh.tcc"

#endif  // MPM_MESH_H_
