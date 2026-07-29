#ifndef MPM_ASSEMBLER_BASE_H_
#define MPM_ASSEMBLER_BASE_H_

#include <algorithm>
#include <array>
#include <memory>
#include <utility>
#include <vector>

#include <Eigen/Sparse>

#include "mesh.h"
#include "node_base.h"

namespace mpm {

// Matrix assembler base class
//! \brief Assemble matrixs (stiffness matrix)
//! \details Get local stiffness matrix and node ids from cell
//! \tparam Tdim Dimension
template <unsigned Tdim>
class AssemblerBase {
 public:
  AssemblerBase() {
    //! Global degrees of freedom
    active_dof_ = 0;
  }
  // Virtual destructor
  virtual ~AssemblerBase() = default;

  //! Copy constructor
  // AssemblerBase(const AssemblerBase<Tdim>&) = default;

  //! Assignment operator
  // AssemblerBase& operator=(const AssemblerBase<Tdim>&) = default;

  //! Move constructor
  // AssemblerBase(AssemblerBase<Tdim>&&) = default;

  //! Assign mesh pointer
  void assign_mesh_pointer(std::shared_ptr<mpm::Mesh<Tdim>>& mesh) {
    mesh_ = mesh;
  }

  //! Create a pair between nodes and index in Matrix / Vector
  virtual bool assign_global_node_indices(unsigned active_dof) = 0;

  //! Assemble stiffness matrix (semi-implicit)
  virtual bool assemble_stiffness_matrix(unsigned dir, 
                                          double dt,
                                          bool implicit_drag_force, 
                                          int entries_number) = 0;

  //! Assemble force vector (semi-implicit)
  virtual bool assemble_force_vector(double dt) = 0;

  // //! Assemble force vector (semi-implicit)
  // virtual bool assemble_inter_force_vector(double dt) = 0;  

  //! Assemble displacement vector
  // virtual void assemble_disp_vector() = 0;

  //! Apply displacement to nodes
  // virtual void apply_displacement_nodes() = 0;

  //! Apply forces to nodes
  // virtual void apply_forces_nodes() = 0;

  //! Apply restraints
  // virtual Eigen::VectorXd global_restraints() = 0;

  //! Initialise force vector to zero
  // virtual void initialise_force_zero() = 0;

  virtual bool assemble_K_cor_matrix(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                                     double dt, int entries_number) = 0;

  //! Assemble laplacian matrix
  virtual bool assemble_laplacian_matrix(double dt, int entries_number) = 0;

  //! Assemble laplacian matrix
  virtual bool assemble_stab_matrix(double dt) = 0;  

  //! Assemble poisson right
  virtual bool assemble_poisson_right(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                                      double dt, int entries_number) = 0;
                                  
  //! Assemble poisson right - thermal part
  virtual bool assemble_poisson_right_thermal(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                                      double dt, int entries_number) = 0;                                                                       

  //! Assemble poisson right - thermal part
  virtual bool assemble_poisson_right_pressure(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                                      double dt, int entries_number) = 0;

  //! Return stiffness matrix
  virtual Eigen::SparseMatrix<double>& stiffness_matrix(unsigned dir) = 0;

  //! Return intermediate force
  virtual Eigen::MatrixXd& force_inter() = 0;
  // virtual Eigen::MatrixXd& force_cor() = 0;  

  virtual void assign_intermediate_acceleration(
      unsigned dim, Eigen::VectorXd acceleration_inter) = 0;     

  virtual void assign_pore_pressure_increment(
      Eigen::VectorXd pore_pressure_increment) = 0;

  //! Return intermediate velocity
  virtual Eigen::MatrixXd& acceleration_inter() = 0;

  virtual Eigen::SparseMatrix<double>& K_cor_matrix() = 0;

  virtual Eigen::SparseMatrix<double>& laplacian_matrix() = 0;

  virtual Eigen::SparseMatrix<double>& stab_matrix() = 0;  

  virtual Eigen::VectorXd& force_laplacian_matrix() = 0;

  virtual Eigen::VectorXd& pore_pressure_increment() = 0;

  virtual void apply_pressure_constraints() = 0;

  virtual bool apply_velocity_constraints() = 0; 

  virtual bool assign_pressure_constraints(double beta,
                                           const double current_time) = 0;

  virtual bool assign_velocity_constraints() = 0;

  virtual std::set<mpm::Index> free_surface() = 0;

  virtual void assign_free_surface(
      const std::set<mpm::Index>& free_surface_id) = 0;

  virtual unsigned active_dof() { return active_dof_; };

  //! Assemble_KTT_matrix
  virtual bool assemble_KTT_matrix(double dt) = 0; 

  //! Assemble_MTT_matrix
  virtual bool assemble_MTT_matrix(double dt) = 0; 

  //! Assemble_FT_vector
  virtual bool assemble_FT_vector(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_, double dt) = 0;  

  //! Assign_temperature_constraint
  virtual bool assign_temperature_constraints(const double current_time) = 0; 

  //! Apply_temperature_constraint
  virtual void apply_temperature_constraints() = 0; 

  virtual Eigen::SparseMatrix<double>& MTT_matrix() = 0;  

  virtual Eigen::SparseMatrix<double>& KTT_matrix() = 0;    

  virtual void assign_temperature_rate(Eigen::VectorXd temperature_rate) = 0;  

  virtual Eigen::VectorXd& FT_vector() = 0;

  virtual Eigen::VectorXd& temperature_rate() = 0;

  virtual Eigen::VectorXd& temperature() = 0;

 protected:
  //! Active node number
  unsigned active_dof_;
  //! Mesh object
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
  //! Sparse Stiffness Matrix
  std::shared_ptr<Eigen::SparseMatrix<double>> stiffness_matrix_;
  //! Force vector
  std::shared_ptr<Eigen::MatrixXd> force_inter_;
  // std::shared_ptr<Eigen::MatrixXd> force_cor_;  
  //! Intermediate velocity vector
  std::shared_ptr<Eigen::VectorXd> velocity_inter_vector_;
  //! Displacement vector
  std::shared_ptr<Eigen::MatrixXd> displacement_vector_;
};
}  // namespace mpm

#endif  // MPM_ASSEMBLER_BASE_H_