#ifndef MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_TWOPHASE_H_
#define MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_TWOPHASE_H_

#include <Eigen/Sparse>
#include <string>
#include <bits/stdc++.h>

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// Speed log
#include "assembler_base.h"
#include "spdlog/spdlog.h"

#include "cg_eigen.h"
#include "mesh.h"

namespace mpm {
template <unsigned Tdim>
class AssemblerEigenSemiImplicitTwoPhase : public AssemblerBase<Tdim> {
 public:
  //! Constructor
  AssemblerEigenSemiImplicitTwoPhase();

  //! Create a pair between nodes and index in Matrix / Vector
  bool assign_global_node_indices(unsigned active_dof) override;

  //! Assemble stiffness matrix (semi-implicit)
  bool assemble_stiffness_matrix(unsigned dir, double dt, bool implicit_drag_force, int entries_number) override;

  //! Assemble force vector (semi-implicit)
  bool assemble_force_vector(double dt) override;

  //! Resize containers of matrix
  bool resize_semi_implicit_matrix();

  //! Assemble K_cor matrix (used in correcting nodal velocity)
  bool assemble_K_cor_matrix(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                             double dt, int entries_number) override;

  //! Assemble_laplacian_matrix
  bool assemble_laplacian_matrix(double dt, int entries_number) override;

  //! Assemble_laplacian_matrix
  bool assemble_stab_matrix(double dt) override;   

  //! Assemble_poisson_right
  bool assemble_poisson_right(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                              double dt, int entries_number) override;                             
                              
  //! Assemble_poisson_right - thermal part
  bool assemble_poisson_right_thermal(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                              double dt, int entries_number) override;                           

  //! Assemble_poisson_right - thermal part
  bool assemble_poisson_right_pressure(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                              double dt, int entries_number) override;                                

  //! Return stiffness matrix
  Eigen::SparseMatrix<double>& stiffness_matrix(unsigned dir) override {
    return stiffness_matrix_.at(dir);
  }
  //! Return intermediate force
  Eigen::MatrixXd& force_inter() override { return force_inter_; } 

  //! Assign inter
  void assign_intermediate_acceleration(
      unsigned dir, Eigen::VectorXd acceleration_inter) override {
    acceleration_inter_.col(dir) = acceleration_inter;
  } 

  void assign_pore_pressure_increment(
      Eigen::VectorXd pore_pressure_increment) override {
    pore_pressure_increment_ = pore_pressure_increment;
  }

  //! Return intermediate acceleration
  Eigen::MatrixXd& acceleration_inter() override { return acceleration_inter_; }

  Eigen::SparseMatrix<double>& K_cor_matrix() override { return K_cor_matrix_; }

  Eigen::SparseMatrix<double>& laplacian_matrix() override {
    return laplacian_matrix_;
  }

  Eigen::SparseMatrix<double>& stab_matrix() override {
    return stab_matrix_;
  }  

  Eigen::VectorXd& force_laplacian_matrix() override {
    return force_laplacian_matrix_;
  }

  Eigen::VectorXd& pore_pressure_increment() override {
    return pore_pressure_increment_;
  }

  std::set<mpm::Index> free_surface() override { return free_surface_; }

  void assign_free_surface(
      const std::set<mpm::Index>& free_surface_id) override {
    free_surface_ = free_surface_id;
  }

  bool assign_velocity_constraints() override;

  bool assign_pressure_constraints(double beta,
                                   const double current_time) override;

  bool apply_velocity_constraints() override;

  void apply_pressure_constraints();

  //! Assemble_KTT_matrix
  bool assemble_KTT_matrix(double dt) override;

  //! Assemble_MTT_matrix
  bool assemble_MTT_matrix(double dt) override;

  //! Assemble_FT_vector
  bool assemble_FT_vector(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_, double dt) override;  

  //! Assemble_temperature_constraint
  bool assign_temperature_constraints(const double current_time) override;  

  //! Apply_temperature_constraint
  void apply_temperature_constraints() override;  

  Eigen::SparseMatrix<double>& MTT_matrix() override {
    return MTT_matrix_;
  }

  Eigen::SparseMatrix<double>& KTT_matrix() override {
    return KTT_matrix_;
  }

  Eigen::VectorXd& FT_vector() override {
    return FT_vector_;
  }

  void assign_temperature_rate(
      Eigen::VectorXd temperature_rate) override {
    temperature_rate_ = temperature_rate;
  }

  Eigen::VectorXd& temperature_rate() override {
    return temperature_rate_;
  }

  Eigen::VectorXd& temperature() override {
    return temperature_;
  }

 protected:
  //! Logger
  std::shared_ptr<spdlog::logger> console_;

 private:
  //! number of nodes
  using AssemblerBase<Tdim>::active_dof_;
  //! Mesh object
  using AssemblerBase<Tdim>::mesh_;
  //! Stiffness_matrix
  std::map<unsigned, Eigen::SparseMatrix<double>> stiffness_matrix_;
  //! Force vector
  Eigen::MatrixXd force_inter_;
  //! Intermediate acceleration vector (each column represent one direction)
  Eigen::MatrixXd acceleration_inter_;  
  //! poisson equation RHS (F31 and F32)
  Eigen::VectorXd force_laplacian_matrix_;
  //! Laplacian matrix
  Eigen::SparseMatrix<double> laplacian_matrix_;
  //! Stab matrix
  Eigen::SparseMatrix<double> stab_matrix_;  
  //! K_cor_matrix
  Eigen::SparseMatrix<double> K_cor_matrix_;
  // Eigen::SparseMatrix<double> K_cor_matrix_res_;  
  //! p^(t+1) - beta * p^(t)
  Eigen::VectorXd pore_pressure_increment_;
  //! Laplacian coefficient
  Eigen::VectorXd poisson_right_vector_;
  //! Solver base
  std::shared_ptr<mpm::SolverBase<Tdim>> solver_;
  //! Global node indices
  std::vector<Eigen::VectorXi> global_node_indices_;
  //! Velocity constraints
  Eigen::SparseMatrix<double> velocity_constraints_;
  //! Pressure constraints
  Eigen::SparseVector<double> pressure_constraints_;
  //! Free surface
  std::set<mpm::Index> free_surface_;

  //! MTT matrix
  Eigen::SparseMatrix<double> MTT_matrix_; 
  //! KTT matrix
  Eigen::SparseMatrix<double> KTT_matrix_;   
  //! FT vector
  Eigen::VectorXd FT_vector_;
  //! T_rate vector
  Eigen::VectorXd temperature_rate_;  
  //! T vector
  Eigen::VectorXd temperature_; 
  //! Temperature constraints
  Eigen::SparseVector<double> temperature_constraints_;

};
}  // namespace mpm

#include "assembler_eigen_semi_implicit_twophase.tcc"
#endif  // MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_TWOPHASE_H_
