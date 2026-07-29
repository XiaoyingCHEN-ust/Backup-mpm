//! Construct a semi-implicit eigen matrix assembler
template <unsigned Tdim>
mpm::AssemblerEigenSemiImplicitTwoPhase<
    Tdim>::AssemblerEigenSemiImplicitTwoPhase()
    : mpm::AssemblerBase<Tdim>() {
  //! Logger
  console_ = spdlog::stdout_color_mt("AssemblerEigenSemiImplicitTwoPhase");
}

//! Assign global node indices
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assign_global_node_indices(
    unsigned active_dof) {
  bool status = true;
  try {
    active_dof_ = active_dof;

    global_node_indices_ = mesh_->global_node_indices();

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble stiffness matrix (semi-implicit)
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_stiffness_matrix(
    unsigned dir, double dt, bool implicit_drag_force, int entries_number) {
  bool status = true;
  try {

    stiffness_matrix_.erase(dir);
    // Initialise stiffness_matrix
    Eigen::SparseMatrix<double> stiffness_matrix;
    // Resize stiffness matrix
    stiffness_matrix.resize(2 * active_dof_, 2 * active_dof_);
    stiffness_matrix.setZero();
    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        stiffness_matrix.reserve(
            Eigen::VectorXi::Constant(2 * active_dof_, entries_number));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        stiffness_matrix.reserve(
            Eigen::VectorXi::Constant(2 * active_dof_, entries_number));
        break;
      }
    }

    const unsigned nnodes_per_cell = global_node_indices_.at(0).size();

    const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
    const unsigned pore_liquid = mpm::ParticlePhase::Liquid;

    const auto& cells = mesh_->cells();
    const auto& nodes = mesh_->active_nodes();

    if (implicit_drag_force) {
      Eigen::MatrixXd k_inter_element;
      unsigned cell_index = 0;
      for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
        if ((*cell_itr)->status()) {
          k_inter_element = (*cell_itr)->K_inter_element();
          for (unsigned i = 0; i < nnodes_per_cell; ++i) {
            for (unsigned j = 0; j < nnodes_per_cell; ++j) {

              auto row = global_node_indices_.at(cell_index)(i);
              auto col = global_node_indices_.at(cell_index)(j);
              stiffness_matrix.coeffRef(row + active_dof_, col) +=
                  -k_inter_element(i, j) * dt;
              stiffness_matrix.coeffRef(row + active_dof_, col + active_dof_) +=
                  k_inter_element(i, j) * dt;
            }
          }
          cell_index++;
        }
      }
    }
    // mass is not a lumped
    for (auto node_itr = nodes.cbegin(); node_itr != nodes.cend(); ++node_itr) {
      auto active_id = (*node_itr)->active_id();
      stiffness_matrix.coeffRef(active_id, active_id) +=
          (*node_itr)->mass(soil_skeleton);

      stiffness_matrix.coeffRef(active_id + active_dof_,
                                active_id + active_dof_) +=
          (*node_itr)->mass(pore_liquid);

      stiffness_matrix.coeffRef(active_id, active_id + active_dof_) +=
          (*node_itr)->mass(pore_liquid);
    }

    // Add stiffness matrix to map
    stiffness_matrix_.insert(
        std::make_pair<unsigned, Eigen::SparseMatrix<double>>(
            static_cast<unsigned>(dir),
            static_cast<Eigen::SparseMatrix<double>>(stiffness_matrix)));

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble force vector (semi-implicit)
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_force_vector(
    double dt) {
  bool status = true;
  try {

    // Resize force vector
    force_inter_.resize(active_dof_ * 2, Tdim);
    force_inter_.setZero();

    // Resize intermediate velocity vector
    acceleration_inter_.resize(active_dof_ * 2, Tdim);
    acceleration_inter_.setZero();

    Eigen::MatrixXd relative_velocity;
    Eigen::MatrixXd drag_force;

    const auto& cells = mesh_->cells();
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        auto k_inter_element = (*cell_itr)->K_inter_element();
        auto nodes = (*cell_itr)->nodes();
        relative_velocity.resize(nodes.size(), Tdim);
        drag_force.resize(nodes.size(), Tdim);
        relative_velocity.setZero();
        drag_force.setZero();

        for (unsigned i = 0; i < nodes.size(); ++i) {
          relative_velocity.row(i) =
              (nodes.at(i)->velocity(1) - nodes.at(i)->velocity(0)).transpose();
        }

        drag_force = k_inter_element * relative_velocity;

        for (unsigned i = 0; i < nodes.size(); ++i) {
          nodes.at(i)->update_drag_force(drag_force.row(i).transpose());
        }
      }
    }

    // Pointer to active nodes
    const auto& active_nodes = mesh_->active_nodes();

    // Active node index
    unsigned node_index = 0;
    for (auto node_itr = active_nodes.cbegin(); node_itr != active_nodes.cend();
         ++node_itr) {
      // Compute nodal intermediate force
      (*node_itr)->compute_intermediate_force(dt);
      // Assemble intermediate force vector
      force_inter_.row(node_index) =
          (*node_itr)->force_total_inter().transpose();
      force_inter_.row(node_index + active_dof_) =
          (*node_itr)->force_fluid_inter().transpose();
      node_index++;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble Laplacian matrix
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_laplacian_matrix(
    double dt, int entries_number) {
  bool status = true;
  try {
    // Initialise Laplacian matrix
    laplacian_matrix_.resize(active_dof_, active_dof_);
    laplacian_matrix_.setZero();
    // TODO: Make the storage being able to adaptively resize
    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        laplacian_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_, entries_number));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        laplacian_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_, entries_number));
        break;
      }
    }
    // Cell pointer
    const auto& cells = mesh_->cells();
    // Iterate over cells
    // active cell id
    mpm::Index cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        // Node ids in each cell
        const auto nids = global_node_indices_.at(cid);
        // Laplacian element of cell
        const auto L_element = (*cell_itr)->L_element();
        const auto P_element = (*cell_itr)->P_element();
        const auto S_element = (*cell_itr)->S_element();        
        // Compute Laplacian elements
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
            laplacian_matrix_.coeffRef(global_node_indices_.at(cid)(i),
                                       global_node_indices_.at(cid)(j)) +=
                L_element(i, j) * dt + P_element(i, j) / dt + S_element(i, j);                
          }
        }
        ++cid;
      }
    }
    laplacian_matrix_.makeCompressed(); 
    // laplacian_matrix_ *= dt;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble Laplacian matrix
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_stab_matrix(
    double dt) {
  bool status = true;
  try {
    // Initialise Laplacian matrix
    stab_matrix_.resize(active_dof_, active_dof_);
    stab_matrix_.setZero();
    // TODO: Make the storage being able to adaptively resize
    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 1200 entries /column
      case (2): {
        stab_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_, 100));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        stab_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_, 30));
        break;
      }
    }
    // Cell pointer
    const auto& cells = mesh_->cells();
    // Iterate over cells
    // active cell id
    mpm::Index cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        // Node ids in each cell
        const auto nids = global_node_indices_.at(cid);
        // Laplacian element of cell
        const auto S_element = (*cell_itr)->S_element();
        // Compute Laplacian elements
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
            stab_matrix_.coeffRef(global_node_indices_.at(cid)(i),
                                       global_node_indices_.at(cid)(j)) +=
                S_element(i, j);               
          }
        }
        ++cid;
      }
    }
    // stab_matrix_.makeCompressed();    
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_poisson_right(
    std::shared_ptr<mpm::Mesh<Tdim>>& mesh_, double dt, int entries_number) {
  bool status = true;
  try {
    // Initialise Fs & Fm matrix
    Eigen::SparseMatrix<double> F_s_matrix, F_m_matrix;
    // Resize Fs matrix
    F_s_matrix.resize(active_dof_, active_dof_ * Tdim);
    F_s_matrix.setZero();
    // Resize Fm matrix
    F_m_matrix.resize(active_dof_, active_dof_ * Tdim);
    F_m_matrix.setZero();
    // TODO: Make the storage being able to adaptively resize
    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        F_s_matrix.reserve(Eigen::VectorXi::Constant(active_dof_ * Tdim, entries_number));
        F_m_matrix.reserve(Eigen::VectorXi::Constant(active_dof_ * Tdim, entries_number));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        F_s_matrix.reserve(Eigen::VectorXi::Constant(active_dof_ * Tdim, entries_number));
        F_m_matrix.reserve(Eigen::VectorXi::Constant(active_dof_ * Tdim, entries_number));
        break;
      }
    }
    // Cell pointer
    const auto& cells = mesh_->cells();
    // Iterate over cells
    // active cell id
    mpm::Index cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      // Cell id
      if ((*cell_itr)->status()) {
        // Node ids in each cell
        const auto nids = global_node_indices_.at(cid);
        // Fs element of cell
        auto F_s_element = (*cell_itr)->F_s_element();
        // Fm element of cell
        auto F_m_element = (*cell_itr)->F_m_element();
        // Compute Fs & Fm
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
            for (unsigned k = 0; k < Tdim; ++k) {
              F_s_matrix.coeffRef(
                  nids(i), nids(j) + k * active_dof_) +=
                  F_s_element(i, j + k * nids.size());
              F_m_matrix.coeffRef(
                  nids(i), nids(j) + k * active_dof_) +=
                  F_m_element(i, j + k * nids.size());

            }
          }
        }
        cid++;
      }          
    }

    // Resize poisson right matrix
    force_laplacian_matrix_.resize(active_dof_);
    force_laplacian_matrix_.setZero();
    // Compute velocity
    Eigen::MatrixXd solid_velocity, relative_velocity;
    solid_velocity.resize(active_dof_, Tdim);
    solid_velocity.setZero();

    relative_velocity.resize(active_dof_, Tdim);
    relative_velocity.setZero();

    // Active nodes
    const auto& active_nodes = mesh_->active_nodes();
    const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
    const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
    unsigned node_index = 0;

    for (auto node_itr = active_nodes.cbegin(); node_itr != active_nodes.cend();
         ++node_itr) {
      // Compute nodal intermediate force
      solid_velocity.row(node_index) =
          (*node_itr)->intermediate_velocity(soil_skeleton).transpose();      
      relative_velocity.row(node_index) =
          ((*node_itr)->intermediate_velocity(pore_liquid) -
           (*node_itr)->intermediate_velocity(soil_skeleton))
              .transpose();
      node_index++;   
    }

    solid_velocity.resize(active_dof_ * Tdim, 1);
    relative_velocity.resize(active_dof_ * Tdim, 1);

    force_laplacian_matrix_ =
        -F_s_matrix * solid_velocity + F_m_matrix * relative_velocity;
                   
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_poisson_right_thermal(
    std::shared_ptr<mpm::Mesh<Tdim>>& mesh_, double dt, int entries_number) {
  bool status = true;
  try {
    // Initialise Ts & Tw matrix
    Eigen::SparseMatrix<double> T_matrix;
    // Resize Ts matrix
    T_matrix.resize(active_dof_, active_dof_);
    T_matrix.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
    // For 2d: 10 entries /column
      case (2): {
        T_matrix.reserve(Eigen::VectorXi::Constant(active_dof_, entries_number));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        T_matrix.reserve(Eigen::VectorXi::Constant(active_dof_, entries_number));
        break;
      }
    }
    // Cell pointer
    const auto& cells = mesh_->cells();
    // Iterate over cells
    // active cell id
    mpm::Index cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      // Cell id
      if ((*cell_itr)->status()) {
        // Node ids in each cell
        const auto nids = global_node_indices_.at(cid);
        // Ts element of cell
        auto T_element = (*cell_itr)->T_element();
        // Compute Ts & Tw
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
              T_matrix.coeffRef(global_node_indices_.at(cid)(i), 
                                  global_node_indices_.at(cid)(j)) += T_element(i, j);
          }
        }
        cid++;
      }
    }

    // Compute temperature acceleration
    Eigen::VectorXd temperature_change_rate;
    temperature_change_rate.resize(active_dof_);
    temperature_change_rate.setZero();

    // Active nodes
    const auto& active_nodes = mesh_->active_nodes();
    const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
    unsigned node_index = 0;

    for (auto node_itr = active_nodes.cbegin(); node_itr != active_nodes.cend();
         ++node_itr) {
      // Compute nodal temperature acceleration
     temperature_change_rate(node_index) =
          (*node_itr)->temperature_acceleration(soil_skeleton);
      node_index++;
    }

    temperature_change_rate.resize(active_dof_, 1);
    force_laplacian_matrix_ += T_matrix * temperature_change_rate;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_poisson_right_pressure(
    std::shared_ptr<mpm::Mesh<Tdim>>& mesh_, double dt, int entries_number) {
  bool status = true;
  try {
    // Initialise Ts & Tw matrix
    Eigen::SparseMatrix<double> P_w_matrix;
    // Resize Ts matrix
    P_w_matrix.resize(active_dof_, active_dof_);
    P_w_matrix.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
    // For 2d: 10 entries /column
      case (2): {
        P_w_matrix.reserve(Eigen::VectorXi::Constant(active_dof_, entries_number));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        P_w_matrix.reserve(Eigen::VectorXi::Constant(active_dof_, entries_number));
        break;
      }
    }
    // Cell pointer
    const auto& cells = mesh_->cells();
    // Iterate over cells
    // active cell id
    mpm::Index cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      // Cell id
      if ((*cell_itr)->status()) {
        // Node ids in each cell
        const auto nids = global_node_indices_.at(cid);
        // Ts element of cell
        auto P_w_element = (*cell_itr)->P_w_element();
        // Compute Ts & Tw
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
              P_w_matrix.coeffRef(global_node_indices_.at(cid)(i), 
                  global_node_indices_.at(cid)(j)) = P_w_element(i, j) / dt;
          }
        }
        cid++;
      }
    }

    force_laplacian_matrix_ += P_w_matrix;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble K_cor_matrix
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_K_cor_matrix(
    std::shared_ptr<mpm::Mesh<Tdim>>& mesh_, double dt, int entries_number) {
  bool status = true;
  try {
    K_cor_matrix_.resize(2 * active_dof_, active_dof_ * Tdim);
    K_cor_matrix_.setZero();
    // Reserve storage for sparse matrix
    switch (Tdim) {
    // For 2d: 10 entries /column
      case (2): {
        K_cor_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_ * Tdim, entries_number));     
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        K_cor_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_ * Tdim, entries_number));
        break;
      }
    }

    unsigned nnodes_per_cell = global_node_indices_.at(0).size();

    const auto& cell = mesh_->cells();

    unsigned cid = 0;
    for (auto cell_itr = cell.cbegin(); cell_itr != cell.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        auto k_cor_element_solid = (*cell_itr)->K_cor_s_element();
        auto k_cor_element_water = (*cell_itr)->K_cor_w_element();

        for (unsigned k = 0; k < Tdim; k++) {
          for (unsigned i = 0; i < nnodes_per_cell; i++) {
            for (unsigned j = 0; j < nnodes_per_cell; j++) {
              // Solid phase
              K_cor_matrix_.coeffRef(
                  global_node_indices_.at(cid)(i),
                  k * active_dof_ + global_node_indices_.at(cid)(j)) +=
                  k_cor_element_solid(i, j + k * nnodes_per_cell);
              // Water phase
              K_cor_matrix_.coeffRef(
                  global_node_indices_.at(cid)(i) + active_dof_,
                  k * active_dof_ + global_node_indices_.at(cid)(j)) +=
                  k_cor_element_water(i, j + k * nnodes_per_cell);                  
            }
          }
        }
        cid++;
      }
    }
    K_cor_matrix_.makeCompressed();
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply pressure constraints vector
template <unsigned Tdim>
void mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::apply_pressure_constraints() {
  try {
    // Modify force_laplacian_matrix_
    force_laplacian_matrix_ -= laplacian_matrix_ * pressure_constraints_;
    // Apply free surface
    if (!free_surface_.empty()) {
      const auto& nodes = mesh_->nodes();
      for (const auto& free_node : free_surface_) {
        const auto column_index = nodes[free_node]->active_id();
        // Modify force_laplacian_matrix        
        force_laplacian_matrix_(column_index) = 0;        
        // Modify laplacian_matrix
        laplacian_matrix_.row(column_index) *= 0;
        laplacian_matrix_.col(column_index) *= 0;
        laplacian_matrix_.coeffRef(column_index, column_index) = 1;
      }
      // Clear the vector
      free_surface_.clear();
    }

    // Iterate over pressure constraints
    for (Eigen::SparseVector<double>::InnerIterator it(pressure_constraints_);
         it; ++it) {
      // Modify force_laplacian_matrix
      force_laplacian_matrix_(it.index()) = it.value();
      // Modify laplacian_matrix
      laplacian_matrix_.row(it.index()) *= 0;
      laplacian_matrix_.col(it.index()) *= 0;
      laplacian_matrix_.coeffRef(it.index(), it.index()) = 1;
      std::cout << it.index() << "\n";
    }  
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
}

//! Apply velocity constraints to force vector
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<
    Tdim>::apply_velocity_constraints() {
  bool status = false;
  try {
    // Modify the force vector(b = b - A * bc)
    for (unsigned i = 0; i < Tdim; i++) {
      force_inter_.col(i) -=
          stiffness_matrix_.at(i) * velocity_constraints_.col(i);

      // Iterate over velocity constraints (non-zero elements)
      for (unsigned j = 0; j < velocity_constraints_.outerSize(); ++j) {
        for (Eigen::SparseMatrix<double>::InnerIterator itr(
                 velocity_constraints_, j);
             itr; ++itr) {
          // Check direction
          if (itr.col() == i) {
            // Assign 0 to specified column
            stiffness_matrix_.at(i).col(itr.row()) *= 0;
            // Assign 0 to specified row
            stiffness_matrix_.at(i).row(itr.row()) *= 0;
            // Assign 1  to diagnal element
            stiffness_matrix_.at(i).coeffRef(itr.row(), itr.row()) = 1;

            force_inter_(itr.row(), itr.col()) = 0;
          }
        }
      }
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

//! Assign pressure constraints
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assign_pressure_constraints(
    const double beta, const double current_time) {
  bool status = false;
  try {
    // Phase index
    const unsigned pore_liquid = mpm::ParticlePhase::Liquid;
    // Total number of phases
    const unsigned Tnphases = 2;
    // Resize pressure constraints vector
    pressure_constraints_.resize(active_dof_);
    pressure_constraints_.reserve(int(0.5 * active_dof_));

    // Nodes container
    const auto& nodes = mesh_->active_nodes();
    // Iterate over nodes to get pressure constraints
    for (auto node = nodes.cbegin(); node != nodes.cend(); ++node) {
      // Assign total pressure constraint
      const double pressure_constraint =
          (*node)->pressure_constraint(pore_liquid, current_time);
      // Assign pressure increment constraint
      // const double pressure_increment_constraint =
      //    (*node)->pressure_constraint(pore_liquid + Tnphases, current_time);

      // Check if there is a pressure constraint
      if (pressure_constraint != std::numeric_limits<double>::max()) {
        // Insert the pressure constraints
        pressure_constraints_.insert((*node)->active_id()) =
            (1 - beta) * pressure_constraint;
      }
      // // Check if there is a pressure constraint
      // else if (pressure_increment_constraint !=
      //          std::numeric_limits<double>::max()) {
      //   // Insert the pressure constraints
      //   // TODO: If beta != 1
      //   pressure_constraints_.insert((*node)->active_id()) =
      //       pressure_increment_constraint;
      // }
    }
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

//! Assign velocity constraints
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<
    Tdim>::assign_velocity_constraints() {
  bool status = false;
  try {
    // Initialise constraints matrix from triplet
    std::vector<Eigen::Triplet<double>> triplet_list;
    // Nodes container
    const auto& nodes = mesh_->active_nodes();
    // Iterate over nodes
    for (auto node = nodes.cbegin(); node != nodes.cend(); ++node) {
      // Get velocity constraints
      const auto& velocity_constraints = (*node)->velocity_constraints();
      // Assign constraints matrix
      for (const auto constraint : velocity_constraints) {
        // Insert constraint to the matrix
        triplet_list.push_back(Eigen::Triplet<double>(
            (constraint).first / Tdim * active_dof_ + (*node)->active_id(),
            (constraint).first % Tdim, (constraint).second));
      }
    }
    // Reserve the storage for the velocity constraints matrix
    velocity_constraints_.setZero();
    velocity_constraints_.data().squeeze();
    velocity_constraints_.resize(active_dof_ * 2, Tdim);
    velocity_constraints_.reserve(
        Eigen::VectorXi::Constant(Tdim, triplet_list.size() + 1200));
    // Assemble the velocity constraints matrix
    velocity_constraints_.setFromTriplets(triplet_list.begin(),
                                          triplet_list.end());

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}


//------------------------------------------------------------
// Implict mpm

//! Assemble KTT matrix
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_KTT_matrix(
    double dt) {
  bool status = true;
  try {
    // Initialise KTT matrix
    KTT_matrix_.resize(active_dof_, active_dof_);
    KTT_matrix_.setZero();
    KTT_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_, 100));

    // Cell pointer
    const auto& cells = mesh_->cells();
    // Iterate over cells
    // active cell id
    mpm::Index cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        // Node ids in each cell
        const auto nids = global_node_indices_.at(cid);
        // KTT element of cell
        const auto KTT_element = (*cell_itr)->KTT_element();
        // Compute KTT elements
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
            KTT_matrix_.coeffRef(global_node_indices_.at(cid)(i),
                                       global_node_indices_.at(cid)(j)) +=
                KTT_element(i, j);
          }
        }
        ++cid;
      }
    }
    KTT_matrix_.makeCompressed(); 
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble MTT matrix
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_MTT_matrix(double dt) {
  bool status = true;
  try {

    // Resize stiffness matrix
    MTT_matrix_.resize(active_dof_, active_dof_);
    MTT_matrix_.resize(active_dof_, active_dof_);
    MTT_matrix_.setZero();
    // Reserve storage for sparse matrix
    MTT_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_, 100));

    const unsigned nnodes_per_cell = global_node_indices_.at(0).size();

    const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
    const unsigned pore_liquid = mpm::ParticlePhase::Liquid;

    const auto& cells = mesh_->cells();
    const auto& nodes = mesh_->active_nodes();

    // mass is not a lumped
    for (auto node_itr = nodes.cbegin(); node_itr != nodes.cend(); ++node_itr) {
      auto active_id = (*node_itr)->active_id();
      MTT_matrix_.coeffRef(active_id, active_id) +=
          (*node_itr)->heat_capacity(soil_skeleton);
    }

    MTT_matrix_ += dt * KTT_matrix_;
    MTT_matrix_.makeCompressed(); 

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_FT_vector(
    std::shared_ptr<mpm::Mesh<Tdim>>& mesh_, double dt) {
  bool status = true;
  try {
   
    // Assemble temperature
    FT_vector_.resize(active_dof_);
    FT_vector_.setZero();
   
    // Assemble temperature
    temperature_.resize(active_dof_);
    temperature_.setZero();

    // Active nodes
    const auto& active_nodes = mesh_->active_nodes();
    const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
    unsigned node_index = 0;

    for (auto node_itr = active_nodes.cbegin(); node_itr != active_nodes.cend();
         ++node_itr) {
      // Compute nodal temperature acceleration
     temperature_(node_index) =
          (*node_itr)->temperature(soil_skeleton);
      node_index++;
    }

    temperature_.resize(active_dof_, 1);
    FT_vector_ += -KTT_matrix_ * temperature_;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign temperature constraints
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assign_temperature_constraints(
  const double current_time) {
  bool status = false;
  try {
    // Phase index
    const unsigned soil_skeleton = mpm::ParticlePhase::Solid;
    // Resize temperature constraints vector
    temperature_constraints_.resize(active_dof_);
    temperature_constraints_.reserve(int(0.5 * active_dof_));

    // Nodes container
    const auto& nodes = mesh_->active_nodes();
    // Iterate over nodes to get temperature constraints
    for (auto node = nodes.cbegin(); node != nodes.cend(); ++node) {
      // Assign total temperature constraint
      const double temperature_constraint =
          (*node)->temperature_constraint(soil_skeleton, current_time);

      // Check if there is a temperature constraint
      if (temperature_constraint != std::numeric_limits<double>::max()) {
        // Insert the temperature constraints
        temperature_constraints_.insert((*node)->active_id()) =
            temperature_constraint;
      }
    }
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

//! Apply temperature constraints vector
template <unsigned Tdim>
void mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::apply_temperature_constraints() {
  try {
    // // Modify FT_vector_
    // FT_vector_ -= MTT_matrix_ * temperature_constraints_;

    // Iterate over temperature constraints
    for (Eigen::SparseVector<double>::InnerIterator it(temperature_constraints_);
         it; ++it) {
      // Modify FT_vector
      FT_vector_(it.index()) = it.value();
      // Modify laplacian_matrix
      MTT_matrix_.row(it.index()) *= 0;
      MTT_matrix_.col(it.index()) *= 0;
      MTT_matrix_.coeffRef(it.index(), it.index()) = 1;
      std::cout << it.index() << "\n";
    }  
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
}