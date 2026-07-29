
// Compute acceleration and velocity for two phase
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_acc_vel_twophase_explicit(
    unsigned soil_skeleton, unsigned pore_fluid, unsigned mixture, double dt) {
  bool status = true;
  const double tolerance = 1.0E-15;
  try {
    // Compute drag force
    VectorDim drag_force = drag_force_coefficient_(pore_fluid) * (
        velocity_.col(pore_fluid) - velocity_.col(soil_skeleton));

    // if (mass_(soil_skeleton) > tolerance && mass_(pore_fluid) > tolerance) {
    //   // Acceleration of pore fluid (momentume balance of fluid phase)
    //   this->acceleration_.col(pore_fluid) =
    //       (this->external_force_.col(pore_fluid) +
    //        this->internal_force_.col(pore_fluid) - drag_force) /
    //       this->mass_(pore_fluid);

    //   // Acceleration of solid skeleton (momentume balance of mixture)
    //   this->acceleration_.col(soil_skeleton) =
    //       (this->external_force_.col(mixture) +
    //        this->internal_force_.col(mixture) -
    //        this->mass_(pore_fluid) * this->acceleration_.col(pore_fluid)) /
    //       this->mass_(soil_skeleton);


    VectorDim mixture_force = this->external_force_.col(mixture) +
            this->internal_force_.col(mixture);
    VectorDim liquid_force = this->external_force_.col(pore_fluid) +
            this->internal_force_.col(pore_fluid);
    double mass_ratio = 1 / this->mass_(soil_skeleton) +
                        1 / this->mass_(pore_fluid);
    double drag_force_coefficient = drag_force_coefficient_(pore_fluid);

    //! get reaction force
    reaction_force_.setZero();

    if (this->contact_) {
        reaction_force_ = -mixture_force;
    }

    if (mass_(soil_skeleton) > tolerance && mass_(pore_fluid) > tolerance) {
      //! Acceleration of pore fluid (momentume balance of fluid phase)
      //! (1) Explicit - (vl_star - vs_star)
      // this->acceleration_inter_.col(pore_fluid) = 
      //     (liquid_force - drag_force) / this->mass_(pore_fluid);      

      //! Yuan et al. (2023 CG) - (vl_star - vs_k)
      // this->acceleration_inter_.col(pore_fluid) = 
      //     (liquid_force - drag_force) / 
      //     (this->mass_(pore_fluid) + dt * drag_force_coefficient);

      //! Impliict - (vl_star - vs_star)  
      this->acceleration_.col(pore_fluid) = 
          (dt * drag_force_coefficient / this->mass_(soil_skeleton) * 
          mixture_force + liquid_force - drag_force) / this->mass_(pore_fluid) /
          (1 + dt * drag_force_coefficient * mass_ratio);  

      // Acceleration of solid skeleton (momentume balance of mixture)
      this->acceleration_.col(soil_skeleton) =
          (mixture_force -
          this->mass_(pore_fluid) * this->acceleration_.col(pore_fluid)) /
          this->mass_(soil_skeleton);

      // Apply friction constraints
      this->apply_friction_constraints(dt);

      // Velocity += acceleration * dt
      this->velocity_ += this->acceleration_ * dt;

      // Apply velocity constraints, which also sets acceleration to 0,
      // when velocity is set.
      this->apply_velocity_constraints();

      if (this->contact_) {
        this->velocity_ = rigid_velocity_;
        this->acceleration_ = rigid_acceleration_;
      }

      if (this->contact_) {
        // set zero total force for rigid particle influence node
        // get reaction force
        reaction_force_.setZero();
        reaction_force_ = -(this->external_force_.col(mixture) +
                            this->internal_force_.col(mixture));
      }
      
      // Set a threshold
      for (unsigned i = 0; i < Tdim; ++i) {
        if (!(std::abs(velocity_.col(soil_skeleton)(i))) > tolerance)
          velocity_.col(soil_skeleton)(i) = 0.;
        if (!(std::abs(acceleration_.col(soil_skeleton)(i))) > tolerance)
          acceleration_.col(soil_skeleton)(i) = 0.;
        if (!(std::abs(velocity_.col(pore_fluid)(i))) > tolerance)
          velocity_.col(pore_fluid)(i) = 0.;
        if (!(std::abs(acceleration_.col(pore_fluid)(i))) > tolerance)
          acceleration_.col(pore_fluid)(i) = 0.;
      }
    }

  // // Debug
  // if (id_ <= 100 && id_ >= 90) {
  //   std::cout << "nodes_.id_: " << id_ << std::endl;
  //   std::cout << "    drag_force_coefficient_(pore_fluid): " << drag_force_coefficient_(pore_fluid) << std::endl;
  //   std::cout << "    external_force_.col(mixture): " << this->external_force_.col(mixture) << std::endl;
  //   std::cout << "    external_force_.col(pore_fluid): " << this->external_force_.col(pore_fluid) << std::endl;
  //   std::cout << "    internal_force_.col(mixture): " << this->internal_force_.col(mixture) << std::endl;
  //   std::cout << "    internal_force_.col(pore_fluid): " << this->internal_force_.col(pore_fluid) << std::endl;
  // }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute acceleration and velocity for two Phase
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_inter_acc_vel_twophase_semi(
                                                    unsigned soil_skeleton,
                                                    unsigned pore_fluid,
                                                    unsigned mixture,
                                                    double dt) {
  bool status = true;
  const double tolerance = 1.0E-15;
  try {
    // Compute drag force
    VectorDim drag_force = drag_force_coefficient_(pore_fluid) * (
        velocity_.col(pore_fluid) - velocity_.col(soil_skeleton));
    VectorDim mixture_force = this->external_force_.col(mixture) +
           this->internal_force_.col(mixture);
    VectorDim liquid_force = this->external_force_.col(pore_fluid) +
           this->internal_force_.col(pore_fluid);
    double mass_ratio = 1 / this->mass_(soil_skeleton) +
                        1 / this->mass_(pore_fluid);
    double drag_force_coefficient = drag_force_coefficient_(pore_fluid);

    //! set zero total force for rigid particle influence node
    //! get reaction force
    reaction_force_.setZero();

    if (this->contact_) {
        reaction_force_ = -mixture_force;
    }

    if (mass_(soil_skeleton) > tolerance && mass_(pore_fluid) > tolerance) {
      //! Acceleration of pore fluid (momentume balance of fluid phase)
      //! (1) Explicit - (vl_star - vs_star)
      // this->acceleration_inter_.col(pore_fluid) = 
      //     (liquid_force - drag_force) / this->mass_(pore_fluid);      

      //! Yuan et al. (2023 CG) - (vl_star - vs_k)
      // this->acceleration_inter_.col(pore_fluid) = 
      //     (liquid_force - drag_force) / 
      //     (this->mass_(pore_fluid) + dt * drag_force_coefficient);

      //! Impliict - (vl_star - vs_star)  
      this->acceleration_inter_.col(pore_fluid) = 
          (dt * drag_force_coefficient / this->mass_(soil_skeleton) * 
          mixture_force + liquid_force - drag_force) / this->mass_(pore_fluid) /
          (1 + dt * drag_force_coefficient * mass_ratio);  

      // Acceleration of solid skeleton (momentume balance of mixture)
      this->acceleration_inter_.col(soil_skeleton) =
          (mixture_force -
          this->mass_(pore_fluid) * this->acceleration_inter_.col(pore_fluid)) /
          this->mass_(soil_skeleton);

      // Apply friction constraints
      this->apply_friction_constraints(dt);

      // Velocity += acceleration * dt
      this->velocity_inter_ = this->velocity_ + this->acceleration_inter_ * dt;
      }
      // Apply velocity constraints, which also sets acceleration to 0,
      // when velocity is set.
      this->apply_velocity_constraints();

      // Set a threshold
      for (unsigned i = 0; i < Tdim; ++i) {
        if (std::abs(velocity_inter_.col(soil_skeleton)(i)) < tolerance)
          velocity_inter_.col(soil_skeleton)(i) = 0.;
        if ((std::abs(acceleration_inter_.col(soil_skeleton)(i))) < tolerance)
          acceleration_inter_.col(soil_skeleton)(i) = 0.;
        if ((std::abs(velocity_inter_.col(pore_fluid)(i))) < tolerance)
          velocity_inter_.col(pore_fluid)(i) = 0.;
        if ((std::abs(acceleration_inter_.col(pore_fluid)(i))) < tolerance)
          acceleration_inter_.col(pore_fluid)(i) = 0.;

    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute semi-implicit acceleration and velocity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_acc_vel_twophase_semi(
                                                unsigned phase, double dt) {
  bool status = true;
  const double tolerance = std::numeric_limits<double>::min();
  try {

    // Semi-implicit solver
    Eigen::Matrix<double, Tdim, 1> acceleration_corrected = 
        force_cor_.col(phase) / mass_(phase);
    // Acceleration
    this->acceleration_.col(phase) = acceleration_corrected;

    // Update velocity
    velocity_.col(phase) += acceleration_corrected * dt;

    // Apply friction constraints
    this->apply_friction_constraints(dt);

    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints();

    // Set a threshold
    for (unsigned i = 0; i < Tdim; ++i) {
      if (!(std::abs(velocity_.col(phase)(i)) > tolerance))
        velocity_.col(phase)(i) = 0.;
      if (!(std::abs(acceleration_.col(phase)(i)) > tolerance))
        acceleration_.col(phase)(i) = 0.;
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute intermediate force
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_intermediate_force(
    const double dt) {
  bool status = true;
  const unsigned mixture = 0;
  const unsigned pore_liquid = 1;
  try {
    std::lock_guard<std::mutex> guard(node_mutex_);

    auto force_total = internal_force_ + external_force_;
    force_total_inter_ = force_total.col(mixture);
    force_fluid_inter_ = force_total.col(pore_liquid) - drag_force_;

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute nodal corrected force
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_nodal_corrected_force(
    VectorDim& force_cor_part_solid, VectorDim& force_cor_part_water) {
  bool status = true;

  try {
    // Compute corrected force for solid phase
    force_cor_.col(0) =
        mass_(0) * acceleration_inter_.col(0) + force_cor_part_solid;
    // Compute corrected force for water phase
    force_cor_.col(1) =
        mass_(1) * acceleration_inter_.col(1) + force_cor_part_water;
      
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Update pore pressure increment at the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_pore_pressure_increment(
    const Eigen::VectorXd& pore_pressure_increment, double current_time) {
  const unsigned pore_liquid = 1;
  this->pore_pressure_increment_ = pore_pressure_increment(active_id_);

  // If pressure boundary, increment is zero
  if (pressure_constraints_.find(pore_liquid) != pressure_constraints_.end() ||
      this->free_surface()) {
    this->pore_pressure_increment_ = 0; 
  } 
}

// Update temperature rate at the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_temperature_rate(
    const Eigen::VectorXd& temperature_rate, double dt, double current_time) {
  this->temperature_acceleration_(0) = temperature_rate(active_id_);

  if ((temperature_constraints_.find(0) != temperature_constraints_.end())) {
    this->temperature_acceleration_(0) = 0; 
  }

  // If temperature boundary, increment is zero
  if (!(temperature_constraints_.find(0) != temperature_constraints_.end())) {
    this->temperature_ += this->temperature_acceleration_ * dt; 
  }

}
