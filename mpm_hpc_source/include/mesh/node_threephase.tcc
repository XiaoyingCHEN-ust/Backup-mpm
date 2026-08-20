
// // Compute acceleration and velocity for two phase
// template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
// bool mpm::Node<Tdim, Tdof, Tnphases>::compute_acc_vel_threephase_explicit(
//     unsigned soil_skeleton, unsigned pore_fluid, unsigned pore_gas, 
//     unsigned mixture,  double dt) {
//   bool status = true;
//   const double tolerance = 1.0E-15;
//   try {
//     // Compute liquid drag force
//     this->drag_force_liquid_ = drag_force_coefficient_(pore_fluid) * (
//         velocity_.col(pore_fluid) - velocity_.col(soil_skeleton));

//     // Compute gas drag force
//     this->drag_force_gas_ = drag_force_coefficient_(pore_gas) * (
//         velocity_.col(pore_gas) - velocity_.col(soil_skeleton));

//     if (mass_(soil_skeleton) > tolerance && 
//         mass_(pore_fluid) > tolerance && 
//         mass_(pore_gas) > tolerance) {

//       // Acceleration of pore fluid (momentume balance of fluid phase)
//       this->acceleration_.col(pore_fluid) =
//           (this->external_force_.col(pore_fluid) +
//             this->internal_force_.col(pore_fluid) - drag_force_liquid_) /
//             this->mass_(pore_fluid);

//       // Acceleration of pore gas (momentume balance of fluid phase)
//       this->acceleration_.col(pore_gas) =
//           (this->external_force_.col(pore_gas) +
//           this->internal_force_.col(pore_gas) - drag_force_gas_) /
//           this->mass_(pore_gas);

//       // Acceleration of solid skeleton (momentume balance of mixture)
//       this->acceleration_.col(soil_skeleton) =
//           (this->external_force_.col(mixture) +
//           this->internal_force_.col(mixture) -
//           this->mass_(pore_fluid) * this->acceleration_.col(pore_fluid) -
//           this->mass_(pore_gas) * this->acceleration_.col(pore_gas)) /
//           this->mass_(soil_skeleton);

//       // Apply friction constraints
//       this->apply_friction_constraints(dt);

//       // Velocity += acceleration * dt
//       this->velocity_ += this->acceleration_ * dt;

//       // Apply velocity constraints, which also sets acceleration to 0,
//       // when velocity is set.
//       this->apply_velocity_constraints();

//       // Set a threshold
//       for (unsigned i = 0; i < Tdim; ++i) {
//         if (!(std::abs(velocity_.col(soil_skeleton)(i))) > tolerance)
//           velocity_.col(soil_skeleton)(i) = 0.;
//         if (!(std::abs(acceleration_.col(soil_skeleton)(i))) > tolerance)
//           acceleration_.col(soil_skeleton)(i) = 0.;
//         if (!(std::abs(velocity_.col(pore_fluid)(i))) > tolerance)
//           velocity_.col(pore_fluid)(i) = 0.;
//         if (!(std::abs(acceleration_.col(pore_fluid)(i))) > tolerance)
//           acceleration_.col(pore_fluid)(i) = 0.;
//         if (!(std::abs(velocity_.col(pore_gas)(i))) > tolerance)
//           velocity_.col(pore_gas)(i) = 0.;
//         if (!(std::abs(acceleration_.col(pore_gas)(i))) > tolerance)
//           acceleration_.col(pore_gas)(i) = 0.;
//       }
//     }
//   } catch (std::exception& exception) {
//     console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
//     status = false;
//   }
//   return status;
// }

// Compute acceleration and velocity for three phase
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_acc_vel_threephase_explicit(
    unsigned soil_skeleton, unsigned pore_fluid, unsigned pore_gas, 
    unsigned mixture, double dt) {
  bool status = true;
  double tolerance = 1.0E-15;
  try {

    if (minimum_nodal_density_ > 0. && volume_(soil_skeleton) > tolerance &&
        mass_(soil_skeleton) / volume_(soil_skeleton) <
            minimum_nodal_density_) {
      velocity_.setZero();
      acceleration_.setZero();
      this->apply_velocity_constraints();
      return true;
    }

    // Compute liquid drag force
    this->drag_force_liquid_ = drag_force_coefficient_(pore_fluid) * (
        velocity_.col(pore_fluid) - velocity_.col(soil_skeleton));

    // Compute gas drag force
    this->drag_force_gas_ = drag_force_coefficient_(pore_gas) * (
        velocity_.col(pore_gas) - velocity_.col(soil_skeleton));

    if (mass_(soil_skeleton) > tolerance) {

      // Mixture force vector
      Eigen::Matrix<double, Tdim, Tnphases> F_matrix;
      F_matrix.setZero();
      
      F_matrix.col(mixture) = this->external_force_.col(mixture) +
            this->internal_force_.col(mixture);

      F_matrix.col(pore_fluid) = this->external_force_.col(pore_fluid) +
              this->internal_force_.col(pore_fluid) - drag_force_liquid_;

      F_matrix.col(pore_gas) = this->external_force_.col(pore_gas) +
              this->internal_force_.col(pore_gas) - drag_force_gas_; 

      // if (mass_(pore_fluid) > tolerance) 
      //   // Acceleration of pore fluid (momentume balance of fluid phase)
      //   this->acceleration_.col(pore_fluid) = F_matrix.col(pore_fluid) / 
      //       (this->mass_(pore_fluid) + dt * drag_force_coefficient_(pore_fluid));

      // if (mass_(pore_gas) > tolerance) 
      // // Acceleration of pore gas
      // this->acceleration_.col(pore_gas) = F_matrix.col(pore_gas) / 
      //           (this->mass_(pore_gas) + dt * drag_force_coefficient_(pore_gas));

      // if (mass_(pore_fluid) > tolerance) 
      //   // Acceleration of pore fluid (momentume balance of fluid phase)
      //   this->acceleration_.col(pore_fluid) = F_matrix.col(pore_fluid) / this->mass_(pore_fluid);

      // if (mass_(pore_gas) > tolerance) 
      // // Acceleration of pore gas
      // this->acceleration_.col(pore_gas) = F_matrix.col(pore_gas) / this->mass_(pore_gas);                                            

      // // Acceleration of solid skeleton (momentume balance of mixture)
      // this->acceleration_.col(soil_skeleton) = (F_matrix.col(mixture) -
      //      this->mass_(pore_fluid) * this->acceleration_.col(pore_fluid) -
      //      this->mass_(pore_gas) * this->acceleration_.col(pore_gas)) / this->mass_(soil_skeleton);

      Eigen::Matrix<double, Tnphases, Tnphases> M_matrix;
      M_matrix.setZero();
      if ((mass_(pore_fluid) > tolerance) & (mass_(pore_gas) > tolerance)) {

        M_matrix(0, 0) = this->mass_(soil_skeleton);
        M_matrix(0, 1) = this->mass_(pore_fluid);
        M_matrix(0, 2) = this->mass_(pore_gas);
        M_matrix(1, 0) = -dt * drag_force_coefficient_(pore_fluid);
        M_matrix(1, 1) = this->mass_(pore_fluid) + dt * drag_force_coefficient_(pore_fluid);
        M_matrix(1, 2) = 0;
        M_matrix(2, 0) = -dt * drag_force_coefficient_(pore_gas);
        M_matrix(2, 1) = 0;
        M_matrix(2, 2) = this->mass_(pore_gas) + dt * drag_force_coefficient_(pore_gas);

        const auto mass_solver = M_matrix.fullPivLu();
        const double reciprocal_condition = mass_solver.rcond();
        if (!mass_solver.isInvertible() ||
            !std::isfinite(reciprocal_condition))
          throw std::runtime_error(
              "Three-phase nodal mass matrix is singular or non-finite");
        this->acceleration_ =
            mass_solver.solve(F_matrix.transpose()).transpose();

      }
      else if ((mass_(pore_fluid) > tolerance)) {

        M_matrix(0, 0) = this->mass_(soil_skeleton);
        M_matrix(0, 1) = this->mass_(pore_fluid);
        M_matrix(0, 2) = 0;
        M_matrix(1, 0) = -dt * drag_force_coefficient_(pore_fluid);
        M_matrix(1, 1) = this->mass_(pore_fluid) + dt * drag_force_coefficient_(pore_fluid);
        M_matrix(1, 2) = 0;
        M_matrix(2, 0) = 0;
        M_matrix(2, 1) = 0;
        M_matrix(2, 2) = 1;
        F_matrix.col(pore_gas).setZero();

        const auto mass_solver = M_matrix.fullPivLu();
        const double reciprocal_condition = mass_solver.rcond();
        if (!mass_solver.isInvertible() ||
            !std::isfinite(reciprocal_condition))
          throw std::runtime_error(
              "Liquid-solid nodal mass matrix is singular or non-finite");
        this->acceleration_ =
            mass_solver.solve(F_matrix.transpose()).transpose();
      }
      else {
        M_matrix(0, 0) = this->mass_(soil_skeleton);
        M_matrix(0, 1) = 0;
        M_matrix(0, 2) = 0;
        M_matrix(1, 0) = 0;
        M_matrix(1, 1) = 1;
        M_matrix(1, 2) = 0;
        M_matrix(2, 0) = 0;
        M_matrix(2, 1) = 0;
        M_matrix(2, 2) = 1;
        F_matrix.col(pore_fluid).setZero();
        F_matrix.col(pore_gas).setZero();

        const auto mass_solver = M_matrix.fullPivLu();
        const double reciprocal_condition = mass_solver.rcond();
        if (!mass_solver.isInvertible() ||
            !std::isfinite(reciprocal_condition))
          throw std::runtime_error(
              "Solid nodal mass matrix is singular or non-finite");
        this->acceleration_ =
            mass_solver.solve(F_matrix.transpose()).transpose();
      }

      // Preserve the free-surface kinematic condition through the force
      // update. Particle PIC updates interpolate velocity, while FLIP updates
      // interpolate acceleration, so both quantities must be phase consistent
      // before boundary-condition corrections are applied.
      if (this->phase_kinematic_boundary_ && Tnphases == 3) {
        this->velocity_.col(pore_fluid) =
            this->velocity_.col(soil_skeleton);
        this->velocity_.col(pore_gas) =
            this->velocity_.col(soil_skeleton);
        this->acceleration_.col(pore_fluid) =
            this->acceleration_.col(soil_skeleton);
        this->acceleration_.col(pore_gas) =
            this->acceleration_.col(soil_skeleton);
      }

      // Friction modifies acceleration, so apply it before integrating the
      // corrected nodal velocity. Because all free-surface phases now have
      // identical inputs, the shared friction law preserves their equality.
      this->apply_friction_constraints(dt);

      // Velocity += acceleration * dt
      this->velocity_ += this->acceleration_ * dt;

      if (!this->acceleration_.allFinite()) {
        std::ostringstream message;
        message << "Non-finite three-phase nodal acceleration at node "
                << id_ << ", coordinates=(" << coordinates_.transpose()
                << "), mass=" << mass_ << ", drag="
                << drag_force_coefficient_ << ", force=" << F_matrix
                << ", matrix=" << M_matrix << ", acceleration="
                << acceleration_;
        throw std::runtime_error(message.str());
      }

      tolerance = 1E-15;
      // Set a threshold
      for (unsigned i = 0; i < Tdim; ++i) {
        if ((std::abs(velocity_.col(soil_skeleton)(i))) < tolerance)
          velocity_.col(soil_skeleton)(i) = 0.;
        if ((std::abs(acceleration_.col(soil_skeleton)(i))) < tolerance)
          acceleration_.col(soil_skeleton)(i) = 0.;
        if ((std::abs(velocity_.col(pore_fluid)(i))) < tolerance)
          velocity_.col(pore_fluid)(i) = 0.;
        if ((std::abs(acceleration_.col(pore_fluid)(i))) < tolerance)
          acceleration_.col(pore_fluid)(i) = 0.;
        if ((std::abs(velocity_.col(pore_gas)(i))) < tolerance)
          velocity_.col(pore_gas)(i) = 0.;
        if ((std::abs(acceleration_.col(pore_gas)(i))) < tolerance)
          acceleration_.col(pore_gas)(i) = 0.;
      }

      // Rigid contact is authoritative over the unconstrained free-surface
      // update and retains the existing mixture-force reaction definition.
      if (this->contact_) {
        this->velocity_ = rigid_velocity_;
        this->acceleration_ = rigid_acceleration_;
        reaction_force_.setZero();
        reaction_force_ = -(this->external_force_.col(mixture) +
                            this->internal_force_.col(mixture));
      }

      // Prescribed velocities are the final authority when a free-surface
      // node is also constrained or belongs to a rigid contact. Applying the
      // constraint also zeros the corresponding acceleration for FLIP.
      this->apply_velocity_constraints();

      if (!this->acceleration_.allFinite())
        throw std::runtime_error(
            "Non-finite three-phase nodal acceleration after constraints at "
            "node " +
            std::to_string(id_));
    }

  // // Debug
  // if (id_ <= 10 && id_ >= 0) {
  //   std::cout << "nodes_.id_: " << id_ << std::endl;
  //   // std::cout << "    velocity_.col(mixture): " << velocity_.col(mixture) << std::endl;
  //   std::cout << "    velocity_.col(pore_fluid): " << velocity_.col(pore_fluid) << std::endl;
  //   std::cout << "    drag_force_liquid_: " << drag_force_liquid_ << std::endl;
  //   // std::cout << "    drag_force_gas_: " << drag_force_gas_ << std::endl;
  //   // std::cout << "    external_force_.col(mixture): " << this->external_force_.col(mixture) << std::endl;
  //   // std::cout << "    external_force_.col(pore_fluid): " << this->external_force_.col(pore_fluid) << std::endl;
  //   // std::cout << "    external_force_.col(pore_gas): " << this->external_force_.col(pore_gas) << std::endl;
  //   // std::cout << "    internal_force_.col(mixture): " << this->internal_force_.col(mixture) << std::endl;
  //   // std::cout << "    internal_force_.col(pore_fluid): " << this->internal_force_.col(pore_fluid) << std::endl;
  //   // std::cout << "    internal_force_.col(pore_gas): " << this->internal_force_.col(pore_gas) << std::endl;
  // }

  } catch (std::exception& exception) {
    throw std::runtime_error("Three-phase nodal update failed at node " +
                             std::to_string(id_) + ": " + exception.what());
  }
  return status;
}

// Compute nodal pressure acceleration and update nodal pressure
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_acceleration_pressure_threephase(
    unsigned pore_liquid, unsigned pore_gas, double dt) noexcept {
  bool status = false;
  const double tolerance = 1.0E-15;

  if (this->mass_(pore_liquid) > tolerance || this->mass_(pore_gas) > tolerance) {

    Eigen::Matrix<double, 2, 2> K_matrix, K_matrix_inverse;
    K_matrix.setZero();
    K_matrix(0,0) = K_coeff_(pore_liquid) + mass_(pore_liquid);
    K_matrix(0,1) = -K_coeff_(pore_liquid);
    K_matrix(1,0) = -K_coeff_(pore_gas);
    K_matrix(1,1) = K_coeff_(pore_gas) + mass_(pore_gas);

    Eigen::Matrix<double, 1, 2> F_matrix;
    F_matrix.setZero();
    F_matrix[0] = this->hydraulic_conduction_(pore_liquid) + this->mass_source_(pore_liquid);// - this->volumetric_strain_(pore_liquid);// + this->mass_convection_(pore_liquid);
    F_matrix[1] = this->hydraulic_conduction_(pore_gas) + this->mass_source_(pore_gas);// - this->volumetric_strain_(pore_gas);// + this->mass_convection_(pore_gas);

    if (this->mass_(pore_liquid) > tolerance & this->mass_(pore_gas) > tolerance) {
      K_matrix_inverse = K_matrix.inverse();
      this->pressure_acceleration_(pore_liquid) = K_matrix_inverse(0,0) * F_matrix[0] + K_matrix_inverse(0,1) * F_matrix[1];
      this->pressure_acceleration_(pore_gas) = K_matrix_inverse(1,0) * F_matrix[0] + K_matrix_inverse(1,1) * F_matrix[1];
    } 
    else if (this->mass_(pore_liquid) > tolerance) {
      this->pressure_acceleration_(pore_liquid) = F_matrix[0] / K_matrix(0,0);
      this->pressure_acceleration_(pore_gas) = 0;
    }
    else if (this->mass_(pore_gas) > tolerance) {
      this->pressure_acceleration_(pore_gas) = F_matrix[1] / K_matrix(1,1);
      this->pressure_acceleration_(pore_liquid) = 0;
    }

    this->pressure_ += this->pressure_acceleration_ * dt;

    if (this->free_surface_) {
      this->pressure_acceleration_.setZero();
      this->pressure_.setZero();
    }

    status = true;
  }
  return status;
}
