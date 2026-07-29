#ifndef MPM_RIGID_PIPELINE2D_H_
#define MPM_RIGID_PIPELINE2D_H_

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <Eigen/Dense>

namespace mpm {
namespace pipeline {

//! Circular pipe section represented per unit out-of-plane length.
struct RigidPipelineSection2D {
  double outer_radius{0.};
  double wall_thickness{0.};
  double density{0.};
  double contents_mass_per_length{0.};
  double contents_inertia_per_length{0.};
};

//! Result of one soil-particle/contact-point interaction with the pipe.
struct RigidPipelineContact2D {
  Eigen::Vector2d soil_force{Eigen::Vector2d::Zero()};
  Eigen::Vector2d pipe_force{Eigen::Vector2d::Zero()};
  Eigen::Vector2d application_point{Eigen::Vector2d::Zero()};
  double pipe_moment{0.};
  double penetration{0.};
  bool active{false};
};

//! Penalty contact between a circular rigid pipe and one soil particle.
//! Forces are expressed per unit out-of-plane pipe length.
inline RigidPipelineContact2D compute_rigid_pipeline_contact(
    const Eigen::Vector2d& particle_position,
    const Eigen::Vector2d& particle_velocity,
    const Eigen::Vector2d& pipe_centre,
    const Eigen::Vector2d& pipe_velocity, double pipe_angular_velocity,
    double pipe_radius, double particle_radius, double normal_penalty,
    double normal_damping, double tangential_damping,
    double friction_coefficient) {
  RigidPipelineContact2D result;
  const Eigen::Vector2d offset = particle_position - pipe_centre;
  const double distance = offset.norm();
  const double contact_distance = pipe_radius + particle_radius;
  if (!particle_position.allFinite() || !particle_velocity.allFinite() ||
      !pipe_centre.allFinite() || !pipe_velocity.allFinite() ||
      !std::isfinite(pipe_angular_velocity) || !std::isfinite(pipe_radius) ||
      !std::isfinite(particle_radius) || !std::isfinite(normal_penalty) ||
      !std::isfinite(normal_damping) || !std::isfinite(tangential_damping) ||
      !std::isfinite(friction_coefficient) || pipe_radius <= 0. ||
      particle_radius < 0. || normal_penalty <= 0. || normal_damping < 0. ||
      tangential_damping < 0. || friction_coefficient < 0.)
    throw std::invalid_argument("Rigid pipeline contact input is invalid");
  if (distance <= 1.E-14 || distance >= contact_distance) return result;

  const Eigen::Vector2d normal = offset / distance;
  const Eigen::Vector2d tangent(-normal.y(), normal.x());
  result.application_point = pipe_centre + pipe_radius * normal;
  const Eigen::Vector2d surface_radius = result.application_point - pipe_centre;
  const Eigen::Vector2d surface_velocity(
      pipe_velocity.x() - pipe_angular_velocity * surface_radius.y(),
      pipe_velocity.y() + pipe_angular_velocity * surface_radius.x());
  const Eigen::Vector2d relative_velocity =
      particle_velocity - surface_velocity;
  const double normal_velocity = relative_velocity.dot(normal);
  result.penetration = contact_distance - distance;
  const double normal_force = normal_penalty * result.penetration +
                              normal_damping * std::max(0., -normal_velocity);
  if (normal_force <= 0.) return result;

  const double tangential_velocity = relative_velocity.dot(tangent);
  const double friction_limit = friction_coefficient * normal_force;
  const double tangential_force =
      std::max(-friction_limit,
               std::min(friction_limit,
                        -tangential_damping * tangential_velocity));
  result.soil_force = normal_force * normal + tangential_force * tangent;
  result.pipe_force = -result.soil_force;
  result.pipe_moment = surface_radius.x() * result.pipe_force.y() -
                       surface_radius.y() * result.pipe_force.x();
  result.active = true;
  return result;
}

//! Free rigid circular pipeline with translation and in-plane rotation.
//! \brief The pipe owns only three kinematic degrees of freedom (ux, uy,
//! theta). Surface markers follow the rigid motion and cannot ovalise.
class RigidPipeline2D {
 public:
  RigidPipeline2D(const Eigen::Vector2d& centre,
                  const RigidPipelineSection2D& section,
                  unsigned surface_markers)
      : reference_centre_(centre), centre_(centre), section_(section) {
    this->validate_input(surface_markers);
    constexpr double pi = 3.14159265358979323846;
    local_markers_.resize(surface_markers, 2);
    for (unsigned marker = 0; marker < surface_markers; ++marker) {
      const double angle = 2. * pi * marker / surface_markers;
      local_markers_(marker, 0) = section_.outer_radius * std::cos(angle);
      local_markers_(marker, 1) = section_.outer_radius * std::sin(angle);
    }
  }

  double mass_per_length() const {
    constexpr double pi = 3.14159265358979323846;
    const double inner_radius = section_.outer_radius - section_.wall_thickness;
    const double wall_mass = section_.density * pi *
                             (section_.outer_radius * section_.outer_radius -
                              inner_radius * inner_radius);
    return wall_mass + section_.contents_mass_per_length;
  }

  double inertia_per_length() const {
    constexpr double pi = 3.14159265358979323846;
    const double inner_radius = section_.outer_radius - section_.wall_thickness;
    const double wall_inertia =
        0.5 * section_.density * pi *
        (std::pow(section_.outer_radius, 4) - std::pow(inner_radius, 4));
    return wall_inertia + section_.contents_inertia_per_length;
  }

  unsigned nmarkers() const {
    return static_cast<unsigned>(local_markers_.rows());
  }
  const Eigen::Vector2d& reference_centre() const { return reference_centre_; }
  const Eigen::Vector2d& centre() const { return centre_; }
  double outer_radius() const { return section_.outer_radius; }
  Eigen::Vector2d displacement() const { return centre_ - reference_centre_; }
  const Eigen::Vector2d& velocity() const { return velocity_; }
  const Eigen::Vector2d& acceleration() const { return acceleration_; }
  double angle() const { return angle_; }
  double angular_velocity() const { return angular_velocity_; }
  double angular_acceleration() const { return angular_acceleration_; }
  const Eigen::Vector2d& resultant_force() const { return resultant_force_; }
  double resultant_moment() const { return resultant_moment_; }

  Eigen::Matrix<double, Eigen::Dynamic, 2> marker_positions() const {
    const Eigen::Matrix2d rotation = this->rotation_matrix();
    Eigen::Matrix<double, Eigen::Dynamic, 2> positions(local_markers_.rows(), 2);
    for (Eigen::Index marker = 0; marker < local_markers_.rows(); ++marker)
      positions.row(marker) =
          (centre_ + rotation * local_markers_.row(marker).transpose()).transpose();
    return positions;
  }

  Eigen::Matrix<double, Eigen::Dynamic, 2> marker_velocities() const {
    const auto positions = this->marker_positions();
    Eigen::Matrix<double, Eigen::Dynamic, 2> velocities(positions.rows(), 2);
    for (Eigen::Index marker = 0; marker < positions.rows(); ++marker) {
      const Eigen::Vector2d radius = positions.row(marker).transpose() - centre_;
      velocities(marker, 0) = velocity_.x() - angular_velocity_ * radius.y();
      velocities(marker, 1) = velocity_.y() + angular_velocity_ * radius.x();
    }
    return velocities;
  }

  //! Add a force per unit pipe length at a world-coordinate point.
  void apply_force(const Eigen::Vector2d& force,
                   const Eigen::Vector2d& application_point) {
    if (!force.allFinite() || !application_point.allFinite())
      throw std::invalid_argument("Rigid pipeline force must be finite");
    const Eigen::Vector2d radius = application_point - centre_;
    resultant_force_ += force;
    resultant_moment_ += radius.x() * force.y() - radius.y() * force.x();
  }

  //! Add a uniform body acceleration, for example gravity.
  void apply_body_acceleration(const Eigen::Vector2d& body_acceleration) {
    if (!body_acceleration.allFinite())
      throw std::invalid_argument("Rigid pipeline body acceleration must be finite");
    resultant_force_ += this->mass_per_length() * body_acceleration;
  }

  //! Add a free moment per unit pipe length about the pipe centre.
  void apply_moment(double moment) {
    if (!std::isfinite(moment))
      throw std::invalid_argument("Rigid pipeline moment must be finite");
    resultant_moment_ += moment;
  }

  void clear_forces() {
    resultant_force_.setZero();
    resultant_moment_ = 0.;
  }

  //! Advance the three rigid-body degrees of freedom with symplectic Euler.
  void advance(double dt, double translational_damping = 0.,
               double rotational_damping = 0.) {
    if (!std::isfinite(dt) || dt <= 0.)
      throw std::invalid_argument("Rigid pipeline time step must be positive");
    if (!std::isfinite(translational_damping) || translational_damping < 0. ||
        !std::isfinite(rotational_damping) || rotational_damping < 0.)
      throw std::invalid_argument("Rigid pipeline damping must be nonnegative");

    acceleration_ = resultant_force_ / this->mass_per_length() -
                    translational_damping * velocity_;
    angular_acceleration_ =
        resultant_moment_ / this->inertia_per_length() -
        rotational_damping * angular_velocity_;
    velocity_ += dt * acceleration_;
    angular_velocity_ += dt * angular_acceleration_;
    centre_ += dt * velocity_;
    angle_ += dt * angular_velocity_;
    this->clear_forces();
  }

  void set_state(const Eigen::Vector2d& centre, double angle,
                 const Eigen::Vector2d& velocity, double angular_velocity) {
    if (!centre.allFinite() || !velocity.allFinite() || !std::isfinite(angle) ||
        !std::isfinite(angular_velocity))
      throw std::invalid_argument("Rigid pipeline state must be finite");
    centre_ = centre;
    angle_ = angle;
    velocity_ = velocity;
    angular_velocity_ = angular_velocity;
    acceleration_.setZero();
    angular_acceleration_ = 0.;
    this->clear_forces();
  }

 private:
  void validate_input(unsigned surface_markers) const {
    if (!reference_centre_.allFinite())
      throw std::invalid_argument("Rigid pipeline centre must be finite");
    if (surface_markers < 3)
      throw std::invalid_argument("Rigid pipeline needs at least three markers");
    if (!std::isfinite(section_.outer_radius) || section_.outer_radius <= 0. ||
        !std::isfinite(section_.wall_thickness) ||
        section_.wall_thickness <= 0. ||
        section_.wall_thickness >= section_.outer_radius ||
        !std::isfinite(section_.density) || section_.density <= 0. ||
        !std::isfinite(section_.contents_mass_per_length) ||
        section_.contents_mass_per_length < 0. ||
        !std::isfinite(section_.contents_inertia_per_length) ||
        section_.contents_inertia_per_length < 0.)
      throw std::invalid_argument("Rigid pipeline section is invalid");
  }

  Eigen::Matrix2d rotation_matrix() const {
    Eigen::Matrix2d rotation;
    rotation << std::cos(angle_), -std::sin(angle_), std::sin(angle_),
        std::cos(angle_);
    return rotation;
  }

 private:
  Eigen::Vector2d reference_centre_;
  Eigen::Vector2d centre_;
  RigidPipelineSection2D section_;
  Eigen::Matrix<double, Eigen::Dynamic, 2> local_markers_;
  Eigen::Vector2d velocity_{Eigen::Vector2d::Zero()};
  Eigen::Vector2d acceleration_{Eigen::Vector2d::Zero()};
  double angle_{0.};
  double angular_velocity_{0.};
  double angular_acceleration_{0.};
  Eigen::Vector2d resultant_force_{Eigen::Vector2d::Zero()};
  double resultant_moment_{0.};
};

}  // namespace pipeline
}  // namespace mpm

#endif  // MPM_RIGID_PIPELINE2D_H_
