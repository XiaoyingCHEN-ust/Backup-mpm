#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <cmath>

#include "pipeline/rigid_pipeline2d.h"

namespace {

mpm::pipeline::RigidPipeline2D make_pipeline(unsigned markers = 64) {
  mpm::pipeline::RigidPipelineSection2D section;
  section.outer_radius = 0.06;
  section.wall_thickness = 0.12 / 17.;
  section.density = 959.;
  return mpm::pipeline::RigidPipeline2D(Eigen::Vector2d(0.7, 0.38), section,
                                        markers);
}

}  // namespace

TEST_CASE("Rigid pipeline section has analytical mass and inertia") {
  constexpr double pi = 3.14159265358979323846;
  const auto pipeline = make_pipeline();
  const double inner_radius = 0.06 - 0.12 / 17.;
  const double expected_mass =
      959. * pi * (0.06 * 0.06 - inner_radius * inner_radius);
  const double expected_inertia =
      0.5 * expected_mass * (0.06 * 0.06 + inner_radius * inner_radius);
  REQUIRE(pipeline.mass_per_length() == Approx(expected_mass));
  REQUIRE(pipeline.inertia_per_length() == Approx(expected_inertia));
}

TEST_CASE("Rigid pipeline contact obeys action-reaction and Coulomb limit") {
  const Eigen::Vector2d particle_position(0.764, 0.38);
  const Eigen::Vector2d particle_velocity(0., 1.);
  const Eigen::Vector2d pipe_centre(0.7, 0.38);
  const auto contact = mpm::pipeline::compute_rigid_pipeline_contact(
      particle_position, particle_velocity, pipe_centre,
      Eigen::Vector2d::Zero(), 0., 0.06, 0.005, 238000., 75., 75., 0.15);

  REQUIRE(contact.active);
  REQUIRE(contact.penetration == Approx(0.001));
  REQUIRE((contact.soil_force + contact.pipe_force).isZero(1.E-14));
  REQUIRE(contact.soil_force.x() == Approx(238.));
  REQUIRE(std::abs(contact.soil_force.y()) <=
          0.15 * std::abs(contact.soil_force.x()) + 1.E-14);
  REQUIRE(contact.pipe_moment > 0.);
}

TEST_CASE("Rigid pipeline translates without changing its circular shape") {
  auto pipeline = make_pipeline();
  const Eigen::Vector2d initial_centre = pipeline.centre();
  const auto initial_markers = pipeline.marker_positions();
  pipeline.set_state(initial_centre, 0., Eigen::Vector2d(0.2, -0.1), 0.);
  pipeline.advance(0.01);

  REQUIRE(pipeline.displacement().x() == Approx(0.002));
  REQUIRE(pipeline.displacement().y() == Approx(-0.001));
  const auto markers = pipeline.marker_positions();
  for (Eigen::Index marker = 0; marker < markers.rows(); ++marker) {
    REQUIRE((markers.row(marker).transpose() - pipeline.centre()).norm() ==
            Approx(0.06));
    REQUIRE((markers.row(marker) - initial_markers.row(marker)).transpose()
                .isApprox(pipeline.displacement(), 1.E-14));
  }
}

TEST_CASE("Rigid pipeline force updates centre-of-mass motion") {
  auto pipeline = make_pipeline();
  const Eigen::Vector2d force(10., -4.);
  pipeline.apply_force(force, pipeline.centre());
  pipeline.advance(0.02);

  REQUIRE(pipeline.acceleration().x() ==
          Approx(force.x() / pipeline.mass_per_length()));
  REQUIRE(pipeline.acceleration().y() ==
          Approx(force.y() / pipeline.mass_per_length()));
  REQUIRE(pipeline.angular_acceleration() == Approx(0.));
  REQUIRE(pipeline.resultant_force().isZero());
  REQUIRE(pipeline.resultant_moment() == Approx(0.));
}

TEST_CASE("Rigid pipeline force couple produces rotation only") {
  auto pipeline = make_pipeline();
  const Eigen::Vector2d centre = pipeline.centre();
  pipeline.apply_force(Eigen::Vector2d(0., 5.),
                       centre + Eigen::Vector2d(0.06, 0.));
  pipeline.apply_force(Eigen::Vector2d(0., -5.),
                       centre - Eigen::Vector2d(0.06, 0.));
  pipeline.advance(0.001);

  REQUIRE(pipeline.acceleration().isZero(1.E-14));
  REQUIRE(pipeline.angular_acceleration() ==
          Approx(0.6 / pipeline.inertia_per_length()));
  REQUIRE(pipeline.angle() > 0.);
  for (Eigen::Index marker = 0; marker < pipeline.marker_positions().rows();
       ++marker)
    REQUIRE((pipeline.marker_positions().row(marker).transpose() -
             pipeline.centre())
                .norm() == Approx(0.06));
}
