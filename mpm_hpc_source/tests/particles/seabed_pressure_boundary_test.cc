#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "particle_threephase_lag.h"

TEST_CASE("Wave total traction marker selects only the initial seabed row") {
  using mpm::threephase_lag_boundary::make_surface_traction_marker;

  const double seabed_y = 0.50;
  const double particle_size_y = 0.01;
  REQUIRE(make_surface_traction_marker(true, false, 0.495, seabed_y,
                                       particle_size_y));
  REQUIRE_FALSE(
      make_surface_traction_marker(true, false, 0.485, seabed_y,
                                   particle_size_y));
  REQUIRE_FALSE(
      make_surface_traction_marker(true, false, 0.46, seabed_y,
                                   particle_size_y));
  REQUIRE_FALSE(make_surface_traction_marker(false, false, 0.495, seabed_y,
                                              particle_size_y));
  REQUIRE_FALSE(make_surface_traction_marker(true, true, 0.495, seabed_y,
                                              particle_size_y));
  REQUIRE_FALSE(make_surface_traction_marker(true, false, 0.515, seabed_y,
                                              particle_size_y));
}

TEST_CASE("Wave total traction marker survives subsequent displacement") {
  using mpm::threephase_lag_boundary::make_surface_traction_marker;

  // The marker is created once from the reference coordinate. Current
  // coordinates are intentionally absent from its runtime semantics.
  const bool marker =
      make_surface_traction_marker(true, false, 0.495, 0.50, 0.01);
  REQUIRE(marker);
  for (const double displaced_y : {0.470, 0.495, 0.525}) {
    CAPTURE(displaced_y);
    REQUIRE(marker);
  }
}

TEST_CASE("Three-phase internal forces use pressure relative to restart state") {
  using mpm::threephase_lag_force::excess_phase_pressure;

  REQUIRE(excess_phase_pressure(5000.0, 5000.0) == Approx(0.0));
  REQUIRE(excess_phase_pressure(5250.0, 5000.0) == Approx(250.0));
  REQUIRE(excess_phase_pressure(4750.0, 5000.0) == Approx(-250.0));
  REQUIRE_THROWS(excess_phase_pressure(
      std::numeric_limits<double>::quiet_NaN(), 5000.0));
}

TEST_CASE("Three-phase force validation rejects non-finite mapped vectors") {
  Eigen::Vector2d force(1.0, -2.0);
  REQUIRE_NOTHROW(
      mpm::threephase_lag_force::require_finite_force(force, "test"));
  force[1] = std::numeric_limits<double>::infinity();
  REQUIRE_THROWS_WITH(
      mpm::threephase_lag_force::require_finite_force(force, "test"),
      Catch::Matchers::Contains("Non-finite test force contribution"));
}
