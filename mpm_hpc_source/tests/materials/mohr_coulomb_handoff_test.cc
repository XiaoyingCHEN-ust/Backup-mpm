#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <cmath>
#include <limits>

#include "material/mohr_coulomb.h"

namespace {

Json mohr_coulomb_properties() {
  return Json{{"density", 2650.0},
              {"youngs_modulus", 9.100807964990605E6},
              {"poisson_ratio", -0.00796252868675392},
              {"friction", 31.14738992150898},
              {"cohesion", 0.0},
              {"residual_friction", 31.14738992150898},
              {"residual_cohesion", 0.0},
              {"dilation", 0.0},
              {"residual_dilation", 0.0},
              {"softening", false},
              {"softening_type", "none"},
              {"peak_pdstrain", 0.01},
              {"residual_pdstrain", 0.05},
              {"tension_cutoff", 0.0},
              {"eta", 10.0},
              {"eta_T", 0.2},
              {"eta_p", 1.E-5},
              {"porosity", 0.485},
              {"intrinsic_permeability", 9.79E-12}};
}

}  // namespace

TEST_CASE("Mohr-Coulomb explicitly supports equilibrium-state handoff") {
  mpm::MohrCoulomb<2> material(0, mohr_coulomb_properties());
  const auto state = material.initialise_state_variables_from_particle(0.485);

  REQUIRE(state.size() == 8);
  REQUIRE(state.at("pdstrain") == Approx(0.0));
  REQUIRE(state.at("cohesion") == Approx(0.0));
  REQUIRE(state.at("psi") == Approx(0.0));
  REQUIRE(state.at("phi") ==
          Approx(31.14738992150898 * std::acos(-1.0) / 180.0));
}

TEST_CASE("Mohr-Coulomb handoff rejects invalid restored porosity") {
  mpm::MohrCoulomb<2> material(0, mohr_coulomb_properties());
  REQUIRE_THROWS(material.initialise_state_variables_from_particle(0.0));
  REQUIRE_THROWS(material.initialise_state_variables_from_particle(1.0));
  REQUIRE_THROWS(material.initialise_state_variables_from_particle(
      std::numeric_limits<double>::quiet_NaN()));
}

TEST_CASE("Mohr-Coulomb return-map convergence targets the active surface") {
  Eigen::Matrix<double, 2, 1> residuals;

  residuals << -100., 5.E-7;
  REQUIRE(mpm::mohrcoulomb::active_surface_converged(
      mpm::mohrcoulomb::FailureState::Shear, residuals, 1.E-6));
  REQUIRE_FALSE(mpm::mohrcoulomb::active_surface_converged(
      mpm::mohrcoulomb::FailureState::Tensile, residuals, 1.E-6));

  // The former one-sided test accepted this deep interior overshoot because
  // both residuals were below +tolerance.
  residuals << -100., -10.;
  REQUIRE_FALSE(mpm::mohrcoulomb::active_surface_converged(
      mpm::mohrcoulomb::FailureState::Shear, residuals, 1.E-6));

  residuals << -5.E-7, -100.;
  REQUIRE(mpm::mohrcoulomb::active_surface_converged(
      mpm::mohrcoulomb::FailureState::Tensile, residuals, 1.E-6));
  residuals[0] = std::numeric_limits<double>::quiet_NaN();
  REQUIRE_FALSE(mpm::mohrcoulomb::active_surface_converged(
      mpm::mohrcoulomb::FailureState::Tensile, residuals, 1.E-6));
}

TEST_CASE("Mohr-Coulomb zero-cohesion apex remains finite") {
  mpm::MohrCoulomb<2> material(0, mohr_coulomb_properties());
  auto state = material.initialise_state_variables_from_particle(0.485);
  const Eigen::Matrix<double, 6, 1> stress =
      Eigen::Matrix<double, 6, 1>::Zero();
  const Eigen::Matrix<double, 6, 1> dstrain =
      Eigen::Matrix<double, 6, 1>::Zero();

  const auto updated =
      material.compute_stress(stress, dstrain, nullptr, &state);

  REQUIRE(updated.allFinite());
  REQUIRE(updated.norm() == Approx(0.0).margin(1.E-14));
  REQUIRE(std::isfinite(state.at("pdstrain")));
  REQUIRE(state.at("pdstrain") == Approx(0.0));
}

TEST_CASE("Mohr-Coulomb states adjacent to zero pressure remain finite") {
  mpm::MohrCoulomb<2> material(0, mohr_coulomb_properties());

  SECTION("sub-tolerance strain increment at the apex") {
    auto state = material.initialise_state_variables_from_particle(0.485);
    const Eigen::Matrix<double, 6, 1> stress =
        Eigen::Matrix<double, 6, 1>::Zero();
    Eigen::Matrix<double, 6, 1> dstrain =
        Eigen::Matrix<double, 6, 1>::Zero();
    dstrain(0) = -1.E-14;

    const auto updated =
        material.compute_stress(stress, dstrain, nullptr, &state);
    REQUIRE(updated.allFinite());
    REQUIRE(updated.norm() < 1.E-6);
    REQUIRE(std::isfinite(state.at("pdstrain")));
  }

  SECTION("representative hydrostatic compression") {
    auto state = material.initialise_state_variables_from_particle(0.485);
    Eigen::Matrix<double, 6, 1> stress =
        Eigen::Matrix<double, 6, 1>::Zero();
    stress(0) = stress(1) = stress(2) = -1.;
    const Eigen::Matrix<double, 6, 1> dstrain =
        Eigen::Matrix<double, 6, 1>::Zero();

    const auto updated =
        material.compute_stress(stress, dstrain, nullptr, &state);
    REQUIRE(updated.allFinite());
    REQUIRE((updated - stress).norm() == Approx(0.0).margin(1.E-14));
    REQUIRE(std::isfinite(state.at("pdstrain")));
  }

  SECTION("small hydrostatic tension returns from the apex") {
    auto state = material.initialise_state_variables_from_particle(0.485);
    Eigen::Matrix<double, 6, 1> stress =
        Eigen::Matrix<double, 6, 1>::Zero();
    stress(0) = stress(1) = stress(2) = 1.E-5;
    const Eigen::Matrix<double, 6, 1> dstrain =
        Eigen::Matrix<double, 6, 1>::Zero();

    const auto updated =
        material.compute_stress(stress, dstrain, nullptr, &state);
    REQUIRE(updated.allFinite());
    REQUIRE(updated.cwiseAbs().maxCoeff() <= 1.E-6);
    Eigen::Matrix<double, 2, 1> yield_function;
    material.compute_stress_invariants(updated, &state);
    material.compute_yield_state(&yield_function, state);
    REQUIRE(yield_function(0) < 1.E-6);
    REQUIRE(yield_function(1) < 1.E-6);
    REQUIRE(std::isfinite(state.at("pdstrain")));
  }
}

TEST_CASE("Mohr-Coulomb equilibrium handoff returns to the shear surface") {
  mpm::MohrCoulomb<2> material(0, mohr_coulomb_properties());
  auto state = material.initialise_state_variables_from_particle(0.485);
  Eigen::Matrix<double, 6, 1> stress;
  // Representative restored stress from particle 6454 in the shortened HS
  // equilibrium.  It exercises the near-zero shear residual that the return
  // map must accept without confusing it with elastic classification.
  stress << -9.68000132, -1016.36570135, -375.78115639, -158.56789724, 0.,
      0.;
  const Eigen::Matrix<double, 6, 1> dstrain =
      Eigen::Matrix<double, 6, 1>::Zero();

  const auto updated =
      material.compute_stress(stress, dstrain, nullptr, &state);

  REQUIRE(updated.allFinite());
  REQUIRE(std::isfinite(state.at("pdstrain")));
  Eigen::Matrix<double, 2, 1> yield_function;
  material.compute_stress_invariants(updated, &state);
  material.compute_yield_state(&yield_function, state);
  REQUIRE(yield_function.allFinite());
  REQUIRE(yield_function(0) <= 1.E-6);
  REQUIRE(std::fabs(yield_function(1)) <= 1.E-6);
}
