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
