#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <utility>

#include "material/sanisand.h"

namespace {

Json sanisand_properties() {
  return Json{{"density", 2650.0},
              {"youngs_modulus", 2.38E7},
              {"poisson_ratio", 0.4},
              {"porosity", 0.3},
              {"intrinsic_permeability", 1.E-10},
              {"theta", 0.0},
              {"thermal_expansivity", 1.E-5},
              {"thermal_conductivity", 1.E5},
              {"specific_heat", 2080.0},
              {"density_ratio", 0.8},
              {"initial_temperature", 0.0},
              {"p_ref", 100000.0},
              {"G0", 125.0},
              {"K0", 150.0},
              {"Mc", 1.25},
              {"Lambda", 0.37},
              {"N_c", 240.902666166616},
              {"alpha_c", 3370000.0},
              {"n_b", 1.25},
              {"ch", 0.968},
              {"n_d", 2.3},
              {"h0", 12.0},
              {"A0", 0.4},
              {"Me", 0.89},
              {"cz", 600.0},
              {"zmax", 4.0},
              {"m_iso", 0.01},
              {"Patm", 100000.0},
              {"P_min", 100.0},
              {"FTOL", 1.E-5},
              {"STOL", 1.E-5},
              {"LTOL", 1.E-6}};
}

Json source_kpa_properties() {
  auto properties = sanisand_properties();
  properties["porosity"] = 0.907 / 1.907;
  properties["N_c"] = 18.7;
  properties["alpha_c"] = 3370.0;
  properties["Patm"] = 100.0;
  properties["P_min"] = 0.1;
  return properties;
}

}  // namespace

TEST_CASE("SANISAND initial state has twenty values") {
  mpm::Sanisand<2> material(0, sanisand_properties());
  const auto state = material.initialise_state_variables();

  const std::array<const char*, 20> state_names{{
      "AlphaXX",        "AlphaYY",        "AlphaZZ",
      "AlphaZY",        "AlphaZX",        "AlphaXY",
      "AlphaInitialXX", "AlphaInitialYY", "AlphaInitialZZ",
      "AlphaInitialZY", "AlphaInitialZX", "AlphaInitialXY",
      "ZXX",            "ZYY",            "ZZZ",
      "ZZY",            "ZZX",            "ZXY",
      "eps_p_q",        "void_ratio",
  }};

  REQUIRE(state.size() == 20);
  for (const auto* name : state_names) REQUIRE(state.count(name) == 1);
  REQUIRE(state.at("AlphaXX") == Approx(0.0));
  REQUIRE(state.at("AlphaInitialXX") == Approx(0.0));
  REQUIRE(state.at("ZXX") == Approx(0.0));
  REQUIRE(state.at("eps_p_q") == Approx(0.0));
  REQUIRE(state.at("void_ratio") == Approx(0.3 / 0.7));
}

TEST_CASE("SANISAND eps_p_q removes one third of the plastic strain trace") {
  using Vector6d = mpm::Sanisand<2>::Vector6d;

  Vector6d hydrostatic;
  hydrostatic << 3., 3., 3., 0., 0., 0.;
  REQUIRE(mpm::sanisand::detail::equivalent_plastic_deviatoric_strain(
              hydrostatic) == Approx(0.).margin(1.E-15));

  Vector6d plastic_strain;
  plastic_strain << 4., 1., 1., 0., 0., 0.;
  const double equivalent =
      mpm::sanisand::detail::equivalent_plastic_deviatoric_strain(
          plastic_strain);
  REQUIRE(equivalent == Approx(2.));

  Vector6d shifted = plastic_strain + 7. * hydrostatic;
  REQUIRE(mpm::sanisand::detail::equivalent_plastic_deviatoric_strain(
              shifted) == Approx(equivalent));

  Vector6d pure_engineering_shear = Vector6d::Zero();
  pure_engineering_shear[5] = 3.;
  REQUIRE(mpm::sanisand::detail::equivalent_plastic_deviatoric_strain(
              pure_engineering_shear) == Approx(std::sqrt(3.)));
}

TEST_CASE("SANISAND rejects missing and invalid parameters") {
  SECTION("missing Patm") {
    auto properties = sanisand_properties();
    properties.erase("Patm");
    REQUIRE_THROWS(mpm::Sanisand<2>(0, properties));
  }

  SECTION("nonpositive Lambda") {
    auto properties = sanisand_properties();
    properties["Lambda"] = 0.0;
    REQUIRE_THROWS(mpm::Sanisand<2>(0, properties));
  }

  SECTION("nonpositive P_min") {
    auto properties = sanisand_properties();
    properties["P_min"] = 0.0;
    REQUIRE_THROWS(mpm::Sanisand<2>(0, properties));
  }

  SECTION("nonpositive critical-state slopes") {
    auto properties = sanisand_properties();
    properties["Mc"] = 0.0;
    REQUIRE_THROWS(mpm::Sanisand<2>(0, properties));
    properties = sanisand_properties();
    properties["Me"] = 0.0;
    REQUIRE_THROWS(mpm::Sanisand<2>(0, properties));
  }

  SECTION("nonpositive tolerances") {
    auto properties = sanisand_properties();
    properties["FTOL"] = 0.0;
    REQUIRE_THROWS(mpm::Sanisand<2>(0, properties));
    properties = sanisand_properties();
    properties["LTOL"] = 0.0;
    REQUIRE_THROWS(mpm::Sanisand<2>(0, properties));
  }

  SECTION("zero porosity") {
    auto properties = sanisand_properties();
    properties["porosity"] = 0.0;
    REQUIRE_THROWS(mpm::Sanisand<2>(0, properties));
  }

  SECTION("unit porosity") {
    auto properties = sanisand_properties();
    properties["porosity"] = 1.0;
    REQUIRE_THROWS(mpm::Sanisand<2>(0, properties));
  }

  SECTION("nonfinite parameter") {
    auto properties = sanisand_properties();
    properties["G0"] = std::numeric_limits<double>::infinity();
    REQUIRE_THROWS(mpm::Sanisand<2>(0, properties));
  }
}

TEST_CASE("SANISAND hydrostatic compression updates stress and void ratio") {
  mpm::Sanisand<2> material(0, sanisand_properties());
  auto state = material.initialise_state_variables();

  Eigen::Matrix<double, 6, 1> stress;
  stress << -100000.0, -100000.0, -100000.0, 0.0, 0.0, 0.0;
  Eigen::Matrix<double, 6, 1> dstrain;
  dstrain << -1.E-7, -1.E-7, -1.E-7, 0.0, 0.0, 0.0;

  const double initial_void_ratio = state.at("void_ratio");
  const auto updated =
      material.compute_stress(stress, dstrain, nullptr, &state);

  REQUIRE(updated.allFinite());
  REQUIRE(updated[0] < stress[0]);
  REQUIRE(updated[0] == Approx(updated[1]));
  REQUIRE(updated[1] == Approx(updated[2]));
  REQUIRE(state.at("void_ratio") < initial_void_ratio);
}

TEST_CASE("SANISAND low-pressure branch remains elastic") {
  mpm::Sanisand<2> material(0, sanisand_properties());
  auto state = material.initialise_state_variables();

  Eigen::Matrix<double, 6, 1> stress;
  stress << -80.0, -50.0, -20.0, 0.0, 0.0, 0.0;
  Eigen::Matrix<double, 6, 1> dstrain;
  dstrain << -1.E-7, 0.0, 1.E-7, 0.0, 0.0, 0.0;

  const auto updated =
      material.compute_stress(stress, dstrain, nullptr, &state);

  REQUIRE(updated.allFinite());
  REQUIRE(updated[0] < stress[0]);
  REQUIRE(updated[2] > stress[2]);
  REQUIRE(state.at("AlphaXX") == Approx(0.0).margin(1.E-15));
  REQUIRE(state.at("AlphaYY") == Approx(0.0).margin(1.E-15));
  REQUIRE(state.at("AlphaZZ") == Approx(0.0).margin(1.E-15));
  REQUIRE(state.at("ZXX") == Approx(0.0).margin(1.E-15));
  REQUIRE(state.at("eps_p_q") == Approx(0.0).margin(1.E-15));
}

TEST_CASE("SANISAND maps every MPM engineering shear component") {
  mpm::Sanisand<2> material(0, sanisand_properties());
  Eigen::Matrix<double, 6, 1> stress;
  stress << -100000.0, -100000.0, -100000.0, 0.0, 0.0, 0.0;

  for (unsigned shear = 3; shear < 6; ++shear) {
    auto state = material.initialise_state_variables();
    Eigen::Matrix<double, 6, 1> dstrain =
        Eigen::Matrix<double, 6, 1>::Zero();
    dstrain[shear] = 1.E-7;

    const auto updated =
        material.compute_stress(stress, dstrain, nullptr, &state);

    REQUIRE(updated.allFinite());
    REQUIRE(std::abs(updated[shear]) > 0.0);
    for (unsigned other = 3; other < 6; ++other) {
      if (other != shear) REQUIRE(updated[other] == Approx(0.0).margin(1.E-12));
    }
  }
}

TEST_CASE("SANISAND handoff derives void ratio from restored porosity") {
  mpm::Sanisand<2> material(0, sanisand_properties());
  const auto state = material.initialise_state_variables_from_particle(0.42);

  REQUIRE(state.size() == 20);
  REQUIRE(state.at("void_ratio") == Approx(0.42 / 0.58));
  REQUIRE(state.at("AlphaXX") == Approx(0.0));
  REQUIRE(state.at("AlphaInitialXX") == Approx(0.0));
  REQUIRE(state.at("ZXX") == Approx(0.0));
  REQUIRE(state.at("eps_p_q") == Approx(0.0));
}

TEST_CASE("SANISAND handoff rejects invalid restored porosity") {
  mpm::Sanisand<2> material(0, sanisand_properties());
  REQUIRE_THROWS(material.initialise_state_variables_from_particle(0.0));
  REQUIRE_THROWS(material.initialise_state_variables_from_particle(1.0));
  REQUIRE_THROWS(material.initialise_state_variables_from_particle(
      std::numeric_limits<double>::quiet_NaN()));
}

TEST_CASE("SANISAND2D is available through the material factory") {
  unsigned id = 0;
  const auto material =
      Factory<mpm::Material<2>, unsigned, const Json&>::instance()->create(
          "SANISAND2D", std::move(id), sanisand_properties());

  REQUIRE(material != nullptr);
  REQUIRE(material->id() == 0);
  REQUIRE(material->initialise_state_variables().size() == 20);
}

TEST_CASE("SANISAND monotonic constant-volume loading evolves plastic state") {
  auto properties = sanisand_properties();
  properties["porosity"] = 0.96 / 1.96;
  mpm::Sanisand<2> material(0, properties);
  auto state = material.initialise_state_variables();

  Eigen::Matrix<double, 6, 1> stress;
  stress << -500000.0, -500000.0, -500000.0, 0.0, 0.0, 0.0;
  Eigen::Matrix<double, 6, 1> dstrain;
  dstrain << 5.E-5, -1.E-4, 5.E-5, 0.0, 0.0, 0.0;

  for (unsigned step = 0; step < 25; ++step) {
    stress = material.compute_stress(stress, dstrain, nullptr, &state);
    REQUIRE(stress.allFinite());
  }

  const double alpha_norm =
      std::sqrt(std::pow(state.at("AlphaXX"), 2.) +
                std::pow(state.at("AlphaYY"), 2.) +
                std::pow(state.at("AlphaZZ"), 2.));
  REQUIRE(alpha_norm > 0.0);
  REQUIRE(state.at("eps_p_q") > 0.0);
  REQUIRE(state.at("void_ratio") == Approx(0.96).margin(1.E-10));
}

TEST_CASE("SANISAND SI parameters preserve the source kPa calibration") {
  auto source_properties = source_kpa_properties();
  auto si_properties = sanisand_properties();
  si_properties["porosity"] = 0.907 / 1.907;

  mpm::Sanisand<2> source_material(0, source_properties);
  mpm::Sanisand<2> si_material(0, si_properties);
  auto source_state = source_material.initialise_state_variables();
  auto si_state = si_material.initialise_state_variables();

  Eigen::Matrix<double, 6, 1> source_stress;
  source_stress << -1000.0, -1000.0, -1000.0, 0.0, 0.0, 0.0;
  Eigen::Matrix<double, 6, 1> si_stress = 1000.0 * source_stress;
  Eigen::Matrix<double, 6, 1> dstrain;
  dstrain << 5.E-5, -1.E-4, 5.E-5, 0.0, 0.0, 0.0;

  source_stress = source_material.compute_stress(
      source_stress, dstrain, nullptr, &source_state);
  si_stress =
      si_material.compute_stress(si_stress, dstrain, nullptr, &si_state);

  for (unsigned component = 0; component < 6; ++component)
    REQUIRE(si_stress[component] / 1000.0 ==
            Approx(source_stress[component]).epsilon(1.E-8));
  for (const auto& source_value : source_state)
    REQUIRE(si_state.at(source_value.first) ==
            Approx(source_value.second).epsilon(1.E-5).margin(1.E-12));
}

TEST_CASE("SANISAND replays the standalone undrained material path") {
  auto properties = sanisand_properties();
  properties["porosity"] = 0.907 / 1.907;
  mpm::Sanisand<2> material(0, properties);
  auto state = material.initialise_state_variables();

  Eigen::Matrix<double, 6, 1> stress;
  stress << -1.E6, -1.E6, -1.E6, 0.0, 0.0, 0.0;
  Eigen::Matrix<double, 6, 1> dstrain;
  dstrain << 5.E-5, -1.E-4, 5.E-5, 0.0, 0.0, 0.0;

  struct ReferencePoint {
    unsigned step;
    double mean_stress_kpa;
    double deviatoric_stress_kpa;
  };
  const std::array<ReferencePoint, 4> reference{{
      {1, 999.502, 25.1821},
      {10, 971.845, 183.903},
      {100, 594.607, 477.295},
      {2500, 202.497, 253.065},
  }};

  unsigned reference_index = 0;
  // The standalone path is evaluated in kPa while the MPM path is evaluated
  // in Pa. Adaptive substep boundaries can diverge slightly across MSVC and
  // GCC after thousands of increments; 1 kPa keeps the cross-compiler error
  // below 0.5% while remaining sensitive to calibration or mapping defects.
  constexpr double reference_tolerance = 1000.;
  for (unsigned step = 1; step <= reference.back().step; ++step) {
    stress = material.compute_stress(stress, dstrain, nullptr, &state);
    if (step != reference[reference_index].step) continue;

    const double mean_stress = -(stress[0] + stress[1] + stress[2]) / 3.;
    const double deviatoric_stress = -(stress[1] - stress[0]);
    REQUIRE(mean_stress ==
            Approx(1000. * reference[reference_index].mean_stress_kpa)
                .margin(reference_tolerance));
    REQUIRE(deviatoric_stress ==
            Approx(1000. * reference[reference_index].deviatoric_stress_kpa)
                .margin(reference_tolerance));
    REQUIRE(state.at("void_ratio") + 1. == Approx(1.907).margin(1.E-10));
    ++reference_index;
  }
  REQUIRE(reference_index == reference.size());
  REQUIRE(state.at("eps_p_q") > 0.0);
}

TEST_CASE("SANISAND replays the standalone cyclic reversal path") {
  auto properties = sanisand_properties();
  properties["porosity"] = 0.735 / 1.735;
  mpm::Sanisand<2> material(0, properties);
  auto state = material.initialise_state_variables();

  Eigen::Matrix<double, 6, 1> stress;
  stress << -1.E5, -1.E5, -1.E5, 0.0, 0.0, 0.0;
  Eigen::Matrix<double, 6, 1> dstrain;
  dstrain << 5.E-5, -1.E-4, 5.E-5, 0.0, 0.0, 0.0;

  struct ReferencePoint {
    unsigned step;
    double mean_stress_kpa;
    double deviatoric_stress_kpa;
  };
  const std::array<ReferencePoint, 3> reference{{
      {1, 99.8995, 9.96663},
      {10, 97.631, 33.4882},
      {100, 69.2965, 46.5886},
  }};

  unsigned reference_index = 0;
  constexpr double reversal_limit = 50.E3;
  constexpr double reference_tolerance = 500.;
  for (unsigned step = 1; step <= 2500; ++step) {
    const double deviatoric_stress = -(stress[1] - stress[0]);
    if (std::abs(deviatoric_stress) >= reversal_limit) dstrain *= -1.;
    stress = material.compute_stress(stress, dstrain, nullptr, &state);
    if (reference_index >= reference.size() ||
        step != reference[reference_index].step)
      continue;

    const double mean_stress = -(stress[0] + stress[1] + stress[2]) / 3.;
    const double updated_deviatoric_stress = -(stress[1] - stress[0]);
    REQUIRE(mean_stress ==
            Approx(1000. * reference[reference_index].mean_stress_kpa)
                .margin(reference_tolerance));
    REQUIRE(updated_deviatoric_stress ==
            Approx(1000. * reference[reference_index].deviatoric_stress_kpa)
                .margin(reference_tolerance));
    REQUIRE(state.at("void_ratio") + 1. == Approx(1.735).margin(1.E-10));
    ++reference_index;
  }

  const double fabric_norm =
      std::sqrt(std::pow(state.at("ZXX"), 2.) +
                std::pow(state.at("ZYY"), 2.) +
                std::pow(state.at("ZZZ"), 2.));
  const double alpha_initial_norm =
      std::sqrt(std::pow(state.at("AlphaInitialXX"), 2.) +
                std::pow(state.at("AlphaInitialYY"), 2.) +
                std::pow(state.at("AlphaInitialZZ"), 2.));
  const double final_mean_stress = -(stress[0] + stress[1] + stress[2]) / 3.;
  const double final_deviatoric_stress = -(stress[1] - stress[0]);
  REQUIRE(reference_index == reference.size());
  REQUIRE(std::isfinite(final_mean_stress));
  REQUIRE(std::isfinite(final_deviatoric_stress));
  REQUIRE(final_mean_stress > 20.E3);
  REQUIRE(final_mean_stress < 80.E3);
  REQUIRE(std::abs(final_deviatoric_stress) < 60.E3);
  REQUIRE(fabric_norm > 0.0);
  REQUIRE(alpha_initial_norm > 0.0);
}

TEST_CASE("SANISAND failed integration does not commit partial history") {
  mpm::Sanisand<2> material(0, sanisand_properties());
  auto state = material.initialise_state_variables();
  const auto original_state = state;

  Eigen::Matrix<double, 6, 1> stress;
  stress << -100000.0, -100000.0, -100000.0, 0.0, 0.0, 0.0;
  Eigen::Matrix<double, 6, 1> dstrain;
  dstrain << -1.0, -1.0, -1.0, 0.0, 0.0, 0.0;

  REQUIRE_THROWS(material.compute_stress(stress, dstrain, nullptr, &state));
  REQUIRE(state.size() == original_state.size());
  for (const auto& value : original_state)
    REQUIRE(state.at(value.first) == Approx(value.second));
}

TEST_CASE("SANISAND low-confinement path does not overshoot backstress") {
  // Particle 6768 from sanisand_short_v2 at frame 350 entered a tensile
  // principal-stress state and its Alpha norm jumped from 1.36 to 7.01 by
  // frame 400. Replaying the recorded increment in 50 material updates keeps
  // the history bounded when tensile trial states use the elastic safety path.
  constexpr double maximum_admissible_alpha_norm = 2.0;
  mpm::Sanisand<2> material(0, sanisand_properties());
  auto state = material.initialise_state_variables();

  state["AlphaXX"] = 0.4688929820405134;
  state["AlphaYY"] = -1.1024528485115528;
  state["AlphaZZ"] = 0.6335598664710393;
  state["AlphaXY"] = -0.06295600669178116;
  state["AlphaInitialXX"] = 0.5368358497671355;
  state["AlphaInitialYY"] = -1.11137272875813;
  state["AlphaInitialZZ"] = 0.5745368789947805;
  state["AlphaInitialXY"] = 0.27693588937390384;
  state["ZXX"] = -0.05968086935416454;
  state["ZYY"] = 0.16410643628639807;
  state["ZZZ"] = -0.1044255669322318;
  state["ZXY"] = 0.010386314890552842;
  state["eps_p_q"] = 0.0002511120970443136;
  state["void_ratio"] = 0.42924575481360283;

  Eigen::Matrix<double, 6, 1> stress;
  stress << -419.00177520926985, 122.07838074701844,
      -475.50032784659055, 21.945702377874138, 0., 0.;
  Eigen::Matrix<double, 6, 1> strain_increment;
  strain_increment << 3.55760224235077E-6, 1.1486167257227956E-4, 0.,
      -1.3161600936367215E-5, 0., 0.;

  constexpr unsigned replay_steps = 50;
  strain_increment /= replay_steps;
  for (unsigned step = 0; step < replay_steps; ++step) {
    stress =
        material.compute_stress(stress, strain_increment, nullptr, &state);
    const double step_alpha_norm =
        std::sqrt(state.at("AlphaXX") * state.at("AlphaXX") +
                  state.at("AlphaYY") * state.at("AlphaYY") +
                  state.at("AlphaZZ") * state.at("AlphaZZ") +
                  2. * state.at("AlphaXY") * state.at("AlphaXY"));
    const double mean_stress = -(stress[0] + stress[1] + stress[2]) / 3.;
    CAPTURE(step);
    CAPTURE(mean_stress);
    CAPTURE(stress[0]);
    CAPTURE(stress[1]);
    CAPTURE(stress[2]);
    CAPTURE(stress[3]);
    CAPTURE(state.at("AlphaXX"));
    CAPTURE(state.at("AlphaYY"));
    CAPTURE(state.at("AlphaZZ"));
    CAPTURE(state.at("AlphaXY"));
    REQUIRE(step_alpha_norm < maximum_admissible_alpha_norm);
  }
  REQUIRE(stress.allFinite());
  REQUIRE(state.at("AlphaXX") == Approx(0.4688929820405134));
  REQUIRE(state.at("AlphaYY") == Approx(-1.1024528485115528));
  REQUIRE(state.at("AlphaZZ") == Approx(0.6335598664710393));
  REQUIRE(state.at("AlphaXY") == Approx(-0.06295600669178116));
  REQUIRE(state.at("ZXX") == Approx(-0.05968086935416454));
  REQUIRE(state.at("ZYY") == Approx(0.16410643628639807));
  REQUIRE(state.at("ZZZ") == Approx(-0.1044255669322318));
  REQUIRE(state.at("ZXY") == Approx(0.010386314890552842));
  REQUIRE(state.at("eps_p_q") == Approx(0.0002511120970443136));
}
