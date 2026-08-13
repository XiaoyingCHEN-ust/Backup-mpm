#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "material/mohr_coulomb.h"
#include "material/sanisand.h"

// Rate-independent, single-material-point comparison for constitutive
// mechanism diagnostics.  This is deliberately not a replacement for a
// field-scale boundary-value problem.  Both materials start from the same
// admissible isotropic compression and receive the exact same prescribed,
// constant-volume cyclic simple-shear history.
namespace {

using Vector6d = Eigen::Matrix<double, 6, 1>;

constexpr double kPorosity = 0.485;
constexpr double kVoidRatio = kPorosity / (1. - kPorosity);
constexpr double kDenseVoidRatio = 0.735;
constexpr double kDensePorosity = kDenseVoidRatio / (1. + kDenseVoidRatio);
constexpr double kInitialPressure = 3000.;
constexpr double kGammaAmplitude = 4.E-4;
constexpr unsigned kCycles = 6;
constexpr unsigned kStepsPerQuarterCycle = 100;
constexpr unsigned kStepsPerCycle = 4 * kStepsPerQuarterCycle;
constexpr unsigned kTotalSteps = kCycles * kStepsPerCycle;

Json sanisand_properties(double porosity = kPorosity) {
  return Json{{"density", 2650.0},
              {"porosity", porosity},
              {"intrinsic_permeability", 9.79E-12},
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

Json mohr_coulomb_properties() {
  // E and nu reproduce the SANISAND initial elastic G and K at p'=3 kPa,
  // n=0.485.  The calibration is defined in prepare_study.py too; duplicating
  // the numeric result here makes this executable independent and auditable.
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
              {"porosity", kPorosity},
              {"intrinsic_permeability", 9.79E-12}};
}

double target_gamma(unsigned step) {
  const unsigned within_cycle = step % kStepsPerCycle;
  const unsigned quarter = within_cycle / kStepsPerQuarterCycle;
  const double fraction =
      static_cast<double>(within_cycle % kStepsPerQuarterCycle) /
      kStepsPerQuarterCycle;
  switch (quarter) {
    case 0:
      return kGammaAmplitude * fraction;
    case 1:
      return kGammaAmplitude * (1. - fraction);
    case 2:
      return -kGammaAmplitude * fraction;
    default:
      return -kGammaAmplitude * (1. - fraction);
  }
}

double pressure(const Vector6d& stress) {
  return -(stress[0] + stress[1] + stress[2]) / 3.;
}

double q_invariant(const Vector6d& stress) {
  const double mean = (stress[0] + stress[1] + stress[2]) / 3.;
  Vector6d deviator = stress;
  deviator.head<3>().array() -= mean;
  const double j2 =
      0.5 * (deviator.head<3>().squaredNorm() +
             2. * deviator.tail<3>().squaredNorm());
  return std::sqrt(std::max(0., 3. * j2));
}

double tensor_norm(const mpm::dense_map& state,
                   const std::array<const char*, 6>& names) {
  return std::sqrt(state.at(names[0]) * state.at(names[0]) +
                   state.at(names[1]) * state.at(names[1]) +
                   state.at(names[2]) * state.at(names[2]) +
                   2. * (state.at(names[3]) * state.at(names[3]) +
                         state.at(names[4]) * state.at(names[4]) +
                         state.at(names[5]) * state.at(names[5])));
}

const std::array<const char*, 6> kAlphaNames{
    {"AlphaXX", "AlphaYY", "AlphaZZ", "AlphaZY", "AlphaZX", "AlphaXY"}};
const std::array<const char*, 6> kAlphaInitialNames{{
    "AlphaInitialXX", "AlphaInitialYY", "AlphaInitialZZ", "AlphaInitialZY",
    "AlphaInitialZX", "AlphaInitialXY"}};
const std::array<const char*, 6> kFabricNames{
    {"ZXX", "ZYY", "ZZZ", "ZZY", "ZZX", "ZXY"}};

struct Audit {
  double sanisand_min_pressure{std::numeric_limits<double>::max()};
  double sanisand_max_pressure{std::numeric_limits<double>::lowest()};
  double mc_min_pressure{std::numeric_limits<double>::max()};
  double mc_max_pressure{std::numeric_limits<double>::lowest()};
  double maximum_sanisand_void_ratio_drift{0.};
  double maximum_mc_positive_yield_residual{0.};
  double maximum_alpha_norm{0.};
  double maximum_fabric_norm{0.};
  double final_sanisand_plastic_strain{0.};
  double final_mc_plastic_strain{0.};
  double dense_minimum_pressure{std::numeric_limits<double>::max()};
  double dense_final_pressure{0.};
  double dense_maximum_fabric_norm{0.};
  double dense_maximum_void_ratio_drift{0.};
};

void emit_header() {
  std::cout
      << "record_type,material,path_id,step,cycle_time,cycle_index,"
         "cycle_coordinate,reversal_count,loading_direction,gamma_xy,"
         "dgamma_xy,volumetric_strain,deviatoric_strain,"
         "sigma_xx_pa,sigma_yy_pa,sigma_zz_pa,"
         "tau_xy_pa,p_effective_pa,q_pa,signed_q_pa,q_over_p,"
         "plastic_strain,void_ratio,alpha_norm,alpha_initial_norm,"
         "fabric_norm,alpha_xy,alpha_initial_xy,fabric_xy,"
         "mc_tension_residual_pa,mc_shear_residual_pa,"
         "reference_pressure_pa,measured_tangent_shear_pa\n";
}

void emit_path_row(const char* material, unsigned step, unsigned reversals,
                   int loading_direction, double gamma, double dgamma,
                   const Vector6d& stress, double plastic_strain,
                   double void_ratio, double alpha_norm,
                   double alpha_initial_norm, double fabric_norm,
                   double alpha_xy, double alpha_initial_xy, double fabric_xy,
                   double mc_tension, double mc_shear) {
  const double p = pressure(stress);
  const double q = q_invariant(stress);
  const double signed_q = std::copysign(q, stress[3]);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const unsigned cycle_index =
      step == 0 ? 0 : std::min((step - 1) / kStepsPerCycle, kCycles - 1);
  const double cycle_coordinate =
      step == 0
          ? 0.
          : static_cast<double>((step - 1) % kStepsPerCycle + 1) /
                kStepsPerCycle;
  std::cout << "path," << material << ",loose_cyclic_simple_shear," << step
            << ','
            << static_cast<double>(step) / kStepsPerCycle << ','
            << cycle_index << ',' << cycle_coordinate << ',' << reversals
            << ',' << loading_direction << ',' << gamma << ',' << dgamma
            << ",0,0," << stress[0] << ',' << stress[1] << ',' << stress[2]
            << ',' << stress[3] << ',' << p << ',' << q << ',' << signed_q
            << ',' << q / p << ',' << plastic_strain << ',' << void_ratio
            << ',' << alpha_norm << ',' << alpha_initial_norm << ','
            << fabric_norm << ',' << alpha_xy << ',' << alpha_initial_xy
            << ',' << fabric_xy << ',' << mc_tension << ',' << mc_shear << ','
            << nan << ',' << nan << '\n';
}

void emit_dense_probe_row(unsigned step, double deviatoric_strain,
                          const Vector6d& stress,
                          const mpm::dense_map& state) {
  const double p = pressure(stress);
  const double q = q_invariant(stress);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::cout << "state_probe,SANISAND,dense_monotonic_constant_volume," << step
            << ',' << static_cast<double>(step) / 200. << ",-1," << nan
            << ",0,1," << nan << ',' << nan << ",0," << deviatoric_strain
            << ',' << stress[0] << ',' << stress[1] << ',' << stress[2] << ','
            << stress[3] << ',' << p << ',' << q << ',' << q << ',' << q / p
            << ',' << state.at("eps_p_q") << ',' << state.at("void_ratio")
            << ',' << tensor_norm(state, kAlphaNames) << ','
            << tensor_norm(state, kAlphaInitialNames) << ','
            << tensor_norm(state, kFabricNames) << ',' << state.at("AlphaXY")
            << ',' << state.at("AlphaInitialXY") << ',' << state.at("ZXY")
            << ',' << nan << ',' << nan << ',' << nan << ',' << nan << '\n';
}

template <typename Material>
double elastic_shear_tangent(Material* material, double reference_pressure) {
  Vector6d stress = Vector6d::Zero();
  stress.head<3>().setConstant(-reference_pressure);
  mpm::dense_map state;
  if constexpr (std::is_same_v<Material, mpm::Sanisand<2>>)
    state = material->initialise_state_variables_from_particle(kPorosity,
                                                               stress);
  else
    state = material->initialise_state_variables_from_particle(kPorosity);
  Vector6d dstrain = Vector6d::Zero();
  constexpr double probe = 1.E-8;
  dstrain[3] = probe;
  const Vector6d updated =
      material->compute_stress(stress, dstrain, nullptr, &state);
  return std::abs(updated[3] - stress[3]) / probe;
}

void emit_tangent_row(const char* material, double reference_pressure,
                      double tangent) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::cout << "tangent_audit," << material << ",pressure_tangent_probe";
  // Columns 3--30 are path-only.  Preserve a rectangular CSV so parsers can
  // distinguish the two record types without special tokenisation.
  for (unsigned column = 3; column <= 30; ++column)
    std::cout << ',' << ((column == 3 || column == 5) ? -1. : nan);
  std::cout << ',' << reference_pressure << ',' << tangent << '\n';
}

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error("cyclic comparison audit: " + message);
}

}  // namespace

int main() {
  try {
    std::cout << std::setprecision(17);
    mpm::Sanisand<2> sanisand(0, sanisand_properties());
    mpm::MohrCoulomb<2> mohr_coulomb(1, mohr_coulomb_properties());

    Vector6d sanisand_stress = Vector6d::Zero();
    sanisand_stress.head<3>().setConstant(-kInitialPressure);
    Vector6d mc_stress = sanisand_stress;
    auto sanisand_state = sanisand.initialise_state_variables_from_particle(
        kPorosity, sanisand_stress);
    auto mc_state =
        mohr_coulomb.initialise_state_variables_from_particle(kPorosity);

    // An unchanged zero-increment state and non-positive MC yield residuals
    // are direct, executable checks that the common initial state is accepted.
    const Vector6d zero = Vector6d::Zero();
    require((sanisand.compute_stress(sanisand_stress, zero, nullptr,
                                     &sanisand_state) -
             sanisand_stress)
                .cwiseAbs()
                .maxCoeff() < 1.E-10,
            "SANISAND rejects the common compressed initial state");
    Eigen::Matrix<double, 2, 1> initial_mc_yield;
    mohr_coulomb.compute_stress_invariants(mc_stress, &mc_state);
    mohr_coulomb.compute_yield_state(&initial_mc_yield, mc_state);
    require(initial_mc_yield.maxCoeff() <= 0.,
            "Mohr-Coulomb rejects the common compressed initial state");

    emit_header();
    Audit audit;
    unsigned reversals = 0;
    int previous_direction = 1;
    double previous_gamma = 0.;

    for (unsigned step = 0; step <= kTotalSteps; ++step) {
      const double gamma = step == kTotalSteps ? 0. : target_gamma(step);
      const double dgamma = step == 0 ? 0. : gamma - previous_gamma;
      int direction = previous_direction;
      if (dgamma > 0.) direction = 1;
      if (dgamma < 0.) direction = -1;
      if (step > 1 && direction != previous_direction) ++reversals;

      if (step > 0) {
        Vector6d dstrain = Vector6d::Zero();
        dstrain[3] = dgamma;
        sanisand_stress = sanisand.compute_stress(
            sanisand_stress, dstrain, nullptr, &sanisand_state);
        mc_stress = mohr_coulomb.compute_stress(mc_stress, dstrain, nullptr,
                                                &mc_state);
      }

      Eigen::Matrix<double, 2, 1> mc_yield;
      mohr_coulomb.compute_stress_invariants(mc_stress, &mc_state);
      mohr_coulomb.compute_yield_state(&mc_yield, mc_state);
      const double alpha_norm = tensor_norm(sanisand_state, kAlphaNames);
      const double alpha_initial_norm =
          tensor_norm(sanisand_state, kAlphaInitialNames);
      const double fabric_norm = tensor_norm(sanisand_state, kFabricNames);
      emit_path_row("SANISAND", step, reversals, direction, gamma, dgamma,
                    sanisand_stress, sanisand_state.at("eps_p_q"),
                    sanisand_state.at("void_ratio"), alpha_norm,
                    alpha_initial_norm, fabric_norm,
                    sanisand_state.at("AlphaXY"),
                    sanisand_state.at("AlphaInitialXY"),
                    sanisand_state.at("ZXY"),
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN());
      emit_path_row("MohrCoulomb", step, reversals, direction, gamma, dgamma,
                    mc_stress, mc_state.at("pdstrain"), kVoidRatio, 0., 0.,
                    0., 0., 0., 0., mc_yield[0], mc_yield[1]);

      require(sanisand_stress.allFinite() && mc_stress.allFinite(),
              "non-finite stress");
      require(pressure(sanisand_stress) > 100. && pressure(mc_stress) > 0.,
              "effective compression was lost");
      audit.sanisand_min_pressure =
          std::min(audit.sanisand_min_pressure, pressure(sanisand_stress));
      audit.sanisand_max_pressure =
          std::max(audit.sanisand_max_pressure, pressure(sanisand_stress));
      audit.mc_min_pressure =
          std::min(audit.mc_min_pressure, pressure(mc_stress));
      audit.mc_max_pressure =
          std::max(audit.mc_max_pressure, pressure(mc_stress));
      audit.maximum_sanisand_void_ratio_drift =
          std::max(audit.maximum_sanisand_void_ratio_drift,
                   std::abs(sanisand_state.at("void_ratio") - kVoidRatio));
      audit.maximum_mc_positive_yield_residual =
          std::max(audit.maximum_mc_positive_yield_residual,
                   std::max(0., mc_yield.maxCoeff()));
      audit.maximum_alpha_norm = std::max(audit.maximum_alpha_norm, alpha_norm);
      audit.maximum_fabric_norm =
          std::max(audit.maximum_fabric_norm, fabric_norm);
      audit.final_sanisand_plastic_strain = sanisand_state.at("eps_p_q");
      audit.final_mc_plastic_strain = mc_state.at("pdstrain");

      previous_gamma = gamma;
      previous_direction = direction;
    }

    // A supplemental dense-state probe exercises the opposite dilatancy
    // branch and fabric tensor.  It is labelled separately and is not part of
    // the two-model cyclic comparison above.
    mpm::Sanisand<2> dense_sanisand(2,
                                    sanisand_properties(kDensePorosity));
    Vector6d dense_stress = Vector6d::Zero();
    dense_stress.head<3>().setConstant(-kInitialPressure);
    auto dense_state = dense_sanisand.initialise_state_variables_from_particle(
        kDensePorosity, dense_stress);
    Vector6d dense_increment = Vector6d::Zero();
    dense_increment << 5.E-6, -1.E-5, 5.E-6, 0., 0., 0.;
    emit_dense_probe_row(0, 0., dense_stress, dense_state);
    for (unsigned step = 1; step <= 200; ++step) {
      dense_stress = dense_sanisand.compute_stress(
          dense_stress, dense_increment, nullptr, &dense_state);
      emit_dense_probe_row(step, step * 1.E-5, dense_stress, dense_state);
      audit.dense_minimum_pressure =
          std::min(audit.dense_minimum_pressure, pressure(dense_stress));
      audit.dense_final_pressure = pressure(dense_stress);
      audit.dense_maximum_fabric_norm =
          std::max(audit.dense_maximum_fabric_norm,
                   tensor_norm(dense_state, kFabricNames));
      audit.dense_maximum_void_ratio_drift =
          std::max(audit.dense_maximum_void_ratio_drift,
                   std::abs(dense_state.at("void_ratio") - kDenseVoidRatio));
    }

    std::array<double, 3> sanisand_tangent{};
    std::array<double, 3> mc_tangent{};
    constexpr std::array<double, 3> pressures{{1500., 3000., 6000.}};
    for (unsigned i = 0; i < pressures.size(); ++i) {
      sanisand_tangent[i] =
          elastic_shear_tangent(&sanisand, pressures[i]);
      mc_tangent[i] =
          elastic_shear_tangent(&mohr_coulomb, pressures[i]);
      emit_tangent_row("SANISAND", pressures[i], sanisand_tangent[i]);
      emit_tangent_row("MohrCoulomb", pressures[i], mc_tangent[i]);
    }

    require(reversals >= 2 * kCycles - 1,
            "the imposed history did not contain the registered reversals");
    require(std::abs(previous_gamma) < 1.E-15,
            "the imposed cyclic path did not close at zero strain");
    require(audit.maximum_sanisand_void_ratio_drift < 1.E-10,
            "SANISAND violated the constant-volume path");
    require(audit.maximum_mc_positive_yield_residual <= 1.1E-6,
            "Mohr-Coulomb left an active yield surface");
    require(audit.final_sanisand_plastic_strain > 1.E-5 &&
                audit.final_mc_plastic_strain > 1.E-5,
            "the comparison did not mobilise both plastic models");
    require(audit.maximum_alpha_norm > 1.E-4,
            "SANISAND back-stress did not evolve under reversals");
    require(audit.sanisand_min_pressure < 0.95 * kInitialPressure,
            "the selected SANISAND state did not show contractive p' loss");
    require(audit.dense_minimum_pressure < kInitialPressure &&
                audit.dense_final_pressure > 1.2 * kInitialPressure,
            "the dense SANISAND probe did not cross from contraction to dilation");
    require(audit.dense_maximum_fabric_norm > 0.1,
            "the dense SANISAND probe did not evolve fabric");
    require(audit.dense_maximum_void_ratio_drift < 1.E-10,
            "the dense SANISAND probe violated constant volume");
    require(audit.mc_max_pressure - audit.mc_min_pressure < 1.E-6,
            "zero-dilation Mohr-Coulomb changed p' on constant-volume shear");
    require(sanisand_tangent[0] < sanisand_tangent[1] &&
                sanisand_tangent[1] < sanisand_tangent[2],
            "SANISAND tangent did not increase with confinement");
    require(std::abs(mc_tangent[2] - mc_tangent[0]) / mc_tangent[1] < 1.E-12,
            "Mohr-Coulomb elastic tangent unexpectedly depends on pressure");
    require(std::abs(sanisand_tangent[1] - mc_tangent[1]) /
                    mc_tangent[1] <
                1.E-3,
            "the 3 kPa initial tangent calibration is inconsistent");

    std::cerr << std::setprecision(9)
              << "cyclic material-point audit passed: reversals=" << reversals
              << ", SANISAND p' range=[" << audit.sanisand_min_pressure << ','
              << audit.sanisand_max_pressure << "] Pa, MC p' range=["
              << audit.mc_min_pressure << ',' << audit.mc_max_pressure
              << "] Pa, final plastic strains=("
              << audit.final_sanisand_plastic_strain << ','
              << audit.final_mc_plastic_strain << "), max(|alpha|,|Z|)=("
              << audit.maximum_alpha_norm << ',' << audit.maximum_fabric_norm
              << "), dense p' min/final=(" << audit.dense_minimum_pressure
              << ',' << audit.dense_final_pressure << ") Pa, dense |Z|max="
              << audit.dense_maximum_fabric_norm << "\n";
    return 0;
  } catch (const std::exception& exception) {
    std::cerr << exception.what() << '\n';
    return 2;
  }
}
