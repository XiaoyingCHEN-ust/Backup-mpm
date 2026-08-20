#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "material/sanisand.h"

// Single-material-point validation driver for the in-repository SANISAND
// implementation.  The loading histories are prescribed here and duplicated
// in the independent Pisano UDSM adapter.  This executable does not link to,
// or contain, any code from the GPL-licensed reference implementation.
namespace {

using Vector6d = Eigen::Matrix<double, 6, 1>;

constexpr double kInitialPressurePa = 294000.;
constexpr double kInitialVoidRatio = 0.808;
constexpr double kPorosity = kInitialVoidRatio / (1. + kInitialVoidRatio);
constexpr double kCyclicAxialAmplitude = 0.0025;
constexpr unsigned kStepsPerCycle = 400;
constexpr unsigned kCycles = 8;
constexpr unsigned kLiteratureCycles = 4;
constexpr double kLiteratureQAmplitudePa = 114200.;
constexpr unsigned kMonotonicSteps = 1000;
constexpr double kMonotonicFinalAxialStrain = -0.02;

enum class Parameterization {
  LegacyMpm,
  Liu2019CriticalState,
  Liu2019Elasticity,
  Liu2019Parent
};

const char* implementation_name(Parameterization parameterization) {
  switch (parameterization) {
    case Parameterization::Liu2019CriticalState:
      return "mpm_sanisand04_direct_csl";
    case Parameterization::Liu2019Elasticity:
      return "mpm_sanisand04_direct_elasticity";
    case Parameterization::Liu2019Parent:
      return "mpm_sanisand04_liu2019_parent";
    default:
      return "mpm_sanisand04_style";
  }
}

Json properties(double stress_tolerance, Parameterization parameterization) {
  // Toyoura parameters from Liu et al. (2019), Table A1.  The current model
  // uses a different critical-state-line and bulk-modulus parameterisation.
  // Nc/Lambda/alpha_c are a least-squares mapping of
  // e_c=0.934-0.019(p/Patm)^0.7 over p'=10--500 kPa; K0 matches nu=0.05 at
  // the literature initial state p'=294 kPa, e=0.808.  Those approximations
  // are deliberately exposed in the CSV/report rather than hidden.
  Json result{{"density", 2650.0},
              {"porosity", kPorosity},
              {"intrinsic_permeability", 1.E-12},
              {"theta", 0.0},
              {"thermal_expansivity", 0.0},
              {"thermal_conductivity", 1.0},
              {"specific_heat", 1.0},
              {"density_ratio", 1.0},
              {"initial_temperature", 0.0},
              {"p_ref", 100000.0},
              {"G0", 125.0},
              {"K0", 93.84984137456378},
              {"Mc", 1.25},
              {"Lambda", 0.06624461855271319},
              {"N_c", 2.1625333788559744},
              {"alpha_c", 333424.70764453075},
              {"n_b", 1.1},
              {"ch", 0.968},
              {"n_d", 3.5},
              {"h0", 7.05},
              {"A0", 0.704},
              {"Me", 0.89},
              {"cz", 600.0},
              {"zmax", 5.0},
              {"m_iso", 0.01},
              {"Patm", 100000.0},
              {"P_min", 10.0},
              {"FTOL", stress_tolerance},
              {"STOL", stress_tolerance},
              {"LTOL", 1.E-8}};
  if (parameterization == Parameterization::Liu2019CriticalState ||
      parameterization == Parameterization::Liu2019Parent) {
    result["nu"] = 0.05;
    result["e0"] = 0.934;
    result["lambda_c"] = 0.019;
    result["xi"] = 0.7;
  }
  if (parameterization == Parameterization::Liu2019Elasticity) {
    result["parameterization"] = "liu2019_elasticity";
    result["nu"] = 0.05;
  } else if (parameterization == Parameterization::Liu2019CriticalState) {
    result["parameterization"] = "liu2019_critical_state";
  } else if (parameterization == Parameterization::Liu2019Parent) {
    result["parameterization"] = "liu2019_parent";
  }
  return result;
}

double pressure(const Vector6d& stress) {
  return -(stress[0] + stress[1] + stress[2]) / 3.;
}

double signed_triaxial_q(const Vector6d& stress) {
  return stress[0] - stress[2];
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

double triangle(double cycle_time) {
  const double fraction = cycle_time - std::floor(cycle_time);
  if (fraction < 0.25) return 4. * fraction;
  if (fraction < 0.75) return 2. - 4. * fraction;
  return -4. + 4. * fraction;
}

Vector6d constant_volume_increment(double axial_increment) {
  Vector6d increment = Vector6d::Zero();
  increment[0] = -0.5 * axial_increment;
  increment[1] = -0.5 * axial_increment;
  increment[2] = axial_increment;
  return increment;
}

struct PointState {
  Vector6d stress{Vector6d::Zero()};
  mpm::dense_map state;
  double axial_strain{0.};
};

PointState initialise(mpm::Sanisand<2>* material) {
  PointState point;
  point.stress.head<3>().setConstant(-kInitialPressurePa);
  point.state = material->initialise_state_variables_from_particle(
      kPorosity, point.stress);
  return point;
}

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error("SANISAND validation: " + message);
}

void emit_header() {
  std::cout
      << "implementation,path_id,tolerance,step,cycle_time,axial_strain,"
         "p_effective_pa,q_signed_pa,ru,void_ratio,plastic_strain,"
         "alpha_norm,alpha_in_norm,alpha_memory_norm,memory_radius,"
         "fabric_norm,controller_target_q_pa,controller_residual_pa,status\n";
}

void emit(Parameterization parameterization, const char* path_id,
          double tolerance, unsigned step,
          double cycle_time, const PointState& point,
          double controller_target =
              std::numeric_limits<double>::quiet_NaN(),
          double controller_residual =
              std::numeric_limits<double>::quiet_NaN(),
          const char* status = "ok") {
  const double p = pressure(point.stress);
  const double q = signed_triaxial_q(point.stress);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  require(std::isfinite(p) && std::isfinite(q) && p > 0.,
          "non-finite or non-compressive material-point state");
  std::cout << implementation_name(parameterization) << ',' << path_id << ','
            << tolerance << ','
            << step << ',' << cycle_time << ',' << point.axial_strain << ','
            << p << ',' << q << ',' << 1. - p / kInitialPressurePa << ','
            << point.state.at("void_ratio") << ','
            << point.state.at("eps_p_q") << ','
            << tensor_norm(point.state, kAlphaNames) << ','
            << tensor_norm(point.state, kAlphaInitialNames) << ',' << nan
            << ',' << nan << ','
            << tensor_norm(point.state, kFabricNames) << ','
            << controller_target << ',' << controller_residual << ','
            << status << '\n';
}

void run_monotonic(double tolerance, Parameterization parameterization) {
  mpm::Sanisand<2> material(0, properties(tolerance, parameterization));
  auto point = initialise(&material);
  emit(parameterization, "toyoura_monotonic_constant_volume", tolerance, 0,
       0., point);
  const double increment = kMonotonicFinalAxialStrain / kMonotonicSteps;
  for (unsigned step = 1; step <= kMonotonicSteps; ++step) {
    point.stress = material.compute_stress(
        point.stress, constant_volume_increment(increment), nullptr,
        &point.state);
    point.axial_strain += increment;
    if (step % 5 == 0)
      emit(parameterization, "toyoura_monotonic_constant_volume", tolerance,
           step, 0., point);
  }
}

void run_cyclic(double tolerance, Parameterization parameterization) {
  mpm::Sanisand<2> material(0, properties(tolerance, parameterization));
  auto point = initialise(&material);
  emit(parameterization, "toyoura_cyclic_constant_volume", tolerance, 0, 0.,
       point);
  double previous_target = 0.;
  const unsigned total_steps = kCycles * kStepsPerCycle;
  for (unsigned step = 1; step <= total_steps; ++step) {
    const double cycle_time = static_cast<double>(step) / kStepsPerCycle;
    const double target = -kCyclicAxialAmplitude * triangle(cycle_time);
    const double increment = target - previous_target;
    previous_target = target;
    point.stress = material.compute_stress(
        point.stress, constant_volume_increment(increment), nullptr,
        &point.state);
    point.axial_strain = target;
    if (step % 5 == 0)
      emit(parameterization, "toyoura_cyclic_constant_volume", tolerance,
           step, cycle_time, point);
  }
}

PointState trial_update(mpm::Sanisand<2>* material, const PointState& source,
                        double axial_increment) {
  PointState trial = source;
  trial.stress = material->compute_stress(
      trial.stress, constant_volume_increment(axial_increment), nullptr,
      &trial.state);
  trial.axial_strain += axial_increment;
  return trial;
}

PointState solve_q_target(mpm::Sanisand<2>* material,
                          const PointState& source, double target_q) {
  // The constitutive update contains adaptive substep branch points close to
  // liquefaction.  The controller is limited to 0.5% of the 114.2 kPa test
  // amplitude; the achieved residual is recorded and audited downstream.
  constexpr double q_tolerance_pa = 0.005 * kLiteratureQAmplitudePa;
  double span = 2.E-5;
  PointState lower;
  PointState upper;
  double lower_residual = 0.;
  double upper_residual = 0.;
  bool bracketed = false;
  for (unsigned expansion = 0; expansion < 12; ++expansion) {
    lower = trial_update(material, source, -span);
    upper = trial_update(material, source, span);
    lower_residual = signed_triaxial_q(lower.stress) - target_q;
    upper_residual = signed_triaxial_q(upper.stress) - target_q;
    if (lower_residual == 0. || upper_residual == 0. ||
        std::signbit(lower_residual) != std::signbit(upper_residual)) {
      bracketed = true;
      break;
    }
    span *= 2.;
  }
  require(bracketed, "could not bracket the prescribed q increment");
  if (std::abs(lower_residual) <= q_tolerance_pa) return lower;
  if (std::abs(upper_residual) <= q_tolerance_pa) return upper;
  for (unsigned iteration = 0; iteration < 60; ++iteration) {
    const double midpoint_increment =
        0.5 * ((lower.axial_strain - source.axial_strain) +
               (upper.axial_strain - source.axial_strain));
    auto midpoint = trial_update(material, source, midpoint_increment);
    const double residual = signed_triaxial_q(midpoint.stress) - target_q;
    if (std::abs(residual) <= q_tolerance_pa) return midpoint;
    if (std::signbit(residual) == std::signbit(lower_residual)) {
      lower = std::move(midpoint);
      lower_residual = residual;
    } else {
      upper = std::move(midpoint);
      upper_residual = residual;
    }
  }
  return std::abs(lower_residual) < std::abs(upper_residual) ? lower : upper;
}

void run_literature_q_controlled(double tolerance,
                                 Parameterization parameterization) {
  mpm::Sanisand<2> material(0, properties(tolerance, parameterization));
  auto point = initialise(&material);
  emit(parameterization, "toyoura_literature_q_controlled", tolerance, 0, 0.,
       point, 0., 0.);
  const unsigned total_steps = kLiteratureCycles * kStepsPerCycle;
  for (unsigned step = 1; step <= total_steps; ++step) {
    const double cycle_time = static_cast<double>(step) / kStepsPerCycle;
    const double target_q = kLiteratureQAmplitudePa * triangle(cycle_time);
    point = solve_q_target(&material, point, target_q);
    const double controller_residual =
        signed_triaxial_q(point.stress) - target_q;
    if (std::abs(controller_residual) >
        0.00501 * kLiteratureQAmplitudePa) {
      emit(parameterization, "toyoura_literature_q_controlled", tolerance,
           step, cycle_time, point, target_q, controller_residual,
           "controller_limit");
      break;
    }
    emit(parameterization, "toyoura_literature_q_controlled", tolerance, step,
         cycle_time, point, target_q, controller_residual);
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    Parameterization parameterization = Parameterization::LegacyMpm;
    if (argc == 2) {
      const std::string argument(argv[1]);
      if (argument == "--direct-csl")
        parameterization = Parameterization::Liu2019CriticalState;
      else if (argument == "--direct-elasticity")
        parameterization = Parameterization::Liu2019Elasticity;
      else if (argument == "--liu2019-parent")
        parameterization = Parameterization::Liu2019Parent;
      else
        throw std::invalid_argument("unknown validation parameterization");
    } else if (argc != 1) {
      throw std::invalid_argument(
          "usage: sanisand_validation_driver [--direct-csl|"
          "--direct-elasticity|--liu2019-parent]");
    }
    std::cout << std::setprecision(17);
    emit_header();
    for (const double tolerance : {1.E-5, 1.E-7}) {
      run_monotonic(tolerance, parameterization);
      run_cyclic(tolerance, parameterization);
      run_literature_q_controlled(tolerance, parameterization);
    }
    return 0;
  } catch (const std::exception& exception) {
    std::cerr << exception.what() << '\n';
    return 1;
  }
}
