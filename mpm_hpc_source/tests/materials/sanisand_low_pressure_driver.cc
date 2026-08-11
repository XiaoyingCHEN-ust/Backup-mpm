#include <cmath>
#include <iostream>

#include "material/sanisand.h"

namespace {

Json properties() {
  return Json{{"density", 2650.0},
              {"youngs_modulus", 2.38E7},
              {"poisson_ratio", 0.4},
              {"porosity", 0.485},
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

double norm(const mpm::dense_map& state, const char* xx, const char* yy,
            const char* zz, const char* xy) {
  return std::sqrt(state.at(xx) * state.at(xx) +
                   state.at(yy) * state.at(yy) +
                   state.at(zz) * state.at(zz) +
                   2. * state.at(xy) * state.at(xy));
}

}  // namespace

int main() {
  mpm::Sanisand<2> material(0, properties());
  auto state = material.initialise_state_variables();
  Eigen::Matrix<double, 6, 1> stress;
  stress << -3000., -3000., -3000., 0., 0., 0.;
  Eigen::Matrix<double, 6, 1> base_increment;
  base_increment << 2.5E-6, -5.E-6, 2.5E-6, 0., 0., 0.;

  constexpr unsigned maximum_steps = 5000;
  constexpr double reversal_q = 1800.;
  int direction = 1;
  unsigned reversals = 0;

  std::cout
      << "step,reversals,p_effective_pa,q_pa,eps_p_q,void_ratio,alpha_norm,"
         "fabric_norm\n";
  for (unsigned step = 0; step <= maximum_steps; ++step) {
    const double mean_stress = -(stress[0] + stress[1] + stress[2]) / 3.;
    const double q = -(stress[1] - stress[0]);
    if (step % 10 == 0 || step == maximum_steps)
      std::cout << step << ',' << reversals << ',' << mean_stress << ',' << q
                << ',' << state.at("eps_p_q") << ','
                << state.at("void_ratio") << ','
                << norm(state, "AlphaXX", "AlphaYY", "AlphaZZ", "AlphaXY")
                << ',' << norm(state, "ZXX", "ZYY", "ZZZ", "ZXY") << '\n';
    if (step == maximum_steps) break;

    if (direction > 0 && q >= reversal_q) {
      direction = -1;
      ++reversals;
    } else if (direction < 0 && q <= -reversal_q) {
      direction = 1;
      ++reversals;
    }
    stress = material.compute_stress(
        stress, static_cast<double>(direction) * base_increment, nullptr,
        &state);
    if (!stress.allFinite() || !std::isfinite(state.at("eps_p_q")) ||
        !std::isfinite(state.at("void_ratio"))) {
      std::cerr << "Non-finite state at low confinement step " << step + 1
                << '\n';
      return 2;
    }
  }

  std::cerr << "3 kPa SANISAND diagnostic completed with " << reversals
            << " stress reversals; inspect the CSV before production runs\n";
  return 0;
}
