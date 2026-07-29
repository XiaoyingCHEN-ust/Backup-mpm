#ifndef MPM_MATERIAL_SANISAND_H_
#define MPM_MATERIAL_SANISAND_H_

#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! SANISAND class
//! \brief Dafalias-Manzari style anisotropic critical-state sand model
//! \tparam Tdim Dimension (only Tdim = 2 is registered in this phase)
template <unsigned Tdim>
class Sanisand : public Material<Tdim> {
 public:
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id and material properties
  Sanisand(unsigned id, const Json& material_properties);

  ~Sanisand() override = default;

  Sanisand(const Sanisand&) = delete;
  Sanisand& operator=(const Sanisand&) = delete;

  //! Initialise the twenty in-memory SANISAND state variables
  mpm::dense_map initialise_state_variables() override;

  //! Initialise SANISAND history from restored particle porosity
  mpm::dense_map initialise_state_variables_from_particle(
      double porosity) override;

  //! Compute effective stress from an MPM strain increment
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr,
                          mpm::dense_map* state_vars) override;

 protected:
  using Material<Tdim>::console_;
  using Material<Tdim>::id_;
  using Material<Tdim>::properties_;

 private:
  struct State {
    Vector6d alpha{Vector6d::Zero()};
    Vector6d alpha_initial{Vector6d::Zero()};
    Vector6d fabric{Vector6d::Zero()};
    double eps_p_q{0.};
    double void_ratio{0.};
  };

  struct Invariants {
    double mean_stress{0.};
    double j2{0.};
    double j3{0.};
    double q{0.};
    double cos_three_theta{1.};
    Vector6d deviator{Vector6d::Zero()};
    Vector6d stress_ratio{Vector6d::Zero()};
    Vector6d normal{Vector6d::Zero()};
  };

  struct Increment {
    Vector6d stress{Vector6d::Zero()};
    Vector6d alpha{Vector6d::Zero()};
    Vector6d fabric{Vector6d::Zero()};
    double eps_p_q{0.};
  };

  struct IntegrationResult {
    Vector6d stress{Vector6d::Zero()};
    State state;
    unsigned accepted_substeps{0};
    unsigned rejected_substeps{0};
  };

  static double require_finite(const Json& properties,
                               const std::string& key);
  static double require_positive(const Json& properties,
                                 const std::string& key);

  static Vector6d to_internal(const Vector6d& mpm_tensor);
  static Vector6d to_mpm(const Vector6d& internal_tensor);
  static State unpack_state(const mpm::dense_map& state_vars);
  static void pack_state(const State& state, mpm::dense_map* state_vars);
  static bool is_finite(const State& state);

  Invariants compute_invariants(const Vector6d& stress,
                                const Vector6d& alpha) const;
  static bool is_strictly_compressive(const Vector6d& stress);
  Matrix6x6 compute_moduli(const Vector6d& stress, const Vector6d& alpha,
                           double specific_volume) const;
  Vector6d yield_gradient(const Vector6d& stress,
                          const Vector6d& alpha) const;
  Vector6d plastic_gradient(const Vector6d& stress, const Vector6d& alpha,
                            const Vector6d& fabric,
                            double specific_volume) const;
  double critical_void_ratio(double mean_stress) const;
  std::pair<Vector6d, double> plastic_modulus(
      const Vector6d& stress, const Vector6d& alpha,
      const Vector6d& fabric, const Vector6d& alpha_initial,
      double specific_volume) const;
  Increment compute_increment(const Vector6d& stress,
                              const Vector6d& elastic_stress_increment,
                              double substep, const State& state,
                              double specific_volume) const;
  void update_alpha_initial(Vector6d* alpha_initial, const Vector6d& stress,
                            const Vector6d& alpha) const;
  double relative_error(const Vector6d& stress1, const Vector6d& stress2,
                        const Increment& increment1,
                        const Increment& increment2,
                        const Vector6d& alpha2,
                        const Vector6d& fabric2) const;
  IntegrationResult integrate(const Vector6d& stress,
                              const Vector6d& strain_increment,
                              const State& old_state) const;

  double density_{std::numeric_limits<double>::quiet_NaN()};
  double porosity_{std::numeric_limits<double>::quiet_NaN()};
  double G0_{std::numeric_limits<double>::quiet_NaN()};
  double K0_{std::numeric_limits<double>::quiet_NaN()};
  double Mc_{std::numeric_limits<double>::quiet_NaN()};
  double lambda_{std::numeric_limits<double>::quiet_NaN()};
  double N_c_{std::numeric_limits<double>::quiet_NaN()};
  double alpha_c_{std::numeric_limits<double>::quiet_NaN()};
  double n_b_{std::numeric_limits<double>::quiet_NaN()};
  double ch_{std::numeric_limits<double>::quiet_NaN()};
  double n_d_{std::numeric_limits<double>::quiet_NaN()};
  double h0_{std::numeric_limits<double>::quiet_NaN()};
  double A0_{std::numeric_limits<double>::quiet_NaN()};
  double Me_{std::numeric_limits<double>::quiet_NaN()};
  double cz_{std::numeric_limits<double>::quiet_NaN()};
  double zmax_{std::numeric_limits<double>::quiet_NaN()};
  double m_iso_{std::numeric_limits<double>::quiet_NaN()};
  double patm_{std::numeric_limits<double>::quiet_NaN()};
  double p_min_{std::numeric_limits<double>::quiet_NaN()};
  double stol_{std::numeric_limits<double>::quiet_NaN()};
  double ftol_{std::numeric_limits<double>::quiet_NaN()};
  double ltol_{std::numeric_limits<double>::quiet_NaN()};
};

}  // namespace mpm

#include "sanisand.tcc"

#endif  // MPM_MATERIAL_SANISAND_H_
