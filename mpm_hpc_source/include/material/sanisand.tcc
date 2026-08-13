namespace mpm::sanisand::detail {

//! Double contraction of symmetric stress-like tensors stored as
//! [xx, yy, zz, yz, xz, xy]. The shear entries are physical tensor
//! components, so each appears twice in A:B.
inline double stress_double_contraction(
    const Eigen::Matrix<double, 6, 1>& first,
    const Eigen::Matrix<double, 6, 1>& second) {
  return first.head<3>().dot(second.head<3>()) +
         2. * first.tail<3>().dot(second.tail<3>());
}

inline double stress_norm(const Eigen::Matrix<double, 6, 1>& tensor) {
  return std::sqrt(
      std::max(0., stress_double_contraction(tensor, tensor)));
}

//! Norm of a vector that is dual to a physical stress-like Voigt vector.
//! If g = [n_xx, n_yy, n_zz, 2 n_yz, 2 n_xz, 2 n_xy], this recovers
//! sqrt(n:n), so the loading cosine is independent of the chosen axes.
inline double stress_dual_norm(
    const Eigen::Matrix<double, 6, 1>& dual_tensor) {
  const double norm_squared = dual_tensor.head<3>().squaredNorm() +
                              0.5 * dual_tensor.tail<3>().squaredNorm();
  return std::sqrt(std::max(0., norm_squared));
}

inline Eigen::Matrix3d stress_matrix(
    const Eigen::Matrix<double, 6, 1>& tensor) {
  Eigen::Matrix3d matrix;
  matrix << tensor[0], tensor[5], tensor[4], tensor[5], tensor[1], tensor[3],
      tensor[4], tensor[3], tensor[2];
  return matrix;
}

//! For a unit deviatoric direction n, sqrt(6) tr(n^3) is the Lode-angle
//! invariant used by the SANISAND interpolation functions.
inline double lode_cosine(
    const Eigen::Matrix<double, 6, 1>& unit_deviator) {
  const Eigen::Matrix3d matrix = stress_matrix(unit_deviator);
  return std::clamp(std::sqrt(6.) * (matrix * matrix * matrix).trace(), -1.,
                    1.);
}

//! Convert a physical symmetric tensor direction to the vector gradient that
//! is work-conjugate to [sigma_xx, sigma_yy, sigma_zz, tau_yz, tau_xz,
//! tau_xy]. Its shear entries are doubled because dF = gradient . dSigma.
inline Eigen::Matrix<double, 6, 1> stress_direction_to_dual(
    Eigen::Matrix<double, 6, 1> direction) {
  direction.tail<3>() *= 2.;
  return direction;
}

//! Direction of the pressure-scaled yield gradient for
//! f = ||s/p - alpha|| - sqrt(2/3) m.  The projection is n:r (not n:alpha)
//! when m is finite, and the returned vector is dual to the stress-like
//! Voigt storage.
inline Eigen::Matrix<double, 6, 1> yield_direction(
    const Eigen::Matrix<double, 6, 1>& stress_ratio,
    const Eigen::Matrix<double, 6, 1>& unit_normal) {
  auto direction = unit_normal;
  const double projection =
      stress_double_contraction(unit_normal, stress_ratio);
  direction.head<3>().array() -= projection / 3.;
  return stress_direction_to_dual(direction);
}

inline double equivalent_plastic_deviatoric_strain(
    const Eigen::Matrix<double, 6, 1>& plastic_strain) {
  auto deviatoric = plastic_strain;
  const double mean_normal_strain =
      (plastic_strain[0] + plastic_strain[1] + plastic_strain[2]) / 3.;
  for (unsigned i = 0; i < 3; ++i) deviatoric[i] -= mean_normal_strain;
  // MPM stores engineering shear strains gamma_ij = 2 epsilon_ij in the
  // last three Voigt components.  Hence e_dev:e_dev contains one half of
  // their squared engineering values, not their unweighted squared norm.
  const double tensor_double_contraction =
      deviatoric.head<3>().squaredNorm() +
      0.5 * deviatoric.tail<3>().squaredNorm();
  return std::sqrt(
      std::max(0., 2. / 3. * tensor_double_contraction));
}

}  // namespace mpm::sanisand::detail

template <unsigned Tdim>
double mpm::Sanisand<Tdim>::require_finite(const Json& properties,
                                           const std::string& key) {
  double value;
  try {
    value = properties.at(key).template get<double>();
  } catch (const std::exception& exception) {
    throw std::invalid_argument("SANISAND parameter '" + key +
                                "' is missing or invalid: " +
                                exception.what());
  }
  if (!std::isfinite(value))
    throw std::invalid_argument("SANISAND parameter '" + key +
                                "' must be finite");
  return value;
}

template <unsigned Tdim>
double mpm::Sanisand<Tdim>::require_positive(const Json& properties,
                                             const std::string& key) {
  const double value = require_finite(properties, key);
  if (value <= 0.)
    throw std::invalid_argument("SANISAND parameter '" + key +
                                "' must be greater than zero");
  return value;
}

template <unsigned Tdim>
mpm::Sanisand<Tdim>::Sanisand(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  static_assert(Tdim == 2, "Only SANISAND2D is supported in this phase");

  density_ = require_positive(material_properties, "density");
  porosity_ = require_finite(material_properties, "porosity");
  if (porosity_ <= 0. || porosity_ >= 1.)
    throw std::invalid_argument("SANISAND porosity must be between zero and one");

  G0_ = require_positive(material_properties, "G0");
  K0_ = require_positive(material_properties, "K0");
  Mc_ = require_positive(material_properties, "Mc");
  lambda_ = require_positive(material_properties, "Lambda");
  N_c_ = require_positive(material_properties, "N_c");
  alpha_c_ = require_positive(material_properties, "alpha_c");
  n_b_ = require_positive(material_properties, "n_b");
  ch_ = require_positive(material_properties, "ch");
  n_d_ = require_positive(material_properties, "n_d");
  h0_ = require_positive(material_properties, "h0");
  A0_ = require_positive(material_properties, "A0");
  Me_ = require_positive(material_properties, "Me");
  cz_ = require_positive(material_properties, "cz");
  zmax_ = require_positive(material_properties, "zmax");
  m_iso_ = require_finite(material_properties, "m_iso");
  if (m_iso_ < 0.)
    throw std::invalid_argument(
        "SANISAND parameter 'm_iso' must be nonnegative");
  if (m_iso_ == 0. &&
      material_properties.value("resume_state_policy", std::string()) ==
          "reinitialize_from_particle")
    throw std::invalid_argument(
        "SANISAND cross-model handoff requires a positive m_iso");
  patm_ = require_positive(material_properties, "Patm");
  p_min_ = require_positive(material_properties, "P_min");
  stol_ = require_positive(material_properties, "STOL");
  ftol_ = require_positive(material_properties, "FTOL");
  ltol_ = require_positive(material_properties, "LTOL");

  properties_ = material_properties;
}

template <unsigned Tdim>
mpm::dense_map mpm::Sanisand<Tdim>::initialise_state_variables() {
  const double void_ratio = porosity_ / (1. - porosity_);
  return {{"AlphaXX", 0.},
          {"AlphaYY", 0.},
          {"AlphaZZ", 0.},
          {"AlphaZY", 0.},
          {"AlphaZX", 0.},
          {"AlphaXY", 0.},
          {"AlphaInitialXX", 0.},
          {"AlphaInitialYY", 0.},
          {"AlphaInitialZZ", 0.},
          {"AlphaInitialZY", 0.},
          {"AlphaInitialZX", 0.},
          {"AlphaInitialXY", 0.},
          {"ZXX", 0.},
          {"ZYY", 0.},
          {"ZZZ", 0.},
          {"ZZY", 0.},
          {"ZZX", 0.},
          {"ZXY", 0.},
          {"eps_p_q", 0.},
          {"void_ratio", void_ratio}};
}

template <unsigned Tdim>
mpm::dense_map mpm::Sanisand<Tdim>::initialise_state_variables_from_particle(
    double porosity) {
  if (!std::isfinite(porosity) || porosity <= 0. || porosity >= 1.)
    throw std::invalid_argument(
        "SANISAND restored porosity must be between zero and one");
  auto state_vars = initialise_state_variables();
  state_vars.at("void_ratio") = porosity / (1. - porosity);
  return state_vars;
}

template <unsigned Tdim>
mpm::dense_map mpm::Sanisand<Tdim>::initialise_state_variables_from_particle(
    double porosity, const Vector6d& stress) {
  auto state_vars = this->initialise_state_variables_from_particle(porosity);
  if (!stress.allFinite())
    throw std::invalid_argument("SANISAND restored stress must be finite");

  // A SANISAND stress state lies on
  //   ||s/p - alpha|| = sqrt(2/3) m.
  // Cross-model checkpoints carry stress and porosity but no SANISAND
  // history. Centre the initial yield surface on the restored stress-ratio
  // point, retaining alpha=0 when that point is already inside the initial
  // cone. Starting every anisotropic K0 state with alpha=0 violates the yield
  // equation and creates a spurious plastic impulse during the first handoff
  // increments.
  const Vector6d internal_stress = to_internal(stress);
  const double raw_mean_stress = internal_stress.head<3>().sum() / 3.;
  if (raw_mean_stress < p_min_ ||
      !is_strictly_compressive(internal_stress))
    return state_vars;

  Vector6d deviator = internal_stress;
  deviator.head<3>().array() -= raw_mean_stress;
  const Vector6d stress_ratio = deviator / raw_mean_stress;
  const double ratio_norm = sanisand::detail::stress_norm(stress_ratio);
  const double yield_radius = std::sqrt(2. / 3.) * m_iso_;
  // A zero-radius surface is retained for reproducing legacy SANISAND
  // calibrations.  It has no unique normal at r == alpha, so an anisotropic
  // cross-model handoff cannot be centred safely in this limiting case.
  if (yield_radius == 0.) return state_vars;
  if (ratio_norm <= yield_radius) return state_vars;

  const Vector6d normal = stress_ratio / ratio_norm;
  const Vector6d alpha = stress_ratio - yield_radius * normal;
  const char* alpha_names[6] = {"AlphaXX", "AlphaYY", "AlphaZZ",
                                "AlphaZY", "AlphaZX", "AlphaXY"};
  const char* alpha_initial_names[6] = {
      "AlphaInitialXX", "AlphaInitialYY", "AlphaInitialZZ",
      "AlphaInitialZY", "AlphaInitialZX", "AlphaInitialXY"};
  for (unsigned i = 0; i < 6; ++i) {
    state_vars.at(alpha_names[i]) = alpha[i];
    state_vars.at(alpha_initial_names[i]) = alpha[i];
  }
  return state_vars;
}

template <unsigned Tdim>
typename mpm::Sanisand<Tdim>::Vector6d mpm::Sanisand<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {
  if (state_vars == nullptr || state_vars->size() != 20)
    throw std::invalid_argument("SANISAND requires twenty state variables");

  const State old_state = unpack_state(*state_vars);
  const Vector6d internal_stress = to_internal(stress);
  const Vector6d internal_dstrain = to_internal(dstrain);

  try {
    if (!dstrain.allFinite())
      throw std::runtime_error(
          "SANISAND received a non-finite strain increment");
    const auto result = integrate(internal_stress, internal_dstrain, old_state);
    if (!result.stress.allFinite() || !is_finite(result.state))
      throw std::runtime_error("SANISAND integration produced non-finite state");
    pack_state(result.state, state_vars);
    return to_mpm(result.stress);
  } catch (const std::exception& exception) {
    console_->error(
        "SANISAND integration failed for particle {}: {}; "
        "stress=[{},{},{},{},{},{}], dstrain=[{},{},{},{},{},{}], "
        "void_ratio={}, eps_p_q={}",
        ptr != nullptr ? ptr->id() : std::numeric_limits<Index>::max(),
        exception.what(), stress[0], stress[1], stress[2], stress[3], stress[4],
        stress[5], dstrain[0], dstrain[1], dstrain[2], dstrain[3], dstrain[4],
        dstrain[5], old_state.void_ratio, old_state.eps_p_q);
    throw;
  }
}

template <unsigned Tdim>
typename mpm::Sanisand<Tdim>::Vector6d mpm::Sanisand<Tdim>::to_internal(
    const Vector6d& tensor) {
  Vector6d result;
  result << -tensor[0], -tensor[1], -tensor[2], -tensor[4], -tensor[5],
      -tensor[3];
  return result;
}

template <unsigned Tdim>
typename mpm::Sanisand<Tdim>::Vector6d mpm::Sanisand<Tdim>::to_mpm(
    const Vector6d& tensor) {
  Vector6d result;
  result << -tensor[0], -tensor[1], -tensor[2], -tensor[5], -tensor[3],
      -tensor[4];
  return result;
}

template <unsigned Tdim>
typename mpm::Sanisand<Tdim>::State mpm::Sanisand<Tdim>::unpack_state(
    const mpm::dense_map& state_vars) {
  State state;
  state.alpha << state_vars.at("AlphaXX"), state_vars.at("AlphaYY"),
      state_vars.at("AlphaZZ"), state_vars.at("AlphaZY"),
      state_vars.at("AlphaZX"), state_vars.at("AlphaXY");
  state.alpha_initial << state_vars.at("AlphaInitialXX"),
      state_vars.at("AlphaInitialYY"), state_vars.at("AlphaInitialZZ"),
      state_vars.at("AlphaInitialZY"), state_vars.at("AlphaInitialZX"),
      state_vars.at("AlphaInitialXY");
  state.fabric << state_vars.at("ZXX"), state_vars.at("ZYY"),
      state_vars.at("ZZZ"), state_vars.at("ZZY"), state_vars.at("ZZX"),
      state_vars.at("ZXY");
  state.eps_p_q = state_vars.at("eps_p_q");
  state.void_ratio = state_vars.at("void_ratio");
  if (!is_finite(state) || state.void_ratio <= 0.)
    throw std::invalid_argument("SANISAND state is non-finite or has invalid void ratio");
  return state;
}

template <unsigned Tdim>
void mpm::Sanisand<Tdim>::pack_state(const State& state,
                                     mpm::dense_map* state_vars) {
  const char* alpha_names[6] = {"AlphaXX", "AlphaYY", "AlphaZZ",
                                "AlphaZY", "AlphaZX", "AlphaXY"};
  const char* alpha_initial_names[6] = {
      "AlphaInitialXX", "AlphaInitialYY", "AlphaInitialZZ",
      "AlphaInitialZY", "AlphaInitialZX", "AlphaInitialXY"};
  const char* fabric_names[6] = {"ZXX", "ZYY", "ZZZ", "ZZY", "ZZX",
                                 "ZXY"};
  for (unsigned i = 0; i < 6; ++i) {
    state_vars->at(alpha_names[i]) = state.alpha[i];
    state_vars->at(alpha_initial_names[i]) = state.alpha_initial[i];
    state_vars->at(fabric_names[i]) = state.fabric[i];
  }
  state_vars->at("eps_p_q") = state.eps_p_q;
  state_vars->at("void_ratio") = state.void_ratio;
}

template <unsigned Tdim>
bool mpm::Sanisand<Tdim>::is_finite(const State& state) {
  return state.alpha.allFinite() && state.alpha_initial.allFinite() &&
         state.fabric.allFinite() && std::isfinite(state.eps_p_q) &&
         std::isfinite(state.void_ratio);
}

template <unsigned Tdim>
typename mpm::Sanisand<Tdim>::Invariants
mpm::Sanisand<Tdim>::compute_invariants(const Vector6d& stress,
                                        const Vector6d& alpha) const {
  constexpr double tiny = 1.E-12;
  Invariants result;
  const double raw_mean = (stress[0] + stress[1] + stress[2]) / 3.;
  result.mean_stress = std::max(p_min_, raw_mean);
  result.deviator = stress;
  result.deviator[0] -= raw_mean;
  result.deviator[1] -= raw_mean;
  result.deviator[2] -= raw_mean;

  result.j2 = 0.5 *
              (result.deviator[0] * result.deviator[0] +
               result.deviator[1] * result.deviator[1] +
               result.deviator[2] * result.deviator[2] +
               2. * (result.deviator[3] * result.deviator[3] +
                     result.deviator[4] * result.deviator[4] +
                     result.deviator[5] * result.deviator[5]));

  const Eigen::Matrix3d deviator_matrix =
      sanisand::detail::stress_matrix(result.deviator);
  result.j3 = deviator_matrix.determinant();
  if (std::abs(result.j3) < tiny) result.j3 = 0.;

  if (result.j2 > tiny) {
    const double sbar = std::sqrt(result.j2);
    double sine_three_theta =
        -4.5 * result.j3 / (std::sqrt(3.) * sbar * result.j2);
    sine_three_theta = std::clamp(sine_three_theta, -1., 1.);
    const double theta = std::asin(sine_three_theta) / 3.;
    result.q = std::sqrt(3. * result.j2);
    static_cast<void>(theta);
  } else {
    result.j2 = tiny * tiny;
    result.q = tiny;
  }

  result.stress_ratio = result.deviator / result.mean_stress;
  result.normal = result.stress_ratio - alpha;
  const double normal_norm = sanisand::detail::stress_norm(result.normal);
  if (result.j2 <= tiny && normal_norm < tiny) {
    result.normal << 1. / 3., 1. / 3., 1. / 3., 0., 0., 0.;
  } else if (normal_norm > tiny) {
    result.normal /= normal_norm;
  } else {
    throw std::runtime_error("SANISAND yield normal is undefined");
  }

  result.cos_three_theta = sanisand::detail::lode_cosine(result.normal);
  return result;
}

template <unsigned Tdim>
bool mpm::Sanisand<Tdim>::is_strictly_compressive(
    const Vector6d& stress) {
  // A symmetric stress tensor is compression-positive in every principal
  // direction exactly when it is positive definite. Sylvester's criterion
  // avoids a per-particle eigenvalue solve in the explicit time loop.
  const double leading_minor_1 = stress[0];
  const double leading_minor_2 =
      stress[0] * stress[1] - stress[5] * stress[5];
  const double determinant =
      stress[0] * stress[1] * stress[2] +
      2. * stress[5] * stress[4] * stress[3] -
      stress[0] * stress[3] * stress[3] -
      stress[1] * stress[4] * stress[4] -
      stress[2] * stress[5] * stress[5];
  return leading_minor_1 > 0. && leading_minor_2 > 0. && determinant > 0.;
}

template <unsigned Tdim>
typename mpm::Sanisand<Tdim>::Matrix6x6
mpm::Sanisand<Tdim>::compute_moduli(const Vector6d& stress,
                                    const Vector6d& alpha,
                                    double specific_volume) const {
  if (!std::isfinite(specific_volume) || specific_volume <= 1.)
    throw std::runtime_error("SANISAND specific volume must be greater than one");
  const auto invariants = compute_invariants(stress, alpha);
  const double pressure_ratio = invariants.mean_stress / patm_;
  const double void_ratio = specific_volume - 1.;
  const double shear = std::pow(pressure_ratio, 0.5) * G0_ * patm_ *
                       std::pow(2.97 - void_ratio, 2.) / specific_volume;
  const double bulk = std::pow(pressure_ratio, 2. / 3.) * K0_ * patm_ *
                      specific_volume / void_ratio;
  if (!std::isfinite(shear) || !std::isfinite(bulk) || shear <= 0. ||
      bulk <= 0.)
    throw std::runtime_error("SANISAND elastic moduli are invalid");

  const double lame = bulk - 2. * shear / 3.;
  Matrix6x6 matrix = Matrix6x6::Zero();
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) matrix(i, j) = lame;
    matrix(i, i) += 2. * shear;
  }
  matrix(3, 3) = shear;
  matrix(4, 4) = shear;
  matrix(5, 5) = shear;
  return matrix;
}

template <unsigned Tdim>
typename mpm::Sanisand<Tdim>::Vector6d
mpm::Sanisand<Tdim>::yield_gradient(const Vector6d& stress,
                                    const Vector6d& alpha) const {
  const auto invariants = compute_invariants(stress, alpha);
  // df/dsigma is proportional to n - (n:r) I/3 for
  // f=||r-alpha||-sqrt(2/3)m. Using n:alpha omits the finite yield radius.
  return sanisand::detail::yield_direction(invariants.stress_ratio,
                                           invariants.normal);
}

template <unsigned Tdim>
double mpm::Sanisand<Tdim>::critical_void_ratio(double mean_stress) const {
  const double argument = std::max(p_min_, mean_stress) + alpha_c_;
  if (argument <= 0.)
    throw std::runtime_error("SANISAND critical-state logarithm is invalid");
  return std::exp(std::log(N_c_) - lambda_ * std::log(argument));
}

template <unsigned Tdim>
typename mpm::Sanisand<Tdim>::Vector6d
mpm::Sanisand<Tdim>::plastic_gradient(const Vector6d& stress,
                                      const Vector6d& alpha,
                                      const Vector6d& fabric,
                                      double specific_volume) const {
  constexpr double tiny = 1.E-14;
  if (std::abs(Mc_) < tiny)
    throw std::runtime_error("SANISAND Mc is too small");
  const auto invariants = compute_invariants(stress, alpha);
  const double c = Me_ / Mc_;
  const double c4 = std::pow(c, 4.);
  const double denominator =
      (1. + c4) - (1. - c4) * invariants.cos_three_theta;
  if (denominator <= tiny)
    throw std::runtime_error("SANISAND Lode interpolation is singular");
  const double g = std::pow(2. * c4 / denominator, 0.25);
  const double state_parameter =
      (specific_volume - 1.) - critical_void_ratio(invariants.mean_stress);
  const double alpha_d =
      g * Mc_ * std::exp(n_d_ * state_parameter) - m_iso_;
  const Vector6d alpha_d_tensor =
      std::sqrt(2. / 3.) * alpha_d * invariants.normal;
  const double fabric_projection =
      sanisand::detail::stress_double_contraction(fabric, invariants.normal);
  const double amplitude =
      A0_ * (1. + std::max(std::sqrt(3. / 2.) * fabric_projection, 0.));
  double dilatancy =
      amplitude * sanisand::detail::stress_double_contraction(
                      invariants.normal, alpha_d_tensor - alpha);
  dilatancy *= std::sqrt(3. / 2.);

  Vector6d gradient = invariants.normal;
  for (unsigned i = 0; i < 3; ++i) gradient[i] += dilatancy / 3.;
  return sanisand::detail::stress_direction_to_dual(gradient);
}

template <unsigned Tdim>
std::pair<typename mpm::Sanisand<Tdim>::Vector6d, double>
mpm::Sanisand<Tdim>::plastic_modulus(
    const Vector6d& stress, const Vector6d& alpha, const Vector6d& fabric,
    const Vector6d& alpha_initial, double specific_volume) const {
  constexpr double tiny = 1.E-15;
  const auto invariants = compute_invariants(stress, alpha);
  const double c = Me_ / Mc_;
  const double c4 = std::pow(c, 4.);
  const double denominator =
      (1. + c4) - (1. - c4) * invariants.cos_three_theta;
  if (denominator <= tiny)
    throw std::runtime_error("SANISAND bounding-surface interpolation is singular");
  const double g = std::pow(2. * c4 / denominator, 0.25);
  const double void_ratio = specific_volume - 1.;
  const double state_parameter =
      void_ratio - critical_void_ratio(invariants.mean_stress);
  const double alpha_b =
      g * Mc_ * std::max(std::exp(-n_b_ * state_parameter), 0.) - m_iso_;
  const Vector6d alpha_b_tensor =
      std::sqrt(2. / 3.) * alpha_b * invariants.normal;
  const Vector6d bounding_direction = alpha_b_tensor - alpha;
  const double reference =
      G0_ * h0_ * std::max(1. - ch_ * void_ratio, 1.E-9) *
      std::pow(invariants.mean_stress / patm_, -0.5);
  const double distance =
      std::max(tiny, sanisand::detail::stress_double_contraction(
                         alpha - alpha_initial, invariants.normal));
  const Vector6d alpha_rate = reference / distance * bounding_direction;
  const double hardening =
      2. / 3. * invariants.mean_stress *
      sanisand::detail::stress_double_contraction(alpha_rate,
                                                   invariants.normal);
  return {2. / 3. * alpha_rate, hardening};
}

template <unsigned Tdim>
typename mpm::Sanisand<Tdim>::Increment
mpm::Sanisand<Tdim>::compute_increment(
    const Vector6d& stress, const Vector6d& elastic_stress_increment,
    double substep, const State& state, double specific_volume) const {
  constexpr double tiny = 1.E-18;
  Increment increment;
  const Vector6d yield = yield_gradient(stress, state.alpha);
  const Vector6d flow =
      plastic_gradient(stress, state.alpha, state.fabric, specific_volume);
  const Matrix6x6 elastic = compute_moduli(stress, state.alpha, specific_volume);
  const auto hardening = plastic_modulus(stress, state.alpha, state.fabric,
                                         state.alpha_initial, specific_volume);
  const double coupling = yield.dot(elastic * flow);
  const double denominator = hardening.second + coupling;
  if (!std::isfinite(denominator) || std::abs(denominator) < tiny)
    throw std::runtime_error("SANISAND consistency denominator is singular");

  const double yield_norm = sanisand::detail::stress_dual_norm(yield);
  const double increment_norm =
      sanisand::detail::stress_norm(elastic_stress_increment);
  double plastic_multiplier = 0.;
  if (yield_norm > tiny && increment_norm > tiny) {
    const double cosine =
        yield.dot(elastic_stress_increment) / (yield_norm * increment_norm);
    if (cosine >= -ltol_) {
      plastic_multiplier =
          substep * yield.dot(elastic_stress_increment) / denominator;
      plastic_multiplier = std::max(plastic_multiplier, 0.);
    }
  }

  increment.stress = substep * elastic_stress_increment -
                     plastic_multiplier * (elastic * flow);
  increment.alpha = plastic_multiplier * hardening.first;

  const Vector6d plastic_strain = plastic_multiplier * flow;
  const double plastic_volumetric =
      plastic_strain[0] + plastic_strain[1] + plastic_strain[2];
  increment.eps_p_q =
      sanisand::detail::equivalent_plastic_deviatoric_strain(plastic_strain);

  const double fabric_factor = -cz_ * std::max(-plastic_volumetric, 0.);
  increment.fabric =
      fabric_factor *
      (std::sqrt(2. / 3.) * zmax_ *
           compute_invariants(stress, state.alpha).normal +
       state.fabric);

  if (!increment.stress.allFinite() || !increment.alpha.allFinite() ||
      !increment.fabric.allFinite() || !std::isfinite(increment.eps_p_q))
    throw std::runtime_error("SANISAND substep increment is non-finite");
  return increment;
}

template <unsigned Tdim>
void mpm::Sanisand<Tdim>::update_alpha_initial(Vector6d* alpha_initial,
                                               const Vector6d& stress,
                                               const Vector6d& alpha) const {
  const auto invariants = compute_invariants(stress, alpha);
  if (sanisand::detail::stress_double_contraction(
          invariants.normal, alpha - *alpha_initial) < 0.)
    *alpha_initial = alpha;
}

template <unsigned Tdim>
double mpm::Sanisand<Tdim>::relative_error(
    const Vector6d& stress1, const Vector6d& stress2,
    const Increment& increment1, const Increment& increment2,
    const Vector6d& alpha2, const Vector6d& fabric2) const {
  constexpr double tiny = 1.E-12;
  const double stress_max =
      (stress2 - stress1).cwiseAbs().maxCoeff() /
      (stress2.cwiseAbs().maxCoeff() + tiny);
  const double alpha_max =
      (increment1.alpha - increment2.alpha).cwiseAbs().maxCoeff() /
      (alpha2.cwiseAbs().maxCoeff() + tiny);
  const double fabric_max =
      (increment1.fabric - increment2.fabric).cwiseAbs().maxCoeff() /
      (fabric2.cwiseAbs().maxCoeff() + tiny);
  const double stress_l2 = sanisand::detail::stress_norm(stress2 - stress1) /
                           (sanisand::detail::stress_norm(stress2) + 1.);
  const double fabric_l2 =
      sanisand::detail::stress_norm(increment1.fabric - increment2.fabric) /
      (sanisand::detail::stress_norm(fabric2) + 1.);
  const double stress_norm = sanisand::detail::stress_norm(stress2);
  const double epsilon =
      tiny * stress_norm * stress_norm / (stress_norm * stress_norm + 1.);
  return std::max({epsilon, stress_max, alpha_max, fabric_max, stress_l2,
                   fabric_l2});
}

template <unsigned Tdim>
typename mpm::Sanisand<Tdim>::IntegrationResult
mpm::Sanisand<Tdim>::integrate(const Vector6d& stress,
                               const Vector6d& strain_increment,
                               const State& old_state) const {
  constexpr double tiny = 1.E-12;
  constexpr double minimum_substep = 0.001;
  constexpr unsigned maximum_attempts = 20000;

  IntegrationResult result;
  result.stress = stress;
  result.state = old_state;

  const double raw_mean_stress =
      (stress[0] + stress[1] + stress[2]) / 3.;
  const double initial_specific_volume = 1. + result.state.void_ratio;
  const Matrix6x6 initial_elastic = compute_moduli(
      result.stress, result.state.alpha, initial_specific_volume);
  const Vector6d elastic_trial_stress =
      result.stress + initial_elastic * strain_increment;
  const double volumetric_increment =
      strain_increment[0] + strain_increment[1] + strain_increment[2];
  const auto elastic_fallback = [&]() {
    IntegrationResult fallback;
    fallback.stress = elastic_trial_stress;
    fallback.state = old_state;
    const double updated_specific_volume =
        initial_specific_volume * (1. - volumetric_increment);
    if (!fallback.stress.allFinite() ||
        !std::isfinite(updated_specific_volume) ||
        updated_specific_volume <= 1.)
      throw std::runtime_error("SANISAND elastic fallback became invalid");
    fallback.state.void_ratio = updated_specific_volume - 1.;
    fallback.accepted_substeps = 1;
    return fallback;
  };
  if (raw_mean_stress < p_min_ ||
      !is_strictly_compressive(result.stress) ||
      !is_strictly_compressive(elastic_trial_stress))
    return elastic_fallback();

  update_alpha_initial(&result.state.alpha_initial, result.stress,
                       result.state.alpha);

  double specific_volume = 1. + result.state.void_ratio;
  double time = 0.;
  double substep = 1.;
  bool previous_rejection = false;

  for (unsigned attempt = 0; attempt < maximum_attempts; ++attempt) {
    const auto reduce_tensile_substep = [&]() {
      if (substep <= minimum_substep + tiny) return false;
      ++result.rejected_substeps;
      previous_rejection = true;
      substep = std::max(0.5 * substep, minimum_substep);
      return true;
    };
    const Matrix6x6 elastic0 =
        compute_moduli(result.stress, result.state.alpha, specific_volume);
    const Vector6d elastic_increment0 = elastic0 * strain_increment;
    const Increment first = compute_increment(
        result.stress, elastic_increment0, substep, result.state,
        specific_volume);
    const Vector6d stress1 = result.stress + first.stress;
    if (!stress1.allFinite())
      throw std::runtime_error("SANISAND predictor stress is non-finite");
    if (!is_strictly_compressive(stress1)) {
      if (reduce_tensile_substep()) continue;
      return elastic_fallback();
    }
    const Vector6d alpha1 = result.state.alpha + first.alpha;
    const Vector6d fabric1 = result.state.fabric + first.fabric;

    State predictor_state = result.state;
    predictor_state.alpha = alpha1;
    predictor_state.fabric = fabric1;
    const Matrix6x6 elastic1 =
        compute_moduli(stress1, alpha1, specific_volume);
    const Increment second = compute_increment(
        stress1, elastic1 * strain_increment, substep, predictor_state,
        specific_volume);

    const Vector6d stress2 =
        result.stress + 0.5 * (first.stress + second.stress);
    if (!stress2.allFinite())
      throw std::runtime_error("SANISAND corrector stress is non-finite");
    if (!is_strictly_compressive(stress2)) {
      if (reduce_tensile_substep()) continue;
      return elastic_fallback();
    }
    const Vector6d alpha2 =
        result.state.alpha + 0.5 * (first.alpha + second.alpha);
    const Vector6d fabric2 =
        result.state.fabric + 0.5 * (first.fabric + second.fabric);
    const double eps_p_q2 =
        result.state.eps_p_q + 0.5 * (first.eps_p_q + second.eps_p_q);
    const double error = relative_error(stress1, stress2, first, second,
                                        alpha2, fabric2);

    if (!std::isfinite(error))
      throw std::runtime_error("SANISAND substep error estimate is non-finite");
    if (error > stol_ && substep > minimum_substep + tiny) {
      ++result.rejected_substeps;
      previous_rejection = true;
      double factor = 0.9 * std::sqrt(stol_ / error);
      factor = std::max({factor, 0.1, minimum_substep / substep});
      substep *= factor;
      continue;
    }

    ++result.accepted_substeps;
    time += substep;
    result.stress = stress2;
    result.state.alpha = alpha2;
    result.state.fabric = fabric2;
    result.state.eps_p_q = eps_p_q2;
    specific_volume *= 1. - substep * volumetric_increment;
    if (!std::isfinite(specific_volume) || specific_volume <= 1.)
      throw std::runtime_error("SANISAND void ratio became invalid");
    result.state.void_ratio = specific_volume - 1.;

    if (std::abs(1. - time) <= tiny) return result;
    if (time > 1. + tiny)
      throw std::runtime_error("SANISAND substepping exceeded unit time");

    double factor =
        error > tiny ? 0.9 * std::sqrt(stol_ / error) : 1.1;
    factor = previous_rejection ? std::min(factor, 1.)
                                : std::min(factor, 1.1);
    previous_rejection = false;
    substep = std::max(substep * factor, minimum_substep);
    substep = std::min(substep, 1. - time);
  }
  throw std::runtime_error("SANISAND exceeded maximum substep attempts");
}
