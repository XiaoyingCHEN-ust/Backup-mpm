#ifndef MPM_THREEPHASE_PRESSURE_STATE_H_
#define MPM_THREEPHASE_PRESSURE_STATE_H_

#include <vector>

#include <Eigen/Core>

#include "data_types.h"

namespace mpm {

//! Frozen particle data used to assemble the semi-implicit pressure system.
template <unsigned Tdim>
struct ThreePhasePressureState {
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  std::vector<Index> active_node_ids;
  Eigen::VectorXd shapefn;
  Eigen::MatrixXd dn_dx;
  double volume{0.};
  double porosity{0.};
  double liquid_saturation{0.};
  double gas_saturation{0.};
  double liquid_pressure{0.};
  double gas_pressure{0.};
  double saturation_pressure_tangent{0.};
  double liquid_compressibility{0.};
  double gas_compressibility{0.};
  double liquid_conductivity{0.};
  double gas_conductivity{0.};
  double liquid_density{0.};
  double gas_density{0.};
  double liquid_mass_source{0.};
  double gas_mass_source{0.};
  bool fixed_gas_pressure{false};
  bool pressure_boundary{false};
  double boundary_liquid_pressure{0.};
  double boundary_gas_pressure{0.};
};

//! Narrow interface implemented by particles that support the coupled solve.
template <unsigned Tdim>
class ThreePhasePressureParticle {
 public:
  virtual ~ThreePhasePressureParticle() = default;
  virtual ThreePhasePressureState<Tdim> semi_implicit_pressure_state() = 0;
  virtual int update_semi_implicit_pressure(
      const Eigen::VectorXd& liquid_pressure_increment,
      const Eigen::VectorXd& gas_pressure_increment,
      const Eigen::VectorXd& nodal_liquid_pressure,
      const Eigen::VectorXd& nodal_gas_pressure,
      bool reconstruct_pressure_gradient, bool bounded_transfer, double dt) = 0;
};

}  // namespace mpm

#endif  // MPM_THREEPHASE_PRESSURE_STATE_H_
