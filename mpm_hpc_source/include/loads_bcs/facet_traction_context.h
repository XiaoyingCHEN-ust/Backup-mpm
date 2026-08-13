#ifndef MPM_FACET_TRACTION_CONTEXT_H_
#define MPM_FACET_TRACTION_CONTEXT_H_

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include "Eigen/Dense"

#include "element.h"

namespace mpm {
namespace facet_traction {

//! Immutable geometry used to map one particle surface traction.
//!
//! The context is intentionally not part of HDF5Particle. A new stage may
//! reconstruct it from the checkpoint geometry before its first load step,
//! which is sufficient for the pipeline handoff workflow. An exact mid-stage
//! restart cannot reproduce the original reference context and is therefore
//! unsupported until the checkpoint schema stores this data explicitly.
struct FacetTractionContext {
  unsigned public_facet{std::numeric_limits<unsigned>::max()};
  unsigned element_face{std::numeric_limits<unsigned>::max()};
  unsigned normal_direction{std::numeric_limits<unsigned>::max()};
  Eigen::VectorXd reference_shape_weights;
  Eigen::VectorXd reference_position;
  double reference_surface_measure{
      std::numeric_limits<double>::quiet_NaN()};

  bool initialised() const noexcept {
    return public_facet != std::numeric_limits<unsigned>::max();
  }

  bool valid(unsigned dimension, Eigen::Index node_count) const {
    return initialised() && public_facet < 2 * dimension &&
           element_face < 2 * dimension && normal_direction < dimension &&
           reference_shape_weights.size() == node_count &&
           reference_shape_weights.allFinite() &&
           std::abs(reference_shape_weights.sum() - 1.) <= 1.e-10 &&
           reference_position.size() == static_cast<Eigen::Index>(dimension) &&
           reference_position.allFinite() &&
           std::isfinite(reference_surface_measure) &&
           reference_surface_measure > 0.;
  }

  void reset() {
    public_facet = std::numeric_limits<unsigned>::max();
    element_face = std::numeric_limits<unsigned>::max();
    normal_direction = std::numeric_limits<unsigned>::max();
    reference_shape_weights.resize(0);
    reference_position.resize(0);
    reference_surface_measure = std::numeric_limits<double>::quiet_NaN();
  }
};

//! Convert the public Cartesian convention
//! 0..5 = +x,-x,+y,-y,+z,-z to the element's tensor-product face id.
inline unsigned public_facet_to_element_face(unsigned dimension,
                                              unsigned public_facet,
                                              unsigned element_nfaces) {
  if (dimension != 2 && dimension != 3)
    throw std::invalid_argument(
        "Facet traction supports only two- and three-dimensional meshes");

  const unsigned expected_nfaces = 2 * dimension;
  if (element_nfaces != expected_nfaces)
    throw std::invalid_argument(
        "Facet traction requires a tensor-product element with " +
        std::to_string(expected_nfaces) + " faces");
  if (public_facet >= expected_nfaces)
    throw std::out_of_range("Public Cartesian traction facet is invalid");

  static constexpr std::array<unsigned, 4> map_2d{{1, 3, 2, 0}};
  static constexpr std::array<unsigned, 6> map_3d{{1, 3, 2, 0, 5, 4}};
  return dimension == 2 ? map_2d.at(public_facet)
                        : map_3d.at(public_facet);
}

namespace detail {

//! Geometry needed to interpolate a public Cartesian facet in the particle's
//! current element. Unlike FacetTractionContext, this data is intentionally
//! rebuilt after particle relocation and carries no physical surface measure.
template <unsigned Tdim>
struct FacetInterpolationData {
  unsigned element_face;
  unsigned normal_direction;
  Eigen::MatrixXd unit_cell;
  Eigen::VectorXi face_indices;
  Eigen::Matrix<double, Tdim, 1> unit_min;
  Eigen::Matrix<double, Tdim, 1> unit_max;
  Eigen::VectorXd shape_weights;
};

template <unsigned Tdim>
FacetInterpolationData<Tdim> make_interpolation_data(
    const std::shared_ptr<const Element<Tdim>>& element,
    const Eigen::Matrix<double, Tdim, 1>& particle_xi,
    unsigned public_facet) {
  static_assert(Tdim == 2 || Tdim == 3,
                "Facet traction supports only 2D and 3D");
  if (element == nullptr)
    throw std::invalid_argument("Facet traction element is null");

  const unsigned expected_nfunctions = 1u << Tdim;
  if (element->shapefn_type() != mpm::ShapefnType::NORMAL_MPM ||
      element->nfunctions() != expected_nfunctions)
    throw std::invalid_argument(
        "Facet traction supports only normal-MPM Q4 or H8 elements");

  FacetInterpolationData<Tdim> data;
  data.element_face = public_facet_to_element_face(
      Tdim, public_facet, element->nfaces());
  data.normal_direction = public_facet / 2;
  data.unit_cell = element->unit_cell_coordinates();
  if (data.unit_cell.rows() !=
          static_cast<Eigen::Index>(expected_nfunctions) ||
      data.unit_cell.cols() != static_cast<Eigen::Index>(Tdim) ||
      !data.unit_cell.allFinite() || !particle_xi.allFinite())
    throw std::invalid_argument(
        "Facet traction interpolation geometry is invalid");

  constexpr double tolerance = 1.e-12;
  for (unsigned dimension = 0; dimension < Tdim; ++dimension) {
    data.unit_min[dimension] = data.unit_cell.col(dimension).minCoeff();
    data.unit_max[dimension] = data.unit_cell.col(dimension).maxCoeff();
    const double unit_span =
        data.unit_max[dimension] - data.unit_min[dimension];
    if (!std::isfinite(unit_span) || unit_span <= tolerance)
      throw std::invalid_argument(
          "Facet traction unit-cell span is invalid");

    // Corner-only Cartesian tensor product: every coordinate lies on one of
    // the two bounds. This deliberately rejects triangles, higher-order
    // elements and GIMP/B-spline variants in the formal ED2Q4 pipeline path.
    for (Eigen::Index node = 0; node < data.unit_cell.rows(); ++node) {
      const double coordinate = data.unit_cell(node, dimension);
      if (std::abs(coordinate - data.unit_min[dimension]) > tolerance &&
          std::abs(coordinate - data.unit_max[dimension]) > tolerance)
        throw std::invalid_argument(
            "Facet traction element is not a corner tensor product");
    }
    if (particle_xi[dimension] < data.unit_min[dimension] - tolerance ||
        particle_xi[dimension] > data.unit_max[dimension] + tolerance)
      throw std::invalid_argument(
          "Facet traction particle coordinate is outside the unit cell");
  }

  data.face_indices = element->face_indices(data.element_face);
  const Eigen::Index expected_face_nodes = 1u << (Tdim - 1);
  if (data.face_indices.size() != expected_face_nodes)
    throw std::invalid_argument(
        "Facet traction element face has an unsupported node count");

  const double expected_normal_coordinate =
      public_facet % 2 == 0 ? data.unit_max[data.normal_direction]
                            : data.unit_min[data.normal_direction];
  std::array<bool, 1u << Tdim> node_on_face{};
  for (Eigen::Index i = 0; i < data.face_indices.size(); ++i) {
    const int node = data.face_indices[i];
    if (node < 0 || node >= data.unit_cell.rows() || node_on_face.at(node))
      throw std::invalid_argument(
          "Facet traction element face indices are invalid");
    node_on_face.at(node) = true;
    if (std::abs(data.unit_cell(node, data.normal_direction) -
                 expected_normal_coordinate) > tolerance)
      throw std::invalid_argument(
          "Element face does not match the public Cartesian facet");
  }

  Eigen::Matrix<double, Tdim, 1> facet_xi = particle_xi;
  facet_xi[data.normal_direction] = expected_normal_coordinate;
  const auto zero = Eigen::Matrix<double, Tdim, 1>::Zero();
  data.shape_weights = element->shapefn_local(facet_xi, zero, zero);
  if (data.shape_weights.size() != data.unit_cell.rows() ||
      !data.shape_weights.allFinite() ||
      std::abs(data.shape_weights.sum() - 1.) > 1.e-10)
    throw std::invalid_argument(
        "Facet traction shape weights are invalid");
  for (Eigen::Index node = 0; node < data.shape_weights.size(); ++node) {
    const double weight = data.shape_weights[node];
    if (weight < -1.e-10 ||
        (!node_on_face.at(node) && std::abs(weight) > 1.e-10))
      throw std::invalid_argument(
          "Facet traction weights leak to off-face nodes");
  }
  return data;
}

}  // namespace detail

//! Recompute only the interpolation weights on the particle's current cell.
//! The caller continues to use the stage-start reference surface measure.
template <unsigned Tdim>
Eigen::VectorXd current_shape_weights(
    const std::shared_ptr<const Element<Tdim>>& element,
    const Eigen::Matrix<double, Tdim, 1>& particle_xi,
    unsigned public_facet) {
  return detail::make_interpolation_data<Tdim>(element, particle_xi,
                                                public_facet)
      .shape_weights;
}

//! Build a frozen traction context for a normal-MPM Q4 or H8 element.
template <unsigned Tdim>
FacetTractionContext make_context(
    const std::shared_ptr<const Element<Tdim>>& element,
    const Eigen::MatrixXd& nodal_coordinates,
    const Eigen::Matrix<double, Tdim, 1>& particle_xi,
    const Eigen::Matrix<double, Tdim, 1>& natural_size,
    unsigned public_facet) {
  const auto interpolation = detail::make_interpolation_data<Tdim>(
      element, particle_xi, public_facet);

  FacetTractionContext context;
  context.public_facet = public_facet;
  context.element_face = interpolation.element_face;
  context.normal_direction = interpolation.normal_direction;
  context.reference_shape_weights = interpolation.shape_weights;

  if (nodal_coordinates.rows() != interpolation.unit_cell.rows() ||
      nodal_coordinates.cols() != interpolation.unit_cell.cols() ||
      !nodal_coordinates.allFinite() || !natural_size.allFinite())
    throw std::invalid_argument(
        "Facet traction reference geometry is invalid");
  constexpr double tolerance = 1.e-12;
  const auto& face_indices = interpolation.face_indices;

  double full_face_measure = 0.;
  if constexpr (Tdim == 2) {
    const Eigen::Vector2d first =
        nodal_coordinates.row(face_indices[0]).transpose();
    const Eigen::Vector2d second =
        nodal_coordinates.row(face_indices[1]).transpose();
    full_face_measure = (second - first).norm();
  } else {
    const Eigen::Vector3d p0 =
        nodal_coordinates.row(face_indices[0]).transpose();
    const Eigen::Vector3d p1 =
        nodal_coordinates.row(face_indices[1]).transpose();
    const Eigen::Vector3d p2 =
        nodal_coordinates.row(face_indices[2]).transpose();
    const Eigen::Vector3d p3 =
        nodal_coordinates.row(face_indices[3]).transpose();
    // Half the sum of all four three-vertex triangle areas equals the area of
    // a planar convex quadrilateral. Unlike a fixed two-triangle split, this
    // expression is independent of the face-index starting point and order;
    // it also gives a symmetric measure for a mildly warped H8 face.
    const auto triangle_area = [](const Eigen::Vector3d& a,
                                  const Eigen::Vector3d& b,
                                  const Eigen::Vector3d& c) {
      return 0.5 * (b - a).cross(c - a).norm();
    };
    full_face_measure =
        0.5 * (triangle_area(p0, p1, p2) +
               triangle_area(p0, p1, p3) +
               triangle_area(p0, p2, p3) +
               triangle_area(p1, p2, p3));
  }
  if (!std::isfinite(full_face_measure) || full_face_measure <= 0.)
    throw std::invalid_argument(
        "Facet traction physical face measure is invalid");

  context.reference_surface_measure = full_face_measure;
  for (unsigned dimension = 0; dimension < Tdim; ++dimension) {
    if (dimension == context.normal_direction) continue;
    const double unit_span = interpolation.unit_max[dimension] -
                             interpolation.unit_min[dimension];
    if (natural_size[dimension] <= 0. ||
        natural_size[dimension] > unit_span + tolerance)
      throw std::invalid_argument(
          "Facet traction natural tangential size is invalid");
    context.reference_surface_measure *=
        natural_size[dimension] / unit_span;
  }
  if (!std::isfinite(context.reference_surface_measure) ||
      context.reference_surface_measure <= 0.)
    throw std::invalid_argument(
        "Facet traction reference surface measure is invalid");

  context.reference_position =
      nodal_coordinates.transpose() * context.reference_shape_weights;
  if (!context.valid(Tdim, nodal_coordinates.rows()))
    throw std::invalid_argument("Facet traction context is invalid");
  return context;
}

}  // namespace facet_traction
}  // namespace mpm

#endif  // MPM_FACET_TRACTION_CONTEXT_H_
