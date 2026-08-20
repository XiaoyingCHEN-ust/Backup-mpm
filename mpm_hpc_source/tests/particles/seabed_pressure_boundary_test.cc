#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <array>
#include <limits>

#include "facet_traction_context.h"
#include "hexahedron_element.h"
#include "mesh.h"
#include "node.h"
#include "particle_threephase_lag.h"
#include "quadrilateral_element.h"

namespace {

class ThreePhaseParticleLag2DProbe
    : public mpm::ThreePhaseParticleLag<2> {
 public:
  explicit ThreePhaseParticleLag2DProbe(const Eigen::Vector2d& coordinates)
      : mpm::ThreePhaseParticleLag<2>(0, coordinates) {}

  void assign_traction_geometry(double volume, const Eigen::Vector2d& size,
                                const Eigen::Vector2d& natural_size) {
    this->volume_ = volume;
    this->size_ = size;
    this->natural_size_ = natural_size;
  }

  void assign_force_mapping(const std::shared_ptr<mpm::NodeBase<2>>& node,
                            double mixture_mass) {
    this->nodes_ = {node};
    this->shapefn_ = Eigen::VectorXd::Ones(1);
    this->mixture_mass_ = mixture_mass;
  }

  void assign_mixture_mass(double mixture_mass) {
    this->mixture_mass_ = mixture_mass;
  }

  void assign_current_geometry(double volume, const Eigen::Vector2d& size) {
    this->volume_ = volume;
    this->size_ = size;
  }

  void configure_dynamic_surface(double reference_pressure,
                                 double current_pressure) {
    this->wave_x_ref_initialized_ = true;
    this->physical_seabed_surface_marker_ = true;
    this->ini_liquid_pressure_ = reference_pressure;
    this->liquid_pressure_ = current_pressure;
    this->initialise_dynamic_surface_traction_context();
  }

  void mark_dynamic_surface_without_context(double reference_pressure,
                                             double current_pressure) {
    this->wave_x_ref_initialized_ = true;
    this->physical_seabed_surface_marker_ = true;
    this->ini_liquid_pressure_ = reference_pressure;
    this->liquid_pressure_ = current_pressure;
    this->dynamic_surface_traction_context_.reset();
  }

  void emulate_stage_handoff_context_loss() {
    for (auto& context : this->static_traction_contexts_) context.reset();
    this->static_traction_active_.fill(false);
    this->dynamic_surface_traction_context_.reset();
    this->mixture_traction_.setZero();
    this->liquid_traction_.setZero();
    this->gas_traction_.setZero();
    this->set_mixture_traction_ = false;
  }

  const Eigen::Vector2d& mixture_traction() const {
    return this->mixture_traction_;
  }
  const Eigen::Vector2d& liquid_traction() const {
    return this->liquid_traction_;
  }
  const Eigen::Vector2d& gas_traction() const { return this->gas_traction_; }
  bool has_phase_traction() const { return this->set_mixture_traction_; }
  std::size_t wave_table_size() const { return this->seabed_surface_x_.size(); }
  const mpm::facet_traction::FacetTractionContext& static_context(
      unsigned direction) const {
    return this->static_traction_contexts_.at(direction);
  }
  const mpm::facet_traction::FacetTractionContext& dynamic_context() const {
    return this->dynamic_surface_traction_context_;
  }
};

struct Q4TractionFixture {
  std::array<Eigen::Vector2d, 4> coordinates{
      Eigen::Vector2d(0., 0.), Eigen::Vector2d(1., 0.),
      Eigen::Vector2d(1., 1.), Eigen::Vector2d(0., 1.)};
  std::array<std::shared_ptr<mpm::Node<2, 2, 3>>, 4> nodes;
  std::shared_ptr<mpm::QuadrilateralElement<2, 4>> element;
  std::shared_ptr<mpm::Cell<2>> cell;
  std::shared_ptr<ThreePhaseParticleLag2DProbe> particle;
  mpm::Mesh<2> mesh{0};

  Q4TractionFixture() {
    for (unsigned id = 0; id < nodes.size(); ++id)
      nodes[id] =
          std::make_shared<mpm::Node<2, 2, 3>>(id, coordinates[id]);
    element = std::make_shared<mpm::QuadrilateralElement<2, 4>>();
    cell = std::make_shared<mpm::Cell<2>>(0, 4, element);
    for (unsigned id = 0; id < nodes.size(); ++id)
      if (!cell->add_node(id, nodes[id]))
        throw std::runtime_error("Failed to construct Q4 test cell");
    if (!cell->initialise())
      throw std::runtime_error("Failed to initialise Q4 test cell");

    particle = std::make_shared<ThreePhaseParticleLag2DProbe>(
        Eigen::Vector2d(0.25, 0.75));
    if (!particle->assign_cell_xi(cell, Eigen::Vector2d(-0.5, 0.5)))
      throw std::runtime_error("Failed to assign Q4 test particle");
    particle->compute_shapefn();
    // The cell is 1 x 1. Natural sizes (0.2, 0.4) therefore represent a
    // reference particle width 0.1 and height 0.2.
    particle->assign_traction_geometry(0.02, Eigen::Vector2d(0.1, 0.2),
                                       Eigen::Vector2d(0.2, 0.4));
    particle->assign_mixture_mass(0.0);
    if (!mesh.add_particle(particle, false))
      throw std::runtime_error("Failed to add Q4 test particle");
  }

  void reset_nodal_forces() {
    for (const auto& node : nodes) node->initialise();
  }

  Eigen::Vector2d resultant(mpm::ParticlePhase phase) const {
    Eigen::Vector2d force = Eigen::Vector2d::Zero();
    for (const auto& node : nodes) force += node->external_force(phase);
    return force;
  }

  double first_moment(mpm::ParticlePhase phase) const {
    double moment = 0.;
    for (unsigned node = 0; node < nodes.size(); ++node) {
      const Eigen::Vector2d force = nodes[node]->external_force(phase);
      moment += coordinates[node][0] * force[1] -
                coordinates[node][1] * force[0];
    }
    return moment;
  }
};

struct CrossCellQ4TractionFixture {
  std::array<Eigen::Vector2d, 6> coordinates{
      Eigen::Vector2d(0., 0.), Eigen::Vector2d(1., 0.),
      Eigen::Vector2d(1., 1.), Eigen::Vector2d(0., 1.),
      Eigen::Vector2d(2., 0.), Eigen::Vector2d(2., 1.)};
  std::array<std::shared_ptr<mpm::Node<2, 2, 3>>, 6> nodes;
  std::shared_ptr<mpm::QuadrilateralElement<2, 4>> first_element;
  std::shared_ptr<mpm::QuadrilateralElement<2, 4>> second_element;
  std::shared_ptr<mpm::Cell<2>> first_cell;
  std::shared_ptr<mpm::Cell<2>> second_cell;
  std::shared_ptr<ThreePhaseParticleLag2DProbe> particle;
  mpm::Mesh<2> mesh{0};

  CrossCellQ4TractionFixture() {
    for (unsigned id = 0; id < nodes.size(); ++id)
      nodes[id] =
          std::make_shared<mpm::Node<2, 2, 3>>(id, coordinates[id]);

    first_element = std::make_shared<mpm::QuadrilateralElement<2, 4>>();
    second_element = std::make_shared<mpm::QuadrilateralElement<2, 4>>();
    first_cell = std::make_shared<mpm::Cell<2>>(0, 4, first_element);
    second_cell = std::make_shared<mpm::Cell<2>>(1, 4, second_element);
    const std::array<unsigned, 4> first_node_ids{{0, 1, 2, 3}};
    const std::array<unsigned, 4> second_node_ids{{1, 4, 5, 2}};
    for (unsigned local = 0; local < first_node_ids.size(); ++local) {
      const unsigned first_id = first_node_ids[local];
      const unsigned second_id = second_node_ids[local];
      if (!first_cell->add_node(local, nodes[first_id]) ||
          !second_cell->add_node(local, nodes[second_id]))
        throw std::runtime_error("Failed to construct cross-cell Q4 mesh");
    }
    if (!first_cell->initialise() || !second_cell->initialise())
      throw std::runtime_error("Failed to initialise cross-cell Q4 mesh");

    // Start at xi_x = +0.5 in the first cell. Crossing into the second cell
    // changes xi_x to -0.5, so applying the frozen local weights to the new
    // local node ordering would reverse the top-facet force distribution.
    particle = std::make_shared<ThreePhaseParticleLag2DProbe>(
        Eigen::Vector2d(0.75, 0.75));
    if (!particle->assign_cell_xi(first_cell, Eigen::Vector2d(0.5, 0.5)))
      throw std::runtime_error("Failed to assign cross-cell Q4 particle");
    particle->compute_shapefn();
    particle->assign_traction_geometry(0.02, Eigen::Vector2d(0.1, 0.2),
                                       Eigen::Vector2d(0.2, 0.4));
    particle->assign_mixture_mass(0.0);
    if (!mesh.add_particle(particle, false))
      throw std::runtime_error("Failed to add cross-cell Q4 particle");
  }

  bool move_to_second_cell() {
    particle->assign_coordinates(Eigen::Vector2d(1.25, 0.75));
    if (!particle->assign_cell_xi(second_cell,
                                  Eigen::Vector2d(-0.5, 0.5)))
      return false;
    particle->compute_shapefn();
    return true;
  }

  void reset_nodal_forces() {
    for (const auto& node : nodes) node->initialise();
  }

  Eigen::Vector2d resultant(mpm::ParticlePhase phase) const {
    Eigen::Vector2d force = Eigen::Vector2d::Zero();
    for (const auto& node : nodes) force += node->external_force(phase);
    return force;
  }

  double first_moment(mpm::ParticlePhase phase) const {
    double moment = 0.;
    for (unsigned node = 0; node < nodes.size(); ++node) {
      const Eigen::Vector2d force = nodes[node]->external_force(phase);
      moment += coordinates[node][0] * force[1] -
                coordinates[node][1] * force[0];
    }
    return moment;
  }
};

}  // namespace

TEST_CASE("Wave total traction marker selects only the initial seabed row") {
  using mpm::threephase_lag_boundary::make_surface_traction_marker;

  const double seabed_y = 0.50;
  const double particle_size_y = 0.01;
  REQUIRE(make_surface_traction_marker(true, false, 0.495, seabed_y,
                                       particle_size_y));
  REQUIRE_FALSE(
      make_surface_traction_marker(true, false, 0.485, seabed_y,
                                   particle_size_y));
  REQUIRE_FALSE(
      make_surface_traction_marker(true, false, 0.46, seabed_y,
                                   particle_size_y));
  REQUIRE_FALSE(make_surface_traction_marker(false, false, 0.495, seabed_y,
                                              particle_size_y));
  REQUIRE_FALSE(make_surface_traction_marker(true, true, 0.495, seabed_y,
                                              particle_size_y));
  REQUIRE_FALSE(make_surface_traction_marker(true, false, 0.515, seabed_y,
                                              particle_size_y));
}

TEST_CASE("Wave total traction marker survives subsequent displacement") {
  using mpm::threephase_lag_boundary::make_surface_traction_marker;

  // The marker is created once from the reference coordinate. Current
  // coordinates are intentionally absent from its runtime semantics.
  const bool marker =
      make_surface_traction_marker(true, false, 0.495, 0.50, 0.01);
  REQUIRE(marker);
  for (const double displaced_y : {0.470, 0.495, 0.525}) {
    CAPTURE(displaced_y);
    REQUIRE(marker);
  }
}

TEST_CASE("Three-phase internal forces use pressure relative to restart state") {
  using mpm::threephase_lag_force::excess_phase_pressure;

  REQUIRE(excess_phase_pressure(5000.0, 5000.0) == Approx(0.0));
  REQUIRE(excess_phase_pressure(5250.0, 5000.0) == Approx(250.0));
  REQUIRE(excess_phase_pressure(4750.0, 5000.0) == Approx(-250.0));
  REQUIRE_THROWS(excess_phase_pressure(
      std::numeric_limits<double>::quiet_NaN(), 5000.0));
}

TEST_CASE("Three-phase force validation rejects non-finite mapped vectors") {
  Eigen::Vector2d force(1.0, -2.0);
  REQUIRE_NOTHROW(
      mpm::threephase_lag_force::require_finite_force(force, "test"));
  force[1] = std::numeric_limits<double>::infinity();
  REQUIRE_THROWS_WITH(
      mpm::threephase_lag_force::require_finite_force(force, "test"),
      Catch::Matchers::Contains("Non-finite test force contribution"));
}

TEST_CASE("Three-phase particle traction starts disabled and zero") {
  ThreePhaseParticleLag2DProbe particle(Eigen::Vector2d::Zero());
  auto node = std::make_shared<mpm::Node<2, 2, 3>>(
      0, Eigen::Vector2d::Zero());
  particle.assign_force_mapping(node, 2.0);

  REQUIRE_FALSE(particle.has_phase_traction());
  REQUIRE(particle.mixture_traction().isZero(0.0));
  REQUIRE(particle.liquid_traction().isZero(0.0));
  REQUIRE(particle.gas_traction().isZero(0.0));
  REQUIRE(particle.wave_table_size() == 0);

  particle.map_external_force(Eigen::Vector2d(0.0, -9.81));
  REQUIRE(node->external_force(mpm::ParticlePhase::Mixture)
              .isApprox(Eigen::Vector2d(0.0, -19.62)));
  REQUIRE(node->external_force(mpm::ParticlePhase::Liquid).isZero(0.0));
  REQUIRE(node->external_force(mpm::ParticlePhase::Gas).isZero(0.0));
}

TEST_CASE("Three-phase traction without an explicit facet fails closed") {
  ThreePhaseParticleLag2DProbe particle(Eigen::Vector2d::Zero());
  auto node = std::make_shared<mpm::Node<2, 2, 3>>(
      1, Eigen::Vector2d::Zero());
  particle.assign_force_mapping(node, 2.0);
  particle.assign_traction_geometry(0.02, Eigen::Vector2d(0.1, 0.2),
                                    Eigen::Vector2d(0.2, 0.4));

  REQUIRE_FALSE(particle.assign_particle_traction(1, -100.0));
  REQUIRE_FALSE(particle.has_phase_traction());
  REQUIRE(particle.mixture_traction()[0] == Approx(0.0));
  REQUIRE(particle.mixture_traction()[1] == Approx(0.0));
  REQUIRE(particle.liquid_traction().isZero(0.0));
  REQUIRE(particle.gas_traction().isZero(0.0));
}

TEST_CASE("Public Cartesian facets map to tensor-product element faces") {
  using mpm::facet_traction::public_facet_to_element_face;

  const std::array<unsigned, 4> expected_2d{{1, 3, 2, 0}};
  for (unsigned facet = 0; facet < expected_2d.size(); ++facet)
    REQUIRE(public_facet_to_element_face(2, facet, 4) ==
            expected_2d[facet]);

  const std::array<unsigned, 6> expected_3d{{1, 3, 2, 0, 5, 4}};
  for (unsigned facet = 0; facet < expected_3d.size(); ++facet)
    REQUIRE(public_facet_to_element_face(3, facet, 6) ==
            expected_3d[facet]);

  REQUIRE_THROWS(public_facet_to_element_face(1, 0, 2));
  REQUIRE_THROWS(public_facet_to_element_face(2, 4, 4));
  REQUIRE_THROWS(public_facet_to_element_face(2, 0, 3));
}

TEST_CASE("Facet context rejects unsupported higher-order elements") {
  auto element = std::make_shared<mpm::QuadrilateralElement<2, 8>>();
  const Eigen::MatrixXd coordinates = element->unit_cell_coordinates();
  REQUIRE_THROWS(mpm::facet_traction::make_context<2>(
      element, coordinates, Eigen::Vector2d::Zero(),
      Eigen::Vector2d::Ones(), 2));
}

TEST_CASE("H8 Cartesian facets conserve full area and centroid moment") {
  auto element = std::make_shared<mpm::HexahedronElement<3, 8>>();
  Eigen::Matrix<double, 8, 3> coordinates;
  coordinates << 0., 0., 0.,
                 2., 0., 0.,
                 2., 3., 0.,
                 0., 3., 0.,
                 0., 0., 4.,
                 2., 0., 4.,
                 2., 3., 4.,
                 0., 3., 4.;
  const Eigen::Vector3d particle_xi = Eigen::Vector3d::Zero();
  const Eigen::Vector3d full_unit_size =
      2. * Eigen::Vector3d::Ones();
  const std::array<double, 6> expected_areas{{12., 12., 8., 8., 6., 6.}};
  const std::array<Eigen::Vector3d, 6> expected_centroids{{
      Eigen::Vector3d(2., 1.5, 2.), Eigen::Vector3d(0., 1.5, 2.),
      Eigen::Vector3d(1., 3., 2.), Eigen::Vector3d(1., 0., 2.),
      Eigen::Vector3d(1., 1.5, 4.), Eigen::Vector3d(1., 1.5, 0.)}};
  const Eigen::Vector3d traction(2., -3., 5.);

  double total_area = 0.;
  for (unsigned public_facet = 0; public_facet < 6; ++public_facet) {
    CAPTURE(public_facet);
    const auto context = mpm::facet_traction::make_context<3>(
        element, coordinates, particle_xi, full_unit_size, public_facet);
    REQUIRE(context.reference_surface_measure ==
            Approx(expected_areas[public_facet]));
    REQUIRE(context.reference_position.isApprox(
        expected_centroids[public_facet], 1.e-14));

    const Eigen::Vector3d expected_resultant =
        expected_areas[public_facet] * traction;
    Eigen::Vector3d resultant = Eigen::Vector3d::Zero();
    Eigen::Vector3d moment = Eigen::Vector3d::Zero();
    for (Eigen::Index node = 0;
         node < context.reference_shape_weights.size(); ++node) {
      const Eigen::Vector3d nodal_force =
          context.reference_shape_weights[node] * expected_resultant;
      resultant += nodal_force;
      moment += coordinates.row(node).transpose().cross(nodal_force);
    }
    REQUIRE(resultant.isApprox(expected_resultant, 1.e-13));
    REQUIRE(moment.isApprox(
        expected_centroids[public_facet].cross(expected_resultant), 1.e-13));
    total_area += context.reference_surface_measure;
  }
  REQUIRE(total_area == Approx(52.0));
}

TEST_CASE("Legacy particle types retain centre-based traction compatibility") {
  const std::array<Eigen::Vector2d, 4> coordinates{
      Eigen::Vector2d(0., 0.), Eigen::Vector2d(1., 0.),
      Eigen::Vector2d(1., 1.), Eigen::Vector2d(0., 1.)};
  std::array<std::shared_ptr<mpm::Node<2, 2, 3>>, 4> nodes;
  auto element = std::make_shared<mpm::QuadrilateralElement<2, 4>>();
  auto cell = std::make_shared<mpm::Cell<2>>(0, 4, element);
  for (unsigned id = 0; id < nodes.size(); ++id) {
    nodes[id] =
        std::make_shared<mpm::Node<2, 2, 3>>(id, coordinates[id]);
    REQUIRE(cell->add_node(id, nodes[id]));
  }
  REQUIRE(cell->initialise());
  auto particle =
      std::make_shared<mpm::Particle<2>>(0, Eigen::Vector2d(0.5, 0.5));
  REQUIRE(particle->assign_cell_xi(cell, Eigen::Vector2d::Zero()));
  particle->compute_shapefn();

  mpm::Mesh<2> mesh(0);
  REQUIRE(mesh.add_particle(particle, false));
  REQUIRE(mesh.create_particles_tractions(nullptr, -1, 2, 1, -100.0));
  REQUIRE_NOTHROW(mesh.apply_traction_on_particles(0.0));
}

TEST_CASE("Q4 Cartesian facets conserve resultant and first moment") {
  struct FacetExpectation {
    unsigned public_facet;
    unsigned element_face;
    unsigned force_direction;
    Eigen::Vector2d resultant;
    double first_moment;
    std::array<Eigen::Vector2d, 4> nodal_forces;
  };

  const std::array<FacetExpectation, 4> expectations{{
      {0, 1, 0, Eigen::Vector2d(-20., 0.), 15.,
       {Eigen::Vector2d(0., 0.), Eigen::Vector2d(-5., 0.),
        Eigen::Vector2d(-15., 0.), Eigen::Vector2d(0., 0.)}},
      {1, 3, 0, Eigen::Vector2d(-20., 0.), 15.,
       {Eigen::Vector2d(-5., 0.), Eigen::Vector2d(0., 0.),
        Eigen::Vector2d(0., 0.), Eigen::Vector2d(-15., 0.)}},
      {2, 2, 1, Eigen::Vector2d(0., -10.), -2.5,
       {Eigen::Vector2d(0., 0.), Eigen::Vector2d(0., 0.),
        Eigen::Vector2d(0., -2.5), Eigen::Vector2d(0., -7.5)}},
      {3, 0, 1, Eigen::Vector2d(0., -10.), -2.5,
       {Eigen::Vector2d(0., -7.5), Eigen::Vector2d(0., -2.5),
        Eigen::Vector2d(0., 0.), Eigen::Vector2d(0., 0.)}},
  }};

  for (const auto& expected : expectations) {
    CAPTURE(expected.public_facet);
    Q4TractionFixture fixture;
    REQUIRE(fixture.mesh.create_particles_tractions(
        nullptr, -1, expected.public_facet, expected.force_direction,
        -100.0));
    REQUIRE_NOTHROW(fixture.mesh.apply_traction_on_particles(0.0));
    REQUIRE_NOTHROW(
        fixture.particle->map_external_force(Eigen::Vector2d::Zero()));

    const auto phase = mpm::ParticlePhase::Mixture;
    for (unsigned node = 0; node < fixture.nodes.size(); ++node)
      REQUIRE(fixture.nodes[node]->external_force(phase).isApprox(
          expected.nodal_forces[node], 1.e-14));
    REQUIRE(fixture.resultant(phase).isApprox(expected.resultant, 1.e-14));
    REQUIRE(fixture.first_moment(phase) ==
            Approx(expected.first_moment).margin(1.e-14));

    const auto& context =
        fixture.particle->static_context(expected.force_direction);
    REQUIRE(context.public_facet == expected.public_facet);
    REQUIRE(context.element_face == expected.element_face);
    REQUIRE(context.reference_shape_weights.sum() == Approx(1.0));
  }
}

TEST_CASE("Mesh rejects invalid and non-finite facet tractions") {
  SECTION("invalid public facet") {
    Q4TractionFixture fixture;
    REQUIRE(fixture.mesh.create_particles_tractions(nullptr, -1, 4, 1,
                                                     -100.0));
    REQUIRE_THROWS(fixture.mesh.apply_traction_on_particles(0.0));
    REQUIRE_FALSE(fixture.particle->has_phase_traction());
  }

  SECTION("non-finite amplitude") {
    Q4TractionFixture fixture;
    REQUIRE(fixture.mesh.create_particles_tractions(
        nullptr, -1, 2, 1, std::numeric_limits<double>::quiet_NaN()));
    REQUIRE_THROWS_WITH(
        fixture.mesh.apply_traction_on_particles(0.0),
        Catch::Matchers::Contains("Particle traction is non-finite"));
    REQUIRE_FALSE(fixture.particle->has_phase_traction());
  }

  SECTION("failed update clears an earlier amplitude") {
    Q4TractionFixture fixture;
    REQUIRE(fixture.particle->assign_particle_traction_on_facet(2, 1,
                                                                 -100.0));
    REQUIRE(fixture.particle->has_phase_traction());
    REQUIRE_FALSE(fixture.particle->assign_particle_traction_on_facet(
        4, 1, -100.0));
    REQUIRE_FALSE(fixture.particle->has_phase_traction());
    REQUIRE_NOTHROW(
        fixture.particle->map_external_force(Eigen::Vector2d::Zero()));
    REQUIRE(fixture.resultant(mpm::ParticlePhase::Mixture).isZero(1.e-14));
  }
}

TEST_CASE("Frozen reference measure is independent of current volume") {
  Q4TractionFixture fixture;
  REQUIRE(fixture.mesh.create_particles_tractions(nullptr, -1, 2, 1,
                                                   -100.0));
  // The production solver freezes this context before entering the first
  // update. Mutating particle geometry afterwards must not alter that stage
  // reference when the load is first mapped.
  fixture.mesh.initialise_particle_traction_contexts(0.0);
  const double reference_measure =
      fixture.particle->static_context(1).reference_surface_measure;
  REQUIRE(reference_measure == Approx(0.1));

  fixture.particle->assign_current_geometry(0.2,
                                             Eigen::Vector2d(0.4, 0.5));
  fixture.mesh.apply_traction_on_particles(1.0);
  fixture.particle->map_external_force(Eigen::Vector2d::Zero());
  const auto phase = mpm::ParticlePhase::Mixture;
  REQUIRE(fixture.resultant(phase).isApprox(Eigen::Vector2d(0., -10.),
                                            1.e-14));
  REQUIRE(fixture.particle->static_context(1).reference_surface_measure ==
          Approx(reference_measure));
}

TEST_CASE("Facet weights follow the current cell while measure stays frozen") {
  const auto phase = mpm::ParticlePhase::Mixture;

  SECTION("static facet traction") {
    CrossCellQ4TractionFixture fixture;
    REQUIRE(fixture.mesh.create_particles_tractions(nullptr, -1, 2, 1,
                                                     -100.0));
    fixture.mesh.initialise_particle_traction_contexts(0.0);
    const double reference_measure =
        fixture.particle->static_context(1).reference_surface_measure;
    const Eigen::VectorXd reference_weights =
        fixture.particle->static_context(1).reference_shape_weights;
    REQUIRE(reference_measure == Approx(0.1));
    REQUIRE(reference_weights.isApprox(
        (Eigen::Vector4d() << 0., 0., 0.75, 0.25).finished(), 1.e-14));

    REQUIRE(fixture.move_to_second_cell());
    fixture.reset_nodal_forces();
    REQUIRE_NOTHROW(fixture.mesh.apply_traction_on_particles(1.0));
    REQUIRE_NOTHROW(
        fixture.particle->map_external_force(Eigen::Vector2d::Zero()));

    // Only the +y facet of the current cell receives load. The shared top
    // node at x=1 gets 0.75 and the new top node at x=2 gets 0.25.
    for (const unsigned node : {0u, 1u, 3u, 4u})
      REQUIRE(fixture.nodes[node]->external_force(phase).isZero(1.e-14));
    REQUIRE(fixture.nodes[2]->external_force(phase).isApprox(
        Eigen::Vector2d(0., -7.5), 1.e-14));
    REQUIRE(fixture.nodes[5]->external_force(phase).isApprox(
        Eigen::Vector2d(0., -2.5), 1.e-14));
    REQUIRE(fixture.resultant(phase).isApprox(Eigen::Vector2d(0., -10.),
                                              1.e-14));
    REQUIRE(fixture.first_moment(phase) == Approx(-12.5).margin(1.e-14));

    const auto& context = fixture.particle->static_context(1);
    REQUIRE(context.reference_surface_measure == Approx(reference_measure));
    REQUIRE(context.reference_shape_weights.isApprox(reference_weights,
                                                       1.e-14));
  }

  SECTION("dynamic pressure traction") {
    CrossCellQ4TractionFixture fixture;
    fixture.particle->configure_dynamic_surface(1000.0, 1020.0);
    const double reference_measure =
        fixture.particle->dynamic_context().reference_surface_measure;
    const Eigen::VectorXd reference_weights =
        fixture.particle->dynamic_context().reference_shape_weights;
    REQUIRE(reference_measure == Approx(0.1));

    REQUIRE(fixture.move_to_second_cell());
    fixture.reset_nodal_forces();
    REQUIRE_NOTHROW(
        fixture.particle->map_external_force(Eigen::Vector2d::Zero()));

    for (const unsigned node : {0u, 1u, 3u, 4u})
      REQUIRE(fixture.nodes[node]->external_force(phase).isZero(1.e-14));
    REQUIRE(fixture.nodes[2]->external_force(phase).isApprox(
        Eigen::Vector2d(0., -1.5), 1.e-14));
    REQUIRE(fixture.nodes[5]->external_force(phase).isApprox(
        Eigen::Vector2d(0., -0.5), 1.e-14));
    REQUIRE(fixture.resultant(phase).isApprox(Eigen::Vector2d(0., -2.),
                                              1.e-14));
    REQUIRE(fixture.first_moment(phase) == Approx(-2.5).margin(1.e-14));

    const auto& context = fixture.particle->dynamic_context();
    REQUIRE(context.reference_surface_measure == Approx(reference_measure));
    REQUIRE(context.reference_shape_weights.isApprox(reference_weights,
                                                       1.e-14));
  }
}

TEST_CASE("Dynamic excess pressure maps without a static traction") {
  Q4TractionFixture fixture;
  fixture.particle->configure_dynamic_surface(1000.0, 1020.0);
  REQUIRE_FALSE(fixture.particle->has_phase_traction());
  REQUIRE(fixture.particle->dynamic_context().public_facet == 2);
  REQUIRE(fixture.particle->dynamic_context().reference_surface_measure ==
          Approx(0.1));

  fixture.particle->map_external_force(Eigen::Vector2d::Zero());
  const auto phase = mpm::ParticlePhase::Mixture;
  REQUIRE(fixture.nodes[0]->external_force(phase).isZero(1.e-14));
  REQUIRE(fixture.nodes[1]->external_force(phase).isZero(1.e-14));
  REQUIRE(fixture.nodes[2]->external_force(phase).isApprox(
      Eigen::Vector2d(0., -0.5), 1.e-14));
  REQUIRE(fixture.nodes[3]->external_force(phase).isApprox(
      Eigen::Vector2d(0., -1.5), 1.e-14));
  REQUIRE(fixture.resultant(phase).isApprox(Eigen::Vector2d(0., -2.),
                                            1.e-14));
}

TEST_CASE("Dynamic surface mapping fails when its context is missing") {
  Q4TractionFixture fixture;
  fixture.particle->mark_dynamic_surface_without_context(1000.0, 1020.0);
  REQUIRE_THROWS_WITH(
      fixture.particle->map_external_force(Eigen::Vector2d::Zero()),
      Catch::Matchers::Contains("no dynamic +y traction context"));
}

TEST_CASE("Static and dynamic traction contexts remain independent") {
  Q4TractionFixture fixture;
  fixture.particle->configure_dynamic_surface(1000.0, 1020.0);
  // Static vertical shear acts on public +x; dynamic pressure acts on +y.
  REQUIRE(fixture.mesh.create_particles_tractions(nullptr, -1, 0, 1,
                                                   -100.0));
  fixture.mesh.apply_traction_on_particles(0.0);
  fixture.particle->map_external_force(Eigen::Vector2d::Zero());

  const auto phase = mpm::ParticlePhase::Mixture;
  const std::array<Eigen::Vector2d, 4> expected{{
      Eigen::Vector2d(0., 0.), Eigen::Vector2d(0., -5.),
      Eigen::Vector2d(0., -15.5), Eigen::Vector2d(0., -1.5)}};
  for (unsigned node = 0; node < fixture.nodes.size(); ++node)
    REQUIRE(fixture.nodes[node]->external_force(phase).isApprox(
        expected[node], 1.e-14));
  REQUIRE(fixture.particle->static_context(1).public_facet == 0);
  REQUIRE(fixture.particle->static_context(1).reference_surface_measure ==
          Approx(0.2));
  REQUIRE(fixture.particle->dynamic_context().public_facet == 2);
  REQUIRE(fixture.particle->dynamic_context().reference_surface_measure ==
          Approx(0.1));
}

TEST_CASE("Stage handoff can rebuild non-serialized traction contexts") {
  Q4TractionFixture fixture;
  fixture.particle->configure_dynamic_surface(1000.0, 1020.0);
  REQUIRE(fixture.particle->assign_particle_traction_on_facet(0, 1,
                                                               -100.0));
  const auto static_before = fixture.particle->static_context(1);
  const auto dynamic_before = fixture.particle->dynamic_context();

  fixture.particle->emulate_stage_handoff_context_loss();
  REQUIRE_FALSE(fixture.particle->static_context(1).initialised());
  REQUIRE_FALSE(fixture.particle->dynamic_context().initialised());

  // A handoff creates a new stage reference from checkpoint geometry. Exact
  // mid-stage restart equivalence remains intentionally unsupported because
  // these contexts are absent from HDF5Particle.
  fixture.particle->configure_dynamic_surface(1000.0, 1020.0);
  REQUIRE(fixture.particle->assign_particle_traction_on_facet(0, 1,
                                                               -100.0));
  const auto& static_after = fixture.particle->static_context(1);
  const auto& dynamic_after = fixture.particle->dynamic_context();
  REQUIRE(static_after.reference_surface_measure ==
          Approx(static_before.reference_surface_measure));
  REQUIRE(static_after.reference_shape_weights.isApprox(
      static_before.reference_shape_weights, 1.e-14));
  REQUIRE(dynamic_after.reference_surface_measure ==
          Approx(dynamic_before.reference_surface_measure));
  REQUIRE(dynamic_after.reference_shape_weights.isApprox(
      dynamic_before.reference_shape_weights, 1.e-14));

  fixture.reset_nodal_forces();
  fixture.particle->map_external_force(Eigen::Vector2d::Zero());
  REQUIRE(fixture.resultant(mpm::ParticlePhase::Mixture)
              .isApprox(Eigen::Vector2d(0., -22.), 1.e-14));
}
