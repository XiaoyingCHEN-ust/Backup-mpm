#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <array>
#include <memory>
#include <stdexcept>

#include "mesh.h"
#include "node.h"
#include "particle.h"
#include "quadrilateral_element.h"

namespace {

class ThrowingDensityNode final : public mpm::Node<2, 2, 3> {
 public:
  using mpm::Node<2, 2, 3>::Node;

  void compute_density() override {
    throw std::runtime_error("intentional density failure");
  }
};

class FreeSurfaceParticleProbe final : public mpm::Particle<2> {
 public:
  using mpm::Particle<2>::Particle;

  bool compute_particle_free_surface() override { return true; }
};

}  // namespace

TEST_CASE("Three-phase support threshold is consistent across one time step") {
  Eigen::Vector2d coordinates(0., 0.);
  mpm::Node<2, 2, 3> node(0, coordinates);
  node.assign_minimum_nodal_density(1.0);

  // Solver order: nodal initialise, one cell-volume map, particle mass and
  // momentum map, initial velocity, force map, then acceleration update.
  node.initialise();
  node.update_volume(true, 0, 1.E-4, 0.02);
  node.update_mass_momentum(true, 0, 5.E-5, Eigen::Vector2d(2., -4.));
  node.update_mass_momentum(true, 1, 2.E-5, Eigen::Vector2d(10., 20.));
  node.update_mass_momentum(true, 2, 1.E-8, Eigen::Vector2d(-21., 12.));

  node.compute_velocity(1.E-4);

  REQUIRE(node.velocity(0).isZero());
  REQUIRE(node.velocity(1).isZero());
  REQUIRE(node.velocity(2).isZero());

  node.update_external_force(true, 0, Eigen::Vector2d(1.E8, -1.E8));
  REQUIRE(node.compute_acc_vel_threephase_explicit(0, 1, 2, 0, 1.E-4));
  REQUIRE(node.acceleration(0).isZero());
  REQUIRE(node.acceleration(1).isZero());
  REQUIRE(node.acceleration(2).isZero());
}

TEST_CASE("Three-phase support just above the threshold remains active") {
  Eigen::Vector2d coordinates(0., 0.);
  mpm::Node<2, 2, 3> node(0, coordinates);
  node.assign_minimum_nodal_density(1.0);
  node.initialise();
  node.update_volume(true, 0, 1.E-4, 0.02);
  node.update_mass_momentum(true, 0, 1.5E-4,
                            Eigen::Vector2d(3.E-4, -6.E-4));
  node.update_mass_momentum(true, 1, 2.E-5,
                            Eigen::Vector2d(2.E-5, 4.E-5));
  node.update_mass_momentum(true, 2, 1.E-8,
                            Eigen::Vector2d(-3.E-8, 4.E-8));

  node.compute_velocity(1.E-4);
  REQUIRE(node.velocity(0).isApprox(Eigen::Vector2d(2., -4.)));
  REQUIRE(node.velocity(1).isApprox(Eigen::Vector2d(1., 2.)));
  REQUIRE(node.velocity(2).isApprox(Eigen::Vector2d(-3., 4.)));

  node.update_external_force(true, 0, Eigen::Vector2d(1.5E-4, 0.));
  REQUIRE(node.compute_acc_vel_threephase_explicit(0, 1, 2, 0, 1.E-4));
  REQUIRE(node.acceleration(0)[0] == Approx(1.));
}

TEST_CASE(
    "Three-phase free-surface constraints remain authoritative over contact") {
  Eigen::Vector2d coordinates(0., 0.);
  mpm::Node<2, 2, 3> node(0, coordinates);
  node.initialise();
  node.assign_free_surface(true);
  node.assign_phase_kinematic_boundary(true);

  node.update_mass_momentum(true, 0, 2., Eigen::Vector2d(2., -4.));
  node.update_mass_momentum(true, 1, 1., Eigen::Vector2d(5., 6.));
  node.update_mass_momentum(true, 2, 0.5, Eigen::Vector2d(-3.5, 4.));

  // Exercise the same precedence used by the explicit solver: free-surface
  // synchronisation first, rigid contact second, and prescribed velocity
  // constraints last. The three x constraints are deliberately different.
  const double dt = 0.1;
  node.assign_velocity_from_rigid(0, 0.9, dt);
  node.assign_velocity_from_rigid(1, 0.8, dt);
  node.assign_velocity_from_rigid(2, 0.7, dt);
  node.assign_velocity_from_rigid(3, 0.6, dt);
  node.assign_velocity_from_rigid(4, 0.5, dt);
  node.assign_velocity_from_rigid(5, 0.4, dt);
  REQUIRE(node.assign_velocity_constraint(0, 0.1));
  REQUIRE(node.assign_velocity_constraint(2, -0.2));
  REQUIRE(node.assign_velocity_constraint(4, 0.3));

  node.compute_velocity(dt);
  REQUIRE(node.velocity(0)[0] == Approx(0.1));
  REQUIRE(node.velocity(1)[0] == Approx(-0.2));
  REQUIRE(node.velocity(2)[0] == Approx(0.3));
  REQUIRE(node.acceleration(0)[0] == Approx(0.));
  REQUIRE(node.acceleration(1)[0] == Approx(0.));
  REQUIRE(node.acceleration(2)[0] == Approx(0.));

  const Eigen::Vector2d mixture_force(7., 11.);
  node.update_external_force(true, 0, mixture_force);
  node.update_external_force(true, 1, Eigen::Vector2d(3., 5.));
  node.update_external_force(true, 2, Eigen::Vector2d(-2., 4.));
  REQUIRE(node.compute_acc_vel_threephase_explicit(0, 1, 2, 0, dt));

  REQUIRE(node.velocity(0)[0] == Approx(0.1));
  REQUIRE(node.velocity(1)[0] == Approx(-0.2));
  REQUIRE(node.velocity(2)[0] == Approx(0.3));
  REQUIRE(node.acceleration(0)[0] == Approx(0.));
  REQUIRE(node.acceleration(1)[0] == Approx(0.));
  REQUIRE(node.acceleration(2)[0] == Approx(0.));
  REQUIRE(node.reaction_force().isApprox(-mixture_force));
}

TEST_CASE("Free-surface mode and detection failures propagate") {
  SECTION("unknown mode") {
    mpm::Mesh<2> mesh(0);
    REQUIRE_THROWS_AS(mesh.compute_free_surface("typo", 0.5, false),
                      std::invalid_argument);
  }

  SECTION("node density exception") {
    mpm::Mesh<2> mesh(0);
    const std::array<Eigen::Vector2d, 4> coordinates{
        Eigen::Vector2d(0., 0.), Eigen::Vector2d(1., 0.),
        Eigen::Vector2d(1., 1.), Eigen::Vector2d(0., 1.)};
    std::array<std::shared_ptr<mpm::NodeBase<2>>, 4> nodes;
    nodes[0] = std::make_shared<ThrowingDensityNode>(0, coordinates[0]);
    for (unsigned id = 1; id < nodes.size(); ++id)
      nodes[id] =
          std::make_shared<mpm::Node<2, 2, 3>>(id, coordinates[id]);
    for (const auto& node : nodes) REQUIRE(mesh.add_node(node));

    auto element = std::make_shared<mpm::QuadrilateralElement<2, 4>>();
    auto cell = std::make_shared<mpm::Cell<2>>(0, 4, element);
    for (unsigned id = 0; id < nodes.size(); ++id)
      REQUIRE(cell->add_node(id, nodes[id]));
    REQUIRE(cell->initialise());
    REQUIRE(mesh.add_cell(cell));

    auto particle =
        std::make_shared<mpm::Particle<2>>(0, Eigen::Vector2d(0.5, 0.5));
    REQUIRE(particle->assign_initial_volume(0.25));
    REQUIRE(particle->assign_cell_xi(cell, Eigen::Vector2d::Zero()));
    REQUIRE(mesh.add_particle(particle, false));
    cell->activate_nodes();

    REQUIRE_THROWS(mesh.compute_free_surface("detect", 0.5, false));
  }

  SECTION("assigned kinematics reject higher-order surface support") {
    mpm::Mesh<2> mesh(0);
    const std::array<Eigen::Vector2d, 8> coordinates{
        Eigen::Vector2d(0., 0.),   Eigen::Vector2d(1., 0.),
        Eigen::Vector2d(1., 1.),   Eigen::Vector2d(0., 1.),
        Eigen::Vector2d(0.5, 0.),  Eigen::Vector2d(1., 0.5),
        Eigen::Vector2d(0.5, 1.),  Eigen::Vector2d(0., 0.5)};
    std::array<std::shared_ptr<mpm::Node<2, 2, 3>>, 8> nodes;
    for (unsigned id = 0; id < nodes.size(); ++id) {
      nodes[id] =
          std::make_shared<mpm::Node<2, 2, 3>>(id, coordinates[id]);
      REQUIRE(mesh.add_node(nodes[id]));
    }

    auto element = std::make_shared<mpm::QuadrilateralElement<2, 8>>();
    auto cell = std::make_shared<mpm::Cell<2>>(0, 8, element);
    for (unsigned id = 0; id < nodes.size(); ++id)
      REQUIRE(cell->add_node(id, nodes[id]));
    REQUIRE(cell->initialise());
    REQUIRE(mesh.add_cell(cell));

    auto particle =
        std::make_shared<mpm::Particle<2>>(0, Eigen::Vector2d(0.5, 0.5));
    REQUIRE(particle->assign_initial_volume(0.25));
    REQUIRE(particle->assign_cell_xi(cell, Eigen::Vector2d::Zero()));
    particle->assign_particle_free_surface(true);
    REQUIRE(mesh.add_particle(particle, false));

    REQUIRE_THROWS_AS(mesh.compute_free_surface("assign", 0.5, false),
                      std::invalid_argument);
  }
}

TEST_CASE("Assigned and detected interfaces own distinct phase kinematics") {
  mpm::Mesh<2> mesh(0);
  auto element = std::make_shared<mpm::QuadrilateralElement<2, 4>>();

  const std::array<Eigen::Vector2d, 16> coordinates{
      Eigen::Vector2d(0., 0.), Eigen::Vector2d(1., 0.),
      Eigen::Vector2d(1., 1.), Eigen::Vector2d(0., 1.),
      Eigen::Vector2d(1., 2.), Eigen::Vector2d(0., 2.),
      Eigen::Vector2d(2., 0.), Eigen::Vector2d(3., 0.),
      Eigen::Vector2d(3., 1.), Eigen::Vector2d(2., 1.),
      Eigen::Vector2d(3., 2.), Eigen::Vector2d(2., 2.),
      Eigen::Vector2d(4., 0.), Eigen::Vector2d(5., 0.),
      Eigen::Vector2d(5., 1.), Eigen::Vector2d(4., 1.)};
  std::array<std::shared_ptr<mpm::Node<2, 2, 3>>, 16> nodes;
  for (unsigned id = 0; id < nodes.size(); ++id) {
    nodes[id] = std::make_shared<mpm::Node<2, 2, 3>>(id, coordinates[id]);
    REQUIRE(mesh.add_node(nodes[id]));
  }

  auto assigned_surface_cell =
      std::make_shared<mpm::Cell<2>>(0, 4, element);
  auto assigned_empty_neighbour =
      std::make_shared<mpm::Cell<2>>(1, 4, element);
  auto geometric_cavity_cell =
      std::make_shared<mpm::Cell<2>>(2, 4, element);
  auto cavity_empty_neighbour =
      std::make_shared<mpm::Cell<2>>(3, 4, element);
  auto interior_cell = std::make_shared<mpm::Cell<2>>(4, 4, element);
  const std::array<std::shared_ptr<mpm::Cell<2>>, 5> cells{
      assigned_surface_cell, assigned_empty_neighbour,
      geometric_cavity_cell, cavity_empty_neighbour, interior_cell};
  const std::array<std::array<unsigned, 4>, 5> cell_node_ids{{
      {{0, 1, 2, 3}}, {{3, 2, 4, 5}}, {{6, 7, 8, 9}},
      {{9, 8, 10, 11}}, {{12, 13, 14, 15}}}};
  for (unsigned cell_id = 0; cell_id < cells.size(); ++cell_id) {
    for (unsigned local_id = 0; local_id < 4; ++local_id) {
      const unsigned node_id = cell_node_ids[cell_id][local_id];
      REQUIRE(cells[cell_id]->add_node(local_id, nodes[node_id]));
    }
    REQUIRE(cells[cell_id]->initialise());
    REQUIRE(mesh.add_cell(cells[cell_id]));
  }
  REQUIRE(assigned_surface_cell->add_neighbour(
      assigned_empty_neighbour->id()));
  REQUIRE(geometric_cavity_cell->add_neighbour(
      cavity_empty_neighbour->id()));

  auto surface_particle = std::make_shared<FreeSurfaceParticleProbe>(
      0, Eigen::Vector2d(0.5, 0.5));
  auto cavity_particle = std::make_shared<FreeSurfaceParticleProbe>(
      1, Eigen::Vector2d(2.5, 0.5));
  auto interior_particle =
      std::make_shared<FreeSurfaceParticleProbe>(2,
                                                 Eigen::Vector2d(4.5, 0.5));
  // Deliberately under-fill both boundary cells so every one of their nodes
  // is geometrically free. Assigned kinematic ownership must nevertheless be
  // restricted to the marker cell's positive-y facet.
  REQUIRE(surface_particle->assign_initial_volume(0.25));
  REQUIRE(cavity_particle->assign_initial_volume(0.25));
  REQUIRE(interior_particle->assign_initial_volume(0.75));
  REQUIRE(surface_particle->assign_cell_xi(assigned_surface_cell,
                                            Eigen::Vector2d::Zero()));
  REQUIRE(cavity_particle->assign_cell_xi(geometric_cavity_cell,
                                           Eigen::Vector2d::Zero()));
  REQUIRE(interior_particle->assign_cell_xi(interior_cell,
                                             Eigen::Vector2d::Zero()));
  surface_particle->assign_particle_free_surface(true);
  REQUIRE(mesh.add_particle(surface_particle, false));
  REQUIRE(mesh.add_particle(cavity_particle, false));
  REQUIRE(mesh.add_particle(interior_particle, false));

  // Match one explicit-solver nodal cycle: initialise clears stale geometric
  // and kinematic flags, then active support, volume and momentum are mapped
  // before the same-step interface pass.
  for (const auto& node : nodes) {
    node->assign_free_surface(true);
    node->assign_phase_kinematic_boundary(true);
  }
  mesh.iterate_over_nodes(
      [](const std::shared_ptr<mpm::NodeBase<2>>& node) {
        node->initialise();
        REQUIRE_FALSE(node->free_surface());
        REQUIRE_FALSE(node->phase_kinematic_boundary());
      });
  mesh.iterate_over_cells(
      [](const std::shared_ptr<mpm::Cell<2>>& cell) { cell->activate_nodes(); });
  mesh.iterate_over_particles(
      [](const std::shared_ptr<mpm::ParticleBase<2>>& particle) {
        particle->compute_shapefn();
      });
  mesh.iterate_over_cells([](const std::shared_ptr<mpm::Cell<2>>& cell) {
    REQUIRE(cell->map_cell_volume_to_nodes(0));
  });

  const Eigen::Vector2d solid_momentum(2., -4.);
  const Eigen::Vector2d liquid_momentum(5., 6.);
  const Eigen::Vector2d gas_momentum(-3.5, 4.);
  for (const auto& node : nodes) {
    if (!node->status()) continue;
    node->update_mass_momentum(true, 0, 2., solid_momentum);
    node->update_mass_momentum(true, 1, 1., liquid_momentum);
    node->update_mass_momentum(true, 2, 0.5, gas_momentum);
  }

  const auto verify_phase_kinematics = [&nodes]() {
    constexpr std::array<unsigned, 4> representative_nodes{{2, 0, 8, 12}};
    for (const unsigned id : representative_nodes)
      nodes[id]->compute_velocity(0.1);
    for (const unsigned id : representative_nodes)
      REQUIRE(nodes[id]->velocity(0).isApprox(Eigen::Vector2d(1., -2.)));

    for (const unsigned id : representative_nodes) {
      const bool synchronised = nodes[id]->phase_kinematic_boundary();
      if (synchronised) {
        REQUIRE(nodes[id]->velocity(1).isApprox(nodes[id]->velocity(0)));
        REQUIRE(nodes[id]->velocity(2).isApprox(nodes[id]->velocity(0)));
      } else {
        REQUIRE(nodes[id]->velocity(1).isApprox(Eigen::Vector2d(5., 6.)));
        REQUIRE(nodes[id]->velocity(2).isApprox(Eigen::Vector2d(-7., 8.)));
      }
    }

    const Eigen::Vector2d mixture_force(7., 11.);
    const Eigen::Vector2d liquid_force(3., 5.);
    const Eigen::Vector2d gas_force(-2., 4.);
    for (const unsigned id : representative_nodes) {
      nodes[id]->update_external_force(true, 0, mixture_force);
      nodes[id]->update_external_force(true, 1, liquid_force);
      nodes[id]->update_external_force(true, 2, gas_force);
      REQUIRE(nodes[id]->compute_acc_vel_threephase_explicit(0, 1, 2, 0,
                                                              0.1));
    }

    for (const unsigned id : representative_nodes) {
      const bool synchronised = nodes[id]->phase_kinematic_boundary();
      REQUIRE(nodes[id]->velocity(1).isApprox(nodes[id]->velocity(0)) ==
              synchronised);
      REQUIRE(nodes[id]->velocity(2).isApprox(nodes[id]->velocity(0)) ==
              synchronised);
      REQUIRE(nodes[id]->acceleration(1).isApprox(nodes[id]->acceleration(0)) ==
              synchronised);
      REQUIRE(nodes[id]->acceleration(2).isApprox(nodes[id]->acceleration(0)) ==
              synchronised);
    }
  };

  SECTION("assign synchronises only the marker cell positive-y facet") {
    REQUIRE(mesh.compute_free_surface("assign", 0.5, false));
    // All nodes are geometrically free, including lower/side support nodes,
    // but only Q4 local {2,3} owns the assigned phase condition.
    REQUIRE(nodes[0]->free_surface());
    REQUIRE(nodes[1]->free_surface());
    REQUIRE(nodes[2]->free_surface());
    REQUIRE(nodes[3]->free_surface());
    REQUIRE_FALSE(nodes[0]->phase_kinematic_boundary());
    REQUIRE_FALSE(nodes[1]->phase_kinematic_boundary());
    REQUIRE(nodes[2]->phase_kinematic_boundary());
    REQUIRE(nodes[3]->phase_kinematic_boundary());
    // The cavity is still a geometric interface for pressure and output, but
    // it is not supported by the configured surface particle.
    REQUIRE(nodes[6]->free_surface());
    REQUIRE(nodes[8]->free_surface());
    REQUIRE_FALSE(nodes[6]->phase_kinematic_boundary());
    REQUIRE_FALSE(nodes[8]->phase_kinematic_boundary());
    REQUIRE_FALSE(nodes[12]->free_surface());
    REQUIRE_FALSE(nodes[12]->phase_kinematic_boundary());
    verify_phase_kinematics();
  }

  SECTION("detect synchronises every geometric interface") {
    REQUIRE(mesh.compute_free_surface("detect", 0.5, false));
    REQUIRE(nodes[0]->free_surface());
    REQUIRE(nodes[2]->free_surface());
    REQUIRE(nodes[0]->phase_kinematic_boundary());
    REQUIRE(nodes[2]->phase_kinematic_boundary());
    REQUIRE(nodes[6]->free_surface());
    REQUIRE(nodes[8]->free_surface());
    REQUIRE(nodes[6]->phase_kinematic_boundary());
    REQUIRE(nodes[8]->phase_kinematic_boundary());
    REQUIRE_FALSE(nodes[12]->free_surface());
    REQUIRE_FALSE(nodes[12]->phase_kinematic_boundary());
    verify_phase_kinematics();
  }

  SECTION("assign rebuilds ownership after vertical marker relocation") {
    REQUIRE(mesh.compute_free_surface("assign", 0.5, false));
    REQUIRE(nodes[2]->phase_kinematic_boundary());
    REQUIRE(nodes[3]->phase_kinematic_boundary());

    surface_particle->assign_coordinates(Eigen::Vector2d(4.5, 0.5));
    REQUIRE(surface_particle->assign_cell_xi(interior_cell,
                                              Eigen::Vector2d::Zero()));
    REQUIRE(mesh.compute_free_surface("assign", 0.5, false));

    // The old facet is cleared and the relocated current cell contributes
    // only its positive-y Q4 local nodes {2,3}, even though that cell is not a
    // geometric free-surface cell.
    REQUIRE_FALSE(nodes[2]->phase_kinematic_boundary());
    REQUIRE_FALSE(nodes[3]->phase_kinematic_boundary());
    REQUIRE_FALSE(nodes[12]->phase_kinematic_boundary());
    REQUIRE_FALSE(nodes[13]->phase_kinematic_boundary());
    REQUIRE(nodes[14]->phase_kinematic_boundary());
    REQUIRE(nodes[15]->phase_kinematic_boundary());
    REQUIRE_FALSE(nodes[6]->phase_kinematic_boundary());
    REQUIRE_FALSE(nodes[8]->phase_kinematic_boundary());
  }
}

TEST_CASE("Assigned free-surface pass preserves the solver volume mapping") {
  mpm::Mesh<2> mesh(0);
  const std::array<Eigen::Vector2d, 4> coordinates{
      Eigen::Vector2d(0., 0.), Eigen::Vector2d(1., 0.),
      Eigen::Vector2d(1., 1.), Eigen::Vector2d(0., 1.)};
  std::array<std::shared_ptr<mpm::Node<2, 2, 3>>, 4> nodes;
  for (unsigned id = 0; id < nodes.size(); ++id) {
    nodes[id] = std::make_shared<mpm::Node<2, 2, 3>>(id, coordinates[id]);
    REQUIRE(mesh.add_node(nodes[id]));
  }
  auto element = std::make_shared<mpm::QuadrilateralElement<2, 4>>();
  auto cell = std::make_shared<mpm::Cell<2>>(0, 4, element);
  for (unsigned id = 0; id < nodes.size(); ++id)
    REQUIRE(cell->add_node(id, nodes[id]));
  REQUIRE(cell->initialise());
  REQUIRE(cell->is_initialised());
  REQUIRE(mesh.add_cell(cell));

  auto particle =
      std::make_shared<mpm::Particle<2>>(0, Eigen::Vector2d(0.5, 0.5));
  REQUIRE(particle->assign_initial_volume(0.25));
  REQUIRE(particle->assign_cell_xi(cell, Eigen::Vector2d::Zero()));
  REQUIRE(mesh.add_particle(particle, false));

  mesh.iterate_over_cells([](const std::shared_ptr<mpm::Cell<2>>& cell) {
    REQUIRE(cell->map_cell_volume_to_nodes(0));
  });
  for (const auto& node : nodes) REQUIRE(node->volume(0) == Approx(0.25));

  // This is the exact explicit-solver call after its earlier volume map.
  REQUIRE(mesh.compute_free_surface("assign", 0.5, false));
  for (const auto& node : nodes) REQUIRE(node->volume(0) == Approx(0.25));
}
