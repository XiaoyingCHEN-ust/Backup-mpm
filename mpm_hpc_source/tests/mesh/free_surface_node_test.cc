#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <array>
#include <memory>

#include "mesh.h"
#include "node.h"
#include "particle.h"
#include "quadrilateral_element.h"

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
