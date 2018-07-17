#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "sim_classes.h"
#include "gtest/gtest.h"

// This velocity_verlet test assumes the following is set up as the potential:
// typedef pot::sho potential_type;
// static double omega(0.001);   // omega
// static double atm_mass(2000); // au
// static double params[] = {omega, atm_mass};

TEST(velocity_verlet, make_one_step) {
  double pos[] = {0.1};
  double vel[] = {-0.01};

  const int nbead(1);
  const int ndim(1);
  const int natom(1);

  double beta(1024);
  double dt(1.0);

  parts::vv sim;
  sim.set_up(nbead, ndim, natom, beta, dt, pos, vel);

  EXPECT_EQ(sim.ndofs(), 1);

  EXPECT_DOUBLE_EQ(sim.invariant(), 0.10001);
  EXPECT_DOUBLE_EQ(sim.Ep(), 1e-05);
  EXPECT_DOUBLE_EQ(sim.Ek(), 0.1000000000000000);
  EXPECT_NEAR(sim.Espring(), 0.0, 1.0E-8);
  EXPECT_NEAR(sim.temp_kT(), 0.0, 1.0E-8);
  EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.10000000000000001);
  EXPECT_NEAR(sim.L1_cart_pos(), 0.0, 1.0E-8);
  EXPECT_NEAR(sim.L2_cart_pos(), 0.0, 1.0E-8);
  EXPECT_NEAR(sim.Linf_cart_pos(), 0.0, 1.0E-8);

  // Make one step
  sim.step(dt);

  EXPECT_DOUBLE_EQ(sim.invariant(), 0.10000999999952499);
  EXPECT_DOUBLE_EQ(sim.Ep(), 8.0999910000025029e-06);
  EXPECT_DOUBLE_EQ(sim.Ek(), 0.10000190000852499);
  EXPECT_NEAR(sim.Espring(), 0.0, 1.0E-8);
  EXPECT_DOUBLE_EQ(sim.temp_kT(), 0.20000380001704998);
  EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.089999950000000009);
  EXPECT_NEAR(sim.L1_cart_pos(), 0.0, 1.0E-8);
  EXPECT_NEAR(sim.L2_cart_pos(), 0.0, 1.0E-8);
  EXPECT_NEAR(sim.Linf_cart_pos(), 0.0, 1.0E-8);
}
