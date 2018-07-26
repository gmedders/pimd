#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "sim_classes.h"
#include "gtest/gtest.h"

TEST(pimd_base, make_one_step) {
  double pos[] = {0.1, 0.13};
  double vel[] = {-0.01, -0.015};

  const size_t nbead(2);
  const size_t ndim(1);
  const size_t natom(1);

  double beta(1024);
  double dt(1.0);

  parts::pimd sim;
  sim.set_up(ndim, natom, nbead, beta, dt, pos, vel);

  const size_t ndof(1);
  EXPECT_EQ(sim.ndofs(), ndof);
  EXPECT_EQ(sim.nbeads(), nbead);

  EXPECT_NEAR(sim.invariant(), 1.9025476256457383, 1.0E-8);
  EXPECT_NEAR(sim.Ep(), 1.3450000000000002e-05, 1.0E-8);
  EXPECT_NEAR(sim.Ek(), 1.8999999999999999, 1.0E-8);
  EXPECT_NEAR(sim.Espring(), 3.4332275390625005e-06, 1.0E-8);
  EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.115);
  EXPECT_NEAR(sim.L1_cart_pos(), 0.014999999999999999, 1.0E-8);
  EXPECT_NEAR(sim.L2_cart_pos(), 0.010606601717798213, 1.0E-8);
  EXPECT_NEAR(sim.Linf_cart_pos(), 0.014999999999999999, 1.0E-8);

  // Make one step
  sim.step(dt);

  EXPECT_NEAR(sim.invariant(), 1.902547632980558, 1.0E-8);
  EXPECT_NEAR(sim.Ep(), 1.1925246007206764e-05, 1.0E-8);
  EXPECT_NEAR(sim.Ek(), 1.9001695691296427, 1.0E-8);
  EXPECT_NEAR(sim.Espring(), 1.3733619862773376e-05, 1.0E-8);
  EXPECT_NEAR(sim.temp_kT(), 1.9001695691296427, 1.0E-8);
  EXPECT_NEAR(sim.avg_cart_pos(), 0.10500094997608378, 1.0E-8);
  EXPECT_NEAR(sim.L1_cart_pos(), 0.03000077517876356, 1.0E-8);
  EXPECT_NEAR(sim.L2_cart_pos(), 0.021213751569756776, 1.0E-8);
  EXPECT_NEAR(sim.Linf_cart_pos(), 0.030000775178763567, 1.0E-8);
}

// TEST(pimd_base, no_nhc_make_one_step) {
//   double pos[] = {0.1, 0.13};
//   double vel[] = {-0.01, -0.015};
//
//   const size_t nbead(2);
//   const size_t ndim(1);
//   const size_t natom(1);
//
//   double beta(1024);
//   double dt(1.0);
//
//   parts::pimd sim;
//   sim.set_up(ndim, natom, nbead, beta, dt, pos, vel);
//
//   EXPECT_NEAR(sim.invariant(), 1.900016883227539, 1.0E-8);
//   EXPECT_NEAR(sim.Ep(), 1.3450000000000002e-05, 1.0E-8);
//   EXPECT_NEAR(sim.Ek(), 1.8999999999999999, 1.0E-8);
//   EXPECT_NEAR(sim.Espring(), 3.4332275390625005e-06, 1.0E-8);
//   EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.115);
//   EXPECT_NEAR(sim.L1_cart_pos(), 0.014999999999999999, 1.0E-8);
//   EXPECT_NEAR(sim.L2_cart_pos(), 0.010606601717798213, 1.0E-8);
//   EXPECT_NEAR(sim.Linf_cart_pos(), 0.014999999999999999, 1.0E-8);
//
//   // Make one step
//   sim.step(dt);
//
//   EXPECT_NEAR(sim.invariant(), 1.9000168832325652, 1.0E-8);
//   EXPECT_NEAR(sim.Ep(), 1.1924987010446653e-05, 1.0E-8);
//   EXPECT_NEAR(sim.Ek(), 1.8999912253493534, 1.0E-8);
//   EXPECT_NEAR(sim.Espring(), 1.3732896201222957e-05, 1.0E-8);
//   EXPECT_NEAR(sim.temp_kT(), 1.8999912253493534, 1.0E-8);
//   EXPECT_NEAR(sim.avg_cart_pos(), 0.1049999425, 1.0E-8);
//   EXPECT_NEAR(sim.L1_cart_pos(), 0.029999984757385245, 1.0E-8);
//   EXPECT_NEAR(sim.L2_cart_pos(), 0.021213192657440171, 1.0E-8);
//   EXPECT_NEAR(sim.Linf_cart_pos(), 0.029999984757385245, 1.0E-8);
// }
