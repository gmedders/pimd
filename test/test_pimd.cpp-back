#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "sim_classes.h"
#include "gtest/gtest.h"

TEST(pimd_base, make_one_step) {
  double pos[] = {0.1, 0.13};
  double vel[] = {-0.01, -0.015};

  const int nbead(2);
  const int ndim(1);
  const int natom(1);

  double beta(1024);
  double dt(1.0);

  parts::pimd sim;
  sim.set_up(nbead, ndim, natom, beta);
  std::copy(pos, pos + 2, sim.m_pos_cart);
  sim.pos_c2n();
  std::copy(vel, vel + 2, sim.m_vel_cart);
  sim.vel_c2n();

  EXPECT_EQ(sim.ndofs(), 1);

  EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.115);
  EXPECT_NEAR(sim.L1_cart_pos(), 0.014999999999999999, 1.0E-8);
  EXPECT_NEAR(sim.L2_cart_pos(), 0.010606601717798213, 1.0E-8);
  EXPECT_NEAR(sim.Linf_cart_pos(), 0.014999999999999999, 1.0E-8);

  // Make one step
  sim.step(dt);

  EXPECT_DOUBLE_EQ(sim.invariant(), 0.20879557840065055);
  EXPECT_DOUBLE_EQ(sim.Ep(), 1.0662755303760603e-05);
  EXPECT_DOUBLE_EQ(sim.Ek(), 0.20620193665119149);
  EXPECT_NEAR(sim.Espring(), 2.3841148950758333e-06, 1.0E-8);
  EXPECT_DOUBLE_EQ(sim.temp_kT(), 0.20620193665119149);
  EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.10250126804092188);
  EXPECT_NEAR(sim.L1_cart_pos(), 0.012499814149165969, 1.0E-8);
  EXPECT_NEAR(sim.L2_cart_pos(), 0.0088387033484468112, 1.0E-8);
  EXPECT_NEAR(sim.Linf_cart_pos(), 0.012499814149165969, 1.0E-8);
}
