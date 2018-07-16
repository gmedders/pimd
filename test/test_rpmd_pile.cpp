#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "sim_classes.h"
#include "gtest/gtest.h"

// TEST(rpmd_base, make_one_step) {
//   double pos[] = {0.1, 0.13};
//   double vel[] = {-0.01, -0.015};
//
//   const int nbead(2);
//   const int ndim(1);
//   const int natom(1);
//
//   double beta(1024);
//   double dt(1.0);
//
//   parts::rpmd sim;
//   sim.set_up(nbead, ndim, natom, beta, dt, pos, vel);
//
//   EXPECT_EQ(sim.ndofs(), 1);
//
//   EXPECT_DOUBLE_EQ(sim.invariant(), 0.32503376645507814);
//   EXPECT_DOUBLE_EQ(sim.Ep(), 2.6900000000000003e-05);
//   EXPECT_DOUBLE_EQ(sim.Ek(), 0.32500000000000001);
//   EXPECT_NEAR(sim.Espring(), 6.8664550781250002e-06, 1.0E-8);
//   EXPECT_NEAR(sim.temp_kT(), 0.32500000000000001, 1.0E-8);
//   EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.115);
//   EXPECT_NEAR(sim.L1_cart_pos(), 0.014999999999999999, 1.0E-8);
//   EXPECT_NEAR(sim.L2_cart_pos(), 0.010606601717798213, 1.0E-8);
//   EXPECT_NEAR(sim.Linf_cart_pos(), 0.014999999999999999, 1.0E-8);
//
//   // Make one step
//   sim.step(dt);
//
//   EXPECT_DOUBLE_EQ(sim.invariant(), 0.32503376645298498);
//   EXPECT_DOUBLE_EQ(sim.Ep(), 2.1324970645886869e-05);
//   EXPECT_DOUBLE_EQ(sim.Ek(), 0.32500767319893942);
//   EXPECT_NEAR(sim.Espring(), 4.7682833996619607e-06, 1.0E-8);
//   EXPECT_DOUBLE_EQ(sim.temp_kT(), 0.32500767319893942);
//   EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.10249994250000002);
//   EXPECT_NEAR(sim.L1_cart_pos(), 0.012499884417070548, 1.0E-8);
//   EXPECT_NEAR(sim.L2_cart_pos(), 0.0088387530353586392, 1.0E-8);
//   EXPECT_NEAR(sim.Linf_cart_pos(), 0.012499884417070548, 1.0E-8);
// }

TEST(rpmd_pile, make_one_step_no_centroid_thermostatting) {
  double pos[] = {0.1, 0.13};
  double vel[] = {-0.01, -0.015};

  const int nbead(2);
  const int ndim(1);
  const int natom(1);

  double beta(1024);
  double dt(1.0);

  parts::rpmd sim;
  sim.set_up(nbead, ndim, natom, beta, dt, pos, vel);
  sim.seed_pile_prng(19107);

  EXPECT_EQ(sim.ndofs(), 1);

  EXPECT_DOUBLE_EQ(sim.invariant(), 0.32503376645507814);
  EXPECT_DOUBLE_EQ(sim.Ep(), 2.6900000000000003e-05);
  EXPECT_DOUBLE_EQ(sim.Ek(), 0.32500000000000001);
  EXPECT_NEAR(sim.Espring(), 6.8664550781250002e-06, 1.0E-8);
  EXPECT_NEAR(sim.temp_kT(), 0.32500000000000001, 1.0E-8);
  EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.115);
  EXPECT_NEAR(sim.L1_cart_pos(), 0.014999999999999999, 1.0E-8);
  EXPECT_NEAR(sim.L2_cart_pos(), 0.010606601717798213, 1.0E-8);
  EXPECT_NEAR(sim.Linf_cart_pos(), 0.014999999999999999, 1.0E-8);

  // Make one step
  sim.step(dt);

  EXPECT_DOUBLE_EQ(sim.invariant(), 0.32456647989316095);
  EXPECT_DOUBLE_EQ(sim.Ep(), 2.1328426367135509e-05);
  EXPECT_DOUBLE_EQ(sim.Ek(), 0.32454033045327257);
  EXPECT_NEAR(sim.Espring(), 4.821013521253757e-06, 1.0E-8);
  EXPECT_DOUBLE_EQ(sim.temp_kT(), 0.32454033045327257);
  EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.10249994250000002);
  EXPECT_NEAR(sim.L1_cart_pos(), 0.01256880945294514, 1.0E-8);
  EXPECT_NEAR(sim.L2_cart_pos(), 0.0088874903956190892, 1.0E-8);
  EXPECT_NEAR(sim.Linf_cart_pos(), 0.01256880945294514, 1.0E-8);
}

// TEST(rpmd_pile, make_one_step_with_langevin) {
//   double pos[] = {0.1, 0.13};
//   double vel[] = {-0.01, -0.015};
//
//   const int nbead(2);
//   const int ndim(1);
//   const int natom(1);
//
//   double beta(1024);
//   double dt(1.0);
//
//   parts::rpmd sim;
//   sim.set_up(nbead, ndim, natom, beta, dt, pos, vel);
//   sim.set_gammaTh(0.002);
//
//   EXPECT_EQ(sim.ndofs(), 1);
//
//   EXPECT_DOUBLE_EQ(sim.invariant(), 0.32503376645507814);
//   EXPECT_DOUBLE_EQ(sim.Ep(), 2.6900000000000003e-05);
//   EXPECT_DOUBLE_EQ(sim.Ek(), 0.32500000000000001);
//   EXPECT_NEAR(sim.Espring(), 6.8664550781250002e-06, 1.0E-8);
//   EXPECT_NEAR(sim.temp_kT(), 0.0, 1.0E-8);
//   EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.115);
//   EXPECT_NEAR(sim.L1_cart_pos(), 0.014999999999999999, 1.0E-8);
//   EXPECT_NEAR(sim.L2_cart_pos(), 0.010606601717798213, 1.0E-8);
//   EXPECT_NEAR(sim.Linf_cart_pos(), 0.014999999999999999, 1.0E-8);
//
//   // Make one step
//   sim.step(dt);
//
//   EXPECT_DOUBLE_EQ(sim.invariant(), 0.32378867137049588);
//   EXPECT_DOUBLE_EQ(sim.Ep(), 2.1326257844156981e-05);
//   EXPECT_DOUBLE_EQ(sim.Ek(), 0.32376255718816516);
//   EXPECT_NEAR(sim.Espring(), 4.7879244865471981e-06, 1.0E-8);
//   EXPECT_DOUBLE_EQ(sim.temp_kT(), 0.32376255718816516);
//   EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.10249994250000002);
//   EXPECT_NEAR(sim.L1_cart_pos(), 0.012525602164174726, 1.0E-8);
//   EXPECT_NEAR(sim.L2_cart_pos(), 0.0088569382287328451, 1.0E-8);
//   EXPECT_NEAR(sim.Linf_cart_pos(), 0.012525602164174726, 1.0E-8);
// }
