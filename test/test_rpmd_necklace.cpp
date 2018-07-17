#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "rpmd_necklace.h"
#include "gtest/gtest.h"

const double tol(1.0E-8);

TEST(rpmd_necklace, setup) {
  const int nbead(2);
  const int ndim(1);
  const int natom(1);
  const int ndof(ndim * natom);

  double beta(1024);
  double dt(1.0);
  double mass(2000);
  parts::rpmd_necklace necklace;
  necklace.setup(ndim, natom, nbead, beta, dt, mass);

  EXPECT_NEAR(necklace.m_cart_to_nm(0, 0), 0.70710678118654757, tol);
  EXPECT_NEAR(necklace.m_cart_to_nm(0, 1), 0.70710678118654757, tol);
  EXPECT_NEAR(necklace.m_cart_to_nm(1, 0), 0.70710678118654757, tol);
  EXPECT_NEAR(necklace.m_cart_to_nm(1, 1), -0.70710678118654757, tol);

  EXPECT_NEAR(necklace.m_omega_k(0), 0.0, tol);
  EXPECT_NEAR(necklace.m_omega_k(1), 0.00390625, tol);

  EXPECT_NEAR(necklace.m_freerp_propagator(0, 0, 0), 1.0, tol);
  EXPECT_NEAR(necklace.m_freerp_propagator(0, 1, 0), 0.00050000000000000001,
              tol);
  EXPECT_NEAR(necklace.m_freerp_propagator(1, 0, 0), 0.0, tol);
  EXPECT_NEAR(necklace.m_freerp_propagator(1, 1, 0), 1.0, tol);

  EXPECT_NEAR(necklace.m_freerp_propagator(0, 0, 1), 0.99999237061516999, tol);
  EXPECT_NEAR(necklace.m_freerp_propagator(0, 1, 1), 0.0004999987284352149,
              tol);
  EXPECT_NEAR(necklace.m_freerp_propagator(1, 0, 1), -0.030517500514844659,
              tol);
  EXPECT_NEAR(necklace.m_freerp_propagator(1, 1, 1), 0.99999237061516999, tol);
}

// This function will test the transformation of the positions and momenta
// from cartesian to normal-mode representation and back to cartesian.
TEST(rpmd_necklace, normal_mode_transformations) {
  const int nbead(2);
  const int ndim(2);
  const int natom(1);
  const int ndof(ndim * natom);

  double beta(1024);
  double dt(1.0);
  double mass(2000);
  parts::rpmd_necklace necklace;
  necklace.setup(ndim, natom, nbead, beta, dt, mass);

  // ndof is first index, nbead is second
  necklace.m_pos_cart(0, 0) = 1.1;
  necklace.m_pos_cart(0, 1) = 1.15;
  necklace.m_pos_cart(1, 0) = 0.02;
  necklace.m_pos_cart(1, 1) = 1.01;

  necklace.m_mom_cart(0, 0) = 0.1;
  necklace.m_mom_cart(0, 1) = 0.15;
  necklace.m_mom_cart(1, 0) = -0.08;
  necklace.m_mom_cart(1, 1) = 0.01;

  arma::mat orig_pos_cart(necklace.m_pos_cart);
  arma::mat orig_mom_cart(necklace.m_mom_cart);

  necklace.pos_c2n();
  EXPECT_NEAR(necklace.m_pos_nmode(0, 0), 1.5909902576697319, tol);
  EXPECT_NEAR(necklace.m_pos_nmode(0, 1), -0.035355339059327195, tol);
  EXPECT_NEAR(necklace.m_pos_nmode(1, 0), 0.72831998462214398, tol);
  EXPECT_NEAR(necklace.m_pos_nmode(1, 1), -0.70003571337468218, tol);

  // Make sure that pos_n2c() is actually doing something by zeroing out the
  // m_pos_cart matrix and testing that it is not equal to orig_pos_cart
  necklace.m_pos_cart.zeros();
  EXPECT_FALSE(arma::norm(orig_pos_cart - necklace.m_pos_cart, "inf") < tol);
  necklace.pos_n2c();
  EXPECT_TRUE(arma::norm(orig_pos_cart - necklace.m_pos_cart, "inf") < tol);

  necklace.mom_c2n();
  EXPECT_NEAR(necklace.m_mom_nmode(0, 0), 0.17677669529663689, tol);
  EXPECT_NEAR(necklace.m_mom_nmode(0, 1), -0.035355339059327362, tol);
  EXPECT_NEAR(necklace.m_mom_nmode(1, 0), -0.049497474683058332, tol);
  EXPECT_NEAR(necklace.m_mom_nmode(1, 1), -0.063639610306789288, tol);

  // Make sure that mom_n2c() is actually doing something by zeroing out the
  // m_mom_cart matrix and testing that it is not equal to orig_mom_cart
  necklace.m_mom_cart.zeros();
  EXPECT_FALSE(arma::norm(orig_mom_cart - necklace.m_mom_cart, "inf") < tol);
  necklace.mom_n2c();
  EXPECT_TRUE(arma::norm(orig_mom_cart - necklace.m_mom_cart, "inf") < tol);
}
