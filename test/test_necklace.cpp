#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "necklace.h"
#include "gtest/gtest.h"

const double tol(1.0E-8);

TEST(necklace, setup) {
  const int nbead(4);
  const int ndim(1);
  const int natom(1);
  const int ndof(ndim * natom);

  parts::necklace necklace;
  necklace.setup(ndof, nbead);

  EXPECT_NEAR(necklace.m_lambda[0], 0.0, tol);
  EXPECT_NEAR(necklace.m_lambda[1], 15.999999999999998, tol);
  EXPECT_NEAR(necklace.m_lambda[2], 16, tol);
  EXPECT_NEAR(necklace.m_lambda[3], 15.999999999999998, tol);
}

// This function will test the transformation of the positions and velenta
// from cartesian to normal-mode representation and back to cartesian.
TEST(necklace, normal_mode_transformations) {
  const int nbead(2);
  const int ndim(2);
  const int natom(1);
  const int ndof(ndim * natom);

  double beta(1024);
  double dt(1.0);
  double mass(2000);
  parts::necklace necklace;
  necklace.setup(ndof, nbead);

  necklace.m_pos_cart[0 + 0 * ndof] = 1.1;
  necklace.m_pos_cart[0 + 1 * ndof] = 1.15;
  necklace.m_pos_cart[1 + 0 * ndof] = 0.02;
  necklace.m_pos_cart[1 + 1 * ndof] = 1.01;
  double orig_pos_cart[4];
  std::copy(necklace.m_pos_cart, necklace.m_pos_cart + 4, orig_pos_cart);

  necklace.m_vel_cart[0 + 0 * ndof] = 0.1;
  necklace.m_vel_cart[0 + 1 * ndof] = 0.15;
  necklace.m_vel_cart[1 + 0 * ndof] = -0.08;
  necklace.m_vel_cart[1 + 1 * ndof] = 0.01;
  double orig_vel_cart[4];
  std::copy(necklace.m_vel_cart, necklace.m_vel_cart + 4, orig_vel_cart);

  necklace.pos_c2n();
  EXPECT_NEAR(necklace.m_pos_nmode[0 + 0 * ndof], 1.125, tol);
  EXPECT_NEAR(necklace.m_pos_nmode[0 + 1 * ndof], -0.024999999999999911, tol);
  EXPECT_NEAR(necklace.m_pos_nmode[1 + 0 * ndof], 0.51500000000000001, tol);
  EXPECT_NEAR(necklace.m_pos_nmode[1 + 1 * ndof], -0.495, tol);

  necklace.pos_n2c();
  EXPECT_NEAR(necklace.m_pos_cart[0 + 0 * ndof], orig_pos_cart[0 + 0 * ndof],
              tol);
  EXPECT_NEAR(necklace.m_pos_cart[0 + 1 * ndof], orig_pos_cart[0 + 1 * ndof],
              tol);
  EXPECT_NEAR(necklace.m_pos_cart[1 + 0 * ndof], orig_pos_cart[1 + 0 * ndof],
              tol);
  EXPECT_NEAR(necklace.m_pos_cart[1 + 1 * ndof], orig_pos_cart[1 + 1 * ndof],
              tol);

  necklace.vel_c2n();
  EXPECT_NEAR(necklace.m_vel_nmode[0 + 0 * ndof], 0.125, tol);
  EXPECT_NEAR(necklace.m_vel_nmode[0 + 1 * ndof], -0.024999999999999994, tol);
  EXPECT_NEAR(necklace.m_vel_nmode[1 + 0 * ndof], -0.035000000000000003, tol);
  EXPECT_NEAR(necklace.m_vel_nmode[1 + 1 * ndof], -0.044999999999999998, tol);

  necklace.vel_n2c();
  EXPECT_NEAR(necklace.m_vel_cart[0 + 0 * ndof], orig_vel_cart[0 + 0 * ndof],
              tol);
  EXPECT_NEAR(necklace.m_vel_cart[0 + 1 * ndof], orig_vel_cart[0 + 1 * ndof],
              tol);
  EXPECT_NEAR(necklace.m_vel_cart[1 + 0 * ndof], orig_vel_cart[1 + 0 * ndof],
              tol);
  EXPECT_NEAR(necklace.m_vel_cart[1 + 1 * ndof], orig_vel_cart[1 + 1 * ndof],
              tol);

  necklace.m_frc_cart[0 + 0 * ndof] = -1.1;
  necklace.m_frc_cart[0 + 1 * ndof] = -1.15;
  necklace.m_frc_cart[1 + 0 * ndof] = 1.08;
  necklace.m_frc_cart[1 + 1 * ndof] = -1.01;

  necklace.frc_c2n();
  EXPECT_NEAR(necklace.m_frc_nmode[0 + 0 * ndof], -1.125, tol);
  EXPECT_NEAR(necklace.m_frc_nmode[0 + 1 * ndof], 0.024999999999999911, tol);
  EXPECT_NEAR(necklace.m_frc_nmode[1 + 0 * ndof], 0.035000000000000031, tol);
  EXPECT_NEAR(necklace.m_frc_nmode[1 + 1 * ndof], 1.0449999999999999, tol);
}
