#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "rpmd_necklace.h"
#include "gtest/gtest.h"

TEST(rpmd_necklace, setup) {
  const int nbead(2);
  const int ndim(1);
  const int natom(1);
  const int ndof(ndim * natom);

  double beta(1024);
  double dt(1.0);
  double mass(2000);
  parts::rpmd_necklace necklace;
  necklace.setup(ndof, nbead, beta, dt, mass);

  EXPECT_NEAR(necklace.m_cart_to_nm(0, 0), 0.70710678118654757, 1.0E-8);
  EXPECT_NEAR(necklace.m_cart_to_nm(0, 1), 0.70710678118654757, 1.0E-8);
  EXPECT_NEAR(necklace.m_cart_to_nm(1, 0), 0.70710678118654757, 1.0E-8);
  EXPECT_NEAR(necklace.m_cart_to_nm(1, 1), -0.70710678118654757, 1.0E-8);

  EXPECT_NEAR(necklace.m_omega_k(0), 0.0, 1.0E-8);
  EXPECT_NEAR(necklace.m_omega_k(1), 0.00390625, 1.0E-8);

  EXPECT_NEAR(necklace.m_freerp_propagator(0, 0, 0), 1.0, 1.0E-8);
  EXPECT_NEAR(necklace.m_freerp_propagator(0, 1, 0), 0.00050000000000000001,
              1.0E-8);
  EXPECT_NEAR(necklace.m_freerp_propagator(1, 0, 0), 0.0, 1.0E-8);
  EXPECT_NEAR(necklace.m_freerp_propagator(1, 1, 0), 1.0, 1.0E-8);

  EXPECT_NEAR(necklace.m_freerp_propagator(0, 0, 1), 0.99999237061516999,
              1.0E-8);
  EXPECT_NEAR(necklace.m_freerp_propagator(0, 1, 1), 0.0004999987284352149,
              1.0E-8);
  EXPECT_NEAR(necklace.m_freerp_propagator(1, 0, 1), -0.030517500514844659,
              1.0E-8);
  EXPECT_NEAR(necklace.m_freerp_propagator(1, 1, 1), 0.99999237061516999,
              1.0E-8);
}
