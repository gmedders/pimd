#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "mt19937.h"
#include "nhc.h"
#include "gtest/gtest.h"

TEST(nhc, test_initialization) {
  const int nbead(2);
  const int ndim(1);
  const int natom(1);
  const int ndof(ndim * natom);

  double beta(1024);

  int nchain(4);
  const size_t th_size = parts::nhc::size(nchain);
  double *m_thermostats = new double[ndof * nbead * th_size];

  std::mt19937 prng;
  prng.seed(19107);

  double kT = 1.0 / beta;
  double hbar(1.0);
  double m_omega_M = kT / hbar;
  double m_tau = 2 * M_PI / m_omega_M;

  for (size_t b = 0; b < nbead; ++b) {
    for (size_t i = 0; i < ndof; ++i) {
      const size_t j = i + b * ndof;
      parts::nhc::initialize(nchain, m_thermostats + th_size * j, m_tau, prng);
    }
  }

  for (size_t i = 0; i < ndof * nbead * th_size; ++i) {
    std::cout << "m_thermostats[" << i << "] = " << m_thermostats[i]
              << std::endl;
  }

  EXPECT_NEAR(m_thermostats[4], 0.00014073773847819479, 1.0E-8);
  EXPECT_NEAR(m_thermostats[5], -0.0001492116584795932, 1.0E-8);
  EXPECT_NEAR(m_thermostats[9], -4.3497416049379006e-09, 1.0E-8);
  EXPECT_DOUBLE_EQ(m_thermostats[12], 0.0);
  EXPECT_NEAR(m_thermostats[16], -0.0001059985446239687, 1.0E-8);
  EXPECT_NEAR(m_thermostats[19], 0.00016685709141472772, 1.0E-8);

  delete[] m_thermostats;
}
