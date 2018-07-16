#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

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
  double m_omega_M = kT;
  double m_tau = 2 * M_PI / m_omega_M;

  for (size_t b = 0; b < nbead; ++b) {
    for (size_t i = 0; i < ndof; ++i) {
      const size_t j = i + b * ndof;
      parts::nhc::initialize(nchain, m_thermostats + th_size * j, m_tau, prng);
    }
  }

  EXPECT_NEAR(m_thermostats[4], 0.00014073773847819479, 1.0E-8);
  EXPECT_NEAR(m_thermostats[5], -0.0001492116584795932, 1.0E-8);
  EXPECT_NEAR(m_thermostats[9], -4.3497416049379006e-09, 1.0E-8);
  EXPECT_DOUBLE_EQ(m_thermostats[12], 0.0);
  EXPECT_NEAR(m_thermostats[16], -0.0001059985446239687, 1.0E-8);
  EXPECT_NEAR(m_thermostats[19], 0.00016685709141472772, 1.0E-8);

  delete[] m_thermostats;
}

TEST(nhc, test_advance) {
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
  double m_omega_M = kT;
  double m_tau = 2 * M_PI / m_omega_M;

  for (size_t b = 0; b < nbead; ++b) {
    for (size_t i = 0; i < ndof; ++i) {
      const size_t j = i + b * ndof;
      parts::nhc::initialize(nchain, m_thermostats + th_size * j, m_tau, prng);
    }
  }

  EXPECT_NEAR(m_thermostats[4], 0.00014073773847819479, 1.0E-8);
  EXPECT_NEAR(m_thermostats[5], -0.0001492116584795932, 1.0E-8);
  EXPECT_NEAR(m_thermostats[9], -4.3497416049379006e-09, 1.0E-8);
  EXPECT_DOUBLE_EQ(m_thermostats[12], 0.0);
  EXPECT_NEAR(m_thermostats[16], -0.0001059985446239687, 1.0E-8);
  EXPECT_NEAR(m_thermostats[19], 0.00016685709141472772, 1.0E-8);

  double dt(1.0);

  double beta_step2(1200);
  double Ek2kT = 2.0 / beta_step2;

  double accu(0.0);
  for (size_t n = 0; n < nbead * ndof; ++n)
    accu += parts::nhc::invariant(nchain, m_thermostats + n * th_size, m_tau);
  EXPECT_NEAR(accu, 2.5914802362361549, 1.0E-8);

  std::vector<double> all_aa;
  for (size_t b = 0; b < nbead; ++b) {
    for (size_t i = 0; i < ndof; ++i) {
      const size_t j = i + b * ndof;
      double aa = parts::nhc::advance(nchain, m_thermostats + th_size * j,
                                      m_tau, Ek2kT, dt);
      all_aa.push_back(aa);
    }
  }

  EXPECT_NEAR(all_aa[0], 0.99985927371848149, 1.0E-8);
  EXPECT_NEAR(all_aa[1], 1.0001060115130245, 1.0E-8);

  for (size_t i = 0; i < all_aa.size(); ++i) {
    std::cout << "all_aa[" << i << "] = " << all_aa[i] << std::endl;
  }

  EXPECT_NEAR(m_thermostats[4], 0.00014073773847819479, 1.0E-8);
  EXPECT_NEAR(m_thermostats[5], -0.00014923867549715387, 1.0E-8);
  EXPECT_NEAR(m_thermostats[9], -4.3497416049379006e-09, 1.0E-8);
  EXPECT_NEAR(m_thermostats[12], -0.00010600589420104945, 1.0E-8);
  EXPECT_NEAR(m_thermostats[16], -0.00010601324344690028, 1.0E-8);
  EXPECT_NEAR(m_thermostats[19], 0.00016683983088406498, 1.0E-8);

  accu = 0.0;
  for (size_t n = 0; n < nbead * ndof; ++n)
    accu += parts::nhc::invariant(nchain, m_thermostats + n * th_size, m_tau);
  EXPECT_NEAR(accu, 2.5914802940690054, 1.0E-8);

  delete[] m_thermostats;
}
