#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "ah.h"
#include "finite_difference.h"
#include "gtest/gtest.h"

typedef pot::ah potential;
double omega(0.003);
double atm_mass(2000);
double param_g(0.02);
double param_Ed_bar(0.1333333);
double params[] = {omega, atm_mass, param_g, param_Ed_bar};

TEST(ah, set_params) {

  potential my_pot;
  my_pot.set_params(params);

  EXPECT_DOUBLE_EQ(my_pot.get_w(), omega);
  EXPECT_DOUBLE_EQ(my_pot.get_m(), atm_mass);
  EXPECT_DOUBLE_EQ(my_pot.get_g(), param_g);
  EXPECT_DOUBLE_EQ(my_pot.get_a(), 0.00900000000000000);
  EXPECT_DOUBLE_EQ(my_pot.get_b(), 0.069282032302755092);
  EXPECT_DOUBLE_EQ(my_pot.get_Ed(), 0.26666663333333329);
  EXPECT_DOUBLE_EQ(my_pot.get_Ed_bar(), param_Ed_bar);
}

TEST(ah, VAA_energy) {

  potential my_pot;
  my_pot.set_params(params);

  std::vector<double> crd = {0.5};
  double f[] = {0.0};

  EXPECT_DOUBLE_EQ(my_pot.VAA(&crd[0], f), 0.002250000000000000);
}

TEST(ah, VAA_grad) {

  potential my_pot;
  my_pot.set_params(params);

  std::vector<double> crd = {0.5};
  std::vector<double> fd_f = {0.0};

  auto fp = std::bind(&potential::VAA, my_pot, _1, _2);
  run_finite_differences(fp, crd, fd_f);

  double f[] = {0.0};
  my_pot.VAA(&crd[0], f);
  EXPECT_NEAR(f[0], fd_f[0], 1.0E-8);
}

TEST(ah, VBB_energy) {

  potential my_pot;
  my_pot.set_params(params);

  std::vector<double> crd = {0.5};
  double f[] = {0.0};

  EXPECT_DOUBLE_EQ(my_pot.VBB(&crd[0], f), 0.3035576494847108);
}

TEST(ah, VBB_grad) {

  potential my_pot;
  my_pot.set_params(params);

  std::vector<double> crd = {0.5};
  std::vector<double> fd_f = {0.0};

  auto fp = std::bind(&potential::VBB, my_pot, _1, _2);
  run_finite_differences(fp, crd, fd_f);

  double f[] = {0.0};
  my_pot.VBB(&crd[0], f);
  EXPECT_NEAR(f[0], fd_f[0], 1.0E-8);
}
