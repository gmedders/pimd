#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "anharmonic.h"
#include "finite-difference.h"
#include "gtest/gtest.h"

double atm_mass(1); // au
double param_a(1.0 / 2.0);
double param_b(1.0 / 10.0);
double param_c(1.0 / 100.0);
double params[] = {param_a, param_b, param_c, atm_mass};

TEST(anharmonic, set_params) {

  pot::anharmonic my_pot;
  my_pot.set_params(params);

  EXPECT_DOUBLE_EQ(my_pot.get_w(), std::sqrt(2.0 * param_a / atm_mass));
  EXPECT_DOUBLE_EQ(my_pot.get_m(), atm_mass);
  EXPECT_DOUBLE_EQ(my_pot.get_a(), param_a);
  EXPECT_DOUBLE_EQ(my_pot.get_b(), param_b);
  EXPECT_DOUBLE_EQ(my_pot.get_c(), param_c);
}

TEST(anharmonic, VAA_energy) {

  pot::anharmonic my_pot;
  my_pot.set_params(params);

  std::vector<double> crd = {0.5};
  double f[] = {0.0};

  EXPECT_DOUBLE_EQ(my_pot.VAA(&crd[0], f), 0.138125);
}

TEST(anharmonic, VAA_grad) {

  pot::anharmonic my_pot;
  my_pot.set_params(params);

  std::vector<double> crd = {0.5};
  std::vector<double> fd_f = {0.0};

  run_finite_differences(my_pot, crd, fd_f);

  double f[] = {0.0};
  my_pot.VAA(&crd[0], f);
  EXPECT_NEAR(f[0], fd_f[0], 1.0E-8);
}
