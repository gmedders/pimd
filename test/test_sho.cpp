#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "finite-difference.h"
#include "sho.h"
#include "gtest/gtest.h"

double omega(0.001);   // omega
double atm_mass(2000); // au
double params[] = {omega, atm_mass};

TEST(sho, set_params) {

  pot::sho my_pot;
  my_pot.set_params(params);

  EXPECT_DOUBLE_EQ(my_pot.get_w(), omega);
  EXPECT_DOUBLE_EQ(my_pot.get_m(), atm_mass);
}

TEST(sho, VAA_energy) {

  pot::sho my_pot;
  my_pot.set_params(params);

  std::vector<double> crd = {0.5};
  double f[] = {0.0};

  EXPECT_DOUBLE_EQ(my_pot.VAA(&crd[0], f), 0.0002500000000000000);
}

TEST(sho, VAA_grad) {

  pot::sho my_pot;
  my_pot.set_params(params);

  std::vector<double> crd = {0.5};
  std::vector<double> fd_f = {0.0};

  run_finite_differences(my_pot, crd, fd_f);

  double f[] = {0.0};
  my_pot.VAA(&crd[0], f);
  EXPECT_NEAR(f[0], fd_f[0], 1.0E-8);
}
