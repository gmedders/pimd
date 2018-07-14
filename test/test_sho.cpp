#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "sho.h"
#include "gtest/gtest.h"

// single_state is an abstract class, so there are only a few things to test
// in it. We'll create an instance of the "sho" class, which is the simplest
// child of single_state, to test.

template <class potential>
void run_finite_differences(potential my_pot, std::vector<double> &crd,
                            std::vector<double> &f) {
  double dx = 5.0E-7;
  std::vector<double> weight = {1.0 / 12.0, -2.0 / 3.0, 2.0 / 3.0, -1.0 / 12.0};
  std::vector<double> displacement = {-2.0 * dx, -1.0 * dx, 1.0 * dx, 2.0 * dx};
  std::vector<double> dummy_grad = f;

  for (std::vector<double>::size_type i(0); i < crd.size(); ++i) {
    double grad(0.0);
    double orig(crd[i]);
    for (std::vector<double>::size_type j(0); j < displacement.size(); ++j) {
      crd[i] = orig + displacement[j];
      grad += weight[j] * my_pot.VAA(&crd[0], &dummy_grad[0]);
    }
    crd[i] = orig;
    f[i] = -grad / dx;
  }
}

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
