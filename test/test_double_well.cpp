#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "double_well.h"
#include "finite_difference.h"
#include "gtest/gtest.h"

typedef pot::double_well potential;
double omega(0.0009765625);
double atm_mass(2000);
double param_g(3.9);
double dG(-0.003906252);
double params[] = {omega, atm_mass, param_g, dG};

TEST(double_well, set_params) {

  potential my_pot;
  my_pot.set_params(params);

  EXPECT_DOUBLE_EQ(my_pot.get_w(), omega);
  EXPECT_DOUBLE_EQ(my_pot.get_m(), atm_mass);
  EXPECT_DOUBLE_EQ(my_pot.get_g(), param_g);
  EXPECT_DOUBLE_EQ(my_pot.get_dG(), dG);
  EXPECT_DOUBLE_EQ(my_pot.get_a(), 0.00095367431640625);
}

TEST(double_well, VAA_energy) {

  potential my_pot;
  my_pot.set_params(params);

  std::vector<double> crd = {0.5};
  double f[] = {0.0};

  EXPECT_DOUBLE_EQ(my_pot.VAA(&crd[0], f), 0.0002384185791015625);
}

TEST(double_well, VAA_grad) {

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

TEST(double_well, VBB_energy) {

  potential my_pot;
  my_pot.set_params(params);

  std::vector<double> crd = {0.5};
  double f[] = {0.0};

  EXPECT_DOUBLE_EQ(my_pot.VBB(&crd[0], f), 0.0071182230976562502);
}

TEST(double_well, VBB_grad) {

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
