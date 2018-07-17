#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "double_well.h"
#include "gtest/gtest.h"

// surface_hopping is an abstract class, so there are only a few things to test
// in it. We'll create instances of the "double_well" class, which is the
// simplest child of surface_hopping, to test.

double GammaEl(1.0);
double dt(2.0);
double beta(1024.0);
double voltage(0.0001);
double hop_params[] = {GammaEl, dt, beta, voltage};

TEST(surface_hopping, hopping_params) {
  pot::double_well my_pot;
  my_pot.set_hopping_params(hop_params);

  EXPECT_EQ(my_pot.Gamma, GammaEl);
  EXPECT_EQ(my_pot.dt, dt);
  EXPECT_EQ(my_pot.beta, beta);
  EXPECT_EQ(my_pot.voltage, voltage);
}

TEST(surface_hopping, set_all_bead_states) {
  const int nbead(2);
  const int init_state(0);

  pot::double_well my_pot;
  my_pot.set_all_bead_states(init_state, nbead);

  EXPECT_EQ(my_pot.state_id.size(), nbead);
  for (auto &id : my_pot.state_id) {
    EXPECT_EQ(id, init_state);
  }
}

TEST(surface_hopping, set_individual_bead_states) {
  std::vector<int> state_ids = {0, 1};
  const int nbead = state_ids.size();

  pot::double_well my_pot;
  my_pot.set_individual_bead_states(state_ids);

  EXPECT_EQ(my_pot.state_id.size(), nbead);
  for (int i = 0; i < nbead; ++i) {
    EXPECT_EQ(my_pot.state_id[i], state_ids[i]);
  }
}

TEST(surface_hopping, fermi_function) {
  double mu(-0.02);
  double dE(0.001);

  pot::double_well my_pot;
  my_pot.set_hopping_params(hop_params);

  double f_mu0 = my_pot.fermi_function(dE, 0.0);
  EXPECT_DOUBLE_EQ(f_mu0, 0.26424898168976962);
  double f_mu = my_pot.fermi_function(dE, mu);
  EXPECT_DOUBLE_EQ(f_mu, 4.5806958984748794e-10);
  double f_huge_dE = my_pot.fermi_function(1000, mu);
  EXPECT_DOUBLE_EQ(f_huge_dE, 0.0);
}

TEST(surface_hopping, hop_probability) {
  double dE(0.001);
  int active_state(0);

  pot::double_well my_pot;
  my_pot.set_hopping_params(hop_params);

  double p_hop_0 = my_pot.hop_probability(dE, 0);
  EXPECT_DOUBLE_EQ(p_hop_0, 0.52873820111413861);
  double p_hop_1 = my_pot.hop_probability(dE, 1);
  EXPECT_DOUBLE_EQ(p_hop_1, 1.4712617988858614);
}

TEST(surface_hopping, sum_active_state) {
  std::vector<int> state_ids = {0, 1, 0, 0, 0, 1, 0, 0};
  const int nbead = state_ids.size();

  pot::double_well my_pot;
  my_pot.set_individual_bead_states(state_ids);

  EXPECT_EQ(my_pot.state_id.size(), nbead);
  int sum(0);
  for (int i = 0; i < nbead; ++i) {
    sum += my_pot.state_id[i];
  }
  EXPECT_EQ(my_pot.sum_active_state(), sum);
  EXPECT_EQ(my_pot.avg_active_state(), double(sum) / double(nbead));
}
