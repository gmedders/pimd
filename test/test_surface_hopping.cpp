#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "double-well.h"
#include "gtest/gtest.h"

// surface_hopping is an abstract class, so there are only a few things to test
// in it. We'll create instances of the "double_well" class, which is the
// simplest child of surface_hopping, to test.

// bool init_pot = false;
// bool init_hop = false;
// double Gamma;
// double dt;
// double beta;
// double voltage;
//
// arma::ivec state_id;
//
// int nhops;
//
// double force(size_t, size_t, size_t, const double *, double *);
// virtual double VAA(const double *x, double *f) = 0;
// virtual double VBB(const double *x, double *f) = 0;
// virtual void set_params(double *) = 0;
//
// void check_allocation(size_t, arma::ivec &);
//
// double fermi_function(const double, const double);
// double hop_probability(const double, int);
// void hop();
//
// void set_all_bead_states(const int, int);
// void set_individual_bead_states(std::vector<int> &);
// void set_hopping_params(double *);
//
// double avg_active_state();
// double sum_active_state();
//
// void print_state_params();

double omega(0.0009765625); // omega
double atm_mass(2000);      // au
double param_g(3.9);        // barrier = 2w
double dG(-0.003906252);
double params[] = {omega, atm_mass, param_g, dG};

TEST(surface_hopping, hopping_params) {
  pot::double_well my_pot;

  double GammaEl(1.0);
  double dt(2.0);
  double beta(1024.0);
  double voltage(3.0);
  double hop_params[] = {GammaEl, dt, beta, voltage};

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
