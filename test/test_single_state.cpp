#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "sho.h"
#include "gtest/gtest.h"

// single_state is an abstract class, so there are only a few things to test
// in it. We'll create an instance of the "sho" class, which is the simplest
// child of single_state, to test.

TEST(single_state, set_all_bead_states) {
  const size_t nbead(2);
  const int init_state(0);

  pot::sho my_pot;
  my_pot.set_all_bead_states(init_state, nbead);

  EXPECT_EQ(my_pot.state_id.size(), nbead);
  for (auto &id : my_pot.state_id) {
    EXPECT_EQ(id, init_state);
  }
}

TEST(single_state, set_individual_bead_states) {
  std::vector<int> state_ids = {0, 1};
  const size_t nbead = state_ids.size();

  pot::sho my_pot;
  my_pot.set_individual_bead_states(state_ids);

  EXPECT_EQ(my_pot.state_id.size(), nbead);
  for (size_t i = 0; i < nbead; ++i) {
    EXPECT_EQ(my_pot.state_id[i], state_ids[i]);
  }
}

TEST(single_state, avg_active_state) {
  const int nbead(2);
  const int init_state(1);

  pot::sho my_pot;
  my_pot.set_all_bead_states(init_state, nbead);

  EXPECT_EQ(my_pot.avg_active_state(), 1);
}

TEST(single_state, sum_active_state) {
  const int nbead(2);
  const int init_state(1);

  pot::sho my_pot;
  my_pot.set_all_bead_states(init_state, nbead);

  EXPECT_EQ(my_pot.sum_active_state(), 2);
}
