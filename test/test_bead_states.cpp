#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "bead_states.h"
#include "gtest/gtest.h"

TEST(bead_states, set_all_bead_states) {
  const size_t nbead(2);
  const int init_state(0);

  pot::bead_states my_beads;
  my_beads.set_all_bead_states(init_state, nbead);

  EXPECT_EQ(my_beads.state_id.size(), nbead);
  for (auto &id : my_beads.state_id) {
    EXPECT_EQ(id, init_state);
  }
}

TEST(bead_states, set_individual_bead_states) {
  std::vector<int> state_ids = {0, 1};
  const size_t nbead = state_ids.size();

  pot::bead_states my_beads;
  my_beads.set_individual_bead_states(state_ids);

  EXPECT_EQ(my_beads.state_id.size(), nbead);
  for (size_t i = 0; i < nbead; ++i) {
    EXPECT_EQ(my_beads.state_id[i], state_ids[i]);
  }
}

TEST(bead_states, avg_active_state) {
  const int nbead(2);
  const int init_state(1);

  pot::bead_states my_beads;
  my_beads.set_all_bead_states(init_state, nbead);

  EXPECT_EQ(my_beads.avg_active_state(), 1);
}

TEST(bead_states, sum_active_state) {
  const int nbead(2);
  const int init_state(1);

  pot::bead_states my_beads;
  my_beads.set_all_bead_states(init_state, nbead);

  EXPECT_EQ(my_beads.sum_active_state(), 2);
}
