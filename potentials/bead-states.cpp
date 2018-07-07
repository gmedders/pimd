#include "bead-states.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

// Both determines whether the vector is the correct size
// and reallocates if necessary
void bead_states::check_allocation(size_t n_elements, arma::ivec &myvec) {
  if (myvec.n_elem != n_elements)
    myvec.set_size(n_elements);
}

//----------------------------------------------------------------------------//

void bead_states::set_all_bead_states(int initial_state_id, int nbead) {
  assert(initial_state_id == 0 || initial_state_id == 1);

  m_nbead = nbead;
  // Set for all beads
  check_allocation(m_nbead, state_id);
  for (int n = 0; n < m_nbead; ++n) {
    state_id[n] = initial_state_id;
  }
}

//----------------------------------------------------------------------------//

void bead_states::set_individual_bead_states(
    std::vector<int> &initial_state_ids) {
  m_nbead = initial_state_ids.size();

  // Set for all beads
  check_allocation(m_nbead, state_id);
  for (int n = 0; n < m_nbead; ++n) {
    int id = initial_state_ids[n];
    assert(id == 0 || id == 1);
    state_id[n] = id;
  }
}

//----------------------------------------------------------------------------//

double bead_states::avg_active_state() {
  return sum_active_state() / ((double)m_nbead);
}

//----------------------------------------------------------------------------//

double bead_states::sum_active_state() {
  double sum(0);
  for (int i = 0; i < m_nbead; ++i)
    sum += state_id[i];
  return sum;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
