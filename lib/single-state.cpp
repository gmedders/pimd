#include "single-state.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

// Both determines whether the vector is the correct size
// and reallocates if necessary
void single_state::check_allocation(size_t n_elements, arma::ivec &myvec) {
  if (myvec.n_elem != n_elements)
    myvec.set_size(n_elements);
}

//----------------------------------------------------------------------------//

double single_state::force(size_t ndim, size_t natom, size_t nbead,
                           const double *crd, double *f) {
  assert(nbead > 0);
  assert(init_pot == true);

  size_t ndof = ndim * natom;
  assert(ndof > 0);

  // std::cerr << crd[0] << std::endl;

  check_allocation(nbead, state_id);

  double E(0);
  for (size_t n = 0; n < nbead; ++n) {
    for (size_t m = 0; m < ndof; ++m) {
      E += VAA(crd + ndof * n + m, f + ndof * n + m);
    }
  }

  return E;
}

//----------------------------------------------------------------------------//

void single_state::set_all_bead_states(int initial_state_id, int nbead) {
  assert(initial_state_id == 0 || initial_state_id == 1);

  m_nbead = nbead;
  // Set for all beads
  check_allocation(m_nbead, state_id);
  for (int n = 0; n < m_nbead; ++n) {
    state_id[n] = initial_state_id;
  }
}

//----------------------------------------------------------------------------//

void single_state::set_individual_bead_states(
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

void single_state::print_state_params() {
  std::cout << "# single_state_dynamics" << std::endl;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
