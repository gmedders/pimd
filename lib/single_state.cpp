#include "single_state.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

double single_state::force(size_t ndim, size_t natom, size_t nbead,
                           const double *crd, double *f) {
  assert(nbead > 0);
  assert(init_pot == true);

  size_t ndof = ndim * natom;
  assert(ndof > 0);

  // std::function<double(const double*, double*)> V_1body = VAA;

  // Ugly fix to allocate the state_id array the first time that
  // force() is called
  static bool first_call = true;
  if (first_call) {
    check_allocation(nbead, state_id);
    first_call = false;
  }

  // Now calculate the energy and forces
  double E(0);
  std::fill(f, f + ndof * nbead, 0.0);
  for (size_t n = 0; n < nbead; ++n) {
    const size_t ind_bead = n * ndof;

    // Evaluate the influence of the underlying potential on the particle
    // ... like a 1-body term
    for (size_t i = 0; i < natom; ++i) {
      const size_t ind_i = ind_bead + i * ndim;

      E += VAA(crd + ind_i, f + ind_i);

      // FIXME Not implemented
      // double EBath = bath_force(crd + ndof*n + m, f + ndof*n + m);
      // E += EBath;
    }
  }

  return E;
}

//----------------------------------------------------------------------------//

void single_state::print_state_params() {
  std::cout << "# single_state_dynamics" << std::endl;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
