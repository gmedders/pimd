#include "single-state.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

double single_state::force(size_t ndim, size_t natom, size_t nbead,
                           const double* crd, double* f)
{
    assert(nbead > 0);
    assert(init_pot == true);

    size_t ndof = ndim*natom;
    assert(ndof > 0);

    //std::cerr << crd[0] << std::endl;

    check_allocation(nbead, state_id);

    double E(0);
    for(size_t n = 0; n < nbead; ++n) {
        for(size_t m = 0; m < ndof; ++m) {
            E += VAA(crd + ndof*n + m, f + ndof*n + m);

            // FIXME Not implemented
            //double EBath = bath_force(crd + ndof*n + m, f + ndof*n + m);
            //E += EBath;
        }
    }

    return E;
}

//----------------------------------------------------------------------------//

void single_state::print_state_params()
{
    std::cout << "# single_state_dynamics" << std::endl;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
