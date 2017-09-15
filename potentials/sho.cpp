#include "sho.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {


void sho::set_active_state(const int)
{
    active_state = -1;
}
void sho::set_hopping_params(double*){};

//----------------------------------------------------------------------------//

double sho::force(const size_t ndof, const size_t nbead,
                  const double* crd, double* f)
{
    assert(ndof == 1);
    assert(nbead > 0);
    assert(init == true);

    double energy(0);
    for(size_t n = 0; n < nbead; ++n)
        energy += VAA(crd + ndof*n, f + ndof*n);

    return energy;
}

//----------------------------------------------------------------------------//

double sho::VAA(const double* crd, double* f)
{
    double x = crd[0];

    double e = a * x * x;
    double dedx = 2.0 * a * x;
    f[0] = -dedx;

    return e;
}

//----------------------------------------------------------------------------//

void sho::set_params(double* params)
{
    w = params[0];
    m = params[1];

    // a = 1/2 m w^2
    a = 0.5 * m * w * w;

    init = true;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
