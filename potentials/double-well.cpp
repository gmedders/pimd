#include "double-well.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

double double_well::operator()(size_t ndof, size_t nbead,
                               const double* crd, double* f) const
{
    assert(ndof == 1);
    assert(nbead > 0);
    assert(init == true);

    double energy(0);
    if(active_state == 0)
        for(size_t n = 0; n < nbead, ++n)
            energy += VAA(crd + ndof*n, f + ndof*n);
    else
        for(size_t n = 0; n < nbead, ++n)
            energy += VBB(crd + ndof*n, f + ndof*n);

    double p = m_potential.hop_probability(nbeads(), ndofs(),
                                           crd);

    if(p > m_prng.random_double()){
        m_potential.hop();
    }

    return energy;
}

//----------------------------------------------------------------------------//

double double_well::VAA(const double* crd, double* f)
{
    double x = crd[0];

    double e = a * x * x;
    double dedx = 2.0 * a * x;
    f[0] = -dedx;

    return e;
}

//----------------------------------------------------------------------------//

double double_well::VBB(const double* crd, double* f)
{
    double xg = crd[0] - g;

    double e = a * xg * xg;
    double dedx = 2.0 * a * xg;
    f[0] = -dedx;

    return e;
}

//----------------------------------------------------------------------------//

void double_well::set_params(double* params)
{
    w = params[0];
    m = params[1];
    g = params[2];

    // a = 1/2 m w^2
    a = 0.5 * m * w * w;

    init = true;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
