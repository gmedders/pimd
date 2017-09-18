#include "double-well.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

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

    double e = a * xg * xg + dG;
    double dedx = 2.0 * a * xg;
    f[0] = -dedx;

    return e;
}

//----------------------------------------------------------------------------//

// These parameters define the diabatic potentials
void double_well::set_params(double* params)
{
    w = params[0];
    m = params[1];
    g = params[2];
    dG = params[3];

    // a = 1/2 m w^2
    a = 0.5 * m * w * w;

    init_pot = true;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
