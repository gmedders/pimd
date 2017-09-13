#include "ah.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

double ah::VAA(const double* crd, double* f)
{
    double x = crd[0];

    double e = a * x * x;
    double dedx = 2.0 * a * x;
    f[0] = -dedx;

    return e;
}

//----------------------------------------------------------------------------//

double ah::VBB(const double* crd, double* f)
{
    double x = crd[0];

    double e = a*x*x + b*x + Ed;
    double dedx = 2.0*a*x + b;
    f[0] = -dedx;

    return e;
}

//----------------------------------------------------------------------------//

// These parameters define the diabatic potentials
void ah::set_params(double* params)
{
    w = params[0];
    m = params[1];
    g = params[2];
    double Ed_bar = params[3];

    // a = 1/2 m w^2
    a = 0.5 * m * w * w;

    // b = sqrt(2mw) * g
    b = std::sqrt(2.0 * m * w) * g;

    //Ed = Ed_bar + g^2/w
    double Er = g*g/w;
    Ed = Ed_bar + Er;

    init_pot = true;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
