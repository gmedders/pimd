#include "anharmonic.h"

#include <cassert>

// equation 28 of dx.doi.org/10.1063/1.1777575

namespace pot {

double anharmonic::operator()(size_t natom, const double* crd, double* f) const
{
    assert(natom == 1);
    assert(init == true);

    double x = crd[0];
    double x2 = x*x;
    double x3 = x*x2;
    double x4 = x*x3;

    double e    = a * 1.0/2.0 * x2 + b * 1.0/10.0 * x3 + c * 1.0/100.0 * x4;
    double dedx = a *           x  + b * 3.0/10.0 * x2 + c * 4.0/100.0 * x3;

    f[0] = -dedx;

    return e;
}

void anharmonic::set_params(double* params)
{
    a = params[0];
    b = params[1];
    c = params[2];

    init = true;
}

} // namespace pot
