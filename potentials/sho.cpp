#include "sho.h"

#include <cassert>

namespace pot {

double sho::operator()(size_t natom, const double* x, double* f) const
{
    assert(natom == 1);
    assert(init == true);

    double e = a * x[0]*x[0];
    double dedx = 2.0 * a * x[0];
    f[0] = -dedx;

    return e;
}

void sho::set_params(double* params)
{
    w = params[0];
    m = params[1];

    // a = 1/2 m w^2
    a = 0.5 * m * w * w;

    init = true;
}

} // namespace pot
