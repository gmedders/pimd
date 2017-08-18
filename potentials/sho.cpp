#include "sho.h"

#include <cassert>

namespace pot {

double sho::operator()(size_t natom, const double* x, double* f) const
{
    assert(natom == 1);

    double e = a * x[0]*x[0];
    double dedx = 2.0 * a * x[0];
    f[0] = -dedx;

    return e;
}

} // namespace pot
