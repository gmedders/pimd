#include "sho.h"

namespace pot {

double sho::operator()(const double x, double& f) const
{

    double e = a * x*x;
    double dedx = 2.0 * a * x;
    f = -dedx;

    return e;

}

} // namespace pot
