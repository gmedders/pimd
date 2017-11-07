#include "anharmonic.h"

#include <cassert>

// equation 28 of dx.doi.org/10.1063/1.1777575

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

void anharmonic::assert_ndim(int ndim)
{
    assert(ndim == 1);
}

//----------------------------------------------------------------------------//

double anharmonic::VAA(const double* crd, double* f)
{
    double x = crd[0];
    double x2 = x*x;
    double x3 = x*x2;
    double x4 = x*x3;

    double e    =       a * x2 +       b * x3 +       c * x4;
    double dedx = 2.0 * a * x  + 3.0 * b * x2 + 4.0 * c * x3;

    f[0] = -dedx;

    return e;
}

//----------------------------------------------------------------------------//

void anharmonic::set_params(double* params)
{
    a = params[0];
    b = params[1];
    c = params[2];

    m = params[2];

    init_pot = true;
}

//----------------------------------------------------------------------------//

void anharmonic::print_params()
{
    std::cout << "# pot::anharmonic" << std::endl;
    std::cout << "# a = " << a << std::endl;
    std::cout << "# b = " << b << std::endl;
    std::cout << "# c = " << c << std::endl;
    std::cout << "# m = " << m << std::endl;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
