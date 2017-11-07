#include "sho.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

void sho::assert_ndim(int ndim)
{
    assert(ndim == 1);
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

    init_pot = true;
}

//----------------------------------------------------------------------------//

void sho::print_params()
{
    std::cout << "# pot::sho" << std::endl;
    std::cout << "# w     = " << w << std::endl;
    std::cout << "# m     = " << m << std::endl;
    print_state_params();
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
