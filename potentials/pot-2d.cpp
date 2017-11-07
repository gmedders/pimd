#include "pot-2d.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

void pot_2d::assert_ndim(int ndim)
{
    assert(ndim == 2);
}

//----------------------------------------------------------------------------//

double pot_2d::VAA(const double* crd, double* f)
{
    double x = crd[0];
    double z = crd[1];

    double xg = (x - x0);
    double zg = (z - z0);
    double a = 0.5 * mx * w*w;

    double exp_z = std::exp(C0*zg);

    double R = a * xg * xg;
    double S = A0*exp_z*(exp_z - 2);

    double e = R + S;

    double dedx = 2.0 * a * xg;
    double dedz = 2.0*A0*C0 * exp_z * (exp_z - 1);

    f[0] += -dedx;
    f[1] += -dedz;

    return e;
}

//----------------------------------------------------------------------------//

double pot_2d::VBB(const double* crd, double* f)
{
    double x = crd[0];
    double z = crd[1];

    double xg = (x - x1);
    double a = 0.5 * mx * w*w;

    double R = a * xg * xg;
    double S = 1.0/(4*z) + 1.0/std::pow(z-C1, 6) + B1;

    double e = R + S;

    double dedx = 2.0 * a * xg;
    double dedz = -1.0/(4*z*z) - 6.0/std::pow(z-C1, 7);

    f[0] += -dedx;
    f[1] += -dedz;

    return e;
}

//----------------------------------------------------------------------------//

void pot_2d::set_params(double* params)
{
    A0 = params[0];
    B1 = params[1];
    x0 = params[2];
    x1 = params[3];
    C0 = params[4];
    C1 = params[5];
    w  = params[6];
    mx = params[7];
    z0 = params[8];

    init_pot = true;
}

//----------------------------------------------------------------------------//

void pot_2d::print_params()
{
    std::cout << "# pot::pot_2d" << std::endl;
    std::cout << "# A0     = " << A0 << std::endl;
    std::cout << "# B1     = " << B1 << std::endl;
    std::cout << "# x0     = " << x0 << std::endl;
    std::cout << "# x1     = " << x1 << std::endl;
    std::cout << "# C0     = " << C0 << std::endl;
    std::cout << "# C1     = " << C1 << std::endl;
    std::cout << "# w      = " << w  << std::endl;
    std::cout << "# mx     = " << mx << std::endl;
    std::cout << "# z0     = " << z0 << std::endl;

    print_state_params();
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
