#include "lj.h"

#include <cassert>
#include <iostream>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

double lj::V2B(const size_t ndim, const double* crd_a, const double* crd_b,
           double* f_a, double* f_b)
{
    static double A = 4.0*eps*std::pow(sig, 12);
    static double B = 4.0*eps*std::pow(sig, 6);

    double R(0.0);
    double dR[ndim];
    for (size_t i = 0; i < ndim; i++) {
        dR[i] = crd_b[i] - crd_a[i];
        R += dR[i]*dR[i];
    }

    R = std::sqrt(R);

    double R2 = R*R;
    double R6 = R2*R2*R2;
    double R12 = R6*R6;

    double e = A/R12 - B/R6;
    double dedr = -12*A/(R12*R) + 6*B/(R6/R);

    for (size_t i = 0; i < ndim; i++) {
        f_a[i] -= -dR[i]/R * dedr;
        f_b[i] -= +dR[i]/R * dedr;
    }

    return e;
}

// //----------------------------------------------------------------------------//
//
// void lj::set_params(double* params)
// {
//     eps = params[0];
//     sig = params[1];
//
//
//
//     init_pot = true;
// }
//
//----------------------------------------------------------------------------//

void lj::print_params()
{
    std::cout << "# pot::lj" << std::endl;
    std::cout << "# eps     = " << eps << std::endl;
    std::cout << "# sig     = " << sig << std::endl;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
