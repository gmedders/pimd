#include "anharmonic.h"

#include <cassert>

// equation 28 of dx.doi.org/10.1063/1.1777575

////////////////////////////////////////////////////////////////////////////////

namespace pot {

void anharmonic::set_all_bead_states(const int)
{
    active_state = -1;
}
void anharmonic::set_hopping_params(double*){};

//----------------------------------------------------------------------------//

double anharmonic::force(const size_t ndof, const size_t nbead,
                         const double* crd, double* f)
{
  assert(ndof == 1);
  assert(nbead > 0);
  assert(init == true);

  double energy(0);
  for(size_t n = 0; n < nbead; ++n)
      energy += VAA(crd + ndof*n, f + ndof*n);

  return energy;
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

    init = true;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
