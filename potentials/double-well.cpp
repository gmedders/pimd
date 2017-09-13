#include "double-well.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

// Both determines whether the vector is the correct size
// and reallocates if necessary
void double_well::check_allocation(size_t n_elements, arma::vec& myvec)
{
    if(myvec.n_elem != n_elements)
        myvec.set_size(n_elements);
}

//----------------------------------------------------------------------------//

double double_well::force(size_t ndof, size_t nbead,
                          const double* crd, double* f)
{
    assert(ndof == 1);
    assert(nbead > 0);
    assert(init_pot == true);
    assert(init_hop == true);

    check_allocation(ndof*nbead, fAA);
    check_allocation(ndof*nbead, fBB);

    double EAA(0);
    double EBB(0);
    if(active_state == 0)
        for(size_t n = 0; n < nbead; ++n)
            EAA += VAA(crd + ndof*n, fAA.memptr() + ndof*n);
    else
        for(size_t n = 0; n < nbead; ++n)
            EBB += VBB(crd + ndof*n, fBB.memptr() + ndof*n);

    // Determine if the active state should be changed
    double dE = (EAA - EBB)/nbead;
    if(hop_probability(dE) > m_prng.random_double())
        hop();

    // Set the force and energy to the active state values
    if(active_state == 0){
        std::copy(fAA.memptr(), fAA.memptr() + ndof*nbead, f);
        return EAA;
    } else {
        std::copy(fBB.memptr(), fBB.memptr() + ndof*nbead, f);
        return EBB;
    }

}

//----------------------------------------------------------------------------//

double double_well::fermi_function(const double dE)
{
    double exp_arg = beta * dE;
    if(exp_arg > 100.0)
        return 0.0;
    else
        return 1.0 / (1.0 + std::exp(exp_arg));
}

//----------------------------------------------------------------------------//

double double_well::hop_probability(const double dE)
{
    double fE = fermi_function(dE);

    if(active_state == 0)
        return Gamma * dt * fE;
    else
        return Gamma * dt * (1.0 - fE);
}

//----------------------------------------------------------------------------//

void double_well::hop()
{
  if(active_state == 0)
      active_state = 1;
  else
      active_state = 0;

  ++nhops;
}

//----------------------------------------------------------------------------//

double double_well::VAA(const double* crd, double* f)
{
    double x = crd[0];

    double e = a * x * x;
    double dedx = 2.0 * a * x;
    f[0] = -dedx;

    return e;
}

//----------------------------------------------------------------------------//

double double_well::VBB(const double* crd, double* f)
{
    double xg = crd[0] - g;

    double e = a * xg * xg;
    double dedx = 2.0 * a * xg;
    f[0] = -dedx;

    return e;
}

//----------------------------------------------------------------------------//

// These parameters define the diabatic potentials
void double_well::set_params(double* params)
{
    w = params[0];
    m = params[1];
    g = params[2];

    // a = 1/2 m w^2
    a = 0.5 * m * w * w;

    init_pot = true;
}

//----------------------------------------------------------------------------//

// These parameters define the hopping probability
void double_well::set_hopping_params(double* params)
{
    Gamma = params[0];
    dt = params[1];
    beta = params[2];

    init_hop = true;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
