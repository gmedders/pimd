#include "surface-hopping.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

// Both determines whether the vector is the correct size
// and reallocates if necessary
void surface_hopping::check_allocation(size_t n_elements, arma::vec& myvec)
{
    if(myvec.n_elem != n_elements)
        myvec.set_size(n_elements);
}

//----------------------------------------------------------------------------//

double surface_hopping::force(size_t ndof, size_t nbead,
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
    if(hop_probability(dE) > prng.random_double())
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

double surface_hopping::fermi_function(const double dE)
{
    double exp_arg = beta * dE;
    if(exp_arg > 100.0)
        return 0.0;
    else
        return 1.0 / (1.0 + std::exp(exp_arg));
}

//----------------------------------------------------------------------------//

double surface_hopping::hop_probability(const double dE)
{
    double fE = fermi_function(dE);

    if(active_state == 0)
        return Gamma * dt * fE;
    else
        return Gamma * dt * (1.0 - fE);
}

//----------------------------------------------------------------------------//

void surface_hopping::hop()
{
  if(active_state == 0)
      active_state = 1;
  else
      active_state = 0;

  ++nhops;
}

//----------------------------------------------------------------------------//

// These parameters define the hopping probability
void surface_hopping::set_hopping_params(double* params)
{
    Gamma = params[0];
    dt = params[1];
    beta = params[2];

    prng.seed(19104);

    init_hop = true;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
