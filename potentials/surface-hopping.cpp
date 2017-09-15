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

    //std::cerr << crd[0] << std::endl;

    check_allocation(ndof*nbead, fAA);
    double EAA(0);
    for(size_t n = 0; n < nbead; ++n)
        EAA += VAA(crd + ndof*n, fAA.memptr() + ndof*n);

    check_allocation(ndof*nbead, fBB);
    double EBB(0);
    for(size_t n = 0; n < nbead; ++n)
        EBB += VBB(crd + ndof*n, fBB.memptr() + ndof*n);

    // Determine if the active state should be changed
    double dE = (EBB - EAA) / nbead;
    double random_number = ((double) rand() / (RAND_MAX));
    assert(random_number >= 0.0 && random_number <= 1.0);

    //std::cerr << EBB << ' ' << EAA << ' ' << hop_probability(dE)
    //          << ' ' << random_number << std::endl;

    // Set the force and energy to the active state values
    double Eactive;
    if(active_state == 0){
        std::copy(fAA.memptr(), fAA.memptr() + ndof*nbead, f);
        Eactive = EAA;
    } else {
        std::copy(fBB.memptr(), fBB.memptr() + ndof*nbead, f);
        Eactive = EBB;
    }

    if(hop_probability(dE) > random_number)
        hop();

    return Eactive;
}

//----------------------------------------------------------------------------//

double surface_hopping::fermi_function(const double dE)
{
    double exp_arg = beta * dE;
    //std::cout << " Beta = " << beta << ' ' << dE << ' ' << exp_arg << std::endl;
    //exit(0);
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

void surface_hopping::set_active_state(int state_id)
{
    assert(state_id == 0 || state_id == 1);
    active_state = state_id;
}

//----------------------------------------------------------------------------//

// These parameters define the hopping probability
void surface_hopping::set_hopping_params(double* params)
{
    Gamma = params[0];
    dt = params[1];
    beta = params[2];

    //prng.seed(19104);

    init_hop = true;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
