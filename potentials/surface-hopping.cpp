#include "surface-hopping.h"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////

namespace pot {

//----------------------------------------------------------------------------//

// Both determines whether the vector is the correct size
// and reallocates if necessary
void surface_hopping::check_allocation(size_t n_elements, arma::ivec& myvec)
{
    if(myvec.n_elem != n_elements)
        myvec.set_size(n_elements);
}

//----------------------------------------------------------------------------//

double surface_hopping::force(size_t ndim, size_t nparticles, size_t nbead,
                              const double* crd, double* f)
{
    assert(ndim == 1);
    assert(nparticles > 0);
    assert(nbead > 0);
    assert(init_pot == true);
    assert(init_hop == true);

    size_t ndof = ndim*nparticles;

    //std::cerr << crd[0] << std::endl;

    check_allocation(nbead, state_id);

    double Eactive(0);
    double fAA[ndof];
    double fBB[ndof];
    for(size_t n = 0; n < nbead; ++n) {
        double EAA = VAA(crd + ndof*n, fAA);
        double EBB = VBB(crd + ndof*n, fBB);

        if(state_id[n] == 0) {
            std::copy(fAA, fAA + ndof, f + ndof*n);
            Eactive += EAA;
        } else {
            std::copy(fBB, fBB + ndof, f + ndof*n);
            Eactive += EBB;
        }

        // CAREFUL!
        // Note that we're adding the bath force and energy directly into the
        // active diabatic force/energy
        double EBath = bath_force(crd + ndof*n, f + ndof*n);
        Eactive += EBath;

        // Determine if the active state of this bead should be changed
        double dE = EBB - EAA;
        double random_number = ((double) rand() / (RAND_MAX));
        assert(random_number >= 0.0 && random_number <= 1.0);

        if(hop_probability(dE, state_id[n]) > random_number){
            // Hop successful, change state id of this bead
            if(state_id[n] == 0)
                state_id[n] = 1;
            else
                state_id[n] = 0;
        }
    }

    return Eactive;
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

double surface_hopping::hop_probability(const double dE, int active_state)
{
    double fE = fermi_function(dE);

    if(active_state == 0)
        return Gamma * dt * fE;
    else
        return Gamma * dt * (1.0 - fE);
}

//----------------------------------------------------------------------------//

void surface_hopping::set_all_bead_states(int initial_state_id, int nbead)
{
    assert(initial_state_id == 0 || initial_state_id == 1);

    m_nbead = nbead;
    // Set for all beads
    check_allocation(m_nbead, state_id);
    for(int n = 0; n < m_nbead; ++n){
        state_id[n] = initial_state_id;
    }
}

//----------------------------------------------------------------------------//

void surface_hopping::set_individual_bead_states(
                                            std::vector<int>& initial_state_ids)
{
    m_nbead = initial_state_ids.size();

    // Set for all beads
    check_allocation(m_nbead, state_id);
    for(int n = 0; n < m_nbead; ++n){
        int id = initial_state_ids[n];
        assert(id == 0 || id == 1);
        state_id[n] = id;
    }
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

double surface_hopping::avg_active_state()
{
    return sum_active_state()/((double)m_nbead);
}

//----------------------------------------------------------------------------//

double surface_hopping::sum_active_state()
{
    double sum(0);
    for(int i = 0; i < m_nbead; ++i)
        sum += state_id[i];
    return sum;
}

//----------------------------------------------------------------------------//

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
