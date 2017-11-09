#include <cmath>
#include <cassert>

#include "rpmd-necklace.h"

namespace parts {

//----------------------------------------------------------------------------//

rpmd_necklace::rpmd_necklace()
: m_ndofs(0), m_nbeads(0)
{
}

//----------------------------------------------------------------------------//

void rpmd_necklace::setup(size_t ndim, size_t natom, size_t nbead, double beta,
                          double dt, double mass)
{
    //assert(ndof == 1);
    assert(nbead > 0);
    assert(nbead%2 == 0 || nbead == 1);

    m_ndim = ndim;
    m_natoms = natom;
    size_t ndof = ndim*natom;
    m_ndofs = ndof;
    m_nbeads = nbead;

    m_dt = dt;

    //m_pos_cart = arma::mat(nbead, ndof);
    //m_pos_nmode = arma::mat(nbead, ndof);

    //m_mom_cart = arma::mat(nbead, ndof);
    //m_mom_nmode = arma::mat(nbead, ndof);

    //m_frc_cart = arma::mat(nbead, ndof);

    //m_cart_to_nm = arma::mat(nbead, nbead);
    //m_freerp_propagator = arma::cube(2, 2, nbead);

    //m_omega_k = arma::vec(nbead);

    m_pos_cart = arma::mat(ndof, nbead);
    m_pos_nmode = arma::mat(ndof, nbead);

    m_mom_cart = arma::mat(ndof, nbead);
    m_mom_nmode = arma::mat(ndof, nbead);

    m_frc_cart = arma::mat(ndof, nbead);

    m_cart_to_nm = arma::mat(nbead, nbead);
    m_freerp_propagator = arma::cube(2, 2, nbead);

    m_omega_k = arma::vec(nbead);

    m_mass = arma::vec(ndof);
    m_sqrt_mass = arma::vec(ndof);

    // Catesian <-> normal mode transformation matrix
    for(size_t k = 0; k < nbead; ++k) {
        if(k == 0) {
            for(size_t j = 0; j < nbead; ++j)
                m_cart_to_nm(j, k) = std::sqrt(1.0/nbead);

        } else if(k <= nbead/2 - 1) { // 1 <= k <= nbead/2 - 1
            for(size_t j = 0; j < nbead; ++j)
                m_cart_to_nm(j, k)
                    = std::sqrt(2.0/nbead) * std::cos(2.0*M_PI*j*k/nbead);

        } else if(k == nbead/2) {
            for(size_t j = 0; j < nbead; ++j)
                m_cart_to_nm(j, k)
                    = std::sqrt(1.0/nbead) * std::pow(-1.0, j);

        } else { // nbead/2 + 1 <= k <= nbead - 1
            for(size_t j = 0; j < nbead; ++j)
                m_cart_to_nm(j, k)
                    = std::sqrt(2.0/nbead) * std::sin(2.0*M_PI*j*k/nbead);

        }
    }

    // omega_k

    m_beta_n = beta/nbead;
    m_omega_n = 1.0/(m_beta_n * hbar);

    for (size_t k = 0; k < nbead; ++k)
        m_omega_k(k) = 2.0 * m_omega_n * std::sin(M_PI * k / nbead);

    // Propagator for the free ring polymer hamiltonian
    { // k = 0
        size_t k = 0;
        double mw = mass * m_omega_k(k);

        m_freerp_propagator(0,0,k) = 1.0;
        m_freerp_propagator(1,0,k) = 0.0;
        // lim(x->0) [ 1/(m*x) * sin(x * dt) ] = dt/m
        m_freerp_propagator(0,1,k) = m_dt / mass;
        m_freerp_propagator(1,1,k) = 1.0;
    }

    for(size_t k = 1; k < nbead; ++k) {
        double mw = mass * m_omega_k(k);
        double c = std::cos(m_omega_k(k) * m_dt);
        double s = std::sin(m_omega_k(k) * m_dt);

        m_freerp_propagator(0,0,k) = c;
        m_freerp_propagator(1,0,k) = -mw * s;
        m_freerp_propagator(0,1,k) = 1.0/mw * s;
        m_freerp_propagator(1,1,k) = c;
    }
}

//----------------------------------------------------------------------------//

void rpmd_necklace::pos_c2n()
{
    m_pos_nmode = m_pos_cart * m_cart_to_nm;
}

//----------------------------------------------------------------------------//

void rpmd_necklace::pos_n2c()
{
    m_pos_cart = m_pos_nmode * m_cart_to_nm.t();
}

//----------------------------------------------------------------------------//

void rpmd_necklace::mom_c2n()
{
    m_mom_nmode = m_mom_cart * m_cart_to_nm;
}

//----------------------------------------------------------------------------//

void rpmd_necklace::mom_n2c()
{
    m_mom_cart = m_mom_nmode * m_cart_to_nm.t();
}

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
