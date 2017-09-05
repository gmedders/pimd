#include <cmath>
#include <cassert>

#include <armadillo>

#include "rpmd-necklace.h"

namespace parts {

//----------------------------------------------------------------------------//

rpmd_necklace::rpmd_necklace()
: m_ndofs(0), m_nbeads(0)
{
}

//----------------------------------------------------------------------------//

rpmd_necklace::~rpmd_necklace()
{
    setup(0, 0);
}

//----------------------------------------------------------------------------//

void rpmd_necklace::setup(size_t ndof, size_t nbead,
                          double dt = 1.0, double mass = 2000.0)
{
    if (ndof == m_ndofs && nbead == m_nbeads)
        return;

    if (m_nbeads > 0 || m_ndofs > 0) {
        assert(m_nbeads > 0 && m_ndofs > 0);

        delete[] m_pos_cart;
        delete[] m_pos_nmode;

        delete[] m_vel_cart;
        delete[] m_vel_nmode;

        delete[] m_frc_cart;

        delete[] m_cart_to_nm;
        delete[] m_freerp_propagator;

        delete[] m_lambda;
    }

    if (ndof > 0 || nbead > 0) {
        assert(ndof > 0 && nbead > 0);
        assert(nbead%2 == 0 || nbead == 1);

        m_ndofs = ndof;
        m_nbeads = nbead;

        m_pos_cart = new double[nbead*ndof];
        m_pos_nmode = new double[nbead*ndof];

        m_vel_cart = new double[nbead*ndof];
        m_vel_nmode = new double[nbead*ndof];

        m_frc_cart = new double[nbead*ndof];

        m_cart_to_nm = new double[nbead*nbead];
        m_freerp_propagator = new double[2*2*nbead];

        m_lambda = new double[nbead];

        // Catesian <-> normal mode transformation matrix
        for(size_t k = 0; k < nbead; ++k) {
            if(k == 0) {
                for(size_t j = 0; j < nbead; ++j)
                    m_cart_to_nm[k*nbead + j] = std::sqrt(1.0/nbead);

            } else if(k <= nbead/2 - 1) { // 1 <= k <= nbead/2 - 1
                for(size_t j = 0; j < nbead; ++j)
                    m_cart_to_nm[k*nbead + j] 
                        = std::sqrt(2.0/nbead) * std::cos(2.0*M_PI*j*k/nbead);

            } else if(k == nbead/2) {
                for(size_t j = 0; j < nbead; ++j)
                    m_cart_to_nm[k*nbead + j]
                        = std::sqrt(1.0/nbead) * std::pow(-1.0, j);

            } else { // nbead/2 + 1 <= k <= nbead - 1
                for(size_t j = 0; j < nbead; ++j)
                    m_cart_to_nm[k*nbead + j] 
                        = std::sqrt(2.0/nbead) * std::sin(2.0*M_PI*j*k/nbead);

            }
        }

        // omega_k

        m_beta_n = beta/nbead;
        m_omega_n = 1.0/(m_beta_n * hbar);

        for (size_t k = 0; k < nbead; ++k)
            m_omega_k[k] = 2.0 * m_omega_n * std::sin(M_PI * k / nbead);

        // Propagator for the free ring polymer hamiltonian
        { // k = 0
            double mw = mass * m_omega_k[k];
            m_freerp_propagator[0] = 1.0;
            m_freerp_propagator[1] = 0.0;
            m_freerp_propagator[2] = 0.0;
            m_freerp_propagator[3] = 1.0;
        }

        for(size_t k = 1; k < nbead; ++k) {
            double mw = mass * m_omega_k[k];
            double c = std::cos(m_omega_k[k] * dt);
            double s = std::sin(m_omega_k[k] * dt);
            m_freerp_propagator[4*k + 0] = c;
            m_freerp_propagator[4*k + 1] = -mw * s;
            m_freerp_propagator[4*k + 2] = 1.0/mw * s;
            m_freerp_propagator[4*k + 3] = c;
        }

    }
}

//----------------------------------------------------------------------------//

void rpmd_necklace::pos_c2n()
{
    assert(m_nbeads > 0 && m_ndofs > 0);

    fftw_execute_r2r(m_plan_cart2nmode, m_pos_cart, m_pos_nmode);

    const double factor = 1.0/m_nbeads;
    for (size_t n = 0; n < m_nbeads*m_ndofs; ++n)
        m_pos_nmode[n] *= factor;
}

//----------------------------------------------------------------------------//

void rpmd_necklace::pos_n2c()
{
    assert(m_nbeads > 0 && m_ndofs > 0);

    fftw_execute_r2r(m_plan_nmode2cart, m_pos_nmode, m_pos_cart);
}

//----------------------------------------------------------------------------//

void rpmd_necklace::vel_n2c()
{
    assert(m_nbeads > 0 && m_ndofs > 0);

    fftw_execute_r2r(m_plan_nmode2cart, m_vel_nmode, m_vel_cart);
}

//----------------------------------------------------------------------------//

void rpmd_necklace::vel_c2n()
{
    assert(m_nbeads > 0 && m_ndofs > 0);

    fftw_execute_r2r(m_plan_cart2nmode, m_vel_cart, m_vel_nmode);

    const double factor = 1.0/m_nbeads;
    for (size_t n = 0; n < m_nbeads*m_ndofs; ++n)
        m_vel_nmode[n] *= factor;
}

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
