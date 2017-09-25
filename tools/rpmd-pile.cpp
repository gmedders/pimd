#include <cmath>
#include <cassert>

#include <algorithm>

#include "rpmd-pile.h"

////////////////////////////////////////////////////////////////////////////////

namespace parts {

//----------------------------------------------------------------------------//

void rpmd_pile::init(size_t ndof, size_t nbead,
                     const double& kT, const double& dt,
                     const double* mass, const double* cartpos,
                     const double* cartvel, double tau)
{
    // Intialize the RPMD base class
    rpmd_base::init(ndof, nbead, kT, dt,
                    mass, cartpos, cartvel, tau);

    // Now do all the RPMD-NHC specific stuff

    c1 = arma::vec(nbead);
    c2 = arma::vec(nbead);

    //rand_gaussian = arma::mat(ndof, nbead, arma::zeros);

    m_tau0 = tau;

    for(size_t k = 0; k < nbead; ++k) {
        double gamma;
        if(k == 0)
            gamma = 1.0 / m_tau0;
        else
            gamma = 2.0 * m_omega_k(k);

        c1(k) = std::exp(-0.5*dt * gamma);
        c2(k) = std::sqrt(1.0 - c1(k)*c1(k));
    }

    //std::cerr << "<<< Thermostatting ( tau = " << m_tau0 << " ) >>>"<<std::endl;
}

//----------------------------------------------------------------------------//

void rpmd_pile::step(const double& dt)
{
    const double dt2 = dt/2;
    const double sqrt_beta_n = std::sqrt(m_beta_n);

    // 1. Advance thermostats dt2
    mom_c2n();
    for (size_t i = 0; i < ndofs(); ++i) {
        double fac = m_sqrt_mass(i)/sqrt_beta_n;
        for (size_t k = 0; k < nbeads(); ++k) {
            double rand_gaussian = randn(0, 1);
            m_mom_nmode(i,k) = m_mom_nmode(i,k)*c1(k)
                             + fac * c2(k) * rand_gaussian;
        }
    }
    mom_n2c();

    // 2. Do un-thermostatted RPMD evolution
    rpmd_base::step(dt);

    // 3. Advance thermostats final dt2 and recalc KE
    mom_c2n();
    for (size_t i = 0; i < ndofs(); ++i) {
        double fac = m_sqrt_mass(i)/sqrt_beta_n;
        for (size_t k = 0; k < nbeads(); ++k) {
            double rand_gaussian = randn(0, 1);
            m_mom_nmode(i,k) = m_mom_nmode(i,k)*c1(k)
                             + fac * c2(k) * rand_gaussian;
        }
    }
    mom_n2c();

    m_Ekin = 0.0;
    for (size_t n = 0; n < nbeads(); ++n) {
        for (size_t i = 0; i < ndofs(); ++i) {
            double mass = m_mass(i);
            const double Ekin2 = m_mom_cart(i,n)*m_mom_cart(i,n)/mass;
            m_Ekin += Ekin2;
        }
    }

    m_Ekin /= 2;
    m_temp_kT = m_Ekin*2.0/ndofs()/nbeads(); // not actual temperature, kT
}

//----------------------------------------------------------------------------//

double rpmd_pile::invariant() const
{
    return m_Ekin + m_Espring + m_Epot_sum;
}

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
