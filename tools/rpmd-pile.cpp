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
                     const double* cartvel, double gamma_centroid)
{
    // Intialize the RPMD base class
    rpmd_base::init(ndof, nbead, kT, dt,
                    mass, cartpos, cartvel, gamma_centroid);

    // Now do all the RPMD-NHC specific stuff

    c1 = arma::vec(nbead);
    c2 = arma::vec(nbead);

    //rand_gaussian = arma::mat(ndof, nbead, arma::zeros);

    for(size_t k = 0; k < nbead; ++k) {
        double gamma;
        if(k == 0)
            gamma = gamma_centroid;
        else
            gamma = 2.0 * m_omega_k(k);

        c1(k) = std::exp(-0.5*dt * gamma);
        c2(k) = std::sqrt(1.0 - c1(k)*c1(k));
    }

    //std::cerr << "<<< Thermostatting ( tau = " << 1.0/gamma_centroid 
    //          << " ) >>>"<<std::endl;
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

    double m_Ekin_centroid = 0.0;
    double m_Ekin_higherNM = 0.0;
    for (size_t n = 0; n < nbeads(); ++n) {
        for (size_t i = 0; i < ndofs(); ++i) {
            double mass = m_mass(i);
            if(n == 0)
                m_Ekin_centroid += m_mom_nmode(i,n)*m_mom_nmode(i,n)/mass;
            else
                m_Ekin_higherNM += m_mom_nmode(i,n)*m_mom_nmode(i,n)/mass;
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
    m_temp_kT_centroid = m_Ekin_centroid*2.0/ndofs();
    m_temp_kT_higherNM = m_Ekin_higherNM*2.0/ndofs()/(nbeads()-1);
}

//----------------------------------------------------------------------------//

double rpmd_pile::invariant() const
{
    return m_Ekin + m_Espring + m_Epot_sum;
}

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
