#include <cmath>
#include <cassert>

#include <algorithm>

#include "rpmd-base.h"

////////////////////////////////////////////////////////////////////////////////

namespace parts {

//----------------------------------------------------------------------------//

rpmd_base::rpmd_base()
: rpmd_necklace()
{
}

//----------------------------------------------------------------------------//

rpmd_base::~rpmd_base()
{
}

//----------------------------------------------------------------------------//

void rpmd_base::init(size_t ndof, size_t nbead,
                     const double& kT, const double& dt,
                     const double* mass, const double* cartpos,
                     const double* cartvel)
{
    assert(ndof > 0 && nbead > 0);
    assert(nbead%2 == 0 || nbead == 1);
    assert(kT > 0.0 && mass != 0 && cartpos != 0 && cartvel != 0);

    m_kT = kT;
    rpmd_necklace::setup(ndof, nbead, 1.0/kT, dt, mass[0]);

    // initialize cartesian positions and momenta
    for(size_t n = 0; n < nbeads(); ++n){
        for(size_t i = 0; i < ndofs(); ++i){
            const size_t j = n*ndofs() + i;

            m_pos_cart(i,n) = cartpos[j];
            m_mom_cart(i,n) = mass[i]*cartvel[j];
        }
    }

    for (size_t i = 0; i < ndofs(); ++i){
        m_mass(i) = mass[i];
        m_sqrt_mass(i) = std::sqrt(mass[i]);
    }

    // calculate the initial kinetic energy

    m_Ekin = 0.0;
    for (size_t n = 0; n < nbeads(); ++n)
        for (size_t i = 0; i < ndofs(); ++i) {
            m_Ekin += m_mom_cart(i,n)*m_mom_cart(i,n)/m_mass(i);
        }

    m_Ekin /= 2;

    pimd_force();
    spring_energy();
}

//----------------------------------------------------------------------------//

void rpmd_base::pimd_force()
{
    m_Epot_sum = force(ndofs(), nbeads(),
                       m_pos_cart.memptr(), m_frc_cart.memptr());
}

//----------------------------------------------------------------------------//

void rpmd_base::spring_energy()
{
    m_Espring = 0.0;
    if(nbeads() == 0)
        return;

    for (size_t n = 0; n < nbeads(); ++n){
        size_t bind_a(n);
        size_t bind_b;
        if(n == 0)
            bind_b = nbeads() - 1;
        else
            bind_b = n - 1;

        arma::vec q_a(m_pos_cart.colptr(bind_a), ndofs(), false, true);
        arma::vec q_b(m_pos_cart.colptr(bind_b), ndofs(), false, true);

        arma::vec diff = m_sqrt_mass*(q_a - q_b);

        m_Espring += arma::dot(diff, diff);
    }

    m_Espring *= 0.5 * m_omega_n * m_omega_n;
}

//----------------------------------------------------------------------------//

void rpmd_base::step(const double& dt)
{
    const double dt2 = dt/2;

    // Following equations 21-25 of dx.doi.org/10.1063/1.3489925
    // 1. Evolution of RP momenta under Hamiltonian V_{n}[q(t0)] by dt/2
    m_mom_cart += dt2*m_frc_cart;


    // 2. Transform positions & momenta from Cartesian to normal-mode
    pos_c2n();
    mom_c2n();


    // 3. Evolve RP coords/momenta by dt under free RP Hamiltonian H_{n}^{0}
    //    So, for each bead (in NM representation)...
    for(size_t b = 0; b < nbeads(); ++b){
        arma::mat evolve = m_freerp_propagator.slice(b); // (2,2)
        //evolve.print("> evolve:");
        arma::mat nm_pos_mom_t0(ndofs(),2);
        arma::mat nm_pos_mom_tdt(ndofs(),2);

        //m_mom_nmode.col(b).print("> m_mom_nmode.col(b):");
        //m_pos_nmode.col(b).print("> m_pos_nmode.col(b):");

        nm_pos_mom_t0.col(0) = m_mom_nmode.col(b);
        nm_pos_mom_t0.col(1) = m_pos_nmode.col(b);

        //nm_pos_mom_t0.print("> nm_pos_mom_t0:");

        //// (2,ndof) = (2,2) * (2,ndof)
        //nm_pos_mom_tdt = evolve * nm_pos_mom_t0;

        // (ndof,2) = (ndof,2) * (2,2)
        //nm_pos_mom_t0.print("nm_pos_mom_t0:");
        //evolve.print("evolve:");
        nm_pos_mom_tdt = nm_pos_mom_t0 * evolve;
        //nm_pos_mom_tdt.print("nm_pos_mom_tdt:");

        m_mom_nmode.col(b) = nm_pos_mom_tdt.col(0);
        m_pos_nmode.col(b) = nm_pos_mom_tdt.col(1);
    }


    // 4. Transform positions & momenta from normal-mode to Cartesian
    pos_n2c();
    mom_n2c();


    // 5. Final evolution of RP momenta under Hamiltonian V_{n}[q(t0+dt)]
    pimd_force();
    spring_energy();

    m_mom_cart += dt2*m_frc_cart;

    m_Ekin = 0.0;

    for (size_t n = 0; n < nbeads(); ++n) {
        for (size_t i = 0; i < ndofs(); ++i) {
            m_Ekin += m_mom_cart(i,n)*m_mom_cart(i,n)/m_mass(i);
        }
    }

    m_Ekin /= 2;
    m_temp_kT = m_Ekin*2.0/ndofs()/nbeads(); // not actual temperature, kT
}

//----------------------------------------------------------------------------//

double rpmd_base::invariant() const
{
    // Epot is in kcal/mol already
    return (m_Ekin + m_Espring)/engunit + m_Epot_sum;
}

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
