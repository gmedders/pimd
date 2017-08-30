#include <cmath>
#include <cassert>

#include <algorithm>

#include "rpmd-base.h"

////////////////////////////////////////////////////////////////////////////////

namespace parts {

//----------------------------------------------------------------------------//

rpmd_base::rpmd_base()
: necklace()
{
    m_phys_mass = 0;
}

//----------------------------------------------------------------------------//

rpmd_base::~rpmd_base()
{
    delete[] m_phys_mass;
}

//----------------------------------------------------------------------------//

void rpmd_base::init(size_t ndof, size_t nbead, const double& kT,
                     const double* mass, const double* cartpos,
                     const double* cartvel)
{
    assert(ndof > 0 && nbead > 0);
    assert(nbead%2 == 0 || nbead == 1);
    assert(kT > 0.0 && mass != 0 && cartpos != 0 && cartvel != 0);

    necklace::setup(ndof, nbead);

    delete[] m_phys_mass;

    m_kT = kT;
    m_omega_M = std::sqrt(double(nbead))*kT/hbar;

    // fictitious masses

    m_phys_mass = new double[ndof*nbead];

    for (size_t b = 0; b < nbead; ++b) {
        //const double factor = (b == 0 ? 1.0 : lambda(b));
        for (size_t i = 0; i < ndof; ++i)
            m_phys_mass[i + b*ndof] = mass[i];
    }

    // initialize cartesian positions and velocities

    std::copy(cartpos, cartpos + nbead*ndof, m_pos_cart);
    std::copy(cartvel, cartvel + nbead*ndof, m_vel_cart);

    // calculate the initial kinetic energy

    m_Ekin = 0.0;
    for (size_t b = 0; b < nbead; ++b)
        for (size_t i = 0; i < ndof; ++i) {
            const size_t j = i + b*ndof;

            m_Ekin += m_phys_mass[j]*m_vel_cart[j]*m_vel_cart[j];
        }

    m_Ekin /= 2;

    pimd_force();
}

//----------------------------------------------------------------------------//

void rpmd_base::pimd_force()
{
    const size_t n_total = ndofs()*nbeads();

    // zero out m_frc_cart
    std::fill(m_frc_cart, m_frc_cart + n_total, 0.0);

    // compute forces for each bead
    m_Epot_sum = 0.0;
    for (size_t b = 0; b < nbeads(); ++b)
        m_Epot_sum += force(m_pos_cart + b*ndofs(), m_frc_cart + b*ndofs());

    m_Epot_sum /= nbeads();

    const double omega2 = m_omega_M*m_omega_M;

    m_Espring = 0.0;

    // add the harmonic part

    for (size_t b = 1; b < nbeads(); ++b) {
        for (size_t i = 0; i < ndofs(); ++i) {
            const size_t j = b*ndofs() + i;
            //FIXME
            const double tmp = m_phys_mass[j]*omega2*m_pos_cart[j];
            m_frc_cart[j] -= tmp;
            m_Espring += tmp*m_pos_cart[j];
        }
    }

    m_Espring /= 2;
}

//----------------------------------------------------------------------------//

void rpmd_base::step(const double& dt)
{

    // 1. advance thermostats, velocities by dt/2, nmode position on dt

    const double dt2 = dt/2;

    for (size_t b = 0; b < nbeads(); ++b) {
        for (size_t i = 0; i < ndofs(); ++i) {
            const size_t j = b*ndofs() + i;
            const double mass = m_phys_mass[j];
            const double Ekin2 = mass*m_vel_cart[j]*m_vel_cart[j];

            // FIXME
            m_vel_cart[j] += dt2*m_frc_cart[j]/mass;
            m_pos_cart[j] += dt*m_vel_cart[j];
        }
    }

    // 2. transform position to cartesian, compute forces

    pimd_force(); // computes normal mode forces

    // 3. advance velocities and thermostats by dt/2

    m_Ekin = 0.0;

    for (size_t b = 0; b < nbeads(); ++b) {
        for (size_t i = 0; i < ndofs(); ++i) {
            const size_t j = b*ndofs() + i;
            const double mass = m_phys_mass[j];

            // FIXME
            m_vel_cart[j] += dt2*m_frc_cart[j]/mass;
            m_Ekin += mass*m_vel_cart[j]*m_vel_cart[j];
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
