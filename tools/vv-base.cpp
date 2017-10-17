#include <cmath>
#include <cassert>

#include <algorithm>

#include "vv-base.h"

////////////////////////////////////////////////////////////////////////////////

namespace parts {

//----------------------------------------------------------------------------//

void vv_base::init(size_t ndof,
                   const double& dt,
                   const double* mass, const double* cartpos,
                   const double* cartvel)
{
    assert(ndof > 0);
    assert(mass != 0 && cartpos != 0 && cartvel != 0);

    assert(ndof == 1);

    m_ndofs = ndof;
    m_dt = dt;

    m_pos = arma::vec(ndof);
    m_mom = arma::vec(ndof);
    m_frc = arma::vec(ndof);

    m_mass = arma::vec(ndof);

    // initialize cartesian positions and momenta
    for(size_t i = 0; i < ndofs(); ++i){
        m_pos(i) = cartpos[i];
        m_mom(i) = mass[i]*cartvel[i];
    }

    for (size_t i = 0; i < ndofs(); ++i){
        m_mass(i) = mass[i];
    }

    // calculate the initial kinetic energy

    m_Ekin = 0.0;
    for (size_t i = 0; i < ndofs(); ++i) {
        m_Ekin += m_mom(i)*m_mom(i)/m_mass(i);
    }

    m_Ekin /= 2;

    m_Epot = force(m_pos.memptr(), m_frc.memptr());
}

//----------------------------------------------------------------------------//

void vv_base::step(const double& dt)
{
    const double dt2 = dt/2;

    // Following equations 21-25 of dx.doi.org/10.1063/1.3489925
    // 1. Evolution of RP momenta under Hamiltonian V_{n}[q(t0)] by dt/2
    m_mom += dt2*m_frc;

    // 2. Advance positions
    m_pos += dt*m_mom/m_mass;

    // 3. Final evolution of RP momenta under Hamiltonian V_{n}[q(t0+dt)]
    m_Epot = force(m_pos.memptr(), m_frc.memptr());

    m_mom += dt2*m_frc;

    m_Ekin = 0.0;
    for (size_t i = 0; i < ndofs(); ++i) {
        m_Ekin += m_mom(i)*m_mom(i)/m_mass(i);
    }

    m_Ekin /= 2;
    m_temp_kT = m_Ekin*2.0/ndofs(); // not actual temperature, kT
}

//----------------------------------------------------------------------------//

double vv_base::invariant() const
{
    // Epot is in kcal/mol already
    return m_Ekin + m_Epot;
}

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
