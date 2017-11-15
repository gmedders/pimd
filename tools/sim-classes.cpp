#include "sim-classes.h"

namespace parts {

////////////////////////////////////////////////////////////////////////////////

// Do this only for the first atom and first dimension
void rpmd::calc_pos_stats(void)
{
    m_avg_cart_pos = 0;
    for(size_t n = 0; n < nbeads(); ++n)
        m_avg_cart_pos += m_pos_cart(0,n);
    m_avg_cart_pos /= nbeads();

    m_Linf_cart_pos = 0;
    m_L2_cart_pos = 0;
    m_L1_cart_pos = 0;
    for(size_t n = 0; n < nbeads(); ++n){
        double diff = m_pos_cart(0,n) - m_avg_cart_pos;
        if(diff > m_Linf_cart_pos)
            m_Linf_cart_pos = diff;
        m_L1_cart_pos += std::abs(diff);
        m_L2_cart_pos += diff*diff;
    }
    m_L2_cart_pos = std::sqrt(m_L2_cart_pos)/nbeads();
    m_L1_cart_pos = m_L1_cart_pos/nbeads();
}

//----------------------------------------------------------------------------//

void rpmd::dump_1D_frame(std::ostream& of_traj)
{
    //of_traj << nbeads() << ' ' << ndofs() << ' ' << beta() << std::endl;
    for(size_t n = 0; n < nbeads(); ++n) {
        of_traj << m_potential.state_id[n] << ' ';
        for(size_t i = 0; i < ndofs(); ++i) {
            of_traj << ' ' << m_pos_cart(i,n)
                    << ' ' << m_mom_cart(i,n)/m_mass(i);
        }
        of_traj << std::endl;
    }
}

//----------------------------------------------------------------------------//

void rpmd::print_params()
{
    m_potential.print_params();
    std::cout << "# gammaTh_fac = " << m_gamma << std::endl;
}
//----------------------------------------------------------------------------//

void rpmd::set_up(const size_t ndim, const size_t natoms, const size_t nbead,
                  const double beta, const double dt,
                  double* pos, double* vel)
{
    size_t ndof = ndim*natoms;

    // Generate initial positions, if needed
    std::vector<double> all_bead_crd;
    if(pos == nullptr) {
        for(size_t i = 0; i < nbead*ndof; ++i){
            all_bead_crd.push_back(0.0);
        }
        pos = &all_bead_crd[0];
    }

    // Generate initial velocities, if needed
    std::vector<double> all_bead_vel;
    if(vel == nullptr) {
        double kT = 1.0/beta;

        for(size_t i = 0; i < nbead*ndof; ++i){
            const double sigma = std::sqrt(kT/atm_mass);
            double v = sigma*randn(0,1);
            all_bead_vel.push_back(v);
        }
        vel = &all_bead_vel[0];
    }

    // Now populate the mass array
    double mass[ndof];
    for (size_t i = 0; i < natoms; ++i) {
        const double Mi = atm_mass;
        for (size_t k = 0; k < ndim; ++k)
            mass[k + ndim*i] = Mi;
    }

    // setup the simulation
    m_potential.set_params(params);

    init(ndim, natoms, nbead, 1.0/beta, dt, mass, pos, vel, m_gamma);
}

//----------------------------------------------------------------------------//

void rpmd::set_gammaTh(const double& dt, double gam_fac)
{
    m_gamma = gam_fac*m_potential.get_w();
    init_langevin(dt, m_gamma);
}

//----------------------------------------------------------------------------//

double rpmd::force(size_t ndim, size_t natoms, size_t nbeads,
                   const double* x, double* f)
{
    return m_potential.force(ndim, natoms, nbeads, x, f);
}

////////////////////////////////////////////////////////////////////////////////

void vv::calc_pos_stats(void)
{
    assert(ndim()*natoms() == 1);

    m_avg_cart_pos = m_pos[0];
}

//----------------------------------------------------------------------------//

void vv::set_up(const size_t ndim, const size_t natoms, const size_t nbead,
                double beta, const double dt, double* pos, double* vel)
{
    assert(nbead == 1);

    size_t ndof = ndim*natoms;

    // Generate initial positions, if needed
    std::vector<double> all_bead_crd;
    if(pos == nullptr) {
        for(size_t i = 0; i < nbead*ndof; ++i){
            all_bead_crd.push_back(0.0);
        }
        pos = &all_bead_crd[0];
    }

    // Generate initial velocities, if needed
    std::vector<double> all_bead_vel;
    if(vel == nullptr) {
        double kT = 1.0/beta;

        for(size_t i = 0; i < nbead*ndof; ++i){
            const double sigma = std::sqrt(kT/atm_mass);
            double v = sigma*randn(0,1);
            all_bead_vel.push_back(v);
        }
        vel = &all_bead_vel[0];
    }

    // prepare masses
    double* mass = new double[ndim*natoms]; // for every degree of freedom
    for (size_t i = 0; i < natoms; ++i) {
        const double Mi = atm_mass;
        for (size_t k = 0; k < ndim; ++k)
            mass[k + ndim*i] = Mi;
    }

    // setup the simulation

    m_potential.set_params(params);

    init(ndim, natoms, dt, mass, pos, vel);

    // clean up

    delete[] mass;
}

//----------------------------------------------------------------------------//

double vv::force(size_t ndim, size_t natoms, size_t nbeads,
                 const double* x, double* f)
{
    return m_potential.force(ndim, natoms, nbeads, x, f);
}

//----------------------------------------------------------------------------//

void vv::print_params()
{
    m_potential.print_params();
    std::cout << "# gammaTh_fac = 0.0" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////
