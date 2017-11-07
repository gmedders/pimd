#include "sim-classes.h"

namespace parts {

////////////////////////////////////////////////////////////////////////////////

pimd::~pimd()
{
    delete[] pos;
}

//----------------------------------------------------------------------------//

void pimd::calc_pos_stats(void)
{
    assert(m_ndofs == 1);

    m_avg_cart_pos = 0;
    for(size_t i = 0; i < m_nbead; ++i)
        m_avg_cart_pos += m_pos_cart[i*m_ndofs];
    m_avg_cart_pos /= m_nbead;

    m_Linf_cart_pos = 0;
    m_L2_cart_pos = 0;
    m_L1_cart_pos = 0;
    for(size_t i = 0; i < m_nbead; ++i){
        double diff = m_pos_cart[i*m_ndofs] - m_avg_cart_pos;
        if(diff > m_Linf_cart_pos)
            m_Linf_cart_pos = diff;
        m_L1_cart_pos += std::abs(diff);
        m_L2_cart_pos += diff*diff;
    }
    m_L2_cart_pos = std::sqrt(m_L2_cart_pos)/m_nbead;
    m_L1_cart_pos = m_L1_cart_pos/m_nbead;
}

//----------------------------------------------------------------------------//

void pimd::dump_1D_frame(std::ofstream& of_traj)
{
    vel_n2c();

    of_traj << m_nbead << ' ' << m_ndofs << ' ' << m_beta << std::endl;
    for(size_t i = 0; i < m_nbead; ++i) {
        of_traj << m_potential.state_id[i] << ' ';
        for(size_t j = 0; j < m_ndofs; ++j) {
            of_traj << ' ' << m_pos_cart[i*m_ndofs + j]
                    << ' ' << m_vel_cart[i*m_ndofs + j];
        }
        of_traj << std::endl;
    }
}

//----------------------------------------------------------------------------//

void pimd::print_params()
{
    m_potential.print_params();
    std::cout << "# gammaTh_fac = 0.0" << std::endl;
}

//----------------------------------------------------------------------------//

void pimd::set_up(const size_t nbead, const size_t ndim, const size_t natom,
                  const double beta)
{
    m_natom = natom;
    m_ndim = ndim;
    m_nbead = nbead;

    m_beta = beta;

    m_ndofs = m_natom*m_ndim;

    pos = new double[m_ndofs];

    pos[0] = 0.0;

    // prepare masses
    double* mass = new double[m_ndofs]; // for every degree of freedom
    for (size_t i = 0; i < natom; ++i) {
        const double Mi = atm_mass;
        for (size_t k = 0; k < ndim; ++k)
            mass[k + ndim*i] = Mi;
    }

    // setup the simulation

    m_potential.set_params(params);

    init(m_ndofs, nbead, 1.0/m_beta, mass, pos);

    // clean up

    delete[] mass;
}

//----------------------------------------------------------------------------//

double pimd::force(const double* x, double* f)
{
    return m_potential.force(m_ndim, m_natom, m_nbead, x, f);
}

////////////////////////////////////////////////////////////////////////////////

// Do this only for the first atom and first dimension
void rpmd::calc_pos_stats(void)
{
    m_avg_cart_pos = 0;
    for(size_t n = 0; n < m_nbead; ++n)
        m_avg_cart_pos += m_pos_cart(0,n);
    m_avg_cart_pos /= m_nbead;

    m_Linf_cart_pos = 0;
    m_L2_cart_pos = 0;
    m_L1_cart_pos = 0;
    for(size_t n = 0; n < m_nbead; ++n){
        double diff = m_pos_cart(0,n) - m_avg_cart_pos;
        if(diff > m_Linf_cart_pos)
            m_Linf_cart_pos = diff;
        m_L1_cart_pos += std::abs(diff);
        m_L2_cart_pos += diff*diff;
    }
    m_L2_cart_pos = std::sqrt(m_L2_cart_pos)/m_nbead;
    m_L1_cart_pos = m_L1_cart_pos/m_nbead;
}

//----------------------------------------------------------------------------//

void rpmd::dump_1D_frame(std::ofstream& of_traj)
{
    //of_traj << m_nbead << ' ' << m_ndofs << ' ' << m_beta << std::endl;
    for(size_t n = 0; n < m_nbead; ++n) {
        of_traj << m_potential.state_id[n] << ' ';
        for(size_t i = 0; i < m_ndofs; ++i) {
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

void rpmd::set_up_new_init_cond(const size_t nbead, const size_t ndim,
                                const size_t natom, const double beta,
                                const double dt)
{
    std::vector<double> all_bead_crd;

    size_t ndofs = ndim*natom;

    for(size_t i = 0; i < nbead*ndofs; ++i){
        all_bead_crd.push_back(0.0);
    }
    set_up_new_init_cond(nbead, ndim, natom, beta, dt,
                         &all_bead_crd[0]);
}

//----------------------------------------------------------------------------//

void rpmd::set_up_new_init_cond(const size_t nbead, const size_t ndim,
                                const size_t natom, const double beta,
                                const double dt, double* pos)
{
    std::vector<double> all_bead_vel;

    size_t ndofs = ndim*natom;
    double kT = 1.0/beta;

    for(size_t i = 0; i < nbead*ndofs; ++i){
        const double sigma = std::sqrt(kT/atm_mass);
        double v = sigma*randn(0,1);
        all_bead_vel.push_back(v);
    }
    set_up(nbead, ndim, natom, beta, dt,
           pos, &all_bead_vel[0]);
}

//----------------------------------------------------------------------------//

void rpmd::set_up(const size_t nbead, const size_t ndim, const size_t natom,
                  const double beta, const double dt, double* pos, double* vel)
{
    m_nbead = nbead;
    m_ndim = ndim;
    m_natom = natom;

    m_beta = beta;

    m_ndofs = m_natom*m_ndim;

    // prepare masses
    double* mass = new double[m_ndofs]; // for every degree of freedom
    for (size_t i = 0; i < natom; ++i) {
        const double Mi = atm_mass;
        for (size_t k = 0; k < ndim; ++k)
            mass[k + ndim*i] = Mi;
    }

    // setup the simulation

    m_potential.set_params(params);

    init(m_ndofs, nbead, 1.0/beta, dt, mass, pos, vel, m_gamma);

    // clean up

    delete[] mass;
}

//----------------------------------------------------------------------------//

void rpmd::set_gammaTh(const double& dt, double gam_fac)
{
    m_gamma = gam_fac*m_potential.get_w();
    init_langevin(dt, m_gamma);
}

//----------------------------------------------------------------------------//

double rpmd::force(const double* x, double* f)
{
    return m_potential.force(m_ndim, m_natom, m_nbead, x, f);
}

////////////////////////////////////////////////////////////////////////////////

void vv::calc_pos_stats(void)
{
    assert(m_ndofs == 1);

    m_avg_cart_pos = m_pos[0];
    m_Linf_cart_pos = 0;
    m_L2_cart_pos = 0;
    m_L1_cart_pos = 0;
}

//----------------------------------------------------------------------------//

void vv::set_up_new_init_cond(const size_t nbead, const size_t ndim,
                                const size_t natom, const double beta,
                                const double dt)
{
    std::vector<double> all_crd;

    size_t ndofs = ndim*natom;

    for(size_t i = 0; i < ndofs; ++i){
        all_crd.push_back(0.0);
    }
    set_up_new_init_cond(nbead, ndim, natom, beta, dt,
                         &all_crd[0]);
}

//----------------------------------------------------------------------------//

void vv::set_up_new_init_cond(const size_t nbead, const size_t ndim,
                                const size_t natom, const double beta,
                                const double dt, double* pos)
{
    std::vector<double> all_vel;

    size_t ndofs = ndim*natom;
    double kT = 1.0/beta;

    for(size_t i = 0; i < ndofs; ++i){
        const double sigma = std::sqrt(kT/atm_mass);
        double v = sigma*randn(0,1);
        all_vel.push_back(v);
    }
    set_up(nbead, ndim, natom, beta, dt,
           pos, &all_vel[0]);
}

//----------------------------------------------------------------------------//

void vv::set_up(const size_t nbead, const size_t ndim, const size_t natom,
                double beta, const double dt, double* pos, double* vel)
{
    assert(nbead == 1);

    m_ndim = ndim;
    m_natom = natom;

    m_ndofs = m_natom*m_ndim;

    // prepare masses
    double* mass = new double[m_ndofs]; // for every degree of freedom
    for (size_t i = 0; i < natom; ++i) {
        const double Mi = atm_mass;
        for (size_t k = 0; k < ndim; ++k)
            mass[k + ndim*i] = Mi;
    }

    // setup the simulation

    m_potential.set_params(params);

    init(m_ndofs, dt, mass, pos, vel);

    // clean up

    delete[] mass;
}

//----------------------------------------------------------------------------//

double vv::force(const double* x, double* f)
{
    return m_potential.force(m_ndim, m_natom, 1, x, f);
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
