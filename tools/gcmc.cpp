#include "gcmc.h"

//
// units are au
//

namespace parts {

////////////////////////////////////////////////////////////////////////////////

void set_md_ensemble(std::string& ensemble_name)
{
    if(ensemble_name.compare("rpmd"))
        ensemble_type = rpmd;
    else
        ensemble_type = vv;
}

void set_chemical_potential(double mu)
{
    assert(mu >= 0);
    m_chemical_potential = mu;
}

void gcmc::set_up(const size_t ndim, const size_t natoms, const size_t nbead,
                  const double beta, const double dt, double* pos, double* vel)
{
    m_md_ensemble.set_up(ndim, natom, nbead, beta, dt, pos, vel);
}

bool gcmc::calc_insertion_probability()
{
    return true;
}

void gcmc::step(double dt, double beta)
{
    bool insert = calc_insertion_probability();

    if(insert) {
        size_t curr_ndim = m_md_ensemble.ndim();
        size_t curr_natoms = m_md_ensemble.natoms();
        size_t curr_nbeads = m_md_ensemble.nbeads();

        size_t new_natoms = curr_natoms + 1;

        size_t new_total_ndof = curr_ndim * new_natoms * curr_nbeads;
        double new_crd[curr_ndim * new_natoms * curr_nbeads];
        double new_vel[curr_ndim * new_natoms * curr_nbeads];

        double* curr_crd = m_md_ensemble.crd();
        double* curr_vel = m_md_ensemble.vel();

        std::copy(curr_crd, curr_crd + new_total_ndof, new_crd);
        std::copy(curr_vel, curr_vel + new_total_ndof, new_vel);

        ~m_md_ensemble();
        m_md_ensemble();

        m_md_ensemble.set_up(curr_ndim, new_natom, curr_nbeads, beta, dt,
                             new_crd, new_vel);
    }

    m_md_ensemble.step(dt);
}

////////////////////////////////////////////////////////////////////////////////

} // namespace
