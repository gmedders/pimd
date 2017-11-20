#include "gcmc.h"

//
// units are au
//

namespace parts {

////////////////////////////////////////////////////////////////////////////////

// void gcmc::set_md_ensemble(std::string& ensemble_name)
// {
//     if(ensemble_name.compare("rpmd"))
//         ensemble_type = rpmd;
//     else
//         ensemble_type = vv;
// }

//----------------------------------------------------------------------------//

void gcmc::set_chemical_potential(double mu)
{
    assert(mu >= 0);
    m_chemical_potential = mu;
}

//----------------------------------------------------------------------------//

void gcmc::set_up(const size_t ndim, const size_t natoms, const size_t nbead,
                  const double beta, const double dt, double* pos, double* vel)
{
    //m_md_ensemble = std::unique_ptr<rpmd>();
    m_md_ensemble->set_up(ndim, natoms, nbead, beta, dt, pos, vel);
}

//----------------------------------------------------------------------------//

bool gcmc::calc_insertion_probability()
{
    return true;
}

//----------------------------------------------------------------------------//

void gcmc::step(double dt, double beta)
{
    bool insert = calc_insertion_probability();

    if(insert) {
        size_t ndim = m_md_ensemble->ndim();
        size_t old_natoms = m_md_ensemble->natoms();
        size_t nbeads = m_md_ensemble->nbeads();

        size_t new_natoms = old_natoms + 1;

        size_t new_total_ndof = ndim * new_natoms * nbeads;
        double new_crd[ndim * new_natoms * nbeads];
        double new_vel[ndim * new_natoms * nbeads];

        const double* old_crd = m_md_ensemble->get_crd();
        const double* old_vel = m_md_ensemble->get_vel();

        std::fill(new_crd, new_crd + new_total_ndof, 0.0);
        std::fill(new_vel, new_vel + new_total_ndof, 0.0);

        // Copy the old coordinate and velocity
        // Leave a space for the new atom when incrementing the bead offset
        for(size_t n = 0; n < nbeads; ++n){
            size_t old_offset_bead = n * old_natoms * ndim;
            size_t new_offset_bead = n * new_natoms * ndim;
            for(size_t i = 0; i < old_natoms; ++i){
                size_t offset_atom = i * ndim;
                for(size_t j = 0; j < ndim; ++j){
                    size_t old_offset = old_offset_bead + offset_atom + j;
                    size_t new_offset = new_offset_bead + offset_atom + j;
                    new_crd[new_offset] = old_crd[old_offset];
                    new_vel[new_offset] = old_vel[old_offset];
                }
            }
            // Avoid setting the new crd/vel by hand by initializing array to 0
            //for(size_t i = old_natoms; i < new_natoms; ++i){
            //    size_t offset_atom = i * ndim;
            //    for(size_t j = 0; j < ndim; ++j){
            //        size_t old_offset = old_offset_bead + offset_atom + j;
            //        size_t new_offset = new_offset_bead + offset_atom + j;
            //        new_crd[new_offset] = 0.0;
            //        new_vel[new_offset] = 0.0;
            //    }
            //}
        }

        std::vector<int> states;
        for(size_t i = 0; i < nbeads; ++i)
            states.push_back(m_md_ensemble->m_potential->state_id[i]);

        std::unique_ptr<rpmd> new_ensemble(new rpmd);
        m_md_ensemble = std::move(new_ensemble);

        m_md_ensemble->set_up(ndim, new_natoms, nbeads, beta, dt,
                              new_crd, new_vel);
        m_md_ensemble->m_potential->set_individual_bead_states(states);
    }

    m_md_ensemble->step(dt);
}

//----------------------------------------------------------------------------//

void gcmc::dump(std::ostream& os)
{
    os << m_md_ensemble->ndim() << ' '
       << m_md_ensemble->natoms() << ' '
       << m_md_ensemble->nbeads() << std::endl;
    m_md_ensemble->dump_1D_frame(os);
}

////////////////////////////////////////////////////////////////////////////////

} // namespace
