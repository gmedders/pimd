#ifndef GCMC_H
#define GCMC_H

#include "sim-classes.h"

//
// units are au
//

namespace parts {

////////////////////////////////////////////////////////////////////////////////

class gcmc {

    void set_chemical_potential(double);
    void set_md_ensemble(std::string&);

    void set_up(const size_t, const size_t, const size_t,
                const double, const double,
                double*, double*);

    void step(double);

    typedef ensemble_type m_md_ensemble;

private:
    double m_chemical_potential = 0;
};

////////////////////////////////////////////////////////////////////////////////

} // namespace

#endif // GCMC_H
