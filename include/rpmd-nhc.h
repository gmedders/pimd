#ifndef RPMD_NHC_H
#define RPMD_NHC_H

#include "rpmd-base.h"

namespace parts {

struct rpmd_nhc : public rpmd_base {

    rpmd_nhc();
    ~rpmd_nhc();

    void init(size_t ndof, size_t nbead,
              const double& kT, const double& dt,
              const double* mass, const double* cartpos, const double* cartvel,
              double dummy);

    void step(const double&);
    double invariant() const;

    double m_omega_M;

    double m_tau;
    double* m_thermostats;
};

} // namespace parts

#endif // RPMD_NHC_H
