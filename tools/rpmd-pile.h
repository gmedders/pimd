#ifndef RPMD_PILE_H
#define RPMD_PILE_H

#include "rpmd-base.h"

#include "mt19937.h"

namespace parts {

struct rpmd_pile : public rpmd_base {

    void init(size_t ndof, size_t nbead,
              const double& kT, const double& dt,
              const double* mass, const double* cartpos, const double* cartvel,
              double tau);

    void step(const double&);
    double invariant() const;

    double m_tau0;
    arma::vec c1;
    arma::vec c2;

    mt19937 prg;
};

} // namespace parts

#endif // RPMD_PILE_H
