#ifndef RPMD_PILE_H
#define RPMD_PILE_H

#include "rpmd-base.h"

#include "rand_gauss.h"

namespace parts {

struct rpmd_pile : public rpmd_base {

  void init(size_t, size_t, size_t, const double &, const double &,
            const double *, const double *, const double *, double);
  void init_langevin(const double &dt, double gamma_centroid);

  void step(const double &);
  double invariant() const;

  double calc_KE();

  double m_tau0;
  arma::vec c1;
  arma::vec c2;

  arma::mat saved_mom;
};

} // namespace parts

#endif // RPMD_PILE_H
