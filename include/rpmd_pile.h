#ifndef RPMD_PILE_H
#define RPMD_PILE_H

#include <random>

#include "rpmd_base.h"

namespace parts {

struct rpmd_pile : public rpmd_base {

  rpmd_pile() : rand_01(0, 1){};

  void init(size_t, size_t, size_t, const double &, const double &,
            const double *, const double *, const double *, double);
  void init_langevin(const double &dt, double gamma_centroid);

  void seed_pile_prng(int);

  void step(const double &);
  double invariant() const;

  std::mt19937 pile_prng; // Fixed seed of 0
  std::normal_distribution<> rand_01;

  double calc_KE();

  double m_tau0;
  arma::vec c1;
  arma::vec c2;

  arma::mat saved_mom;
};

} // namespace parts

#endif // RPMD_PILE_H
