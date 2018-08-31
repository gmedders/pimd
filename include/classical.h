#ifndef CLASSICAL_H
#define CLASSICAL_H

#include <armadillo>

#include "eom.h"

namespace parts {

struct classical : public eom {

  void init(size_t, size_t, const double &dt, const double beta,
            const double atm_mass, const double *cartpos,
            const double *cartvel);

  inline double *get_crd() { return m_pos.memptr(); };
  inline double *get_vel() {
    m_vel = m_mom / m_mass;
    return m_vel.memptr();
  };

  void calc_pos_stats(void);
  void dump_1D_frame(std::ostream &);

protected:
  arma::vec m_pos;
  arma::vec m_mom;
  arma::vec m_vel;
  arma::vec m_frc;
  arma::vec m_mass;
};

} // namespace parts

#endif // CLASSICAL_H
