#ifndef RPMD_BASE_H
#define RPMD_BASE_H

#include "rpmd-necklace.h"

namespace parts {

struct rpmd_base : public rpmd_necklace {

  virtual double force(size_t, size_t, size_t, const double *, double *) = 0;

  void init(size_t ndim, size_t natom, size_t nbead, const double &kT,
            const double &dt, const double *mass, const double *cartpos,
            const double *cartvel, double dummy);

  void step(const double &);
  void spring_energy();
  double invariant() const;

  constexpr static size_t nchain = 4;

protected:
  double m_kT;

  double m_Espring;
  double m_Epot_sum;
  double m_Ekin;
  double m_temp_kT;
  double m_temp_kT_centroid;
  double m_temp_kT_higherNM;

  void pimd_force();
};

} // namespace parts

#endif // RPMD_BASE_H
