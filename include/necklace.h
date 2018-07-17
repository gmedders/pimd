#ifndef NECKLACE_H
#define NECKLACE_H

#include <fftw3.h>

#include <cassert>
#include <cstddef>

namespace parts {

//
// necklace for PIMD/CMD; manages memory for cartesian/normal mode
// positions, velocities, forces and handles the transforms
//

struct necklace {

  necklace();
  ~necklace();

  inline size_t ndim() const;
  inline size_t natoms() const;
  inline size_t ndofs() const;
  inline size_t nbeads() const;

  inline const double &lambda(size_t) const;

private:
  necklace(const necklace &);
  necklace &operator=(const necklace &);

  size_t m_ndim;
  size_t m_natom;
  size_t m_ndofs;
  size_t m_nbeads;

  fftw_plan m_plan_cart2nmode;
  fftw_plan m_plan_nmode2cart;

public:
  void setup(size_t ndim, size_t natom, size_t nbead);

  void pos_c2n();
  void pos_n2c();

  void vel_c2n();
  void vel_n2c();

  void frc_c2n();

  // layout is bead1, bead2, ..., beadN, where
  // each bead consists of ndofs elements

  double *m_pos_cart;
  double *m_pos_nmode;

  double *m_vel_cart;
  double *m_vel_nmode;

  double *m_frc_cart;
  double *m_frc_nmode;

  double *m_lambda;
};

inline size_t necklace::ndim() const { return m_ndim; }

inline size_t necklace::natoms() const { return m_natom; }

inline size_t necklace::ndofs() const { return m_ndofs; }

inline size_t necklace::nbeads() const { return m_nbeads; }

inline const double &necklace::lambda(size_t b) const {
  assert(b < nbeads());

  return m_lambda[b];
}

} // namespace parts

#endif // NECKLACE_H
