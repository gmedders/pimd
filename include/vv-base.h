#ifndef VV_H
#define VV_H

#include <armadillo>

namespace parts {

struct vv_base {

  inline size_t ndim() const { return m_ndim; };
  inline size_t natoms() const { return m_natoms; };
  inline size_t ndofs() const { return m_ndof; };
  inline size_t nbeads() const { return m_nbeads; };
  // inline size_t beta() const {return m_beta;};

  inline double *cart_ptr() { return m_pos.memptr(); };

  virtual double force(const size_t, const size_t, const size_t, const double *,
                       double *) = 0;

  void init(size_t, size_t, const double &dt, const double *mass,
            const double *cartpos, const double *cartvel);

  void step(const double &);
  double invariant() const;

protected:
  double m_dt;
  size_t m_ndim;
  size_t m_natoms;
  size_t m_ndof;
  size_t m_nbeads;

  arma::vec m_pos;
  arma::vec m_mom;
  arma::vec m_frc;
  arma::vec m_mass;

  double m_Epot;
  double m_Ekin;
  double m_temp_kT;
};

} // namespace parts

#endif // VV_H
