#ifndef EOM_H
#define EOM_H

#include <memory>

#include "pot_definition.h"

//
// units are au
//

namespace parts {

////////////////////////////////////////////////////////////////////////////////

struct eom {

public:
  eom()
      : m_potential(std::make_unique<potential_type>()), m_dt(1.0), m_ndim(1),
        m_natom(1), m_ndof(1), m_Epot(0.0), m_Ekin(0.0), m_temp_kT(0.0),
        m_avg_cart_pos(0.0), m_L1_cart_pos(0.0), m_L2_cart_pos(0.0),
        m_Linf_cart_pos(0.0), m_gamma(0.0){};

  void set_up(const size_t ndim, const size_t natom, const double dt);

  void generate_boltzmann_velocities(const double beta, const double atm_mass,
                                     double *vel);
  double force(size_t, size_t, size_t, const double *, double *);
  virtual void step(const double &) = 0;
  virtual double invariant() const = 0;

  virtual void calc_pos_stats(void) = 0;
  virtual void dump_1D_frame(std::ostream &) = 0;

  inline size_t ndim() const { return m_ndim; };
  inline size_t natoms() const { return m_natom; };
  inline size_t ndofs() const { return m_ndof; };

  inline double Ep() const { return m_Epot; }
  inline double Ek() const { return m_Ekin; }
  inline double temp_kT() const { return m_temp_kT; }
  double avg_cart_pos() {
    calc_pos_stats();
    return m_avg_cart_pos;
  };
  double L1_cart_pos() const { return m_L1_cart_pos; };
  double L2_cart_pos() const { return m_L2_cart_pos; };
  double Linf_cart_pos() const { return m_Linf_cart_pos; };

  virtual double *get_crd() = 0;
  virtual double *get_vel() = 0;

  void set_gammaTh(const double &, double);
  void print_params();

  std::unique_ptr<potential_type> m_potential;

protected:
  double m_dt;
  size_t m_ndim;
  size_t m_natom;
  size_t m_ndof;

  double m_Epot;
  double m_Ekin;
  double m_temp_kT;

  double m_avg_cart_pos;
  double m_L1_cart_pos;
  double m_L2_cart_pos;
  double m_Linf_cart_pos;

  double m_gamma;
};

////////////////////////////////////////////////////////////////////////////////

} // namespace parts

#endif // EOM_H
