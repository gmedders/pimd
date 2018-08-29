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
      : m_potential(std::make_unique<potential_type>()), m_Epot(0.0),
        m_Ekin(0.0), m_temp_kT(0.0), m_avg_cart_pos(0.0), m_L1_cart_pos(0.0),
        m_L2_cart_pos(0.0), m_Linf_cart_pos(0.0){};

  void set_up(const size_t ndim, const size_t natom, const size_t nbeads,
              const double beta, const double dt, double *crd = nullptr,
              double *vel = nullptr);
  double force(size_t, size_t, size_t, const double *, double *);

  virtual double Espring() = 0;
  virtual void calc_pos_stats(void) = 0;
  virtual void dump_1D_frame(std::ostream &) = 0;

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

  virtual const double *get_crd() = 0;
  virtual const double *get_vel() = 0;

  void set_gammaTh(const double &, double);
  void print_params();

  std::unique_ptr<potential_type> m_potential;

protected:
  double m_Epot;
  double m_Ekin;
  double m_temp_kT;

  double m_avg_cart_pos;
  double m_L1_cart_pos;
  double m_L2_cart_pos;
  double m_Linf_cart_pos;

  double m_gamma = 0.0;
};

////////////////////////////////////////////////////////////////////////////////

} // namespace parts

#endif // EOM_H
