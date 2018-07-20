#ifndef SIM_CLASSES_H
#define SIM_CLASSES_H

#include <memory>

#include "pimd_base.h"

#include "rpmd_base.h"
#include "rpmd_nhc.h"
#include "rpmd_pile.h"

#include "vv_base.h"

#include "pot_definition.h"

//
// units are au
//

namespace parts {

////////////////////////////////////////////////////////////////////////////////

// struct rpmd : public parts::rpmd_base {
class rpmd : public parts::rpmd_pile {
  // struct rpmd : public parts::rpmd_nhc {

public:
  rpmd() : m_potential(new potential_type()){};
  ~rpmd(){};

  void set_up(const size_t ndim, const size_t natom, const size_t nbeads,
              const double beta, const double dt, double *crd = nullptr,
              double *vel = nullptr);
  double force(size_t, size_t, size_t, const double *, double *);

  inline double Espring() const { return m_Espring; }
  inline double Ep() const { return m_Epot_sum; }
  inline double Ek() const { return m_Ekin; }
  inline double temp_kT() const { return m_temp_kT; }
  inline double temp_kT_centroid() const { return m_temp_kT_centroid; }
  inline double temp_kT_higherNM() const { return m_temp_kT_higherNM; }
  double avg_cart_pos() {
    calc_pos_stats();
    return m_avg_cart_pos;
  };
  // double avg_cart_pos() const { return m_avg_cart_pos; };
  double L1_cart_pos() const { return m_L1_cart_pos; };
  double L2_cart_pos() const { return m_L2_cart_pos; };
  double Linf_cart_pos() const { return m_Linf_cart_pos; };
  void calc_pos_stats(void);

  const double *get_crd() { return m_pos_cart.memptr(); }
  const double *get_vel() {
    mom2vel();
    return m_vel_cart.memptr();
  }

  void set_gammaTh(const double &, double);

  void dump_1D_frame(std::ostream &);
  void print_params();

  std::unique_ptr<potential_type> m_potential;

  double m_gamma = 0.0;

private:
  double m_avg_cart_pos;
  double m_L1_cart_pos;
  double m_L2_cart_pos;
  double m_Linf_cart_pos;
};

////////////////////////////////////////////////////////////////////////////////

struct vv : public parts::vv_base {

  vv() : m_potential(new potential_type()){};

  void set_up(const size_t ndim, const size_t natom, const size_t nbead,
              const double beta, const double dt, double *crd = nullptr,
              double *vel = nullptr);
  double force(size_t, size_t, size_t, const double *, double *);

  inline double Espring() const { return 0.0; }
  inline double Ep() const { return m_Epot; }
  inline double Ek() const { return m_Ekin; }
  inline double temp_kT() const { return m_temp_kT; }

  double avg_cart_pos() {
    calc_pos_stats();
    return m_avg_cart_pos;
  };
  double L1_cart_pos() const { return 0.0; };
  double L2_cart_pos() const { return 0.0; };
  double Linf_cart_pos() const { return 0.0; };
  void calc_pos_stats(void);

  void print_params();

  std::unique_ptr<potential_type> m_potential;

private:
  double m_avg_cart_pos;
};

////////////////////////////////////////////////////////////////////////////////

} // namespace parts

#endif // SIM_CLASSES_H
