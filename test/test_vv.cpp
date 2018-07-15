#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "sim_classes.h"
#include "gtest/gtest.h"

// surface_hopping is an abstract class, so there are only a few things to test
// in it. We'll create instances of the "double_well" class, which is the
// simplest child of surface_hopping, to test.

TEST(velocity_verlet, make_one_step) {
  double pos[] = {0.1};
  double vel[] = {-0.01};

  const int nbead(1);
  const int ndim(1);
  const int natom(1);

  double beta(1024);
  double dt(1.0);

  parts::vv sim;
  sim.set_up(nbead, ndim, natom, beta, dt, pos, vel);

  EXPECT_DOUBLE_EQ(sim.invariant(), 0.10001);
  EXPECT_DOUBLE_EQ(sim.Ep(), 1e-05);
  EXPECT_DOUBLE_EQ(sim.temp_kT(), 0.0);
  EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.10000000000000001);

  // Make one step
  sim.step(dt);

  EXPECT_DOUBLE_EQ(sim.invariant(), 0.10000999999952499);
  EXPECT_DOUBLE_EQ(sim.Ep(), 8.0999910000025029e-06);
  EXPECT_DOUBLE_EQ(sim.temp_kT(), 0.20000380001704998);
  EXPECT_DOUBLE_EQ(sim.avg_cart_pos(), 0.089999950000000009);
}

// inline size_t ndofs() const;
//
// inline double *cart_ptr() { return m_pos.memptr(); };
//
// virtual double force(const double *, double *) = 0;
//
// void init(size_t ndof, const double &dt, const double *mass,
//           const double *cartpos, const double *cartvel);
//
// void step(const double &);
// double invariant() const;
//
// protected:
// double m_dt;
// size_t m_ndofs;
//
// arma::vec m_pos;
// arma::vec m_mom;
// arma::vec m_frc;
// arma::vec m_mass;
//
// double m_Epot;
// double m_Ekin;
// double m_temp_kT;

// void set_up_new_init_cond(const size_t, const size_t, const size_t,
//                           const double, const double);
// void set_up_new_init_cond(const size_t, const size_t, const size_t,
//                           const double, const double, double *);
// void set_up(const size_t, const size_t, const size_t, const double,
//             const double, double *, double *);
// double force(const double *, double *);
//
// inline double Espring() const { return 0.0; }
// inline double Ep() const { return m_Epot; }
// inline double Ek() const { return m_Ekin; }
// inline double temp_kT() const { return m_temp_kT; }
// double avg_cart_pos() {
//   calc_pos_stats();
//   return m_avg_cart_pos;
// };
// // double avg_cart_pos() const { return m_avg_cart_pos; };
// double L1_cart_pos() const { return m_L1_cart_pos; };
// double L2_cart_pos() const { return m_L2_cart_pos; };
// double Linf_cart_pos() const { return m_Linf_cart_pos; };
// void calc_pos_stats(void);
//
// void print_params();
//
// potential_type m_potential;
//
// private:
// double m_avg_cart_pos;
// double m_L1_cart_pos;
// double m_L2_cart_pos;
// double m_Linf_cart_pos;
//
// size_t m_natom;
// size_t m_ndim;
// size_t m_ndofs;
