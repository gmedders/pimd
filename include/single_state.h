#ifndef SINGLE_STATE_H
#define SINGLE_STATE_H

#include <cstdlib>
#include <functional>

#include <armadillo>

#include "bead_states.h"

#include "lj.h"

namespace pot {

typedef lj pot_2B;

class single_state : public bead_states {

public:
  // single_state() {std::cerr << "single_state::single_state()" << std::endl;
  // }; ~single_state() {std::cerr << "single_state::~single_state()" <<
  // std::endl; };
  bool init_pot = false;

  double force(size_t, size_t, size_t, const double *, double *);
  virtual double VAA(const double *x, double *f) = 0;
  virtual void set_params(double *) = 0;
  virtual void assert_ndim(int) = 0;

  void set_hopping_params(double *){};

  void print_state_params();

private:
  pot_2B m_pot_2B;
};

} // namespace pot

#endif // SINGLE_STATE_H
