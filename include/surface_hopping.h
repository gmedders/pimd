#ifndef SURFACE_HOPPING_H
#define SURFACE_HOPPING_H

#include <cstdlib>

#include <armadillo>

#include "bead_states.h"
#include "explicit_bath.h"

namespace pot {

class surface_hopping : public explicit_bath, public bead_states {

public:
  bool init_pot = false;
  bool init_hop = false;

  double force(size_t, size_t, size_t, const double *, double *);
  virtual double VAA(const double *x, double *f) = 0;
  virtual double VBB(const double *x, double *f) = 0;
  //    virtual double bath_force(const double* x, double* f) = 0;
  virtual void set_params(double *) = 0;

  double fermi_function(const double, const double);
  double hop_probability(const double, int);

  void set_hopping_params(double *);

  void print_state_params();

  double Gamma;
  double dt;
  double beta;
  double voltage;
};

} // namespace pot

#endif // SURFACE_HOPPING_H
