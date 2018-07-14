#ifndef ANHARMONIC_H
#define ANHARMONIC_H

#include <cstdlib>

#include "single-state.h"

namespace pot {

struct anharmonic : public single_state {

  double VAA(const double *crd, double *f);

  void set_params(double *);

  double get_w() { return std::sqrt(2.0 * a / m); };
  double get_a() { return a; };
  double get_b() { return b; };
  double get_c() { return c; };
  double get_m() { return m; };

  void print_params();

private:
  double a;
  double b;
  double c;

  double m;
};

} // namespace pot

#endif // ANHARMONIC_H
