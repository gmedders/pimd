#ifndef ANHARMONIC_H
#define ANHARMONIC_H

#include <cstdlib>

#include "single-state.h"

namespace pot {

class anharmonic : public single_state {

public:
  double VAA(const double *crd, double *f);

  void assert_ndim(int);

  void set_params(double *);

  double get_w() { return std::sqrt(2.0 * a / m); };
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
