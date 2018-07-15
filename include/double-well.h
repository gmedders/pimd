#ifndef DOUBLE_WELL_H
#define DOUBLE_WELL_H

#include <cstdlib>

#include <armadillo>

#include "surface-hopping.h"

namespace pot {

struct double_well : public surface_hopping {
  double VAA(const double *x, double *f);
  double VBB(const double *x, double *f);

  void set_params(double *);

  double get_w() { return w; };
  double get_m() { return m; };
  double get_g() { return g; };
  double get_dG() { return dG; };

  double get_a() { return a; };

  void print_params();

private:
  double w;
  double m;
  double g;
  double dG;

  double a;
};

} // namespace pot

#endif // DOUBLE_WELL_H
