#ifndef DOUBLE_WELL_H
#define DOUBLE_WELL_H

#include <cstdlib>

#include <armadillo>

#include "surface-hopping.h"

namespace pot {

class double_well : public surface_hopping {

public:
  double VAA(const double *x, double *f);
  double VBB(const double *x, double *f);

  void set_params(double *);

  double get_w() { return w; };
  double get_m() { return m; };

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
