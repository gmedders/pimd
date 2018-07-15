#ifndef AH_H
#define AH_H

#include <cstdlib>

#include "surface-hopping.h"

namespace pot {

struct ah : public surface_hopping {

  double VAA(const double *x, double *f);
  double VBB(const double *x, double *f);

  void set_params(double *);

  double get_w() { return w; };
  double get_m() { return m; };
  double get_g() { return g; };

  double get_a() { return a; };
  double get_b() { return b; };
  double get_Ed() { return Ed; };
  double get_Ed_bar() { return Ed_bar; };

  void print_params();

private:
  double w;
  double m;
  double g;

  double a;
  double b;
  double Ed;
  double Ed_bar;
};

} // namespace pot

#endif // AH_H
