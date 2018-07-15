#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H

#include <functional>
using namespace std::placeholders;

void run_finite_differences(std::function<double(const double *, double *)> V,
                            std::vector<double> &crd, std::vector<double> &f) {

  double dx = 5.0E-7;
  std::vector<double> weight = {1.0 / 12.0, -2.0 / 3.0, 2.0 / 3.0, -1.0 / 12.0};
  std::vector<double> displacement = {-2.0 * dx, -1.0 * dx, 1.0 * dx, 2.0 * dx};
  std::vector<double> dummy_grad = f;

  for (std::vector<double>::size_type i(0); i < crd.size(); ++i) {
    double grad(0.0);
    double orig(crd[i]);
    for (std::vector<double>::size_type j(0); j < displacement.size(); ++j) {
      crd[i] = orig + displacement[j];
      grad += weight[j] * V(&crd[0], &dummy_grad[0]);
    }
    crd[i] = orig;
    f[i] = -grad / dx;
  }
}

#endif // FINITE_DIFFERENCE_H
