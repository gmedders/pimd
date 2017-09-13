#ifndef DOUBLE_WELL_H
#define DOUBLE_WELL_H

#include <cstdlib>

#include <armadillo>

#include "surface-hopping.h"

namespace pot {

struct double_well : public surface_hopping {

    double w;
    double m;
    double g;

    double a;

    double VAA(const double* x, double* f);
    double VBB(const double* x, double* f);

    void set_params(double*);

};

} // namespace pot

#endif // DOUBLE_WELL_H
