#ifndef AH_H
#define AH_H

#include <cstdlib>

#include <armadillo>

#include "surface-hopping.h"

namespace pot {

struct ah : public surface_hopping {

    double w;
    double m;
    double g;

    double a;
    double b;
    double Ed;

    double VAA(const double* x, double* f);
    double VBB(const double* x, double* f);

    void set_params(double*);

};

} // namespace pot

#endif // AH_H
