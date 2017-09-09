#ifndef ANHARMONIC_H
#define ANHARMONIC_H

#include <cstdlib>

namespace pot {

struct anharmonic {

    double a;
    double b;
    double c;

    double operator()(size_t natom, const double* crd, double* f) const;

    void set_params(double, double, double);

    bool init = false;

};
        
} // namespace pot

#endif // ANHARMONIC_H
