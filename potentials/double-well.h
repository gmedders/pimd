#ifndef DOUBLE_WELL_H
#define DOUBLE_WELL_H

#include <cstdlib>

namespace pot {

struct double_well {

    double w;
    double m;
    double g;

    double a;

    int active_state;

    double operator()(size_t, size_t, const double*, double*) const;
    double VAA(const double* x, double* f);
    double VBB(const double* x, double* f);

    double hop_probability()(size_t natom, const double* x, double* f) const;

    void set_params(double*);

    bool init = false;

};
        
} // namespace pot

#endif // DOUBLE_WELL_H
