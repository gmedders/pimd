#ifndef SHO_H
#define SHO_H

#include <cstdlib>

namespace pot {

struct sho {

    double w;
    double m;

    double a;

    double operator()(size_t natom, const double* x, double* f) const;

    void set_params(double*);

    bool init = false;

};
        
} // namespace pot

#endif // SHO_H
