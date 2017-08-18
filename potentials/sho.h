#ifndef SHO_H
#define SHO_H

#include <cstdlib>

namespace pot {

struct sho {

    double w = 1.0;
    double m = 1.0;

    // a = 1/2 m w^2
    double a = 0.5 * m * w * w;

    double operator()(size_t natom, const double* x, double* f) const;

};
        
} // namespace pot

#endif // SHO_H
