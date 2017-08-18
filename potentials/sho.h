#ifndef SHO_H
#define SHO_H

namespace pot {

struct sho {

    double w = 1.0;
    double m = 1.0;

    // a = 1/2 m w^2
    double a = 0.5 * m * w * w;

    double operator()(const double x, double& f) const;

};
        
} // namespace pot

#endif // SHO_H
