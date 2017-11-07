#ifndef LJ_H
#define LJ_H

#include <cstdlib>

namespace pot {

class lj
{

public:
    double V2B(const size_t ndim, const double* crd_a, const double* crd_b,
               double* f_a, double* f_b);

    //void set_params(double*);

    void print_params();

private:
    // sig = 0.3345 nm = 6.323 au
    // eps = 125.7 K = 0.000398066 au
    static constexpr double sig = 6.323;
    static constexpr double eps = 0.000398066;
};

} // namespace pot

#endif // LJ_H
