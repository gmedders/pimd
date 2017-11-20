#ifndef LJ_H
#define LJ_H

#include <cstdlib>
#include <iostream>

namespace pot {

class lj
{

public:
    //lj() {std::cerr << "lj::lj()" << std::endl; };
    //~lj() {std::cerr << "lj::~lj()" << std::endl; };
    double V2B(const size_t ndim, const double* crd_a, const double* crd_b,
               double* f_a, double* f_b);

    //void set_params(double*);

    void print_params();

private:
    // sig = 0.2782 nm = 5.25 au
    // eps = 125.7 K = 0.000398066 au
    static constexpr double sig = 5.25;
    static constexpr double eps = 0.0001180938915;
};

} // namespace pot

#endif // LJ_H
