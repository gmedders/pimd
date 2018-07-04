#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>
#include <cstdlib>

#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "rand_gauss.h"

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // 1. load the coordinates

    std::cout.setf(std::ios_base::showpoint);
    std::cout.precision(9);

    int nsamples(10000);
    do{

        std::cout << parts::randn(0,1) << std::endl;

        --nsamples;
    }while(nsamples > 0);

    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
