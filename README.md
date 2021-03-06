[![Build Status](https://travis-ci.org/gmedders/pimd.svg?branch=master)](https://travis-ci.org/gmedders/pimd)

This repository contains codes to do classical and quantum (i.e., RPMD) molecular dynamics simulations. Either of these MD schemes can be combined with a surface hopping-based classical master equation to model the nonadiabatic dynamics of a molecule interacting with a manifold of electronic states.

Note, this work was intended to test some ideas about how nonadiabatic effects (i.e., when the electronic character of the system is strongly coupled with the nuclear motion) might be impacted by an (approximate) quantum treatment of the nuclei themselves. Ultimately, we decided to invest energy in other projects, so these codes are provided with the hope that they might be helpful others but with no guarantee of their accuracy. (And, more importantly, we are not suggesting that combining the CME with RPMD is even a good idea!).

Finally, as the codes were exploratory, they were not optimized for performance considerations.

How to build:
=============

 - git clone https://github.com/gmedders/pimd.git
 - cd pimd
 - mkdir build
 - cd build
 - cmake ..
 - cmake --build .
 - ctest -VV

To override your default C and C++ compilers, you can use the following instead (for example):
```
cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc-8 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-8 ..
```

# Dependencies

These implementations depend on the FFTW3 and armadillo libraries. On MacOS, you can get these with
```
brew install armadillo fftw
```

Usage:
======

The workflow involves generating a set of initial conditions (positions, velocities, and electronic state IDs) consistent with a given temperature. This is done, for example, by `rpmd.cpp`

```
# usage: rpmd nbeads beta dt GammaEl gammaTh_factor voltage nsamples
./rpmd 1 1024 1.0 1.0E-3 2 0.0 2 > log.dat
```

which stores the requested initial conditions in the file `cart_traj-rpmd_1_1024.dat` and the simulation statistics in (e.g., potential and kinetic energy) in `log.dat`. In the command line input:
- `nbeads` is the number of beads in the ring-polymer molecular dynamics (1 bead corresponds to classical dynamics)
- `beta` is the inverse temperature (1024 is roughly 300 Kelvin)
- `dt` is the timestep
- `GammaEl` is the electronic coupling (see reference below)
- `gammaTh_factor` is a scaling factor the thermal friction for Langevin dynamics (the friction is defined as `gamma = gammaTh_factor * omega`, where omega is the harmonic frequency of the example potential)
- `voltage`, if non-zero, corresponds to CME dynamics with multiple electronic baths
- `nsamples` is the number of initial conditions requested

The time-dependent ensemble average is performed using `ensemble_tcf_rpmd.cpp`. This requires the initial conditions generated by `rpmd.cpp` and similar simulation parameters.

```
# usage: ensemble_tcf_rpmd input_file dt GammaEl gammaTh_fac voltage time
./ensemble_tcf_rpmd cart_traj-rpmd_1_1024.dat 1.0 1.0E-3 2.0 0.0 50000 > ens.dat
```

Several example potentials are available and can be selected in `include/sim-classes.h`
