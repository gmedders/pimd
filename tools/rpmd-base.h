#ifndef RPMD_BASE_H
#define RPMD_BASE_H

#include "necklace.h"

namespace parts {

struct rpmd_base : public necklace {

    rpmd_base();
    ~rpmd_base();

    virtual double force(const double*, double*) = 0;

    void init(size_t ndof, size_t nbead, const double& kT,
              const double* mass, const double* cartpos, const double* cartvel);

    void step(const double&);
    double invariant() const;

    constexpr static size_t nchain = 4;

    //constexpr static double engunit = 418.4; // conversion from internal units to kcal/mol
    //constexpr static double kB = 8.31451115e-01; // Boltzmann constant in internal units
    //constexpr static double hbar = 6.350780668;
    constexpr static double engunit = 1.0; // conversion from internal units to kcal/mol
    constexpr static double kB = 1.0; // Boltzmann constant in internal units
    constexpr static double hbar = 1.0;

protected:
    double m_kT;
    double m_omega_M;

    double* m_phys_mass;

    double m_Espring;
    double m_Epot_sum;
    double m_Ekin;
    double m_temp_kT;

    void pimd_force();
};

} // namespace parts

#endif // RPMD_BASE_H
