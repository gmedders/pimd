#ifndef PIMD_BASE_H
#define PIMD_BASE_H

#include "necklace.h"

namespace parts {

struct pimd_base : public necklace {

    pimd_base();
    ~pimd_base();

    virtual double force(const double*, double*) = 0;

    void init(size_t ndof, size_t nbead, const double& kT,
              const double* mass, const double* cartpos); // 1 bead

    void step(const double&);
    double invariant() const;

    constexpr static size_t nchain = 4;

    constexpr static double engunit = 418.4; // conversion from internal units to kcal/mol
    constexpr static double kB = 8.31451115e-01; // Boltzmann constant in internal units
    constexpr static double hbar = 6.350780668;

protected:
    double m_kT;
    double m_omega_M;

    double* m_fict_mass;
    double* m_thermostats;

    double m_Espring;
    double m_Epot_sum;
    double m_Ekin_fict;

    void pimd_force();
};

} // namespace parts

#endif // PIMD_BASE_H
