#ifndef RPMD_BASE_H
#define RPMD_BASE_H

#include "rpmd-necklace.h"

namespace parts {

struct rpmd_base : public rpmd_necklace {

    rpmd_base();
    ~rpmd_base();

    virtual double force(const double*, double*) = 0;

    void init(size_t ndof, size_t nbead,
              const double& kT, const double& dt,
              const double* mass, const double* cartpos, const double* cartvel);

    void step(const double&);
    double invariant() const;

    constexpr static size_t nchain = 4;

protected:
    double m_kT;
    double m_omega_M;

    arma::mat m_phys_mass;

    double m_Espring;
    double m_Epot_sum;
    double m_Ekin;
    double m_temp_kT;

    void pimd_force();
};

} // namespace parts

#endif // RPMD_BASE_H
