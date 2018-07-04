#ifndef RPMD_NECKLACE_H
#define RPMD_NECKLACE_H

#include <cassert>
#include <cstddef>

#include <armadillo>

namespace parts {

//
// rpmd_necklace for RPMD; manages memory for cartesian/normal mode
// positions, momocities, forces and handles the transforms
//

struct rpmd_necklace {

    rpmd_necklace();

    inline size_t ndofs() const;
    inline size_t nbeads() const;

//    inline const double& lambda(size_t) const;

private:

    rpmd_necklace(const rpmd_necklace&);
    rpmd_necklace& operator=(const rpmd_necklace&);

    size_t m_ndofs;
    size_t m_nbeads;

protected:

    double m_dt;

    double m_beta;
    double m_beta_n;
    double m_omega_n;

    void setup(size_t ndof, size_t nbead, double, double, double);

    void pos_c2n();
    void pos_n2c();

    void mom_c2n();
    void mom_n2c();

    void frc_c2n();

    // layout is bead1, bead2, ..., beadN, where
    // each bead consists of ndofs elements

    arma::mat m_pos_cart;
    arma::mat m_pos_nmode;

    arma::mat m_mom_cart;
    arma::mat m_mom_nmode;

    arma::mat m_frc_cart;

    arma::mat m_cart_to_nm;
    arma::cube m_freerp_propagator;

    arma::vec m_omega_k;

    arma::vec m_mass;
    arma::vec m_sqrt_mass;
    
    constexpr static double engunit = 1.0; // conversion from internal units to kcal/mol
    constexpr static double kB = 1.0; // Boltzmann constant in internal units
    constexpr static double hbar = 1.0;

};

inline size_t rpmd_necklace::ndofs() const
{
    return m_ndofs;
}

inline size_t rpmd_necklace::nbeads() const
{
    return m_nbeads;
}

//inline const double& rpmd_necklace::lambda(size_t b) const
//{
//    assert(b < nbeads());
//
//    return m_lambda[b];
//}

} // namespace parts

#endif // RPMD_NECKLACE_H
