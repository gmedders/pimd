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
    ~rpmd_necklace();

    inline size_t ndofs() const;
    inline size_t nbeads() const;

//    inline const double& lambda(size_t) const;

private:

    rpmd_necklace(const rpmd_necklace&);
    rpmd_necklace& operator=(const rpmd_necklace&);

    size_t m_ndofs;
    size_t m_nbeads;

    double m_dt;

protected:

    void setup(size_t ndof, size_t nbead, double dt);

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
