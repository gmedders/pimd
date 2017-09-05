#ifndef NECKLACE_H
#define NECKLACE_H

#include <cassert>
#include <cstddef>

namespace parts {

//
// rpmd_necklace for RPMD; manages memory for cartesian/normal mode
// positions, velocities, forces and handles the transforms
//

struct rpmd_necklace {

    rpmd_necklace();
    ~rpmd_necklace();

    inline size_t ndofs() const;
    inline size_t nbeads() const;

    inline const double& lambda(size_t) const;

private:

    rpmd_necklace(const rpmd_necklace&);
    rpmd_necklace& operator=(const rpmd_necklace&);

    size_t m_ndofs;
    size_t m_nbeads;

protected:

    void setup(size_t ndof, size_t nbead, double dt);

    void pos_c2n();
    void pos_n2c();

    void vel_c2n();
    void vel_n2c();

    void frc_c2n();

    // layout is bead1, bead2, ..., beadN, where
    // each bead consists of ndofs elements

    double* m_pos_cart;
    double* m_pos_nmode;

    double* m_vel_cart;
    double* m_vel_nmode;

    double* m_frc_cart;

    double* m_cart_to_nm;
    double* m_freerp_propagator;

    double* m_lambda;
};

inline size_t rpmd_necklace::ndofs() const
{
    return m_ndofs;
}

inline size_t rpmd_necklace::nbeads() const
{
    return m_nbeads;
}

inline const double& rpmd_necklace::lambda(size_t b) const
{
    assert(b < nbeads());

    return m_lambda[b];
}

} // namespace parts

#endif // NECKLACE_H
