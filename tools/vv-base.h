#ifndef VV_H
#define VV_H

#include <armadillo>

namespace parts {

struct vv_base {

    inline size_t ndofs() const;

    virtual double force(const size_t, const size_t,
                         const double*, double*) = 0;

    void init(size_t ndof,
              const double& kT, const double& dt,
              const double* mass, const double* cartpos, const double* cartvel);

    void step(const double&);
    double invariant() const;

protected:
    double m_dt;
    size_t m_ndofs;

    arma::vec m_pos;
    arma::vec m_mom;
    arma::vec m_frc;
    arma::vec m_mass;

    double m_Epot;
    double m_Ekin;
    double m_temp_kT;

};

inline size_t vv_base::ndofs() const
{
    return m_ndofs;
}

} // namespace parts

#endif // VV_H
