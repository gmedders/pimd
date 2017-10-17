#include <cassert>

#include "explicit-bath.h"

////////////////////////////////////////////////////////////////////////////////

namespace pot {

////////////////////////////////////////////////////////////////////////////////

void explicit_bath::set_bath_params(size_t ndim, size_t nBathModes,
                                    double gamma, double cutoff,
                                    double thermo_mass)
{
    assert(ndim == 1);

    m_nBathModes = nBathModes;    
    m_mass = thermo_mass;

    if(m_nBathModes > 0){
        c_iMode = new double[m_nBathModes];
        omega_iMode = new double[m_nBathModes];

        double omega_max = 5.0*cutoff;

        double domega = omega_max/nBathModes;

        for(int i = 0; i < nBathModes; ++i){
            double omega = (i + 0.5)*domega;
            double J = gamma*omega*std::exp(-omega/cutoff);

            omega_iMode[i] = omega;
            c_iMode[i] = std::sqrt(gamma*omega*omega*m_mass);
        }
    }
};

//----------------------------------------------------------------------------//

double explicit_bath::bath_force(const double* crd, double* frc)
{
    if(m_nBathModes == 0)
        return 0.0;

    double energy(0);

    double q = crd[0];
    for(size_t i = 1; i < m_nBathModes; ++i){
        double x = crd[i];
        double c = c_iMode[i];
        double mw2 = m_mass*omega_iMode[i]*omega_iMode[i];
        double c_over_mw2 = c / mw2;

        // [xi - q * ci/(m*omega^2)]^2
        // Foil and save explicitly
        double poly_xx = x*x;
        double poly_xq = - x*(c_over_mw2)*q;
        double poly_qq = c_over_mw2*c_over_mw2*q*q;

        double e = 0.5*mw2*(poly_xx + 2*poly_xq + poly_qq);

        double dx = mw2*(poly_xx + poly_xq)/x;
        double dq = mw2*(poly_xq + poly_qq)/q;

        energy += e;
        frc[0] -= dq;
        frc[i] -= dx;
    }
    return energy;
}

//----------------------------------------------------------------------------//

explicit_bath::explicit_bath()
{
    m_nBathModes = 0;

    c_iMode = nullptr;
    omega_iMode = nullptr;
}

//----------------------------------------------------------------------------//

explicit_bath::~explicit_bath()
{
    if(c_iMode)
        delete[] c_iMode;

    if(omega_iMode)
        delete[] omega_iMode;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace pot

////////////////////////////////////////////////////////////////////////////////
