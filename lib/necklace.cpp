#include <cassert>
#include <cmath>

#include "necklace.h"

namespace parts {

//----------------------------------------------------------------------------//

necklace::necklace() : m_ndim(0), m_natom(0), m_ndofs(0), m_nbeads(0) {}

//----------------------------------------------------------------------------//

necklace::~necklace() { setup(0, 0, 0); }

//----------------------------------------------------------------------------//

void necklace::setup(size_t ndim, size_t natom, size_t nbead) {
  size_t ndof = ndim * natom;
  if (ndof == m_ndofs && nbead == m_nbeads)
    return;

  if (m_nbeads > 0 || m_ndofs > 0) {
    assert(m_nbeads > 0 && m_ndofs > 0);

    fftw_destroy_plan(m_plan_cart2nmode);
    fftw_destroy_plan(m_plan_nmode2cart);

    fftw_free(m_pos_cart);
    fftw_free(m_pos_nmode);

    fftw_free(m_vel_cart);
    fftw_free(m_vel_nmode);

    fftw_free(m_frc_cart);
    fftw_free(m_frc_nmode);

    delete[] m_lambda;
  }

  if (ndof > 0 || nbead > 0) {
    assert(ndof > 0 && nbead > 0);
    assert(nbead % 2 == 0 || nbead == 1);

    m_ndim = ndim;
    m_natom = natom;
    m_ndofs = ndof;
    m_nbeads = nbead;

    const size_t nbytes = sizeof(double) * nbead * ndof;

    m_pos_cart = (double *)fftw_malloc(nbytes);
    m_pos_nmode = (double *)fftw_malloc(nbytes);

    m_vel_cart = (double *)fftw_malloc(nbytes);
    m_vel_nmode = (double *)fftw_malloc(nbytes);

    m_frc_cart = (double *)fftw_malloc(nbytes);
    m_frc_nmode = (double *)fftw_malloc(nbytes);

    {
      const int n = (int)nbead;
      const int howmany = (int)ndof;

      const fftw_r2r_kind c2n_kind = FFTW_R2HC;
      const fftw_r2r_kind n2c_kind = FFTW_HC2R;

      m_plan_cart2nmode = fftw_plan_many_r2r(
          1, &n, howmany, m_pos_cart, 0, howmany, 1, m_pos_nmode, 0, howmany, 1,
          &c2n_kind, FFTW_MEASURE);

      m_plan_nmode2cart = fftw_plan_many_r2r(
          1, &n, howmany, m_pos_nmode, 0, howmany, 1, m_pos_cart, 0, howmany, 1,
          &n2c_kind, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
    }

    m_lambda = new double[nbead];

    // nmodes are "scrambled" by fftw (comparing to Francesco's codes)

    m_lambda[nbead / 2] = 4.0 * nbead;
    m_lambda[0] = 0.0;

    for (size_t n = 1; n < nbead / 2; ++n)
      m_lambda[nbead - n] = m_lambda[n] =
          4 * nbead * (1.0 - std::cos(2 * n * (M_PI / nbead)));
  }
}

//----------------------------------------------------------------------------//

void necklace::pos_c2n() {
  assert(m_nbeads > 0 && m_ndofs > 0);

  fftw_execute_r2r(m_plan_cart2nmode, m_pos_cart, m_pos_nmode);

  const double factor = 1.0 / m_nbeads;
  for (size_t n = 0; n < m_nbeads * m_ndofs; ++n)
    m_pos_nmode[n] *= factor;
}

//----------------------------------------------------------------------------//

void necklace::pos_n2c() {
  assert(m_nbeads > 0 && m_ndofs > 0);

  fftw_execute_r2r(m_plan_nmode2cart, m_pos_nmode, m_pos_cart);
}

//----------------------------------------------------------------------------//

void necklace::vel_n2c() {
  assert(m_nbeads > 0 && m_ndofs > 0);

  fftw_execute_r2r(m_plan_nmode2cart, m_vel_nmode, m_vel_cart);
}

//----------------------------------------------------------------------------//

void necklace::vel_c2n() {
  assert(m_nbeads > 0 && m_ndofs > 0);

  fftw_execute_r2r(m_plan_cart2nmode, m_vel_cart, m_vel_nmode);

  const double factor = 1.0 / m_nbeads;
  for (size_t n = 0; n < m_nbeads * m_ndofs; ++n)
    m_vel_nmode[n] *= factor;
}

//----------------------------------------------------------------------------//

void necklace::frc_c2n() {
  assert(m_nbeads > 0 && m_ndofs > 0);
  assert(m_nbeads % 2 == 0 || m_nbeads == 1);

  fftw_execute_r2r(m_plan_cart2nmode, m_frc_cart, m_frc_nmode);

  double factor = 1.0 / m_nbeads;

  for (size_t i = 0; i < m_ndofs; ++i) {
    m_frc_nmode[i] *= factor;
    m_frc_nmode[i + m_ndofs * (m_nbeads / 2)] *= factor;
  }

  factor *= 2;

  for (size_t b = 1; b < m_nbeads / 2; ++b)
    for (size_t i = 0; i < m_ndofs; ++i) {
      m_frc_nmode[i + b * m_ndofs] *= factor;
      m_frc_nmode[i + (m_nbeads - b) * m_ndofs] *= factor;
    }
}

//----------------------------------------------------------------------------//

} // namespace parts

////////////////////////////////////////////////////////////////////////////////
