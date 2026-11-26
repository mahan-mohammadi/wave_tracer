/* wr_propagate.h - scalar angular-spectrum propagation (Phase 2).
 *
 * Model: for a monochromatic scalar field U(x,y,0) in a homogeneous medium of
 * refractive index n (wavelength in the medium lambda = lambda0/n, k = 2 pi/lambda),
 *   U(x,y,z) = F^-1{ F{U(.,.,0)} * H },
 *   H(fx,fy) = exp(i 2 pi z sqrt(1/lambda^2 - fx^2 - fy^2))   propagating, lambda^2 (fx^2+fy^2) <= 1
 *   H(fx,fy) = exp(-2 pi z sqrt(fx^2 + fy^2 - 1/lambda^2))    evanescent (z >= 0 only)
 * This solves the scalar Helmholtz equation exactly for the sampled field; no
 * paraxial approximation is made.
 *
 * Assumptions / limits (documented, not hidden):
 *  - Scalar field: no polarization coupling; poor for features ~ wavelength.
 *  - The grid is PERIODIC: light diffracted past the window edge wraps around.
 *    The window must hold the whole pattern, ~ z*tan(theta_max) wider than the
 *    aperture. Check the result near the edges before trusting it.
 *  - Frequencies are sampled at 1/(n*d); H must be well sampled, which limits
 *    z for a given window (same wrap-around condition seen in frequency space).
 *  - Evanescent components with z < 0 would grow exponentially; they are set
 *    to zero instead, so backward propagation is stable but lossy for them.
 *  - Sampling must resolve the field: d well below lambda/ (2 sin theta_max).
 * Reference for the method: Goodman, Introduction to Fourier Optics (angular
 * spectrum chapter); Born & Wolf.
 */
#ifndef WR_PROPAGATE_H
#define WR_PROPAGATE_H

#include "wr_field.h"
#include "wr_fft2d.h"

/* Propagate f in place by z (metres, may be negative). plan must match f's
 * size. lambda0 = vacuum wavelength, n = real refractive index (> 0).
 * Returns 0 on success, nonzero on invalid arguments. */
int wr_propagate_asm(wr_field *f, const wr_fft2_plan *plan,
                     double lambda0, double n, double z);

/* Band-limited variant. Components whose transfer-function phase changes by
 * more than pi between adjacent frequency samples are discarded:
 *   |f_axis| <= 1 / (lambda * sqrt(1 + (2 |z| / (N d))^2)).
 * These are exactly the components that would wrap around the periodic
 * window, so removing them avoids coherent wrap-around error at the price of
 * dropping the (genuine) highest-angle light. Near-axis results are then
 * reliable for much larger z than the plain method. It is NOT energy
 * conserving or reversible; use the plain version when you need that.
 * The idea matches the band-limited angular spectrum method of Matsushima &
 * Shimobaba (Opt. Express 17, 2009; cited from memory, verify). The limit
 * formula here was derived directly from the phase-sampling condition. */
int wr_propagate_asm_bl(wr_field *f, const wr_fft2_plan *plan,
                        double lambda0, double n, double z);

#endif
