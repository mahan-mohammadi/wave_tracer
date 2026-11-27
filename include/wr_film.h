/* wr_film.h - Fresnel coefficients and single-layer (thin-film) reflection (Phase 3).
 *
 * Scalar, NORMAL incidence, lossless real refractive indices, non-dispersive
 * (n does not depend on wavelength here; dispersion arrives with Phase 6).
 * Amplitude (field) coefficients, interface from medium i into medium t:
 *     r = (n_i - n_t) / (n_i + n_t),   t = 2 n_i / (n_i + n_t).
 * The pi phase shift on reflection from a denser medium is the sign of r.
 *
 * Layer: medium n1 | film n2, thickness d | medium n3, light incident from n1.
 * With e = exp(i delta), delta = 2 k0 n2 d (round-trip phase), the partial
 * waves are summed coherently:
 *     r = r12 + t12 t21 r23 e * sum_j (r21 r23 e)^j
 *     t = t12 t23 exp(i delta/2) * sum_j (r21 r23 e)^j
 * and the closed form (Airy) is
 *     r = (r12 + r23 e) / (1 + r12 r23 e),  t = t12 t23 exp(i delta/2) / (1 + r12 r23 e).
 * Power: R = |r|^2, T = (n3/n1) |t|^2, and R + T = 1 for lossless layers.
 * Limits: normal incidence only (no oblique angles / polarization yet), no
 * absorption, flat parallel interfaces (slowly varying thickness is fine).
 */
#ifndef WR_FILM_H
#define WR_FILM_H

#include "wr_complex.h"

typedef struct { wr_cplx r, t; } wr_rt;

double wr_fresnel_r_normal(double n_i, double n_t);
double wr_fresnel_t_normal(double n_i, double n_t);

/* Explicit sum of the first nterms partial waves (nterms >= 1). */
wr_rt wr_layer_series(double n1, double n2, double n3, double d, double lambda0, int nterms);
/* Closed-form (infinite-sum) result. */
wr_rt wr_layer_airy(double n1, double n2, double n3, double d, double lambda0);

double wr_layer_reflectance(wr_rt x);                        /* |r|^2 */
double wr_layer_transmittance(wr_rt x, double n1, double n3); /* (n3/n1) |t|^2 */

#endif
