/* wr_material.h - optical materials (Phase 4c).
 *
 * A material describes what is INSIDE a closed object (or both sides of a
 * sheet, for mirrors/absorbers). Description only: Fresnel coefficients, phase
 * jumps and transport are the job of the optics layers.
 *
 * Complex refractive index  n_c(lambda0) = n(lambda0) + i kappa,
 * lambda0 = vacuum wavelength in metres. Project time convention exp(-i w t):
 * a wave travelling L in the medium is multiplied by exp(i k0 n_c L), so
 * kappa > 0 is absorption and the intensity decays as exp(-alpha L) with
 * alpha = 4 pi kappa / lambda0 (absorption coefficient, 1/m).
 *
 * Index models (lambda in micrometres inside the formulas):
 *   constant : n = n0
 *   Cauchy   : n = A + B / lambda^2 + C / lambda^4              (empirical, visible/near-IR only)
 *   Sellmeier: n^2 = 1 + sum_i B_i lambda^2 / (lambda^2 - C_i)  (C_i in um^2, physically motivated)
 * kappa is constant (no absorption dispersion yet). Every material carries a
 * validity range [lambda_min, lambda_max] in metres; evaluation outside it
 * still returns a number but flags in_range = 0 (roadmap: document the validity
 * range of every approximation). Defaults: constant [0, inf), Cauchy
 * [0.4, 0.8] um, Sellmeier [0.3, 2.5] um; edit the fields to change them.
 *
 * Surface kinds:
 *   DIELECTRIC   : smooth interface to a medium of index n_c (metals are dielectrics with large kappa)
 *   IDEAL_MIRROR : perfect reflector (scalar amplitude reflection coefficient -1), no interior
 *   ABSORBER     : perfectly black, no interior
 * For mirror/absorber wr_material_index returns 1 + 0i, which has no physical meaning.
 */
#ifndef WR_MATERIAL_H
#define WR_MATERIAL_H

#include "wr_complex.h"

#define WR_NAME_LEN 48

typedef enum { WR_SURF_DIELECTRIC, WR_SURF_IDEAL_MIRROR, WR_SURF_ABSORBER } wr_surface_kind;
typedef enum { WR_IOR_CONSTANT, WR_IOR_CAUCHY, WR_IOR_SELLMEIER } wr_ior_model;

typedef struct {
    char            name[WR_NAME_LEN];
    wr_surface_kind surface;
    wr_ior_model    model;
    double          n0;            /* constant */
    double          cauchy[3];     /* A, B (um^2), C (um^4) */
    double          sell_b[3];     /* B1..B3 */
    double          sell_c[3];     /* C1..C3 (um^2) */
    double          kappa;         /* extinction coefficient, >= 0 */
    double          lambda_min, lambda_max;   /* validity range, metres */
} wr_material;

wr_material wr_material_constant(const char *name, double n, double kappa);
wr_material wr_material_cauchy(const char *name, double a, double b_um2, double c_um4);
wr_material wr_material_sellmeier(const char *name, const double b[3], const double c_um2[3]);
wr_material wr_material_mirror(const char *name);
wr_material wr_material_absorber(const char *name);

/* Complex index at vacuum wavelength lambda0 (m). *in_range (optional) is 1 if lambda0 lies inside the
 * validity range and the model is numerically meaningful (n^2 > 0), else 0 (n may then be NaN). */
wr_cplx wr_material_index(const wr_material *m, double lambda0, int *in_range);
/* Intensity absorption coefficient 4 pi kappa / lambda0, 1/m. */
double  wr_material_absorption_coeff(const wr_material *m, double lambda0);
/* Returns 1 if parameters are finite and consistent, else 0 with a message in err (may be NULL). */
int     wr_material_validate(const wr_material *m, char *err, size_t errlen);

#endif /* WR_MATERIAL_H */
