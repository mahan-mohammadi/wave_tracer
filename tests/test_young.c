/* Validation: Young's experiment with an extended incoherent source
 * (roadmap Phase 3, scene 3).
 *
 * Scene: a line source of width s (uniform along y) at distance z0 in front of
 * two narrow slits (width a, separation d); screen at distance z behind them.
 * The source is a set of independent point (line) emitters:
 *   - each emitter is a 1-pixel slit; its field is propagated to the slit plane by
 *     the angular-spectrum engine (source -> slits), masked by the two slits, then
 *     propagated to the screen (slits -> screen);
 *   - AMPLITUDES are added within one emitter (the two slits), INTENSITIES are
 *     added across emitters (mutually incoherent).
 * Theory (van Cittert-Zernike, narrow slits): fringe visibility
 *     V = |gamma| = |sinc(u)|,  u = s d / (lambda z0),  sinc(u) = sin(pi u)/(pi u),
 * and V changes sign (contrast reversal, fringes shifted by half a period) for
 * 1 < u < 2. Emitters are placed on a regular grid, 50 um apart (midpoint rule).
 * Signed visibility is measured as (I(0) - I(lambda z / 2 d)) / (I(0) + I(lambda z / 2 d)).
 * Independent reference: Fresnel integral by numerical quadrature over both
 * slits with the source's spherical-wave illumination, summed over the same emitters.
 * Window: the 1-pixel emitter has a flat spectrum, so the band-limit cutoff matters:
 * with a 32.8 mm window |V_sim - V_ref| is 2.9e-3, with 65.5 mm it is 1.0e-4 and a
 * 131 mm window gives the same, so 65.5 mm is used.
 * Limits: scalar, paraxial reference; finite slit width a makes the exact V differ
 * from sinc(u) by O(a/d)^2-level terms (the reference includes it exactly).
 */
#include <stdio.h>
#include <stdlib.h>
#include "wr_test.h"
#include "wr_propagate.h"
#include "wr_integrate.h"
#include "wr_optics.h"

#ifndef NX
#define NX 65536   /* 65.5 mm window; 32.8 mm leaves 3e-3 band-limit error (see header) */
#endif
#define NSRC 71           /* emitters j = -35 .. 35 */
#define NOBS 9            /* x = k * 0.25 mm, k = -4 .. 4 */

static const double LAM = 500e-9, A = 6e-6, D = 100e-6, Z0 = 0.5, Z = 0.2, DX = 1e-6, DSRC = 50e-6;

typedef struct { double x, xs; } rc;

static wr_cplx integrand(double xp, void *c)
{
    const rc *r = c;
    double t = r->x - xp, u = xp - r->xs;
    return wr_cexp_i(WR_PI * t * t / (LAM * Z) + WR_PI * u * u / (LAM * Z0));
}

static double ref_intensity(double x, double xs)
{
    rc r = {x, xs};
    wr_cplx s1 = wr_simpson_c(integrand, &r, -0.5 * D - 0.5 * A, -0.5 * D + 0.5 * A, 200);
    wr_cplx s2 = wr_simpson_c(integrand, &r,  0.5 * D - 0.5 * A,  0.5 * D + 0.5 * A, 200);
    return wr_cabs2(wr_cadd(s1, s2));
}

static double sinc(double u) { return (u == 0.0) ? 1.0 : sin(WR_PI * u) / (WR_PI * u); }

int main(void)
{
    wr_fft2_plan *plan = wr_fft2_plan_create(NX, 1);
    wr_field *leg1 = wr_field_create(NX, 1, DX, 1e-6), *f = wr_field_create(NX, 1, DX, 1e-6);
    const size_t c0 = NX / 2;                         /* pixel at x = +0.5 um (the unshifted emitter) */
    *wr_field_at(leg1, c0, 0) = wr_c(1.0, 0.0);
    CHECK(wr_propagate_asm_bl(leg1, plan, LAM, 1.0, Z0) == 0);   /* source -> slit plane, once */

    double isim[NSRC][NOBS], iref[NSRC][NOBS], xs[NSRC], xo[NOBS];
    size_t io[NOBS];
    for (int k = 0; k < NOBS; k++) {
        io[k] = (size_t)((double)c0 + (k - 4) * 0.25e-3 / DX);
        xo[k] = wr_field_x(leg1, io[k]);
    }
    for (int j = 0; j < NSRC; j++) {
        long n = 50L * (j - 35);                      /* integer pixel shift: translation invariance is exact */
        xs[j] = wr_field_x(leg1, c0) + (double)n * DX;
        for (size_t ix = 0; ix < NX; ix++) {
            size_t src = (size_t)(((long)ix - n) % (long)NX + (long)NX) % NX;
            double x = wr_field_x(f, ix);
            int in = (fabs(x - 0.5 * D) < 0.5 * A) || (fabs(x + 0.5 * D) < 0.5 * A);
            f->data[ix] = in ? leg1->data[src] : wr_c(0.0, 0.0);
        }
        CHECK(wr_propagate_asm_bl(f, plan, LAM, 1.0, Z) == 0);   /* slits -> screen */
        for (int k = 0; k < NOBS; k++) {
            isim[j][k] = wr_cabs2(f->data[io[k]]);
            iref[j][k] = ref_intensity(xo[k], xs[j]);
        }
    }

    const double fringe = LAM * Z / D;                 /* 1 mm */
    const int widths[] = {1, 11, 25, 41, 51, 71};      /* emitters in use: 1 = point source */
    printf("Young: a = %.0f um, d = %.0f um, lambda = %.0f nm, z0 = %.0f mm, z = %.0f mm, fringe %.2f mm\n",
           A * 1e6, D * 1e6, LAM * 1e9, Z0 * 1e3, Z * 1e3, fringe * 1e3);
    printf("  s (mm)    u     V simulated   V reference   V = sinc(u)\n");
    double worst_ref = 0, worst_sinc = 0, worst_pat = 0, v_rev = 0, v_zero = 1;
    for (int w = 0; w < 6; w++) {
        int n = widths[w], lo = (NSRC - n) / 2;
        double s = n * DSRC, u = s * D / (LAM * Z0);
        double sim[NOBS] = {0}, ref[NOBS] = {0};
        for (int j = lo; j < lo + n; j++)
            for (int k = 0; k < NOBS; k++) { sim[k] += isim[j][k]; ref[k] += iref[j][k]; }
        double vs = (sim[4] - sim[6]) / (sim[4] + sim[6]), vr = (ref[4] - ref[6]) / (ref[4] + ref[6]);
        double an = sinc(u);
        printf("  %5.2f  %5.3f   %+8.4f      %+8.4f      %+8.4f\n", s * 1e3, u, vs, vr, an);
        worst_ref = fmax(worst_ref, fabs(vs - vr));
        if (n > 1) worst_sinc = fmax(worst_sinc, fabs(vs - an));
        for (int k = 0; k < NOBS; k++) worst_pat = fmax(worst_pat, fabs(sim[k] / sim[4] - ref[k] / ref[4]));
        if (n == 71) v_rev = vs;
        if (n == 51) v_zero = fabs(vs);
    }
    printf("  max |V_sim - V_ref| = %.2e; max |V_sim - sinc(u)| = %.2e; pattern (I/I(0)) max difference = %.2e\n",
           worst_ref, worst_sinc, worst_pat);
    CHECK_NEAR(worst_ref, 0.0, 5e-4);
    CHECK_NEAR(worst_sinc, 0.0, 5e-3);   /* ~2e-3 is the finite-slit-width effect (a/d = 6%) */
    CHECK_NEAR(worst_pat, 0.0, 2e-3);
    CHECK(v_zero < 0.03);              /* first zero of sinc: fringes wash out */
    CHECK(v_rev < -0.15);              /* contrast reversal beyond the first zero */

    wr_field_destroy(leg1); wr_field_destroy(f); wr_fft2_plan_destroy(plan);
    TEST_REPORT("young");
}
