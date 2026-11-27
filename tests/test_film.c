/* Validation: thin-film interference and multiple reflections (Phase 3, scene 6).
 * Reference results are closed-form optics formulas, not the code under test:
 *   - Fresnel normal-incidence coefficients and the pi shift on reflection
 *   - Airy closed form vs the explicit partial-wave sum, and R + T = 1
 *   - film in air (r12 = -r23): extrema where delta = m pi exactly, i.e.
 *       maxima  lambda_m = 4 n d / (2m+1),  minima  lambda_m = 2 n d / m
 *   - quarter-wave single-layer coating: R = ((n1 n3 - n2^2)/(n1 n3 + n2^2))^2
 * Limits: normal incidence, lossless, non-dispersive indices.
 */
#include <stdio.h>
#include "wr_test.h"
#include "wr_film.h"
#include "wr_optics.h"

static double refl(double n1, double n2, double n3, double d, double lam)
{
    return wr_layer_reflectance(wr_layer_airy(n1, n2, n3, d, lam));
}

/* extremum of R(lambda) near l0 by scan + parabolic refinement */
static double extremum_lambda(double n1, double n2, double n3, double d, double l0, int want_max)
{
    const double step = 0.1e-9, half = 25e-9;
    double best = l0, bv = want_max ? -1.0 : 2.0;
    for (double l = l0 - half; l <= l0 + half; l += step) {
        double v = refl(n1, n2, n3, d, l);
        if (want_max ? v > bv : v < bv) { bv = v; best = l; }
    }
    double a = refl(n1, n2, n3, d, best - step), c = refl(n1, n2, n3, d, best + step);
    double den = a - 2.0 * bv + c;
    return best + ((den != 0.0) ? 0.5 * (a - c) / den : 0.0) * step;
}

int main(void)
{
    /* Fresnel coefficients and phase shifts */
    CHECK_NEAR(wr_fresnel_r_normal(1.0, 1.5), -0.2, 1e-15);   /* air -> glass: pi shift */
    CHECK_NEAR(wr_fresnel_r_normal(1.5, 1.0), +0.2, 1e-15);   /* glass -> air: none */
    CHECK_NEAR(wr_fresnel_t_normal(1.0, 1.5), 0.8, 1e-15);
    for (double nt = 1.1; nt < 3.0; nt += 0.3) {              /* single interface: R + T = 1 */
        double r = wr_fresnel_r_normal(1.0, nt), t = wr_fresnel_t_normal(1.0, nt);
        CHECK_NEAR(r * r + nt * t * t, 1.0, 1e-14);
    }

    /* series vs closed form, and energy conservation, over a (lambda, d, n) grid */
    double worst_ser = 0, worst_e = 0;
    const double ns[4][3] = {{1.0, 1.33, 1.0}, {1.0, 1.38, 1.52}, {1.5, 1.0, 1.5}, {1.0, 2.4, 1.5}};
    for (int c = 0; c < 4; c++)
        for (int il = 0; il < 40; il++)
            for (int id = 0; id < 40; id++) {
                double lam = (380.0 + 10.0 * il) * 1e-9, d = (20.0 + 37.0 * id) * 1e-9;
                wr_rt s = wr_layer_series(ns[c][0], ns[c][1], ns[c][2], d, lam, 200);
                wr_rt a = wr_layer_airy(ns[c][0], ns[c][1], ns[c][2], d, lam);
                worst_ser = fmax(worst_ser, wr_cabs(wr_csub(s.r, a.r)) + wr_cabs(wr_csub(s.t, a.t)));
                worst_e = fmax(worst_e, fabs(wr_layer_reflectance(a) + wr_layer_transmittance(a, ns[c][0], ns[c][2]) - 1.0));
            }
    printf("series (200 partial waves) vs Airy closed form: max difference %.2e\n", worst_ser);
    printf("energy: max |R + T - 1| over 6400 films = %.2e\n", worst_e);
    CHECK_NEAR(worst_ser, 0.0, 1e-12);
    CHECK_NEAR(worst_e, 0.0, 1e-13);

    /* truncation: few partial waves are only an approximation, converging geometrically */
    {
        wr_rt a = wr_layer_airy(1.0, 2.4, 1.5, 100e-9, 550e-9);
        double e1 = wr_cabs(wr_csub(wr_layer_series(1.0, 2.4, 1.5, 100e-9, 550e-9, 1).r, a.r));
        double e5 = wr_cabs(wr_csub(wr_layer_series(1.0, 2.4, 1.5, 100e-9, 550e-9, 5).r, a.r));
        printf("truncation (|r12 r23| = 0.2): 1 partial wave off by %.2e, 5 by %.2e\n", e1, e5);
        CHECK(e5 < 0.05 * e1);
    }

    /* soap film in air: n = 1.33, d = 400 nm */
    {
        const double n = 1.33, d = 400e-9;
        printf("soap film (n = 1.33, d = 400 nm): extrema vs closed form\n  type   lambda (nm)  expected (nm)\n");
        const struct { int max; double exp; } ex[3] = {
            {1, 4.0 * n * d / 5.0}, {0, 2.0 * n * d / 2.0}, {1, 4.0 * n * d / 3.0}};
        for (int i = 0; i < 3; i++) {
            double got = extremum_lambda(1.0, n, 1.0, d, ex[i].exp, ex[i].max);
            printf("  %s   %9.3f   %9.3f\n", ex[i].max ? "max" : "min", got * 1e9, ex[i].exp * 1e9);
            CHECK_NEAR(got, ex[i].exp, 0.05e-9);
        }
        double r12 = wr_fresnel_r_normal(1.0, n);
        double rmax = 4.0 * r12 * r12 / ((1.0 + r12 * r12) * (1.0 + r12 * r12));
        printf("  peak reflectance %.5f (closed form %.5f, two-beam estimate 4 r^2 = %.5f)\n",
               refl(1.0, n, 1.0, d, 4.0 * n * d / 3.0), rmax, 4.0 * r12 * r12);
        CHECK_NEAR(refl(1.0, n, 1.0, d, 4.0 * n * d / 3.0), rmax, 1e-12);
    }

    /* quarter-wave anti-reflection coating (MgF2 on glass) at the design wavelength */
    {
        const double n1 = 1.0, n2 = 1.38, n3 = 1.52, lam = 550e-9, d = lam / (4.0 * n2);
        double want = pow((n1 * n3 - n2 * n2) / (n1 * n3 + n2 * n2), 2.0);
        double bare = pow(wr_fresnel_r_normal(n1, n3), 2.0);
        printf("quarter-wave coating: R = %.5f%% (formula %.5f%%), uncoated glass %.3f%%\n",
               100 * refl(n1, n2, n3, d, lam), 100 * want, 100 * bare);
        CHECK_NEAR(refl(n1, n2, n3, d, lam), want, 1e-13);
        CHECK(refl(n1, n2, n3, d, lam) < 0.35 * bare);
        CHECK(refl(n1, n2, n3, d, lam) < refl(n1, n2, n3, d, 1.5 * lam)); /* worse away from design wavelength */
    }
    TEST_REPORT("film");
}
