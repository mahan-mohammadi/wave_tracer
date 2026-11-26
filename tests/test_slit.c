/* Validation: single-slit diffraction.
 *
 * Geometry: slit of width a, uniform along y (so the problem is exactly 1D
 * and the 2D code is run with a thin ny). Unit-amplitude plane wave in, field
 * U0 = 1 inside the slit, 0 outside, then propagate a distance z.
 *
 * Reference: Fresnel diffraction integral, evaluated by NUMERICAL quadrature
 * (wr_simpson_c), not by a closed form:
 *   U(x,z) = exp(ikz) / sqrt(i lam z) * Int_{-a/2}^{a/2} exp(i pi (x-x')^2/(lam z)) dx'
 *   I = |U|^2 = |Int ...|^2 / (lam z).
 * The Fresnel integral is the paraxial approximation of the exact scalar
 * result; its error grows like k z theta^4 / 8 (about 5e-3 rad at the edge of
 * the compared region), so differences below that level are not meaningful.
 * Compared region: |x| < 0.4 mm. Window is periodic: wrap-around of the outer
 * lobes contributes a small known error (tested via the window size).
 */
#include <stdio.h>
#include <stdlib.h>
#include "wr_test.h"
#include "wr_propagate.h"
#include "wr_integrate.h"
#include "wr_optics.h"

typedef struct { double x, lam, z; } ref_ctx;

static wr_cplx fresnel_integrand(double xp, void *c)
{
    const ref_ctx *r = c;
    double dx = r->x - xp;
    return wr_cexp_i(WR_PI * dx * dx / (r->lam * r->z));
}

static double ref_intensity(double x, double a, double lam, double z)
{
    ref_ctx r = {x, lam, z};
    wr_cplx s = wr_simpson_c(fresnel_integrand, &r, -0.5 * a, 0.5 * a, 4000);
    return wr_cabs2(s) / (lam * z);
}

static int use_bl = 0;
static int propagate_fn(wr_field *f, const wr_fft2_plan *p, double lam, double n, double z)
{ return use_bl ? wr_propagate_asm_bl(f, p, lam, n, z) : wr_propagate_asm(f, p, lam, n, z); }

/* returns max |I_asm - I_ref| over |x| < xmax (incident intensity = 1) */
static double run(size_t nx, double dx, double a, double lam, double z, double xmax, double *peak)
{
    const size_t ny = 4;
    wr_fft2_plan *plan = wr_fft2_plan_create(nx, ny);
    wr_field *f = wr_field_create(nx, ny, dx, 1e-6);
    for (size_t iy = 0; iy < ny; iy++)
        for (size_t ix = 0; ix < nx; ix++)
            if (fabs(wr_field_x(f, ix)) < 0.5 * a) *wr_field_at(f, ix, iy) = wr_c(1, 0);
    if (propagate_fn(f, plan, lam, 1.0, z)) { fprintf(stderr, "propagate failed\n"); exit(2); }

    double err = 0, pk = 0;
    for (size_t ix = 0; ix < nx; ix++) {
        double x = wr_field_x(f, ix);
        if (fabs(x) >= xmax) continue;
        double ia = wr_cabs2(*wr_field_at(f, ix, 0));
        double ir = ref_intensity(x, a, lam, z);
        err = fmax(err, fabs(ia - ir));
        pk = fmax(pk, ir);
    }
    if (peak) *peak = pk;
    wr_field_destroy(f); wr_fft2_plan_destroy(plan);
    return err;
}

int main(void)
{
    const double lam = 500e-9, a = 40e-6, z = 20e-3, xmax = 0.4e-3;
    const double dxs[] = {2e-6, 1e-6, 0.5e-6, 0.25e-6};
    double e_plain[4], e_bl[4], peak = 0;

    printf("single slit: a = %.0f um, lambda = %.0f nm, z = %.0f mm, compared on |x| < %.1f mm\n",
           a * 1e6, lam * 1e9, z * 1e3, xmax * 1e3);

    /* (1) small window 2.048 mm: plain method suffers coherent wrap-around */
    printf("window 2.048 mm      max |dI| (incident I = 1)\n  dx (um)   plain       band-limited\n");
    for (int i = 0; i < 4; i++) {
        size_t nx = (size_t)(2.048e-3 / dxs[i] + 0.5);
        use_bl = 0; e_plain[i] = run(nx, dxs[i], a, lam, z, xmax, &peak);
        use_bl = 1; e_bl[i]    = run(nx, dxs[i], a, lam, z, xmax, NULL);
        printf("  %7.3f   %.3e   %.3e\n", dxs[i] * 1e6, e_plain[i], e_bl[i]);
    }
    printf("reference peak intensity in region: %.4f (a^2/(lam z) = %.4f)\n", peak, a * a / (lam * z));

    /* (2) plain method, window (16.384 mm) large enough to hold the pattern:
     *     this isolates the discretisation error of the hard slit edges. */
    use_bl = 0;
    const double dbig[] = {4e-6, 2e-6, 1e-6, 0.5e-6};
    double e_big[4];
    printf("plain, window 16.384 mm:\n  dx (um)   max |dI|\n");
    for (int i = 0; i < 4; i++) {
        e_big[i] = run((size_t)(16.384e-3 / dbig[i] + 0.5), dbig[i], a, lam, z, xmax, NULL);
        printf("  %7.3f   %.3e\n", dbig[i] * 1e6, e_big[i]);
    }

    CHECK_NEAR(peak, a * a / (lam * z), 2e-3);          /* reference equals the expected on-axis value */
    CHECK(e_bl[3] < e_plain[3]);                        /* band limiting removes the wrap-around error */
    CHECK(e_bl[3] < 1e-3);                              /* < ~0.6% of the 0.16 peak */
    CHECK(e_big[1] < e_big[0] && e_big[2] < e_big[1] && e_big[3] < e_big[2]); /* monotone convergence */
    CHECK(e_big[3] < 1e-4);                             /* < ~0.06% of the peak */
    TEST_REPORT("slit");
}
