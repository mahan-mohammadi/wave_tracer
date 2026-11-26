/* Validation: double-slit interference (roadmap Phase 3, scene 2).
 *
 * The scene contains ONLY two apertures: the field is 1 inside two slits of
 * width a, centres at +-d/2, and 0 elsewhere. The fringes are not coded
 * anywhere; they must emerge from propagation plus linear superposition.
 *
 * Checks
 *  1. Superposition: propagating the two-slit field equals the coherent sum
 *     of the two separately propagated single-slit fields (linearity).
 *  2. Interference phase: arg(U1 conj(U2)) from the simulated complex fields
 *     has slope 2 pi d / (lambda z) (the textbook fringe period lambda*z/d).
 *     Intensity-maximum spacing is compared with the SAME measurement on the
 *     independent reference: the single-slit envelope slopes across the
 *     fringes and pulls the maxima of the product inward, so raw maxima are
 *     ~1% closer than lambda*z/d (both simulation and reference show this).
 *  3. Intensity vs the Fresnel diffraction integral over both slits, computed
 *     by numerical quadrature (independent code path).
 *  4. Fringe visibility vs the reference, and the incoherent sum
 *     |U1|^2 + |U2|^2 (what you get by wrongly adding intensities) shows
 *     almost no fringes.
 * Window: 32.8 mm. The plain method would need a window >> 2 z tan(theta) to
 * avoid wrap-around (it gives ~16% error here), so the band-limited variant
 * is used; with a small window (8.2 mm) its sharp spectral cutoff falls inside
 * the slit's main lobe and costs ~1.6% (printed below for comparison).
 * Limits: reference is paraxial (compared region |x| < 2.5 mm = 0.0125 rad).
 * Writes build/double_slit.csv (x, I_sim, I_ref, I_incoherent) for plotting.
 */
#include <stdio.h>
#include <stdlib.h>
#include "wr_test.h"
#include "wr_propagate.h"
#include "wr_integrate.h"
#include "wr_optics.h"

typedef struct { double x, lam, z, a, d; } ref_ctx;

static wr_cplx integrand(double xp, void *c)
{
    const ref_ctx *r = c;
    double t = r->x - xp;
    return wr_cexp_i(WR_PI * t * t / (r->lam * r->z));
}

/* coherent: amplitudes of both slits added, then squared */
static double ref_intensity(double x, double a, double d, double lam, double z)
{
    ref_ctx r = {x, lam, z, a, d};
    wr_cplx s1 = wr_simpson_c(integrand, &r, -0.5 * d - 0.5 * a, -0.5 * d + 0.5 * a, 1000);
    wr_cplx s2 = wr_simpson_c(integrand, &r,  0.5 * d - 0.5 * a,  0.5 * d + 0.5 * a, 1000);
    return wr_cabs2(wr_cadd(s1, s2)) / (lam * z);
}

static void fill(wr_field *f, double c1, double c2, double a)
{
    for (size_t iy = 0; iy < f->ny; iy++)
        for (size_t ix = 0; ix < f->nx; ix++) {
            double x = wr_field_x(f, ix);
            int in = (fabs(x - c1) < 0.5 * a) || (fabs(x - c2) < 0.5 * a);
            *wr_field_at(f, ix, iy) = wr_c(in ? 1.0 : 0.0, 0.0);
        }
}

/* index of the pixel nearest x0 */
static size_t pix(const wr_field *f, double x0)
{
    return (size_t)((x0 / f->dx) + (double)(f->nx / 2) + 1e-9);
}

/* position of the intensity maximum near x0 (+-half), parabolic refinement */
static double peak_near(const wr_field *f, double x0, double half)
{
    size_t lo = pix(f, x0 - half), hi = pix(f, x0 + half), best = lo;
    double bi = -1;
    for (size_t i = lo; i <= hi; i++) {
        double v = wr_cabs2(*wr_field_at((wr_field *)f, i, 0));
        if (v > bi) { bi = v; best = i; }
    }
    double l = wr_cabs2(*wr_field_at((wr_field *)f, best - 1, 0));
    double c = bi;
    double r = wr_cabs2(*wr_field_at((wr_field *)f, best + 1, 0));
    double den = l - 2 * c + r;
    double off = (den != 0.0) ? 0.5 * (l - r) / den : 0.0;
    return wr_field_x(f, best) + off * f->dx;
}


/* max |I_sim - I_ref| and reference peak over |x| < xmax, window nx*dx */
static double pattern_error(size_t nx, double dx, double a, double d, double lam, double z,
                            double xmax, double *peak_out)
{
    wr_fft2_plan *plan = wr_fft2_plan_create(nx, 4);
    wr_field *f = wr_field_create(nx, 4, dx, 1e-6);
    fill(f, -0.5 * d, 0.5 * d, a);
    wr_propagate_asm_bl(f, plan, lam, 1.0, z);
    double err = 0, pk = 0;
    for (size_t ix = 0; ix < nx; ix += 4) {
        double x = wr_field_x(f, ix);
        if (fabs(x) >= xmax) continue;
        double ir = ref_intensity(x, a, d, lam, z);
        err = fmax(err, fabs(wr_cabs2(*wr_field_at(f, ix, 0)) - ir));
        pk = fmax(pk, ir);
    }
    if (peak_out) *peak_out = pk;
    wr_field_destroy(f); wr_fft2_plan_destroy(plan);
    return err;
}

/* maximum of the REFERENCE intensity near x0 (pixel grid, parabolic refinement) */
static double ref_peak_near(double x0, double half, double dx, double a, double d, double lam, double z)
{
    double best = x0, bi = -1;
    for (double x = x0 - half; x <= x0 + half; x += dx) {
        double v = ref_intensity(x, a, d, lam, z);
        if (v > bi) { bi = v; best = x; }
    }
    double l = ref_intensity(best - dx, a, d, lam, z), r = ref_intensity(best + dx, a, d, lam, z);
    double den = l - 2 * bi + r;
    return best + ((den != 0.0) ? 0.5 * (l - r) / den : 0.0) * dx;
}

int main(void)
{
    const double lam = 500e-9, a = 20e-6, d = 100e-6, z = 0.2;
    const double dx = 0.5e-6, xmax = 2.5e-3;
    const size_t nx = 65536, ny = 4;
    const double fringe = lam * z / d; /* 1 mm */
    const double i_axis = 4.0 * a * a / (lam * z); /* Fraunhofer on-axis value, both slits in phase */

    wr_fft2_plan *plan = wr_fft2_plan_create(nx, ny);
    wr_field *both = wr_field_create(nx, ny, dx, 1e-6);
    wr_field *u1   = wr_field_create(nx, ny, dx, 1e-6);
    wr_field *u2   = wr_field_create(nx, ny, dx, 1e-6);
    fill(both, -0.5 * d, 0.5 * d, a);
    fill(u1, -0.5 * d, -0.5 * d, a);   /* slit 1 only */
    fill(u2,  0.5 * d,  0.5 * d, a);   /* slit 2 only */
    CHECK(wr_propagate_asm_bl(both, plan, lam, 1.0, z) == 0);
    CHECK(wr_propagate_asm_bl(u1, plan, lam, 1.0, z) == 0);
    CHECK(wr_propagate_asm_bl(u2, plan, lam, 1.0, z) == 0);

    printf("double slit: a = %.0f um, d = %.0f um, lambda = %.0f nm, z = %.0f mm, window %.1f mm\n",
           a * 1e6, d * 1e6, lam * 1e9, z * 1e3, nx * dx * 1e3);

    /* 1. linearity: two-slit field == U1 + U2 */
    double lin = 0, scale = 0;
    for (size_t i = 0; i < nx * ny; i++) {
        lin = fmax(lin, wr_cabs(wr_csub(both->data[i], wr_cadd(u1->data[i], u2->data[i]))));
        scale = fmax(scale, wr_cabs(both->data[i]));
    }
    printf("linearity  max|U_both - (U1+U2)| / max|U| = %.2e\n", lin / scale);
    CHECK_NEAR(lin / scale, 0.0, 1e-12);

    /* 2a. interference-phase slope from the complex fields, |x| < 2 mm */
    {
        double sx = 0, sp = 0, sxx = 0, sxp = 0, prev = 0, off = 0;
        long cnt = 0;
        for (size_t ix = 0; ix < nx; ix += 4) {
            double x = wr_field_x(both, ix);
            if (fabs(x) >= 2e-3) continue;
            double ph = wr_carg(wr_cmul(*wr_field_at(u1, ix, 0), wr_cconj(*wr_field_at(u2, ix, 0))));
            if (cnt > 0) { while (ph + off - prev > WR_PI) off -= WR_TWO_PI; while (ph + off - prev < -WR_PI) off += WR_TWO_PI; }
            ph += off; prev = ph; cnt++;
            sx += x; sp += ph; sxx += x * x; sxp += x * ph;
        }
        double slope = ((double)cnt * sxp - sx * sp) / ((double)cnt * sxx - sx * sx);
        double period = WR_TWO_PI / slope;
        printf("interference phase: period %.5f mm vs lambda z / d = %.5f mm (rel. diff %.2e)\n",
               period * 1e3, fringe * 1e3, fabs(period - fringe) / fringe);
        CHECK_NEAR(period / fringe, 1.0, 2e-3);
    }

    /* 2b. intensity-maximum spacing: simulation vs independent reference (same bias) */
    {
        double xs[5], xr[5];
        for (int m = -2; m <= 2; m++) {
            xs[m + 2] = peak_near(both, m * fringe, 0.3 * fringe);
            xr[m + 2] = ref_peak_near(m * fringe, 0.3 * fringe, dx, a, d, lam, z);
        }
        double sp_s = (xs[4] - xs[0]) / 4.0, sp_r = (xr[4] - xr[0]) / 4.0;
        printf("max-to-max spacing: simulated %.5f mm, reference %.5f mm (envelope pulls both below %.3f mm)\n",
               sp_s * 1e3, sp_r * 1e3, fringe * 1e3);
        CHECK_NEAR(sp_s / sp_r, 1.0, 1e-3);
    }

    /* 3. pointwise comparison with the Fresnel integral, every 4th pixel */
    FILE *csv = fopen("build/double_slit.csv", "w");
    if (csv) fprintf(csv, "x_mm,I_sim,I_ref,I_incoherent\n");
    double err = 0, peak = 0;
    for (size_t ix = 0; ix < nx; ix += 4) {
        double x = wr_field_x(both, ix);
        if (fabs(x) >= xmax) continue;
        double is = wr_cabs2(*wr_field_at(both, ix, 0));
        double ir = ref_intensity(x, a, d, lam, z);
        err = fmax(err, fabs(is - ir));
        peak = fmax(peak, ir);
        if (csv && ix % 32 == 0)
            fprintf(csv, "%.5f,%.6f,%.6f,%.6f\n", x * 1e3, is, ir,
                    wr_cabs2(*wr_field_at(u1, ix, 0)) + wr_cabs2(*wr_field_at(u2, ix, 0)));
    }
    if (csv) fclose(csv);
    printf("max |I_sim - I_ref| = %.3e  (reference peak %.4f = %.2f%% error; Fraunhofer 4a^2/(lam z) = %.4f)\n",
           err, peak, 100.0 * err / peak, i_axis);
    CHECK_NEAR(peak, i_axis, 2e-4);
    CHECK(err < 0.005 * peak);

    /* 4. visibility between the central maximum and the first minimum (x = lambda z / 2d) */
    {
        size_t ic = pix(both, 0.0), im = pix(both, 0.5 * fringe);
        double xc = wr_field_x(both, ic), xn = wr_field_x(both, im);
        double sc = wr_cabs2(*wr_field_at(both, ic, 0)), sn = wr_cabs2(*wr_field_at(both, im, 0));
        double rc = ref_intensity(xc, a, d, lam, z), rn = ref_intensity(xn, a, d, lam, z);
        double kc = wr_cabs2(*wr_field_at(u1, ic, 0)) + wr_cabs2(*wr_field_at(u2, ic, 0));
        double kn = wr_cabs2(*wr_field_at(u1, im, 0)) + wr_cabs2(*wr_field_at(u2, im, 0));
        double v_sim = (sc - sn) / (sc + sn), v_ref = (rc - rn) / (rc + rn), v_inc = (kc - kn) / (kc + kn);
        printf("visibility: simulated %.4f, reference %.4f, incoherent sum %.4f\n", v_sim, v_ref, v_inc);
        CHECK_NEAR(v_sim, v_ref, 5e-3);
        CHECK(v_sim > 0.95);
        CHECK(v_inc < 0.05);               /* adding intensities loses the fringes */
    }

    /* 5. window-size effect for the band-limited method (diagnostic) */
    double e_small = pattern_error(16384, dx, a, d, lam, z, xmax, NULL);
    printf("band-limited error: window 8.2 mm %.2e (%.2f%% of peak), window 32.8 mm %.2e (%.2f%%)\n",
           e_small, 100.0 * e_small / peak, err, 100.0 * err / peak);
    CHECK(err < e_small);

    wr_field_destroy(both); wr_field_destroy(u1); wr_field_destroy(u2); wr_fft2_plan_destroy(plan);
    TEST_REPORT("double_slit");
}
