/* Validation: Michelson interferometer (roadmap Phase 3, scene 4).
 *
 * Model (scalar, unfolded geometry): a Gaussian beam U0 (waist w0) enters a
 * lossless 50:50 beam splitter. Each arm of length L is "propagate L to the
 * mirror, reflect, propagate L back": a plane mirror adds nothing beyond the
 * propagation; a mirror tilted by alpha multiplies the field by exp(i q x),
 * q = 2 k alpha (thin-phase-screen, small-angle). The arms recombine:
 *   E_det = (U1 - U2)/2     (port with the pi shift of the lossless splitter)
 *   E_src = (U1 + U2)/2     (port returning to the source)
 * Which port is dark at zero path difference is a convention of the splitter
 * coating; what is physical is that the two ports are complementary:
 *   P_det + P_src = P_in   (energy conservation).
 *
 * Checks
 *  1. Propagated Gaussian beam vs the closed-form Gaussian beam (complex field).
 *  2. Tilted-mirror arm vs the exact paraxial identity for a tilted, sheared
 *     beam: P_L{ f e^{iqx} }(x) = e^{iqx - i q^2 L/(2k)} (P_L f)(x - qL/k),
 *     evaluated with the analytic Gaussian (so no interpolation).
 *  3. Energy conservation of the two output ports (with and without tilt).
 *  4. Phase-dependent fringes: P_det/P_in = (1 - cos(2 pi Delta / lambda)) / 2
 *     for 16 round-trip path differences Delta.
 *  5. Partial coherence: for a Gaussian spectrum (std sigma_k in wavenumber) the
 *     visibility is exp(-sigma_k^2 Delta^2 / 2). Intensities are added ACROSS
 *     wavelengths, amplitudes WITHIN each wavelength.
 * Limits: scalar and paraxial reference formulas; small-angle tilt model;
 * perfectly reflecting, non-dispersive mirrors and splitter.
 */
#include <stdio.h>
#include <stdlib.h>
#include "wr_test.h"
#include "wr_propagate.h"
#include "wr_optics.h"
#include "wr_sum.h"

#define N 256
static const double DX = 40e-6, W0 = 0.3e-3, L1 = 0.5, LAM0 = 500e-9;

static void init_gauss(wr_field *f)
{
    for (size_t iy = 0; iy < f->ny; iy++)
        for (size_t ix = 0; ix < f->nx; ix++) {
            double x = wr_field_x(f, ix), y = wr_field_y(f, iy);
            *wr_field_at(f, ix, iy) = wr_c(exp(-(x * x + y * y) / (W0 * W0)), 0.0);
        }
}

/* closed-form circular Gaussian beam, waist W0 at z = 0, time dependence exp(-i w t) */
static wr_cplx gauss_an(double x, double y, double z, double lam)
{
    double k = WR_TWO_PI / lam, zr = WR_PI * W0 * W0 / lam;
    double w = W0 * sqrt(1.0 + (z / zr) * (z / zr)), r2 = x * x + y * y;
    double ph = fmod(k * z, WR_TWO_PI) - atan(z / zr) + k * r2 * z / (2.0 * (z * z + zr * zr));
    return wr_cpolar((W0 / w) * exp(-r2 / (w * w)), ph);
}

static void copy(wr_field *dst, const wr_field *src)
{
    for (size_t i = 0; i < src->nx * src->ny; i++) dst->data[i] = src->data[i];
}

static void ramp(wr_field *f, double q)
{
    for (size_t iy = 0; iy < f->ny; iy++)
        for (size_t ix = 0; ix < f->nx; ix++) {
            wr_cplx *u = wr_field_at(f, ix, iy);
            *u = wr_cmul(*u, wr_cexp_i(fmod(q * wr_field_x(f, ix), WR_TWO_PI)));
        }
}

/* (1/4) * integral |a + sign*b|^2 dA : power of an output port */
static double port_power(const wr_field *a, const wr_field *b, double sign)
{
    wr_ksum k = wr_ksum_zero();
    for (size_t i = 0; i < a->nx * a->ny; i++) {
        wr_cplx s = wr_cadd(a->data[i], wr_cscale(b->data[i], sign));
        wr_ksum_add(&k, wr_cabs2(s));
    }
    return 0.25 * wr_ksum_value(&k) * a->dx * a->dy;
}

static int prop(wr_field *f, const wr_fft2_plan *p, double lam, double z)
{
    return wr_propagate_asm(f, p, lam, 1.0, z);
}

int main(void)
{
    wr_fft2_plan *plan = wr_fft2_plan_create(N, N);
    wr_field *u0 = wr_field_create(N, N, DX, DX), *u1 = wr_field_create(N, N, DX, DX);
    wr_field *u2 = wr_field_create(N, N, DX, DX);
    init_gauss(u0);
    const double p_in = wr_field_power(u0);
    printf("Michelson: lambda = %.0f nm, w0 = %.1f mm, arm = %.0f mm, grid %dx%d at %.0f um\n",
           LAM0 * 1e9, W0 * 1e3, L1 * 1e3, N, N, DX * 1e6);

    /* 1. propagated beam vs closed-form Gaussian beam (round trip = 2 L1 = 1 m) */
    copy(u1, u0);
    CHECK(prop(u1, plan, LAM0, 2.0 * L1) == 0);
    double e1 = 0;
    for (size_t iy = 0; iy < N; iy++)
        for (size_t ix = 0; ix < N; ix++) {
            wr_cplx an = gauss_an(wr_field_x(u1, ix), wr_field_y(u1, iy), 2.0 * L1, LAM0);
            e1 = fmax(e1, wr_cabs(wr_csub(*wr_field_at(u1, ix, iy), an)));
        }
    printf("1. Gaussian beam: max |U_sim - U_analytic| = %.2e (peak amplitude 1)\n", e1);
    CHECK_NEAR(e1, 0.0, 1e-5);

    /* 2. tilted mirror (alpha = 50 urad) vs paraxial tilt+shear identity */
    const double alpha = 5e-5, q = 2.0 * WR_TWO_PI / LAM0 * alpha, k0 = WR_TWO_PI / LAM0;
    const double shear = q * L1 / k0;
    copy(u2, u0);
    CHECK(prop(u2, plan, LAM0, L1) == 0);
    ramp(u2, q);
    CHECK(prop(u2, plan, LAM0, L1) == 0);
    double e2 = 0;
    for (size_t iy = 0; iy < N; iy++)
        for (size_t ix = 0; ix < N; ix++) {
            double x = wr_field_x(u2, ix), y = wr_field_y(u2, iy);
            wr_cplx an = wr_cmul(wr_cexp_i(fmod(q * x, WR_TWO_PI) - q * q * L1 / (2.0 * k0)),
                                 gauss_an(x - shear, y, 2.0 * L1, LAM0));
            e2 = fmax(e2, wr_cabs(wr_csub(*wr_field_at(u2, ix, iy), an)));
        }
    printf("2. tilted arm: tilt q = %.0f rad/m, shear %.1f um, max |U_sim - U_identity| = %.2e\n",
           q, shear * 1e6, e2);
    CHECK_NEAR(e2, 0.0, 1e-5);

    /* 3. energy conservation across the two ports, aligned and tilted */
    double ps = port_power(u1, u2, +1.0), pd = port_power(u1, u2, -1.0);
    printf("3. tilted: P_det/P_in = %.6f, P_src/P_in = %.6f, sum = %.15f\n", pd / p_in, ps / p_in, (pd + ps) / p_in);
    CHECK_NEAR((pd + ps) / p_in, 1.0, 1e-12);
    CHECK(pd / p_in > 0.0 && pd / p_in < 1.0);

    /* 4. phase-dependent fringes: arm 2 is longer by delta, round-trip difference Delta = 2 delta */
    double worst = 0;
    printf("4. Delta/lambda   P_det/P_in   analytic (1-cos)/2\n");
    for (int j = 0; j < 16; j++) {
        double Delta = j * LAM0 / 8.0;
        copy(u2, u0);
        CHECK(prop(u2, plan, LAM0, 2.0 * (L1 + 0.5 * Delta)) == 0);
        double got = port_power(u1, u2, -1.0) / p_in;
        double want = 0.5 * (1.0 - cos(WR_TWO_PI * Delta / LAM0));
        worst = fmax(worst, fabs(got - want));
        if (j % 2 == 0) printf("   %6.3f         %.6f      %.6f\n", Delta / LAM0, got, want);
        double chk = port_power(u1, u2, +1.0) / p_in + got;
        CHECK_NEAR(chk, 1.0, 1e-12);
    }
    printf("   max |difference| over 16 path differences = %.2e\n", worst);
    CHECK_NEAR(worst, 0.0, 1e-6);

    /* 5. partial coherence: Gaussian spectrum in wavenumber, sigma_lambda = 5 nm */
    {
        const double sig_k = WR_TWO_PI * 5e-9 / (LAM0 * LAM0);
        const int M = 41;
        const int ms[4] = {0, 8, 16, 32};   /* Delta = m * lambda0, so cos(k0 Delta) = 1 */
        double acc[4] = {0, 0, 0, 0}, wsum = 0;
        for (int j = 0; j < M; j++) {
            double t = -5.0 + 0.25 * j, w = exp(-0.5 * t * t);
            double lam = WR_TWO_PI / (k0 + sig_k * t);
            copy(u1, u0);
            prop(u1, plan, lam, 2.0 * L1);
            for (int c = 0; c < 4; c++) {
                copy(u2, u1);
                prop(u2, plan, lam, ms[c] * LAM0);       /* arm 2 longer by Delta/2 each way */
                acc[c] += w * port_power(u1, u2, -1.0) / p_in;
            }
            wsum += w;
        }
        printf("5. coherence length 1/sigma_k = %.2f um\n   Delta (um)   V simulated   V = exp(-sk^2 D^2 / 2)\n", 1e6 / sig_k);
        for (int c = 0; c < 4; c++) {
            double D = ms[c] * LAM0, v_sim = 1.0 - 2.0 * acc[c] / wsum, v_an = exp(-0.5 * sig_k * sig_k * D * D);
            printf("   %8.2f     %.5f        %.5f\n", D * 1e6, v_sim, v_an);
            CHECK_NEAR(v_sim, v_an, 1e-4);
        }
    }

    wr_field_destroy(u0); wr_field_destroy(u1); wr_field_destroy(u2); wr_fft2_plan_destroy(plan);
    TEST_REPORT("michelson");
}
