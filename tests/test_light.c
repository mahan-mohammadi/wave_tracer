#include <float.h>
#include "wr_test.h"
#include "wr_rng.h"
#include "wr_light.h"
#include "wr_integrate.h"
#include "wr_optics.h"

static wr_vec3 rv(wr_rng *r, double lo, double hi)
{
    return wr_v3(wr_rng_range(r, lo, hi), wr_rng_range(r, lo, hi), wr_rng_range(r, lo, hi));
}

/* 4th-order central finite-difference Laplacian of a complex field (error ~ (k h)^4 / 90 relative) */
static wr_cplx laplacian(const wr_light *l, wr_vec3 p, double h)
{
    wr_cplx f0 = wr_light_field(l, p, 0), s = wr_c(0, 0);
    for (int a = 0; a < 3; a++) {
        wr_vec3 e = wr_v3(a == 0, a == 1, a == 2);
        wr_cplx p1 = wr_light_field(l, wr_v3_madd(p, e, h), 0), m1 = wr_light_field(l, wr_v3_madd(p, e, -h), 0);
        wr_cplx p2 = wr_light_field(l, wr_v3_madd(p, e, 2 * h), 0), m2 = wr_light_field(l, wr_v3_madd(p, e, -2 * h), 0);
        wr_cplx t = wr_csub(wr_cscale(wr_cadd(p1, m1), 16.0), wr_cadd(p2, m2));
        t = wr_csub(t, wr_cscale(f0, 30.0));
        s = wr_cadd(s, wr_cscale(t, 1.0 / (12.0 * h * h)));
    }
    return s;
}

/* |laplacian E + k^2 E| / (k^2 |E|): zero for an exact Helmholtz solution */
static double helmholtz_residual(const wr_light *l, wr_vec3 p, double lambda)
{
    double k = wr_wavenumber(lambda);
    wr_cplx r = wr_cadd(laplacian(l, p, lambda / 40.0), wr_cscale(wr_light_field(l, p, 0), k * k));
    return wr_cabs(r) / (k * k * wr_cabs(wr_light_field(l, p, 0)));
}

static double beam_power_integrand_z; static const wr_light *beam_light; static wr_vec3 beam_origin;
static double radial_power(double rho, void *ctx)   /* |E|^2 * 2 pi rho across the plane z = const */
{
    (void)ctx;
    wr_vec3 p = wr_v3(rho, 0, beam_power_integrand_z);
    return wr_cabs2(wr_light_field(beam_light, wr_v3_add(p, beam_origin), 0)) * WR_TWO_PI * rho;
}

int main(void)
{
    wr_rng rng = {99};
    double lam = 500e-9, k = wr_wavenumber(lam);

    /* plane wave: phase advances by k * (displacement along d), constant amplitude, transverse shift changes nothing */
    wr_light pw = wr_light_plane_wave("pw", wr_v3(0.1, 0.2, 0.3), wr_v3(1, 2, -2), lam, 2.0);   /* dir normalized by constructor */
    pw.phase = 0.7;
    CHECK_NEAR(wr_v3_len(pw.dir), 1.0, 1e-16);
    double worst_a = 0, worst_p = 0, worst_t = 0;
    for (int i = 0; i < 2000; i++) {
        wr_vec3 p = rv(&rng, -1e-3, 1e-3), q = rv(&rng, -1e-3, 1e-3);
        wr_cplx a = wr_light_field(&pw, p, 0), b = wr_light_field(&pw, q, 0);
        worst_a = fmax(worst_a, fabs(wr_cabs(a) - 2.0));
        double dphi = k * wr_v3_dot(wr_v3_sub(q, p), pw.dir);
        wr_cplx ratio = wr_cdiv(b, a);
        worst_p = fmax(worst_p, wr_cabs(wr_csub(ratio, wr_cexp_i(dphi))));
        wr_vec3 t, bb; wr_v3_basis(pw.dir, &t, &bb);
        wr_cplx c = wr_light_field(&pw, wr_v3_madd(p, t, 3e-4), 0);
        worst_t = fmax(worst_t, wr_cabs(wr_csub(c, a)));
    }
    /* Phase floor of the representation: a double at |x| has spacing ~eps |x|, so phase carries error ~ k eps |x|.
       Sources 0.1-0.3 m from the origin cannot do better than ~1e-9 rad at 500 nm; this is not a code error. */
    double floor_pw = k * DBL_EPSILON * 0.4;
    CHECK_NEAR(worst_a, 0, 1e-15); CHECK(worst_p < 10 * floor_pw); CHECK(worst_t < 10 * floor_pw);
    CHECK_NEAR(wr_carg(wr_light_field(&pw, pw.pos, 0)), 0.7, 1e-15);                          /* initial phase at the reference point */
    printf("light: plane wave amplitude %.1e, phase ratio %.1e, transverse invariance %.1e (representation floor k*eps*|x| = %.1e)\n",
           worst_a, worst_p, worst_t, floor_pw);
    /* control: same test with the source at the origin and probes within 1 mm: error drops with |x| */
    {
        wr_light pc = wr_light_plane_wave("pc", wr_v3_zero(), wr_v3(1, 2, -2), lam, 2.0);
        double w = 0;
        for (int i = 0; i < 2000; i++) {
            wr_vec3 p = rv(&rng, -1e-3, 1e-3), q = rv(&rng, -1e-3, 1e-3);
            wr_cplx r = wr_cdiv(wr_light_field(&pc, q, 0), wr_light_field(&pc, p, 0));
            w = fmax(w, wr_cabs(wr_csub(r, wr_cexp_i(k * wr_v3_dot(wr_v3_sub(q, p), pc.dir)))));
        }
        CHECK(w < 10 * k * DBL_EPSILON * 1e-3);
        printf("light: plane wave control, source at origin, probes within 1 mm: phase ratio %.1e\n", w);
    }

    /* point source: |E| = A/r, phase k r, energy flux through every sphere identical (A^2 per steradian) */
    wr_light ps = wr_light_point("ps", wr_v3(0.5, -0.5, 0.25), lam, 3.0);
    double w_flux = 0, w_ph = 0, w_hh = 0;
    for (int i = 0; i < 500; i++) {
        wr_vec3 d = wr_v3_normalize(rv(&rng, -1, 1));
        double r = wr_rng_range(&rng, 3 * lam, 50e-6);
        wr_vec3 p = wr_v3_madd(ps.pos, d, r);
        wr_cplx e = wr_light_field(&ps, p, 0);
        w_flux = fmax(w_flux, fabs(r * r * wr_cabs2(e) - 9.0) / 9.0);
        w_ph = fmax(w_ph, fabs(wr_wrap_phase(wr_carg(e) - k * r)));
        w_hh = fmax(w_hh, helmholtz_residual(&ps, p, lam));
    }
    /* position rounding: relative error in r ~ eps |pos| / r (r >= 3 lambda), phase error ~ k eps |pos| */
    CHECK(w_flux < 4 * DBL_EPSILON * 0.75 / (3 * lam)); CHECK(w_ph < 10 * k * DBL_EPSILON * 0.75); CHECK(w_hh < 2e-5);
    wr_cplx at_source = wr_light_field(&ps, ps.pos, 0);
    CHECK(isnan(at_source.re) && isnan(at_source.im));
    printf("light: point source r^2|E|^2 = A^2 to %.1e, phase %.1e, Helmholtz residual %.1e\n", w_flux, w_ph, w_hh);

    /* Gaussian beam: closed-form checks */
    double w0 = 40e-6, zr = WR_PI * w0 * w0 / lam;
    wr_light gb = wr_light_gaussian("gb", wr_v3(0, 0, 0), wr_v3(0, 0, 1), w0, lam, 1.5);
    wr_cplx e0 = wr_light_field(&gb, wr_v3_zero(), 0);
    CHECK_NEAR(e0.re, 1.5, 1e-15); CHECK_NEAR(e0.im, 0, 1e-15);                              /* waist, on axis */
    CHECK_NEAR(wr_cabs(wr_light_field(&gb, wr_v3(w0, 0, 0), 0)), 1.5 / exp(1.0), 1e-14);   /* 1/e field radius = w0 */
    wr_cplx ez = wr_light_field(&gb, wr_v3(0, 0, zr), 0);
    CHECK_NEAR(wr_cabs(ez), 1.5 / sqrt(2.0), 1e-14);                                          /* amplitude A/sqrt(2) at z_R */
    CHECK_NEAR(wr_wrap_phase(wr_carg(ez) - (k * zr - WR_PI / 4)), 0, 1e-8);                    /* kz - Gouy(pi/4) */
    for (int i = 1; i <= 4; i++) {                                                            /* width law w(z) = w0 sqrt(1 + (z/zr)^2) */
        double z = i * 0.75 * zr, rho = 0.6 * w0;
        double ratio = wr_cabs(wr_light_field(&gb, wr_v3(rho, 0, z), 0)) / wr_cabs(wr_light_field(&gb, wr_v3(0, 0, z), 0));
        double w = rho / sqrt(-log(ratio));
        CHECK_NEAR(w / (w0 * sqrt(1 + (z / zr) * (z / zr))), 1.0, 1e-12);
    }
    /* power through any transverse plane is pi w0^2 A^2 / 2 */
    beam_light = &gb; beam_origin = wr_v3_zero();
    double pexp = WR_PI * w0 * w0 * 1.5 * 1.5 / 2.0, w_pow = 0;
    for (int i = 0; i < 4; i++) {
        beam_power_integrand_z = i * zr;
        double P = wr_simpson(radial_power, NULL, 0.0, 12.0 * w0 * sqrt(1 + i * i), 4000);
        w_pow = fmax(w_pow, fabs(P - pexp) / pexp);
    }
    CHECK_NEAR(w_pow, 0, 1e-9);
    /* a tilted beam: same field in its own frame */
    wr_light gt = wr_light_gaussian("gt", wr_v3(1, 2, 3), wr_v3(2, -1, 2), w0, lam, 1.5);
    wr_vec3 t1, t2; wr_v3_basis(gt.dir, &t1, &t2);
    wr_vec3 pt = wr_v3_add(gt.pos, wr_v3_add(wr_v3_scale(gt.dir, 0.4 * zr), wr_v3_scale(t1, 0.5 * w0)));
    wr_vec3 pl = wr_v3(0.5 * w0, 0, 0.4 * zr);
    CHECK(wr_cabs(wr_csub(wr_light_field(&gt, pt, 0), wr_light_field(&gb, pl, 0))) < 10 * k * DBL_EPSILON * 3.5 * 1.5);   /* tilted source sits ~3.5 m from the origin: phase floor */
    printf("light: Gaussian beam waist/Rayleigh/Gouy/width law exact to 1e-12, power conserved to %.1e\n", w_pow);

    /* The Gaussian beam is a paraxial solution, so its Helmholtz residual is nonzero and shrinks with the waist
       like (lambda / w0)^4 (error term d^2 psi / dz^2 relative to k^2 psi). Finite-difference floor here ~7e-6. */
    double res_prev = 1.0;
    printf("light: Gaussian Helmholtz residual:");
    for (int m = 0; m < 3; m++) {
        double ww = (m == 0 ? 2.0 : (m == 1 ? 3.0 : 20.0)) * lam;
        wr_light g = wr_light_gaussian("g", wr_v3_zero(), wr_v3(0, 0, 1), ww, lam, 1.0);
        double zrr = WR_PI * ww * ww / lam, hh = 0.0;
        for (int i = 0; i < 20; i++) {
            wr_vec3 p = wr_v3(wr_rng_range(&rng, -0.8, 0.8) * ww, wr_rng_range(&rng, -0.8, 0.8) * ww, wr_rng_range(&rng, -1, 1) * zrr);
            hh = fmax(hh, helmholtz_residual(&g, p, lam));
        }
        printf("  w0=%.0f lambda: %.1e", ww / lam, hh);
        CHECK(hh < res_prev);
        res_prev = hh;
        if (m == 0) CHECK(hh < 1e-3);
        if (m == 2) CHECK(hh < 5e-5);
    }
    printf("\n");

    /* spectrum: lines are incoherent, their intensities add to amplitude^2; weights are power fractions */
    wr_light sp = wr_light_plane_wave("sp", wr_v3_zero(), wr_v3(0, 0, 1), 450e-9, 2.0);
    sp.spectrum.n = 3;
    sp.spectrum.lambda[0] = 450e-9; sp.spectrum.weight[0] = 1;
    sp.spectrum.lambda[1] = 550e-9; sp.spectrum.weight[1] = 2;
    sp.spectrum.lambda[2] = 650e-9; sp.spectrum.weight[2] = 1;
    double isum = 0;
    for (int i = 0; i < 3; i++) isum += wr_cabs2(wr_light_field(&sp, wr_v3(1e-4, 2e-4, 3e-4), i));
    CHECK_NEAR(isum, 4.0, 1e-14);
    CHECK_NEAR(wr_light_line_amplitude(&sp, 1) * wr_light_line_amplitude(&sp, 1), 2.0, 1e-14);   /* half of the power */
    CHECK(wr_light_validate(&sp, NULL, 0));
    /* each line has its own k */
    wr_cplx f1 = wr_light_field(&sp, wr_v3(0, 0, 1e-6), 1);
    CHECK_NEAR(wr_carg(f1), wr_wrap_phase(wr_wavenumber(550e-9) * 1e-6), 1e-9);

    /* polarization basis: orthonormal, right-handed with d, e1 = projection of pol_ref */
    wr_light pl2 = wr_light_plane_wave("pl", wr_v3_zero(), wr_v3(0, 0, 1), lam, 1.0);
    pl2.pol_ref = wr_v3(0.3, 0.9, 0.5);
    double w_b = 0;
    for (int i = 0; i < 3000; i++) {
        wr_vec3 d = wr_v3_normalize(rv(&rng, -1, 1)), e1, e2;
        wr_light_pol_basis(&pl2, d, &e1, &e2);
        w_b = fmax(w_b, fabs(wr_v3_len(e1) - 1)); w_b = fmax(w_b, fabs(wr_v3_len(e2) - 1));
        w_b = fmax(w_b, fabs(wr_v3_dot(e1, d))); w_b = fmax(w_b, fabs(wr_v3_dot(e2, d))); w_b = fmax(w_b, fabs(wr_v3_dot(e1, e2)));
        w_b = fmax(w_b, wr_v3_dist(wr_v3_cross(e1, e2), d));
        wr_vec3 proj = wr_v3_normalize(wr_v3_madd(pl2.pol_ref, d, -wr_v3_dot(pl2.pol_ref, d)));
        w_b = fmax(w_b, wr_v3_dist(e1, proj));
    }
    wr_vec3 e1, e2;
    wr_light_pol_basis(&pl2, wr_v3_normalize(pl2.pol_ref), &e1, &e2);                         /* d parallel to pol_ref: fallback */
    w_b = fmax(w_b, fabs(wr_v3_dot(e1, pl2.pol_ref)) / wr_v3_len(pl2.pol_ref));
    w_b = fmax(w_b, fabs(wr_v3_len(e1) - 1));
    CHECK_NEAR(w_b, 0, 2e-14);

    /* validation */
    char err[128] = "";
    CHECK(wr_light_validate(&gb, err, sizeof err));
    wr_light b1 = gb; b1.amplitude = -1;                         CHECK(!wr_light_validate(&b1, err, sizeof err) && err[0]);
    wr_light b2 = gb; b2.spectrum.weight[0] = 0;                 CHECK(!wr_light_validate(&b2, err, sizeof err));
    wr_light b3 = gb; b3.spectrum.lambda[0] = -5e-7;             CHECK(!wr_light_validate(&b3, err, sizeof err));
    wr_light b4 = gb; b4.dir = wr_v3(0, 0, 2);                   CHECK(!wr_light_validate(&b4, err, sizeof err));
    wr_light b5 = gb; b5.waist = 0;                              CHECK(!wr_light_validate(&b5, err, sizeof err));
    wr_light b6 = gb; b6.jones[0] = wr_c(1, 0); b6.jones[1] = wr_c(1, 0);   CHECK(!wr_light_validate(&b6, err, sizeof err));
    wr_light b7 = gb; b7.spectrum.n = 0;                         CHECK(!wr_light_validate(&b7, err, sizeof err));
    wr_light b8 = gb; b8.coherence_length = -1;                  CHECK(!wr_light_validate(&b8, err, sizeof err));
    wr_light b9 = gb; b9.phase = NAN;                            CHECK(!wr_light_validate(&b9, err, sizeof err));
    wr_light ok1 = gb; ok1.coherence_length = 0.0;               CHECK(wr_light_validate(&ok1, err, sizeof err));     /* incoherent is legal */
    wr_light ok2 = gb; ok2.jones[0] = wr_c(sqrt(0.5), 0); ok2.jones[1] = wr_c(0, sqrt(0.5));   /* circular */
    CHECK(wr_light_validate(&ok2, err, sizeof err));

    TEST_REPORT("light");
}
