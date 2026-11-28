#include <string.h>
#include "wr_test.h"
#include "wr_material.h"
#include "wr_optics.h"

static const double BK7_B[3] = {1.03961212, 0.231792344, 1.01046945};
static const double BK7_C[3] = {0.00600069867, 0.0200179144, 103.560653};
static const double FS_B[3]  = {0.6961663, 0.4079426, 0.8974794};            /* fused silica, Malitson 1965 */
static const double FS_C[3]  = {0.0046791483, 0.0135120631, 97.9340025};

static double n_at(const wr_material *m, double lam, int *ok)
{
    return wr_material_index(m, lam, ok).re;
}

int main(void)
{
    int ok;
    /* Sellmeier BK7 against published catalogue values (Schott): n_d(587.56 nm) = 1.5168,
       n_F(486.13 nm) = 1.52238, n_C(656.27 nm) = 1.51432, Abbe number V_d = 64.17 */
    wr_material bk7 = wr_material_sellmeier("bk7", BK7_B, BK7_C);
    double nd = n_at(&bk7, 587.56e-9, &ok), nF = n_at(&bk7, 486.13e-9, &ok), nC = n_at(&bk7, 656.27e-9, &ok);
    double Vd = (nd - 1.0) / (nF - nC);
    CHECK(ok);
    CHECK_NEAR(nd, 1.5168, 2e-5);
    CHECK_NEAR(nF, 1.52238, 3e-5);
    CHECK_NEAR(nC, 1.51432, 3e-5);
    CHECK_NEAR(Vd, 64.17, 0.1);
    printf("material: BK7 n_d %.5f, n_F %.5f, n_C %.5f, Abbe %.2f (catalogue 1.5168 / 1.52238 / 1.51432 / 64.17)\n", nd, nF, nC, Vd);

    wr_material fs = wr_material_sellmeier("silica", FS_B, FS_C);
    double fd = n_at(&fs, 587.56e-9, &ok), fF = n_at(&fs, 486.13e-9, &ok), fC = n_at(&fs, 656.27e-9, &ok);
    CHECK_NEAR(fd, 1.4585, 1e-4);
    CHECK_NEAR((fd - 1.0) / (fF - fC), 67.8, 0.5);                      /* published V_d of fused silica ~ 67.8 */
    printf("material: fused silica n_d %.5f, Abbe %.2f (published 1.4585 / ~67.8)\n", fd, (fd - 1.0) / (fF - fC));

    /* normal dispersion: n decreases monotonically with wavelength across the visible, for both models */
    wr_material cau = wr_material_cauchy("crown", 1.5046, 0.00420, 0.0);
    int mono = 1;
    double prev_s = 1e9, prev_c = 1e9;
    for (double lam = 400e-9; lam <= 700e-9; lam += 5e-9) {
        double s = n_at(&bk7, lam, &ok), c = n_at(&cau, lam, &ok);
        if (!(s < prev_s) || !(c < prev_c)) mono = 0;
        prev_s = s; prev_c = c;
    }
    CHECK(mono);
    CHECK_NEAR(n_at(&cau, 500e-9, &ok), 1.5046 + 0.00420 / 0.25, 1e-15);   /* hand evaluation of A + B/lambda^2 */
    CHECK_NEAR(n_at(&cau, 1.0e-6, &ok), 1.5046 + 0.00420, 1e-15);
    wr_material cau2 = wr_material_cauchy("c2", 1.4, 0.0, 0.01);
    CHECK_NEAR(n_at(&cau2, 0.5e-6, &ok), 1.4 + 0.01 / (0.25 * 0.25), 1e-15);

    /* constant, dispersion-free, with a complex part */
    wr_material water = wr_material_constant("water", 1.333, 0.0);
    CHECK(n_at(&water, 400e-9, &ok) == n_at(&water, 700e-9, &ok) && ok);
    wr_material gold = wr_material_constant("gold", 0.2, 3.2);
    wr_cplx nc = wr_material_index(&gold, 600e-9, &ok);
    CHECK_NEAR(nc.re, 0.2, 0); CHECK_NEAR(nc.im, 3.2, 0);

    /* validity ranges are reported, values still returned */
    wr_material vac = wr_material_constant("vacuum", 1.0, 0.0);
    CHECK(wr_material_index(&vac, 1e-12, &ok).re == 1.0 && ok);          /* constant: any wavelength */
    CHECK(wr_material_index(&vac, 0.0, &ok).re == 1.0 && !ok);           /* lambda = 0 is never valid */
    wr_material_index(&bk7, 100e-9, &ok);  CHECK(!ok);                   /* far UV */
    wr_material_index(&bk7, 5e-6, &ok);    CHECK(!ok);                   /* beyond 2.5 um */
    wr_material_index(&bk7, 550e-9, &ok);  CHECK(ok);
    wr_material_index(&cau, 2e-6, &ok);    CHECK(!ok);                   /* Cauchy is visible-only */
    wr_material nores = bk7; nores.lambda_min = 0.0; nores.lambda_max = INFINITY;
    double nres = n_at(&nores, sqrt(BK7_C[0]) * 1e-6 * 0.1, &ok);        /* deep UV, n^2 may go negative: flagged, no crash */
    CHECK(isnan(nres) ? !ok : 1);

    /* absorption: field exp(i k0 (n + i kappa) L) decays; intensity exp(-alpha L), alpha = 4 pi kappa / lambda */
    wr_material ab = wr_material_constant("tint", 1.5, 2e-6);
    double lam = 500e-9, L = 0.02;
    wr_cplx n_c = wr_material_index(&ab, lam, &ok);
    wr_cplx z = wr_cexp(wr_c(0.0, wr_wavenumber(lam) * L) );                       /* exp(i k0 L) */
    (void)z;
    wr_cplx field = wr_cexp(wr_cmul(wr_c(0, 1), wr_cscale(n_c, wr_wavenumber(lam) * L)));   /* exp(i k0 n_c L) */
    double alpha = wr_material_absorption_coeff(&ab, lam);
    CHECK_NEAR(alpha, 4.0 * WR_PI * 2e-6 / 500e-9, 1e-12);
    CHECK_NEAR(wr_cabs2(field), exp(-alpha * L), 1e-14);
    CHECK_NEAR(wr_carg(field), wr_wrap_phase(wr_wavenumber(lam) * 1.5 * L), 1e-9);  /* phase uses the real index only */
    CHECK(wr_cabs2(field) < 1.0);

    /* surface kinds */
    wr_material mir = wr_material_mirror("m"), blk = wr_material_absorber("b");
    CHECK(mir.surface == WR_SURF_IDEAL_MIRROR && blk.surface == WR_SURF_ABSORBER);
    wr_cplx mi = wr_material_index(&mir, 500e-9, &ok);
    CHECK(mi.re == 1.0 && mi.im == 0.0);
    CHECK(wr_material_absorption_coeff(&mir, 500e-9) == 0.0);

    /* validation */
    char err[128] = "";
    CHECK(wr_material_validate(&bk7, err, sizeof err) && wr_material_validate(&mir, err, sizeof err));
    wr_material m1 = wr_material_constant("x", -1.0, 0.0);   CHECK(!wr_material_validate(&m1, err, sizeof err) && err[0]);
    wr_material m2 = wr_material_constant("x", 1.5, -0.1);   CHECK(!wr_material_validate(&m2, err, sizeof err));
    wr_material m3 = wr_material_constant("x", NAN, 0.0);    CHECK(!wr_material_validate(&m3, err, sizeof err));
    wr_material m4 = bk7; m4.lambda_max = m4.lambda_min;     CHECK(!wr_material_validate(&m4, err, sizeof err));
    wr_material m5 = bk7; m5.sell_c[1] = -1.0;               CHECK(!wr_material_validate(&m5, err, sizeof err));
    wr_material m6 = wr_material_constant("", 1.0, 0.0);     CHECK(!wr_material_validate(&m6, err, sizeof err));
    wr_material m7 = wr_material_constant("a-very-long-name-that-does-not-fit-in-forty-eight-bytes-at-all", 1.0, 0.0);
    CHECK(strlen(m7.name) == WR_NAME_LEN - 1);                           /* truncated, still NUL terminated */

    TEST_REPORT("material");
}
