#include <stdlib.h>
#include "wr_test.h"
#include "wr_propagate.h"
#include "wr_optics.h"

static double rnd(unsigned long long *s)
{
    *s = *s * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(*s >> 11) / 4503599627370496.0 - 1.0;
}

static double max_diff(const wr_field *a, const wr_field *b)
{
    double e = 0;
    for (size_t i = 0; i < a->nx * a->ny; i++) e = fmax(e, wr_cabs(wr_csub(a->data[i], b->data[i])));
    return e;
}

int main(void)
{
    const size_t nx = 64, ny = 32;
    const double d = 1e-6, lam0 = 0.5e-6;
    wr_fft2_plan *plan = wr_fft2_plan_create(nx, ny);
    wr_field *f = wr_field_create(nx, ny, d, d);

    /* 1. Exact analytic case: a tilted plane wave with whole cycles in the
     *    window just acquires exp(i kz z), kz = k sqrt(1 - (lam fx)^2 - (lam fy)^2). */
    {
        const double fx = 5.0 / ((double)nx * d), fy = 3.0 / ((double)ny * d), z = 137.3e-6;
        for (size_t iy = 0; iy < ny; iy++)
            for (size_t ix = 0; ix < nx; ix++)
                *wr_field_at(f, ix, iy) = wr_cexp_i(WR_TWO_PI * (fx * (double)ix * d + fy * (double)iy * d));
        CHECK(wr_propagate_asm(f, plan, lam0, 1.0, z) == 0);
        double kz = WR_TWO_PI / lam0 * sqrt(1.0 - lam0 * lam0 * (fx * fx + fy * fy));
        double err = 0;
        for (size_t iy = 0; iy < ny; iy++)
            for (size_t ix = 0; ix < nx; ix++) {
                wr_cplx want = wr_cexp_i(WR_TWO_PI * (fx * (double)ix * d + fy * (double)iy * d) + kz * z);
                err = fmax(err, wr_cabs(wr_csub(*wr_field_at(f, ix, iy), want)));
            }
        CHECK_NEAR(err, 0.0, 1e-10);
    }

    /* random field containing every grid frequency (all propagating here) */
    unsigned long long seed = 2024;
    wr_field *orig = wr_field_create(nx, ny, d, d);
    for (size_t i = 0; i < nx * ny; i++) f->data[i] = orig->data[i] = wr_c(rnd(&seed), rnd(&seed));

    /* 2. power conservation */
    double p0 = wr_field_power(f);
    wr_propagate_asm(f, plan, lam0, 1.0, 321e-6);
    CHECK_NEAR(wr_field_power(f) / p0, 1.0, 1e-12);

    /* 3. reversibility: +z then -z */
    wr_propagate_asm(f, plan, lam0, 1.0, -321e-6);
    CHECK_NEAR(max_diff(f, orig), 0.0, 1e-9);

    /* 4. semigroup: z1 then z2 equals z1+z2 */
    wr_field *g = wr_field_create(nx, ny, d, d);
    for (size_t i = 0; i < nx * ny; i++) f->data[i] = g->data[i] = orig->data[i];
    wr_propagate_asm(f, plan, lam0, 1.0, 100e-6);
    wr_propagate_asm(f, plan, lam0, 1.0, 250e-6);
    wr_propagate_asm(g, plan, lam0, 1.0, 350e-6);
    CHECK_NEAR(max_diff(f, g), 0.0, 1e-9);
    wr_field_destroy(g);

    /* 5. refractive index: propagating in medium n equals vacuum with lambda0/n */
    {
        wr_field *a = wr_field_create(nx, ny, d, d), *b = wr_field_create(nx, ny, d, d);
        for (size_t i = 0; i < nx * ny; i++) a->data[i] = b->data[i] = orig->data[i];
        wr_propagate_asm(a, plan, lam0, 1.5, 80e-6);
        wr_propagate_asm(b, plan, lam0 / 1.5, 1.0, 80e-6);
        CHECK_NEAR(max_diff(a, b), 0.0, 1e-9);
        wr_field_destroy(a); wr_field_destroy(b);
    }
    wr_field_destroy(orig); wr_field_destroy(f); wr_fft2_plan_destroy(plan);

    /* 6. evanescent wave (lam*fx > 1): decays as exp(-k z sqrt(lam^2 fx^2 - 1)) for z > 0,
     *    and is discarded (not amplified) for z < 0. Needs d < lam/2. */
    {
        const double de = 0.2e-6;
        wr_fft2_plan *pe = wr_fft2_plan_create(64, 1);
        wr_field *e = wr_field_create(64, 1, de, de);
        const double fx = 28.0 / (64.0 * de), z = 1e-6;
        for (size_t ix = 0; ix < 64; ix++) *wr_field_at(e, ix, 0) = wr_cexp_i(WR_TWO_PI * fx * (double)ix * de);
        wr_propagate_asm(e, pe, lam0, 1.0, z);
        double want = exp(-WR_TWO_PI / lam0 * z * sqrt(lam0 * lam0 * fx * fx - 1.0));
        CHECK_NEAR(wr_cabs(*wr_field_at(e, 10, 0)) / want, 1.0, 1e-10);
        for (size_t ix = 0; ix < 64; ix++) *wr_field_at(e, ix, 0) = wr_cexp_i(WR_TWO_PI * fx * (double)ix * de);
        wr_propagate_asm(e, pe, lam0, 1.0, -z);
        CHECK_NEAR(wr_cabs(*wr_field_at(e, 10, 0)), 0.0, 1e-12);
        wr_field_destroy(e); wr_fft2_plan_destroy(pe);
    }

    /* invalid arguments */
    {
        wr_fft2_plan *p2 = wr_fft2_plan_create(8, 8);
        wr_field *h = wr_field_create(16, 8, d, d);
        CHECK(wr_propagate_asm(h, p2, lam0, 1.0, 1e-6) != 0);   /* size mismatch */
        CHECK(wr_propagate_asm(h, NULL, lam0, 1.0, 1e-6) != 0);
        wr_field_destroy(h); wr_fft2_plan_destroy(p2);
    }
    TEST_REPORT("propagate");
}
