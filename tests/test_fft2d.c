#include <stdlib.h>
#include "wr_test.h"
#include "wr_fft2d.h"
#include "wr_optics.h"

static double rnd(unsigned long long *s)
{
    *s = *s * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(*s >> 11) / 4503599627370496.0 - 1.0;
}

int main(void)
{
    CHECK(wr_fft2_plan_create(6, 4) == NULL);
    CHECK(wr_fft2_plan_create(4, 0) == NULL);

    /* non-square (8 x 4) against the O(N^2) 2D DFT: catches x/y mix-ups */
    {
        const size_t nx = 8, ny = 4;
        wr_fft2_plan *p = wr_fft2_plan_create(nx, ny);
        wr_cplx x[32], y[32];
        unsigned long long seed = 7;
        for (size_t i = 0; i < nx * ny; i++) x[i] = y[i] = wr_c(rnd(&seed), rnd(&seed));
        wr_fft2_execute(p, y, 0);
        double err = 0;
        for (size_t my = 0; my < ny; my++)
            for (size_t mx = 0; mx < nx; mx++) {
                wr_cplx acc = wr_c(0, 0);
                for (size_t jy = 0; jy < ny; jy++)
                    for (size_t jx = 0; jx < nx; jx++) {
                        double ang = -WR_TWO_PI * ((double)((mx * jx) % nx) / (double)nx
                                                 + (double)((my * jy) % ny) / (double)ny);
                        acc = wr_cadd(acc, wr_cmul(x[jy * nx + jx], wr_cexp_i(ang)));
                    }
                err = fmax(err, wr_cabs(wr_csub(y[my * nx + mx], acc)));
            }
        CHECK_NEAR(err, 0.0, 1e-13);
        wr_fft2_plan_destroy(p);
    }

    /* ny = 1 degenerates to the 1D transform; round trip + Parseval at 128 x 64 */
    {
        const size_t nx = 128, ny = 64, n = nx * ny;
        wr_fft2_plan *p = wr_fft2_plan_create(nx, ny);
        wr_cplx *x = malloc(n * sizeof *x), *y = malloc(n * sizeof *y);
        unsigned long long seed = 99;
        double et = 0, ef = 0;
        for (size_t i = 0; i < n; i++) { x[i] = y[i] = wr_c(rnd(&seed), rnd(&seed)); et += wr_cabs2(x[i]); }
        wr_fft2_execute(p, y, 0);
        for (size_t i = 0; i < n; i++) ef += wr_cabs2(y[i]);
        CHECK_NEAR(ef / (double)n / et, 1.0, 1e-12);
        wr_fft2_execute(p, y, 1);
        double err = 0;
        for (size_t i = 0; i < n; i++) err = fmax(err, wr_cabs(wr_csub(y[i], x[i])));
        CHECK_NEAR(err, 0.0, 1e-13);
        free(x); free(y); wr_fft2_plan_destroy(p);
    }
    TEST_REPORT("fft2d");
}
