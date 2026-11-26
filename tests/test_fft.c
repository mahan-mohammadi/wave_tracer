#include <stdlib.h>
#include "wr_test.h"
#include "wr_fft.h"
#include "wr_optics.h"

/* deterministic pseudo-random data in [-1, 1) */
static double rnd(unsigned long long *s)
{
    *s = *s * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(*s >> 11) / 4503599627370496.0 - 1.0;
}

static void naive_dft(const wr_cplx *x, wr_cplx *X, size_t n)
{
    for (size_t m = 0; m < n; m++) {
        wr_cplx acc = wr_c(0, 0);
        for (size_t j = 0; j < n; j++) {
            /* reduce m*j mod n first so the angle stays small and accurate */
            double ang = -WR_TWO_PI * (double)((m * j) % n) / (double)n;
            acc = wr_cadd(acc, wr_cmul(x[j], wr_cexp_i(ang)));
        }
        X[m] = acc;
    }
}

int main(void)
{
    /* invalid sizes */
    CHECK(wr_fft_plan_create(0) == NULL);
    CHECK(wr_fft_plan_create(3) == NULL);
    CHECK(wr_fft_plan_create(12) == NULL);

    /* n = 1 is the identity */
    wr_fft_plan *p1 = wr_fft_plan_create(1);
    wr_cplx one = wr_c(2.0, -3.0);
    CHECK(p1 && wr_fft_execute(p1, &one, 0) == 0);
    CHECK_NEAR(one.re, 2.0, 0); CHECK_NEAR(one.im, -3.0, 0);
    wr_fft_plan_destroy(p1);

    /* impulse -> flat spectrum */
    wr_fft_plan *p8 = wr_fft_plan_create(8);
    wr_cplx imp[8] = {{1, 0}};
    wr_fft_execute(p8, imp, 0);
    for (int k = 0; k < 8; k++) { CHECK_NEAR(imp[k].re, 1.0, 1e-15); CHECK_NEAR(imp[k].im, 0.0, 1e-15); }
    wr_fft_plan_destroy(p8);

    /* agreement with the O(n^2) DFT, several sizes */
    unsigned long long seed = 12345;
    for (size_t n = 2; n <= 256; n <<= 1) {
        wr_fft_plan *p = wr_fft_plan_create(n);
        wr_cplx *x = malloc(n * sizeof *x), *y = malloc(n * sizeof *y), *X = malloc(n * sizeof *X);
        for (size_t i = 0; i < n; i++) x[i] = y[i] = wr_c(rnd(&seed), rnd(&seed));
        naive_dft(x, X, n);
        wr_fft_execute(p, y, 0);
        double err = 0;
        for (size_t i = 0; i < n; i++) err = fmax(err, wr_cabs(wr_csub(y[i], X[i])));
        CHECK_NEAR(err, 0.0, 1e-12 * (double)n);
        free(x); free(y); free(X); wr_fft_plan_destroy(p);
    }

    /* single tone x[j] = exp(+2 pi i m j / n)  ->  X[m] = n, others 0 */
    {
        size_t n = 64, m = 5;
        wr_fft_plan *p = wr_fft_plan_create(n);
        wr_cplx x[64];
        for (size_t j = 0; j < n; j++) x[j] = wr_cexp_i(WR_TWO_PI * (double)((m * j) % n) / (double)n);
        wr_fft_execute(p, x, 0);
        for (size_t k = 0; k < n; k++)
            CHECK_NEAR(wr_cabs(x[k]), k == m ? (double)n : 0.0, 1e-12);
        wr_fft_plan_destroy(p);
    }

    /* round trip and Parseval on a large transform (2^16) */
    {
        size_t n = 1u << 16;
        wr_fft_plan *p = wr_fft_plan_create(n);
        wr_cplx *x = malloc(n * sizeof *x), *y = malloc(n * sizeof *y);
        double e_time = 0;
        for (size_t i = 0; i < n; i++) { x[i] = y[i] = wr_c(rnd(&seed), rnd(&seed)); e_time += wr_cabs2(x[i]); }
        wr_fft_execute(p, y, 0);
        double e_freq = 0;
        for (size_t i = 0; i < n; i++) e_freq += wr_cabs2(y[i]);
        CHECK_NEAR(e_freq / (double)n / e_time, 1.0, 1e-12);          /* Parseval */
        wr_fft_execute(p, y, 1);
        double err = 0;
        for (size_t i = 0; i < n; i++) err = fmax(err, wr_cabs(wr_csub(y[i], x[i])));
        CHECK_NEAR(err, 0.0, 1e-12);                                  /* inverse(forward(x)) = x */
        free(x); free(y); wr_fft_plan_destroy(p);
    }

    TEST_REPORT("fft");
}
