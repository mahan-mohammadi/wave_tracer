#include <stdlib.h>
#include "wr_fft.h"
#include "wr_optics.h" /* WR_TWO_PI */

struct wr_fft_plan {
    size_t n;
    wr_cplx *tw; /* tw[k] = exp(-2 pi i k / n), k = 0 .. n/2-1 */
};

wr_fft_plan *wr_fft_plan_create(size_t n)
{
    if (n == 0 || (n & (n - 1)) != 0) return NULL;
    wr_fft_plan *p = malloc(sizeof *p);
    if (!p) return NULL;
    size_t half = n / 2;
    p->n = n;
    p->tw = malloc((half ? half : 1) * sizeof *p->tw);
    if (!p->tw) { free(p); return NULL; }
    for (size_t k = 0; k < half; k++) {
        double ang = WR_TWO_PI * (double)k / (double)n;
        p->tw[k] = wr_c(cos(ang), -sin(ang));
    }
    return p;
}

void wr_fft_plan_destroy(wr_fft_plan *p)
{
    if (!p) return;
    free(p->tw);
    free(p);
}

size_t wr_fft_plan_size(const wr_fft_plan *p) { return p ? p->n : 0; }

int wr_fft_execute(const wr_fft_plan *p, wr_cplx *x, int inverse)
{
    if (!p || !x) return 1;
    const size_t n = p->n;

    /* bit-reversal permutation */
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { wr_cplx t = x[i]; x[i] = x[j]; x[j] = t; }
    }

    /* iterative Cooley-Tukey butterflies */
    for (size_t len = 2; len <= n; len <<= 1) {
        size_t half = len >> 1, step = n / len;
        for (size_t i = 0; i < n; i += len) {
            for (size_t j = 0; j < half; j++) {
                wr_cplx w = p->tw[j * step];
                if (inverse) w = wr_cconj(w);
                wr_cplx u = x[i + j];
                wr_cplx v = wr_cmul(x[i + j + half], w);
                x[i + j]        = wr_cadd(u, v);
                x[i + j + half] = wr_csub(u, v);
            }
        }
    }

    if (inverse) {
        double s = 1.0 / (double)n;
        for (size_t i = 0; i < n; i++) x[i] = wr_cscale(x[i], s);
    }
    return 0;
}
