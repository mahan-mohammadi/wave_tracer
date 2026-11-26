/* wr_sum.h - Neumaier (improved Kahan) compensated summation.
 * Keeps the rounding error of a long sum at O(eps) instead of O(n*eps),
 * and, unlike plain Kahan, is also correct when a term is larger than the
 * running sum. Requires strict IEEE evaluation (no -ffast-math). */
#ifndef WR_SUM_H
#define WR_SUM_H

#include <math.h>

typedef struct {
    double s; /* running sum          */
    double c; /* running compensation */
} wr_ksum;

static inline wr_ksum wr_ksum_zero(void) { wr_ksum k = {0.0, 0.0}; return k; }

static inline void wr_ksum_add(wr_ksum *k, double x)
{
    double t = k->s + x;
    if (fabs(k->s) >= fabs(x))
        k->c += (k->s - t) + x;
    else
        k->c += (x - t) + k->s;
    k->s = t;
}

static inline double wr_ksum_value(const wr_ksum *k) { return k->s + k->c; }

#endif /* WR_SUM_H */
