/* wr_fft.h - 1D radix-2 FFT (Phase 1).
 *
 * Convention (matches the project's exp(-i w t) time convention for the
 * spatial transform used by the angular-spectrum method):
 *   forward: X[m] = sum_j x[j] exp(-2 pi i m j / n)
 *   inverse: x[j] = (1/n) sum_m X[m] exp(+2 pi i m j / n)
 * Length must be a power of two. A plan holds a read-only twiddle table
 * computed directly with cos/sin (no recurrence, so no error growth), so one
 * plan may be shared by many OpenMP threads, each transforming its own data.
 * Complexity O(n log n); in place.
 */
#ifndef WR_FFT_H
#define WR_FFT_H

#include <stddef.h>
#include "wr_complex.h"

typedef struct wr_fft_plan wr_fft_plan;

/* Returns NULL if n is 0, not a power of two, or allocation fails. */
wr_fft_plan *wr_fft_plan_create(size_t n);
void         wr_fft_plan_destroy(wr_fft_plan *plan);
size_t       wr_fft_plan_size(const wr_fft_plan *plan);

/* In-place transform of plan->n elements. inverse != 0 includes the 1/n
 * scaling. Returns 0 on success, nonzero on bad arguments. */
int wr_fft_execute(const wr_fft_plan *plan, wr_cplx *x, int inverse);

#endif /* WR_FFT_H */
