/* wr_fft2d.h - 2D FFT built from the 1D plan (Phase 2).
 * Same sign/scaling convention as wr_fft.h: forward unscaled with exp(-i...),
 * inverse scaled by 1/(nx*ny). Row transforms then column transforms.
 * Parallelised with OpenMP over independent rows / columns: no shared
 * accumulation, so results are bit-identical for any thread count.
 * The plan is read-only and may be shared. Column pass gathers each column
 * into a thread-local buffer (simple, not cache-optimal; revisit in Phase 10).
 */
#ifndef WR_FFT2D_H
#define WR_FFT2D_H

#include <stddef.h>
#include "wr_complex.h"

typedef struct wr_fft2_plan wr_fft2_plan;

wr_fft2_plan *wr_fft2_plan_create(size_t nx, size_t ny); /* NULL if not powers of two */
void          wr_fft2_plan_destroy(wr_fft2_plan *p);
size_t        wr_fft2_plan_nx(const wr_fft2_plan *p);
size_t        wr_fft2_plan_ny(const wr_fft2_plan *p);

/* data: ny rows of nx elements, x fastest. Returns 0 on success. */
int wr_fft2_execute(const wr_fft2_plan *p, wr_cplx *data, int inverse);

#endif
