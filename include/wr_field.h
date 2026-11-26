/* wr_field.h - sampled complex scalar field on a regular 2D grid (Phase 2).
 *
 * Storage: row-major, x fastest: data[iy*nx + ix] (array of wr_cplx).
 * Pixel centres: x = (ix - nx/2 + 0.5)*dx, y = (iy - ny/2 + 0.5)*dy, so the
 * grid is symmetric about the origin and pixel EDGES fall on multiples of dx,
 * dy (an aperture edge at such a position is represented exactly).
 * A 1D problem is a field with ny = 1.
 * SI units. nx, ny must be powers of two (FFT requirement).
 */
#ifndef WR_FIELD_H
#define WR_FIELD_H

#include <stddef.h>
#include "wr_complex.h"

typedef struct {
    size_t   nx, ny;
    double   dx, dy;
    wr_cplx *data;
} wr_field;

wr_field *wr_field_create(size_t nx, size_t ny, double dx, double dy); /* zero-filled; NULL on failure */
void      wr_field_destroy(wr_field *f);

static inline wr_cplx *wr_field_at(wr_field *f, size_t ix, size_t iy) { return &f->data[iy * f->nx + ix]; }
static inline double wr_field_x(const wr_field *f, size_t ix) { return ((double)ix - (double)(f->nx / 2) + 0.5) * f->dx; }
static inline double wr_field_y(const wr_field *f, size_t iy) { return ((double)iy - (double)(f->ny / 2) + 0.5) * f->dy; }

/* Total power: sum |U|^2 * dx * dy (compensated). Conserved by propagation. */
double wr_field_power(const wr_field *f);

#endif
