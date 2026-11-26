#include <stdlib.h>
#include "wr_field.h"
#include "wr_sum.h"

wr_field *wr_field_create(size_t nx, size_t ny, double dx, double dy)
{
    if (nx == 0 || ny == 0 || !(dx > 0.0) || !(dy > 0.0)) return NULL;
    wr_field *f = malloc(sizeof *f);
    if (!f) return NULL;
    f->nx = nx; f->ny = ny; f->dx = dx; f->dy = dy;
    f->data = calloc(nx * ny, sizeof *f->data);
    if (!f->data) { free(f); return NULL; }
    return f;
}

void wr_field_destroy(wr_field *f)
{
    if (!f) return;
    free(f->data);
    free(f);
}

double wr_field_power(const wr_field *f)
{
    wr_ksum k = wr_ksum_zero();
    for (size_t i = 0; i < f->nx * f->ny; i++) wr_ksum_add(&k, wr_cabs2(f->data[i]));
    return wr_ksum_value(&k) * f->dx * f->dy;
}
