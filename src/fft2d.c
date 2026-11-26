#include <stdlib.h>
#include "wr_fft2d.h"
#include "wr_fft.h"

struct wr_fft2_plan {
    size_t nx, ny;
    wr_fft_plan *px, *py;
};

wr_fft2_plan *wr_fft2_plan_create(size_t nx, size_t ny)
{
    wr_fft2_plan *p = malloc(sizeof *p);
    if (!p) return NULL;
    p->nx = nx; p->ny = ny;
    p->px = wr_fft_plan_create(nx);
    p->py = wr_fft_plan_create(ny);
    if (!p->px || !p->py) { wr_fft2_plan_destroy(p); return NULL; }
    return p;
}

void wr_fft2_plan_destroy(wr_fft2_plan *p)
{
    if (!p) return;
    wr_fft_plan_destroy(p->px);
    wr_fft_plan_destroy(p->py);
    free(p);
}

size_t wr_fft2_plan_nx(const wr_fft2_plan *p) { return p->nx; }
size_t wr_fft2_plan_ny(const wr_fft2_plan *p) { return p->ny; }

int wr_fft2_execute(const wr_fft2_plan *p, wr_cplx *data, int inverse)
{
    if (!p || !data) return 1;
    const long nx = (long)p->nx, ny = (long)p->ny;
    int err = 0;

    #pragma omp parallel for schedule(static)
    for (long iy = 0; iy < ny; iy++)
        wr_fft_execute(p->px, data + (size_t)iy * p->nx, inverse);

    if (ny > 1) {
        #pragma omp parallel
        {
            wr_cplx *buf = malloc(p->ny * sizeof *buf);
            if (!buf) {
                #pragma omp atomic write
                err = 1;
            }
            #pragma omp for schedule(static)
            for (long ix = 0; ix < nx; ix++) {
                if (!buf) continue;
                for (long iy = 0; iy < ny; iy++) buf[iy] = data[(size_t)iy * p->nx + (size_t)ix];
                wr_fft_execute(p->py, buf, inverse);
                for (long iy = 0; iy < ny; iy++) data[(size_t)iy * p->nx + (size_t)ix] = buf[iy];
            }
            free(buf);
        }
    }
    return err;
}
