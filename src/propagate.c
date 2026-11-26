#include "wr_propagate.h"
#include "wr_optics.h"

/* FFT-order index -> frequency (cycles/m). Indices i < ceil(n/2) are
 * non-negative; for even n the Nyquist bin is treated as negative.
 * (For n = 1 the single sample is frequency 0.) */
static inline double freq(long i, long n, double d)
{
    long m = (i < (n + 1) / 2) ? i : i - n;
    return (double)m / ((double)n * d);
}

static int propagate(wr_field *f, const wr_fft2_plan *plan,
                     double lambda0, double n, double z, int band_limit)
{
    if (!f || !plan || !(lambda0 > 0.0) || !(n > 0.0)) return 1;
    if (wr_fft2_plan_nx(plan) != f->nx || wr_fft2_plan_ny(plan) != f->ny) return 1;

    const long nx = (long)f->nx, ny = (long)f->ny;
    const double lam = lambda0 / n;
    const double k = WR_TWO_PI / lam;

    /* Band limits, derived from |d(phase)/df| * df <= pi for each axis:
     *   f_lim = 1 / (lam * sqrt(1 + (2 |z| df)^2)),  df = 1/(N d). */
    double fxlim = 1e300, fylim = 1e300;
    if (band_limit) {
        double dfx = 1.0 / ((double)nx * f->dx), dfy = 1.0 / ((double)ny * f->dy);
        double az = fabs(z);
        fxlim = 1.0 / (lam * sqrt(1.0 + 4.0 * az * az * dfx * dfx));
        fylim = 1.0 / (lam * sqrt(1.0 + 4.0 * az * az * dfy * dfy));
    }

    if (wr_fft2_execute(plan, f->data, 0)) return 2;

    #pragma omp parallel for schedule(static)
    for (long iy = 0; iy < ny; iy++) {
        double fy = freq(iy, ny, f->dy);
        for (long ix = 0; ix < nx; ix++) {
            double fx = freq(ix, nx, f->dx);
            double s = 1.0 - lam * lam * (fx * fx + fy * fy);
            wr_cplx h;
            if (fabs(fx) > fxlim || fabs(fy) > fylim) {
                h = wr_c(0.0, 0.0);
            } else if (s >= 0.0) {
                h = wr_cexp_i(fmod(k * z * sqrt(s), WR_TWO_PI));
            } else if (z >= 0.0) {
                h = wr_c(exp(-k * z * sqrt(-s)), 0.0);
            } else {
                h = wr_c(0.0, 0.0);
            }
            wr_cplx *u = &f->data[(size_t)iy * f->nx + (size_t)ix];
            *u = wr_cmul(*u, h);
        }
    }

    return wr_fft2_execute(plan, f->data, 1) ? 2 : 0;
}

int wr_propagate_asm(wr_field *f, const wr_fft2_plan *plan,
                     double lambda0, double n, double z)
{
    return propagate(f, plan, lambda0, n, z, 0);
}

int wr_propagate_asm_bl(wr_field *f, const wr_fft2_plan *plan,
                        double lambda0, double n, double z)
{
    return propagate(f, plan, lambda0, n, z, 1);
}
