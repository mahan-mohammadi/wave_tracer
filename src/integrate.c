#include "wr_integrate.h"
#include "wr_sum.h"

double wr_trapezoid(wr_fn f, void *ctx, double a, double b, size_t n)
{
    if (n < 1) return NAN;
    double h = (b - a) / (double)n;
    wr_ksum k = wr_ksum_zero();
    wr_ksum_add(&k, 0.5 * (f(a, ctx) + f(b, ctx)));
    for (size_t i = 1; i < n; i++) wr_ksum_add(&k, f(a + (double)i * h, ctx));
    return wr_ksum_value(&k) * h;
}

double wr_simpson(wr_fn f, void *ctx, double a, double b, size_t n)
{
    if (n < 2 || (n & 1u)) return NAN;
    double h = (b - a) / (double)n;
    wr_ksum k = wr_ksum_zero();
    wr_ksum_add(&k, f(a, ctx));
    wr_ksum_add(&k, f(b, ctx));
    for (size_t i = 1; i < n; i++)
        wr_ksum_add(&k, ((i & 1u) ? 4.0 : 2.0) * f(a + (double)i * h, ctx));
    return wr_ksum_value(&k) * h / 3.0;
}

wr_cplx wr_simpson_c(wr_fn_c f, void *ctx, double a, double b, size_t n)
{
    if (n < 2 || (n & 1u)) return wr_c(NAN, NAN);
    double h = (b - a) / (double)n;
    wr_ksum re = wr_ksum_zero(), im = wr_ksum_zero();
    wr_cplx fa = f(a, ctx), fb = f(b, ctx);
    wr_ksum_add(&re, fa.re + fb.re);
    wr_ksum_add(&im, fa.im + fb.im);
    for (size_t i = 1; i < n; i++) {
        double w = (i & 1u) ? 4.0 : 2.0;
        wr_cplx v = f(a + (double)i * h, ctx);
        wr_ksum_add(&re, w * v.re);
        wr_ksum_add(&im, w * v.im);
    }
    return wr_c(wr_ksum_value(&re) * h / 3.0, wr_ksum_value(&im) * h / 3.0);
}

static double asr(wr_fn f, void *ctx, double a, double b,
                  double fa, double fm, double fb, double whole,
                  double tol, int depth)
{
    double m = 0.5 * (a + b);
    double lm = 0.5 * (a + m), rm = 0.5 * (m + b);
    double flm = f(lm, ctx), frm = f(rm, ctx);
    double left  = (m - a) / 6.0 * (fa + 4.0 * flm + fm);
    double right = (b - m) / 6.0 * (fm + 4.0 * frm + fb);
    double delta = left + right - whole;
    if (depth <= 0 || fabs(delta) <= 15.0 * tol)
        return left + right + delta / 15.0;
    return asr(f, ctx, a, m, fa, flm, fm, left,  0.5 * tol, depth - 1)
         + asr(f, ctx, m, b, fm, frm, fb, right, 0.5 * tol, depth - 1);
}

double wr_adaptive_simpson(wr_fn f, void *ctx, double a, double b,
                           double tol, int max_depth)
{
    double fa = f(a, ctx), fb = f(b, ctx), m = 0.5 * (a + b), fm = f(m, ctx);
    double whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb);
    return asr(f, ctx, a, b, fa, fm, fb, whole, tol, max_depth);
}
