/* wr_complex.h - complex numbers for the wave renderer (Phase 1).
 *
 * Own struct instead of <complex.h> so layout is explicit (two doubles,
 * array-of-structs, cache friendly) and behaviour does not depend on the
 * compiler's complex-number extensions.
 *
 * Time convention used throughout the project: fields vary as exp(-i w t),
 * so a wave travelling a distance L accumulates phase +k*L.
 *
 * Do NOT compile with -ffast-math: compensated summation relies on strict
 * IEEE-754 evaluation order.
 */
#ifndef WR_COMPLEX_H
#define WR_COMPLEX_H

#include <math.h>

typedef struct {
    double re;
    double im;
} wr_cplx;

static inline wr_cplx wr_c(double re, double im) { wr_cplx z = {re, im}; return z; }

static inline wr_cplx wr_cadd(wr_cplx a, wr_cplx b) { return wr_c(a.re + b.re, a.im + b.im); }
static inline wr_cplx wr_csub(wr_cplx a, wr_cplx b) { return wr_c(a.re - b.re, a.im - b.im); }
static inline wr_cplx wr_cscale(wr_cplx a, double s) { return wr_c(a.re * s, a.im * s); }
static inline wr_cplx wr_cconj(wr_cplx a) { return wr_c(a.re, -a.im); }

static inline wr_cplx wr_cmul(wr_cplx a, wr_cplx b)
{
    return wr_c(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

/* Smith's algorithm: avoids overflow/underflow of b.re^2 + b.im^2.
 * Division by exactly zero gives inf/NaN, as with real division. */
static inline wr_cplx wr_cdiv(wr_cplx a, wr_cplx b)
{
    if (fabs(b.re) >= fabs(b.im)) {
        double r = b.im / b.re;
        double d = b.re + b.im * r;
        return wr_c((a.re + a.im * r) / d, (a.im - a.re * r) / d);
    } else {
        double r = b.re / b.im;
        double d = b.re * r + b.im;
        return wr_c((a.re * r + a.im) / d, (a.im * r - a.re) / d);
    }
}

/* |z|^2: this is the intensity of a field amplitude (see wr_optics.h). */
static inline double wr_cabs2(wr_cplx a) { return a.re * a.re + a.im * a.im; }
/* |z| via hypot: no spurious overflow. */
static inline double wr_cabs(wr_cplx a) { return hypot(a.re, a.im); }
/* Principal argument in (-pi, pi]. */
static inline double wr_carg(wr_cplx a) { return atan2(a.im, a.re); }

static inline wr_cplx wr_cpolar(double r, double theta)
{
    return wr_c(r * cos(theta), r * sin(theta));
}

/* exp(i*theta) */
static inline wr_cplx wr_cexp_i(double theta) { return wr_c(cos(theta), sin(theta)); }

/* exp(z) = exp(re) * (cos(im) + i sin(im)) */
static inline wr_cplx wr_cexp(wr_cplx z)
{
    double e = exp(z.re);
    return wr_c(e * cos(z.im), e * sin(z.im));
}

/* Principal square root (Re >= 0). Needed for complex refractive indices
 * and for the angular-spectrum transfer function. */
static inline wr_cplx wr_csqrt(wr_cplx z)
{
    if (z.re == 0.0 && z.im == 0.0) return wr_c(0.0, z.im);
    double m = hypot(z.re, z.im);
    double t = sqrt(0.5 * (m + fabs(z.re)));
    if (z.re >= 0.0) return wr_c(t, z.im / (2.0 * t));
    return wr_c(fabs(z.im) / (2.0 * t), copysign(t, z.im));
}

#endif /* WR_COMPLEX_H */
