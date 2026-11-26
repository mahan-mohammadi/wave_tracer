/* wr_integrate.h - 1D numerical integration (Phase 1).
 *
 * Composite trapezoid: error ~ h^2 f''   (order 2)
 * Composite Simpson:   error ~ h^4 f'''' (order 4), needs an even panel count
 * Adaptive Simpson:    recursive subdivision with Richardson correction;
 *                      tolerance is halved at each split.
 * All use compensated summation. Assumes f is smooth on each panel; for
 * rapidly oscillating integrands (diffraction integrals) choose n so that
 * there are many panels per wavelength.
 */
#ifndef WR_INTEGRATE_H
#define WR_INTEGRATE_H

#include <stddef.h>
#include "wr_complex.h"

typedef double  (*wr_fn)(double x, void *ctx);
typedef wr_cplx (*wr_fn_c)(double x, void *ctx);

/* n = number of panels (n >= 1). Returns NAN on invalid n. */
double wr_trapezoid(wr_fn f, void *ctx, double a, double b, size_t n);

/* n = number of panels, must be even and >= 2. Returns NAN otherwise. */
double  wr_simpson(wr_fn f, void *ctx, double a, double b, size_t n);
wr_cplx wr_simpson_c(wr_fn_c f, void *ctx, double a, double b, size_t n);

/* Absolute tolerance tol; recursion stops at max_depth (e.g. 40). */
double wr_adaptive_simpson(wr_fn f, void *ctx, double a, double b,
                           double tol, int max_depth);

#endif /* WR_INTEGRATE_H */
