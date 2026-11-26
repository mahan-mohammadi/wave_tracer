/* Minimal test harness: no external dependencies. */
#ifndef WR_TEST_H
#define WR_TEST_H

#include <stdio.h>
#include <math.h>

static int wr_checks = 0, wr_failures = 0;

#define CHECK(cond) do { wr_checks++; if (!(cond)) { wr_failures++; \
    fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } } while (0)

/* Fails on NaN as well. */
#define CHECK_NEAR(actual, expected, tol) do { \
    double a_ = (actual), e_ = (expected), t_ = (tol); wr_checks++; \
    if (!(fabs(a_ - e_) <= t_)) { wr_failures++; \
        fprintf(stderr, "FAIL %s:%d: %s = %.17g, expected %.17g (tol %g)\n", \
                __FILE__, __LINE__, #actual, a_, e_, t_); } } while (0)

#define TEST_REPORT(name) do { \
    printf("%-12s %d checks, %d failed\n", name, wr_checks, wr_failures); \
    return wr_failures ? 1 : 0; } while (0)

#endif
