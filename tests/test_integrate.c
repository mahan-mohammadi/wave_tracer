#include "wr_test.h"
#include "wr_integrate.h"
#include "wr_optics.h"

static double f_sin(double x, void *c)  { (void)c; return sin(x); }
static double f_exp(double x, void *c)  { (void)c; return exp(x); }
static double f_sqrt(double x, void *c) { (void)c; return sqrt(x); }
static double f_gauss(double x, void *c){ (void)c; return exp(-x * x); }
static wr_cplx f_eix(double x, void *c) { (void)c; return wr_cexp_i(x); }
static wr_cplx f_e2pi(double x, void *c){ (void)c; return wr_cexp_i(WR_TWO_PI * x); }

int main(void)
{
    /* integral of sin over [0, pi] = 2 */
    CHECK_NEAR(wr_trapezoid(f_sin, NULL, 0, WR_PI, 1000), 2.0, 2e-6);
    CHECK_NEAR(wr_simpson(f_sin, NULL, 0, WR_PI, 100), 2.0, 2e-8);

    /* convergence order: halve h -> trapezoid error /4, Simpson error /16 */
    const double exact = exp(1.0) - 1.0;
    double t1 = fabs(wr_trapezoid(f_exp, NULL, 0, 1, 64)  - exact);
    double t2 = fabs(wr_trapezoid(f_exp, NULL, 0, 1, 128) - exact);
    double s1 = fabs(wr_simpson(f_exp, NULL, 0, 1, 16) - exact);
    double s2 = fabs(wr_simpson(f_exp, NULL, 0, 1, 32) - exact);
    CHECK_NEAR(t1 / t2, 4.0, 0.1);
    CHECK_NEAR(s1 / s2, 16.0, 1.0);

    /* invalid panel counts */
    CHECK(isnan(wr_simpson(f_sin, NULL, 0, 1, 3)));
    CHECK(isnan(wr_simpson(f_sin, NULL, 0, 1, 0)));
    CHECK(isnan(wr_trapezoid(f_sin, NULL, 0, 1, 0)));

    /* adaptive: smooth, non-smooth (sqrt at 0), and a Gaussian */
    CHECK_NEAR(wr_adaptive_simpson(f_sin, NULL, 0, WR_PI, 1e-10, 40), 2.0, 1e-9);
    CHECK_NEAR(wr_adaptive_simpson(f_sqrt, NULL, 0, 1, 1e-9, 40), 2.0 / 3.0, 1e-8);
    CHECK_NEAR(wr_adaptive_simpson(f_gauss, NULL, -10, 10, 1e-10, 40),
               1.7724538509055160273, 1e-9);                 /* sqrt(pi) */

    /* complex: int_0^{pi/2} e^{ix} dx = 1 + i ; int_0^1 e^{2 pi i x} dx = 0 */
    wr_cplx c1 = wr_simpson_c(f_eix, NULL, 0, WR_PI / 2, 200);
    CHECK_NEAR(c1.re, 1.0, 1e-9); CHECK_NEAR(c1.im, 1.0, 1e-9);
    wr_cplx c2 = wr_simpson_c(f_e2pi, NULL, 0, 1, 200);
    CHECK_NEAR(c2.re, 0.0, 1e-8); CHECK_NEAR(c2.im, 0.0, 1e-8);

    TEST_REPORT("integrate");
}
