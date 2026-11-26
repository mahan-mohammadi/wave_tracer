#include "wr_test.h"
#include "wr_complex.h"
#include "wr_optics.h"

int main(void)
{
    wr_cplx a = wr_c(1, 2), b = wr_c(3, -4);

    wr_cplx p = wr_cmul(a, b); /* (1+2i)(3-4i) = 11 + 2i */
    CHECK_NEAR(p.re, 11, 1e-15); CHECK_NEAR(p.im, 2, 1e-15);

    wr_cplx q = wr_cdiv(p, b); /* back to a */
    CHECK_NEAR(q.re, 1, 1e-15); CHECK_NEAR(q.im, 2, 1e-15);

    /* Smith division survives where the naive formula overflows */
    wr_cplx big = wr_cdiv(wr_c(1e200, 1e200), wr_c(1e200, 1e200));
    CHECK_NEAR(big.re, 1, 1e-15); CHECK_NEAR(big.im, 0, 1e-15);
    wr_cplx tiny = wr_cdiv(wr_c(1e-200, 0), wr_c(1e-200, 1e-200));
    CHECK_NEAR(tiny.re, 0.5, 1e-15); CHECK_NEAR(tiny.im, -0.5, 1e-15);

    CHECK_NEAR(wr_cabs(wr_c(3, 4)), 5, 1e-15);
    CHECK_NEAR(wr_cabs2(wr_c(3, 4)), 25, 1e-15);
    CHECK_NEAR(wr_cabs(wr_c(1e200, 1e200)), 1e200 * sqrt(2.0), 1e186);
    CHECK_NEAR(wr_carg(wr_c(0, 1)), WR_PI / 2, 1e-15);
    CHECK_NEAR(wr_carg(wr_c(-1, 0)), WR_PI, 1e-15);

    wr_cplx z = wr_cpolar(2.5, 0.7);
    CHECK_NEAR(wr_cabs(z), 2.5, 1e-15); CHECK_NEAR(wr_carg(z), 0.7, 1e-15);

    wr_cplx e = wr_cexp_i(WR_PI); /* Euler identity */
    CHECK_NEAR(e.re, -1, 1e-15); CHECK_NEAR(e.im, 0, 2e-16);
    wr_cplx e2 = wr_cexp(wr_c(1.0, WR_PI / 2)); /* e * i */
    CHECK_NEAR(e2.re, 0, 2e-16 * 3); CHECK_NEAR(e2.im, exp(1.0), 1e-15);

    wr_cplx s1 = wr_csqrt(wr_c(-1, 0));
    CHECK_NEAR(s1.re, 0, 1e-15); CHECK_NEAR(s1.im, 1, 1e-15);
    wr_cplx s2 = wr_csqrt(wr_c(3, 4));
    CHECK_NEAR(s2.re, 2, 1e-15); CHECK_NEAR(s2.im, 1, 1e-15);
    wr_cplx s3 = wr_csqrt(wr_c(-3, 4));
    CHECK_NEAR(s3.re, 1, 1e-15); CHECK_NEAR(s3.im, 2, 1e-15);
    wr_cplx s4 = wr_csqrt(wr_c(-3, -4)); /* lower half-plane keeps Im < 0 */
    CHECK_NEAR(s4.re, 1, 1e-15); CHECK_NEAR(s4.im, -2, 1e-15);
    wr_cplx sq = wr_cmul(s2, s2);
    CHECK_NEAR(sq.re, 3, 1e-14); CHECK_NEAR(sq.im, 4, 1e-14);

    wr_cplx cj = wr_cconj(a);
    CHECK_NEAR(cj.re, 1, 0); CHECK_NEAR(cj.im, -2, 0);

    TEST_REPORT("complex");
}
