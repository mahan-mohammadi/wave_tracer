#include "wr_test.h"
#include "wr_optics.h"
#include "wr_sum.h"

int main(void)
{
    const double lam = 500e-9;

    /* wavelength <-> frequency (c = lambda f) */
    double f = wr_frequency_from_wavelength(lam);
    CHECK_NEAR(f, 299792458.0 / 500e-9, 1.0);          /* 5.99584916e14 Hz */
    CHECK_NEAR(wr_wavelength_from_frequency(f), lam, 1e-22);
    CHECK_NEAR(wr_wavenumber(lam) * lam, WR_TWO_PI, 1e-14);
    CHECK_NEAR(wr_angular_frequency(lam), WR_TWO_PI * f, 1e-1);
    CHECK_NEAR(wr_wavelength_in_medium(lam, 1.5), 500e-9 / 1.5, 1e-22);

    /* OPL: air gap + glass slab + air gap */
    double n[3] = {1.0, 1.5, 1.0}, d[3] = {0.010, 0.002, 0.010};
    CHECK_NEAR(wr_opl_path(n, d, 3), 0.010 + 0.003 + 0.010, 1e-17);

    /* phase from OPL, including a path of ~1e6 wavelengths */
    CHECK_NEAR(wr_phase_from_opl(0.25 * lam, lam), WR_PI / 2, 1e-12);
    CHECK_NEAR(wr_phase_from_opl(1000000.25 * lam, lam), WR_PI / 2, 1e-8);
    CHECK_NEAR(wr_phase_from_opl(-0.25 * lam, lam), 1.5 * WR_PI, 1e-12);
    CHECK(wr_phase_from_opl(-1e-30, lam) < WR_TWO_PI);   /* never returns 2 pi */
    CHECK_NEAR(wr_wrap_phase(3.5 * WR_PI), -0.5 * WR_PI, 1e-12);
    CHECK_NEAR(wr_wrap_phase(-3.5 * WR_PI), 0.5 * WR_PI, 1e-12);
    CHECK_NEAR(wr_wrap_phase(WR_PI), WR_PI, 0);

    /* two equal coherent beams: I = 4 cos^2(delta/2); matches analytic formula */
    for (int i = 0; i <= 16; i++) {
        double delta = WR_TWO_PI * i / 16.0;
        wr_cplx e[2] = {wr_phasor(1.0, 0.0, lam), wr_phasor(1.0, delta / WR_TWO_PI * lam, lam)};
        double coh = wr_intensity(wr_coherent_sum(e, 2));
        CHECK_NEAR(coh, 4.0 * cos(delta / 2) * cos(delta / 2), 1e-12);
        CHECK_NEAR(coh, wr_two_beam_intensity(1, 1, delta), 1e-12);
        CHECK_NEAR(wr_incoherent_sum(e, 2), 2.0, 1e-15);   /* phase-independent */
    }

    /* fully destructive: lambda/2 path difference */
    wr_cplx de[2] = {wr_phasor(1.0, 0.0, lam), wr_phasor(1.0, 0.5 * lam, lam)};
    CHECK_NEAR(wr_intensity(wr_coherent_sum(de, 2)), 0.0, 1e-30);

    /* 7 equal phasors evenly spaced around the circle cancel exactly */
    wr_cplx ring[7];
    for (int k = 0; k < 7; k++) ring[k] = wr_cexp_i(WR_TWO_PI * k / 7.0);
    CHECK_NEAR(wr_intensity(wr_coherent_sum(ring, 7)), 0.0, 1e-28);
    CHECK_NEAR(wr_incoherent_sum(ring, 7), 7.0, 1e-14);

    /* compensated summation recovers what naive summation loses */
    wr_cplx cs[3] = {wr_c(1e16, 0), wr_c(1.0, 0), wr_c(-1e16, 0)};
    double naive = (1e16 + 1.0) + -1e16;
    wr_cplx sum = wr_coherent_sum(cs, 3);
    CHECK_NEAR(naive, 0.0, 0);
    CHECK_NEAR(sum.re, 1.0, 0);

    TEST_REPORT("optics");
}
