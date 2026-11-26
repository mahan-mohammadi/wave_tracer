#include "wr_optics.h"
#include "wr_sum.h"

double wr_frequency_from_wavelength(double lambda0) { return WR_C0 / lambda0; }
double wr_wavelength_from_frequency(double f)       { return WR_C0 / f; }
double wr_wavenumber(double lambda0)                { return WR_TWO_PI / lambda0; }
double wr_angular_frequency(double lambda0)         { return WR_TWO_PI * WR_C0 / lambda0; }
double wr_wavelength_in_medium(double lambda0, double n) { return lambda0 / n; }

double wr_opl(double n, double d) { return n * d; }

double wr_opl_path(const double *n, const double *d, size_t count)
{
    wr_ksum k = wr_ksum_zero();
    for (size_t i = 0; i < count; i++) wr_ksum_add(&k, n[i] * d[i]);
    return wr_ksum_value(&k);
}

double wr_phase_from_opl(double opl, double lambda0)
{
    double cycles = opl / lambda0;
    double frac = cycles - floor(cycles);
    if (frac >= 1.0) frac = 0.0; /* tiny negative cycles can round up to 1 */
    return WR_TWO_PI * frac;
}

double wr_wrap_phase(double phi)
{
    phi = fmod(phi, WR_TWO_PI);
    if (phi > WR_PI) phi -= WR_TWO_PI;
    else if (phi <= -WR_PI) phi += WR_TWO_PI;
    return phi;
}

wr_cplx wr_phasor(double amp, double opl, double lambda0)
{
    return wr_cpolar(amp, wr_phase_from_opl(opl, lambda0));
}

double wr_intensity(wr_cplx e) { return wr_cabs2(e); }

wr_cplx wr_coherent_sum(const wr_cplx *e, size_t count)
{
    wr_ksum re = wr_ksum_zero(), im = wr_ksum_zero();
    for (size_t i = 0; i < count; i++) {
        wr_ksum_add(&re, e[i].re);
        wr_ksum_add(&im, e[i].im);
    }
    return wr_c(wr_ksum_value(&re), wr_ksum_value(&im));
}

double wr_incoherent_sum(const wr_cplx *e, size_t count)
{
    wr_ksum k = wr_ksum_zero();
    for (size_t i = 0; i < count; i++) wr_ksum_add(&k, wr_cabs2(e[i]));
    return wr_ksum_value(&k);
}

double wr_two_beam_intensity(double i1, double i2, double delta)
{
    return i1 + i2 + 2.0 * sqrt(i1 * i2) * cos(delta);
}
