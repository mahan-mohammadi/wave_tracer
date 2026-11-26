/* wr_optics.h - phase, optical path length, wavelength/frequency, and
 * complex-amplitude superposition (Phase 1).
 *
 * Model (scalar, monochromatic):
 *   E = A exp(i phi),   phi = k0 * OPL = 2 pi * OPL / lambda0,
 *   OPL = sum_j n_j d_j          (n_j real here; absorption comes later)
 *   lambda_medium = lambda0 / n,  c = lambda0 * f  (vacuum values).
 * Intensity is defined as I = |E|^2 in the project's field units; physical
 * irradiance would be (1/2) n eps0 c |E|^2, a constant factor we do not need
 * until absolute radiometry matters.
 *
 * Coherent light:    E_total = sum E_j,  I = |E_total|^2   (phases matter)
 * Incoherent light:  I_total = sum |E_j|^2                 (phases ignored)
 * Never add intensities when the beams are coherent (roadmap rule 5).
 */
#ifndef WR_OPTICS_H
#define WR_OPTICS_H

#include <stddef.h>
#include "wr_complex.h"

#define WR_C0     299792458.0                 /* speed of light in vacuum, m/s (exact) */
#define WR_PI     3.14159265358979323846
#define WR_TWO_PI 6.28318530717958647692

/* Vacuum conversions. lambda in metres, f in Hz. */
double wr_frequency_from_wavelength(double lambda0);
double wr_wavelength_from_frequency(double f);
double wr_wavenumber(double lambda0);          /* k0 = 2 pi / lambda0, rad/m */
double wr_angular_frequency(double lambda0);   /* omega = 2 pi c / lambda0   */
double wr_wavelength_in_medium(double lambda0, double n);

/* Optical path length. */
double wr_opl(double n, double d);
double wr_opl_path(const double *n, const double *d, size_t count);

/* Phase handling. wr_phase_from_opl reduces OPL/lambda to its fractional
 * part BEFORE multiplying by 2 pi, so precision is kept even when the path
 * contains millions of wavelengths. Result in [0, 2 pi). */
double wr_phase_from_opl(double opl, double lambda0);
double wr_wrap_phase(double phi);              /* result in (-pi, pi] */

/* Complex amplitude after travelling a given OPL: amp * exp(i k0 OPL). */
wr_cplx wr_phasor(double amp, double opl, double lambda0);

/* Superposition and intensity. */
double  wr_intensity(wr_cplx e);                              /* |E|^2 */
wr_cplx wr_coherent_sum(const wr_cplx *e, size_t count);      /* compensated */
double  wr_incoherent_sum(const wr_cplx *e, size_t count);    /* compensated */

/* Analytic two-beam interference: I1 + I2 + 2 sqrt(I1 I2) cos(delta). */
double wr_two_beam_intensity(double i1, double i2, double delta);

#endif /* WR_OPTICS_H */
