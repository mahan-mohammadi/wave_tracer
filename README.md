# Wave-optical renderer - Phases 1-3 (in progress): foundation, propagation, validation scenes

C17 + OpenMP. No dependencies beyond libc/libm.

## Build and test
    make test          # builds build/libwr.a and runs all tests
    make clean

Do not use -ffast-math (compensated summation needs strict IEEE order).

## Modules
| Header | Content |
|---|---|
| wr_complex.h | complex struct, add/sub/mul, Smith division, abs/arg, polar, exp, principal sqrt |
| wr_sum.h | Neumaier compensated summation |
| wr_optics.h | c, lambda<->f, k, omega, medium wavelength, OPL, phase from OPL, coherent/incoherent sums, two-beam formula |
| wr_integrate.h | trapezoid, Simpson (real and complex), adaptive Simpson |
| wr_fft.h | radix-2 1D FFT with reusable, thread-shareable plan |
| wr_field.h | complex scalar field on a regular grid (pixel centres, power) |
| wr_fft2d.h | 2D FFT, OpenMP over rows/columns (bit-identical for any thread count) |
| wr_propagate.h | angular-spectrum propagation: plain (exact on a periodic grid) and band-limited |

## Conventions
- Fields vary as exp(-i w t); propagating a distance L adds phase +k0*L.
- Intensity I = |E|^2 (field units; constant factors deferred).
- SI units (metres, Hz, radians).

## Validity ranges
- FFT: length must be a power of two. Accurate to ~1e-12 relative at n = 2^16 in tests.
- Simpson/trapezoid assume smooth integrands; use many panels per wavelength for oscillatory ones.
- Adaptive Simpson tolerance is absolute and heuristic, not a guaranteed bound.
- Scalar monochromatic model only; polarization, dispersion, absorption come in later phases.

## Phase 2: single-slit validation (tests/test_slit.c)
Slit a = 40 um, lambda = 500 nm, z = 20 mm, compared on |x| < 0.4 mm against the
Fresnel integral computed by numerical quadrature (incident intensity 1, peak 0.16).

| Case | max abs error |
|---|---|
| plain, window 2.048 mm (pattern wraps around) | 2-3e-2 (~15-20% of peak) |
| band-limited, window 2.048 mm | 3.4e-4 (floor, independent of dx) |
| plain, window 16.384 mm, dx 4 / 2 / 1 / 0.5 um | 5.6e-4 / 1.3e-4 / 3.3e-5 / 1.0e-5 (about second order in dx) |

Findings: the plain method is correct but the periodic window wraps high-angle
light, which adds coherently, so the window must hold the whole pattern. The
band-limited variant removes the wrap-around at the cost of dropping the
highest-angle light. The Fresnel reference itself is paraxial (valid for small angles).
Not yet validated: double slit, 2D apertures with non-separable geometry, anything beyond
free-space propagation.

## Phase 3, scene 2: double slit (tests/test_double_slit.c)
Two slits (a = 20 um, d = 100 um), lambda = 500 nm, z = 200 mm, window 32.8 mm, band-limited
propagation. Only the two apertures are specified; the fringes emerge from propagation + superposition.

| Check | Result |
|---|---|
| Linearity: two-slit field vs U1 + U2 | 4.7e-16 relative |
| Interference-phase period vs lambda z / d (1 mm) | 1.00003 mm (3e-5 rel.) |
| Intensity vs Fresnel integral (numerical quadrature), \|x\| < 2.5 mm | 1.2e-5 abs = 0.08% of peak (peak 0.0160 = 4a^2/(lambda z)) |
| Visibility, simulated / reference / incoherent sum | 1.0000 / 1.0000 / 0.0156 |
| Band-limited error, window 8.2 mm vs 32.8 mm | 1.6% vs 0.08% of peak |

Notes: raw intensity maxima are ~1.4% closer together than lambda z / d (0.986 mm) because the single-slit
envelope slopes across the fringes; the independent reference shows the same (0.985 mm). Adding intensities
instead of amplitudes removes the fringes (visibility 0.016). The plain (non band-limited) method needs a
window far larger than 2 z tan(theta) and gives ~16% error here, so use the band-limited variant for large z.
Known limitation: the band-limit cutoff is sharp (ringing); a smooth taper is a possible improvement.
Reference is paraxial. Output: build/double_slit.csv.

## Phase 3, scene 4: Michelson interferometer (tests/test_michelson.c)
Gaussian beam (w0 = 0.3 mm, 500 nm), 0.5 m arms, 256x256 grid at 40 um. Lossless 50:50 splitter:
E_det = (U1 - U2)/2, E_src = (U1 + U2)/2 (which port is dark is a coating convention; the ports are
complementary). Tilted mirror = phase ramp exp(i 2 k alpha x) at the mirror, so beam shear is included.

| Check | Result |
|---|---|
| Propagated beam vs closed-form Gaussian beam (complex field, 1 m) | 3.0e-8 |
| Tilted arm (50 urad, shear 50 um) vs exact paraxial tilt+shear identity | 4.1e-8 |
| P_det + P_src = P_in (aligned and tilted) | 1e-12 |
| P_det/P_in = (1 - cos(2 pi Delta / lambda))/2, 16 path differences | max error 3.9e-7 |
| Visibility vs exp(-sigma_k^2 Delta^2 / 2), sigma_lambda = 5 nm, Delta = 0/4/8/16 um | 1.00000/0.88132/0.60331/0.13248, agree to 1e-5 |

Limits: scalar; paraxial reference formulas; small-angle tilt model; ideal mirrors and splitter
(no dispersion, absorption, or polarization). The tilted-mirror port power (0.0156) has no
independent closed form in this test; the field itself is validated by check 2.
Partial coherence is handled by adding intensities across wavelengths and amplitudes within each.
