# Wave-optical renderer - Phases 1-3 done; Phase 4 (3D scene system) done

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
| wr_vec3.h | 3D vector ops, reflect, Snell refract (direction only, reports TIR), orthonormal basis |
| wr_xform.h | 3x3 matrices, Rodrigues rotation, affine transforms, inverse, normal transform |
| wr_ray.h | ray (non-unit d allowed, t in units of \|d\|), hit record (outward normal + front flag) |
| wr_aabb.h | bounding box, union/expand, slab ray test |
| wr_geom.h | sphere, plane, triangle: intersection + bounds. Geometry only, no optics |
| wr_mesh.h | indexed triangle mesh, OBJ parser/loader, box and UV-sphere generators, transform (mirror flips winding), area, signed volume |
| wr_bvh.h | binned-SAH BVH over a mesh, closest-hit and any-hit, thread-safe queries, brute-force reference |
| wr_camera.h | pinhole, orthographic, thin-lens ray generation; pinhole projection; concentric disc sampling |
| wr_material.h | complex refractive index n + i kappa: constant / Cauchy / Sellmeier, validity ranges, surface kinds (dielectric, ideal mirror, absorber) |
| wr_light.h | point / plane-wave / Gaussian-beam sources: line spectrum, amplitude, phase, coherence length, Jones vector, velocity; exact source fields |
| wr_scene.h | container: materials, spheres, planes, BVH meshes, lights, camera, ambient medium, time; closest-hit; validation |
| wr_sceneio.h | text scene format and loader (format documented in the header), errors with line numbers |

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

## Phase 3, scenes 5 and 6: thin film and Newton's rings (tests/test_film.c, tests/test_newton.c)
New module wr_film.h: normal-incidence Fresnel coefficients (pi shift = sign of r) and single-layer
reflection by explicit coherent sum of partial waves, checked against the closed-form Airy result.
Normal incidence, scalar, lossless, non-dispersive indices only.

Thin film / multiple reflections (reference = closed-form optics, independent of the code):
| Check | Result |
|---|---|
| Partial-wave series (200 terms) vs Airy, 6400 films | 9.5e-16 |
| R + T = 1 over the same films | 8.9e-16 |
| Soap film n=1.33, d=400 nm: maxima 425.600 / 709.333 nm, minimum 532.000 nm | exact to 0.001 nm |
| Peak reflectance 0.07711 vs closed form (two-beam estimate 0.08024 is 4% too high) | 1e-12 |
| Quarter-wave MgF2 coating, R = 1.26008% (formula 1.26008%), bare glass 4.258% | 1e-13 |

Newton's rings: lens R = 1 m on a flat plate, lambda = 500 nm, exact sphere gap h(r), 512x512 at 12.5 um.
No ring formula is in the code. Dark radii r_m = sqrt(m lambda R / n_gap), bright radii sqrt((m - 1/2) lambda R / n_gap).
| Case | Worst ring-radius error (14-16 dark + 14-16 bright rings) | I_max (all partial waves vs two-beam) |
|---|---|---|
| air gap | 0.187 um (ring spacing 90-700 um, pixel 12.5 um) | 0.14793 (two-beam 4r^2 would be 0.16) |
| water gap n=1.33 | 0.227 um | 0.01433 |
Centre is dark (pi shift) to < 1e-4 of the maximum. R + T - 1 < 2e-15 at every pixel.
Limits: thin-element (local film) approximation, valid while gap slope r/R << 1 (here <= 3.2e-3);
errors are bounded by pixel sampling (about 2% of a pixel) and there is a small unexplained sign bias for air.
Output: build/newton_rings.ppm.

## Phase 3, scene 3: Young's experiment with an extended incoherent source (tests/test_young.c)
Line source of width s (71 emitters 50 um apart, each a 1-pixel slit), z0 = 500 mm to two slits
(a = 6 um, d = 100 um), z = 200 mm to the screen. Source->slits and slits->screen both go through the
band-limited propagator. Amplitudes are added within an emitter, intensities across emitters.
Signed visibility V = (I(0) - I(lambda z/2d)) / (I(0) + I(lambda z/2d)); theory V = sinc(s d / (lambda z0)).
| s (mm) | u | V simulated | V reference (quadrature) | sinc(u) |
|---|---|---|---|---|
| 0.05 | 0.020 | +1.0000 | +1.0000 | +0.9993 |
| 0.55 | 0.220 | +0.9230 | +0.9231 | +0.9223 |
| 1.25 | 0.500 | +0.6380 | +0.6380 | +0.6366 |
| 2.05 | 0.820 | +0.2099 | +0.2099 | +0.2080 |
| 2.55 | 1.020 | -0.0175 | -0.0175 | -0.0196 |
| 3.55 | 1.420 | -0.2151 | -0.2150 | -0.2171 |
max |V_sim - V_ref| = 1.0e-4; max |V_sim - sinc| = 2.1e-3 (finite slit width a/d = 6%, included exactly in the
reference); normalized pattern difference 7e-4. Fringes wash out at u = 1 and the contrast reverses beyond it.
Window finding: with a 32.8 mm window the discrepancy was 2.9e-3 (band-limit cutoff acting on the emitter's flat
spectrum); 65.5 mm gives 1.0e-4 and 131 mm the same, so 65.5 mm is used.

## Phase 4a: 3D geometry foundation (tests/test_vec3.c, test_xform.c, test_geom.c)
Geometry is independent of the transport method: nothing in these headers knows about wavelength,
phase or materials. Conventions: hit normal is the OUTWARD geometric normal and `front` says which
side the ray hit (the optics layer decides which side is glass); triangles are CCW-outward;
`wr_ray_make` normalizes d, intersect routines also accept non-unit d (t in units of |d|).

| Check | Result |
|---|---|
| Cross product / Lagrange identity (20000 random pairs) | 1.1e-15 |
| Orthonormal basis, right-handed | 6.7e-16 |
| Reflect, unit length + angle preserved | 1.0e-15 |
| Refract vs Snell (eta sin_i = sin_t, coplanar), 9287 transmitted / 10713 TIR cases | 9.1e-15; TIR only past critical angle |
| Rotation orthonormality / det = 1 / angle addition | 1.6e-15 / 1.1e-15 / 1.4e-15 |
| Affine inverse roundtrip, compose order, normal vs transformed tangents (non-uniform scale) | 5.3e-15 / <1e-13 / 1.5e-15 |
| Sphere: 19832 random hits on surface, on ray, normal | 1.7e-15 rel / 6.7e-15 / 1.9e-15 |
| Sphere r = 1 mm at 1e6 m (naive discriminant fails here) | t error < 1e-8 m (test tolerance) |
| Triangle: 9479 inside rays barycentric error; 8096 outside rays all miss | 4.2e-14; 0 false hits |
| AABB conservativeness: 1851 sphere/triangle hits, all inside the box slab interval | all pass |

Limitations: Moller-Trumbore is not watertight (rays through a shared edge can in rare rounding cases
miss both triangles); default tmin 1e-9 m is scale dependent; slab test has no ulp padding, so BVH
boxes should be expanded slightly; refract assumes real indices (complex indices come with absorption);
grazing triangle rays (|d.n| < 0.05) are excluded from the randomized accuracy tests.

## Phase 4b: meshes, OBJ, BVH (tests/test_mesh.c, test_bvh.c)
OBJ: `v` and `f` only (all token forms, negative indices, forward references, CRLF); other statements ignored;
polygons fan-triangulated (convex faces only); index 0 / out of range / non-finite coordinates are errors.
Triangles are CCW-outward; a closed mesh has positive signed volume.

| Check | Result |
|---|---|
| Closed + consistently oriented (every directed edge once, its reverse once), all triangles outward: box, UV sphere 64x32, smallest sphere | pass |
| Box volume/area, pentagon fan area, UV sphere volume / area vs exact (inscribed polyhedron) | exact to 1e-13; 0.99599 / 0.99799 |
| Transforms: scale volume x24, rotation volume x1.000000000000, mirror keeps outward winding, singular rejected and mesh untouched | pass |
| 7 malformed OBJ inputs rejected with a message, nothing leaked | pass |
| BVH vs brute force: 3000-triangle soup, 40000 rays incl. axis-aligned, slab-parallel, tmin/tmax-clipped | 0 mismatches; t bit-identical |
| 4 threads vs 1 thread traversal, 40000 rays | 0 differences |
| 500 identical triangles (no spatial split possible) | depth 7, max leaf 4 |
| Coplanar 40x40 sheet (flat boxes), 21681 rays incl. exact grid vertices | 0 mismatches; 20000/20000 interior rays hit; 1681/1681 vertex rays hit |
| Sphere mesh: hit radius in [0.9976, 1.0000] R, outward normals, front = 1 | pass |

BVH: 3968-triangle sphere, 1 thread: about 3 Mray/s vs 0.02 for brute force (~145x; informational, machine dependent).
SAH cost on the random soup: 38 triangle tests per ray vs 3000 brute force.
Bug found by the tests and fixed: the box generator wound all faces inward (volume -48).
Limitations: build is serial; no refit for moving geometry yet (rebuild); no instancing; one BVH per mesh
(analytic spheres/planes are not in the BVH yet; the scene layer will combine them).

## Phase 4c (camera part): geometric cameras (tests/test_camera.c, test_render.c)
Geometric ray generation only. The wave-optical sensor (aperture diffraction, PSF, spectral response) belongs to the
transport phase and is deliberately not faked here. Pixel origin top-left, py downward; right = forward x up_hint.

| Check | Result |
|---|---|
| Orientation, centre ray, vfov/hfov edge angles | exact to 1e-14 |
| Uniform pixel grid in tan-space | 1e-15 |
| 5000 random poses: orthonormal right-handed frame; ray -> point -> project round trip | 1.8e-15; 6.3e-13 px |
| Orthographic: parallel rays, window extent, pixel pitch | exact |
| Thin lens: all lens samples converge on the same focus-plane point; origins inside aperture | 1.8e-15 m |
| Thin lens: blur radius at 2F equals A|1 - z/F| (geometric circle of confusion) | 0.05000 vs 0.050 |
| Concentric disc map: inside unit disc, mean r^2 = 0.5, quadrants balanced | pass |
| Invalid cameras (parallel up, vfov 0 or pi, zero size, focus 0, negative aperture) rejected | pass |
| End to end: sphere silhouette radius 71.79 px; analytic 0.04% and mesh+BVH 0.18% area error; 4 vs 1 thread byte-identical | pass; image in build/geometry_check.ppm |

Limitations: square pixels only; thin lens is ray optics (no diffraction, no lens aberrations); projection is pinhole only.

## Phase 4c: materials, lights, scene, loader (tests/test_material.c, test_light.c, test_scene.c, test_sceneio.c)
Description layer only: nothing here performs light transport. Time convention exp(-i w t): phase grows along propagation,
kappa > 0 is absorption (intensity decays as exp(-alpha L), alpha = 4 pi kappa / lambda). Lines of a spectrum have different
frequencies, so their intensities ADD (they never interfere); within one line amplitudes are complex.

| Check | Result |
|---|---|
| Sellmeier BK7 vs Schott catalogue: n_d / n_F / n_C / Abbe number | 1.51680 / 1.52238 / 1.51432 / 64.17 (catalogue 1.5168 / 1.52238 / 1.51432 / 64.17) |
| Sellmeier fused silica (Malitson): n_d, Abbe number | 1.45846, 67.82 (published 1.4585, ~67.8) |
| Normal dispersion (n falls with wavelength, 400-700 nm), Cauchy hand evaluation | pass / exact |
| Absorption: exp(i k0 n_c L) intensity = exp(-alpha L), phase from the real index only | 1e-14 |
| Plane wave: constant amplitude, phase advance k d.dx, transverse invariance | amplitude 2e-16; phase at the representation floor (below) |
| Point source: r^2 |E|^2 = A^2, phase k r, Helmholtz residual (4th-order finite differences) | 5.6e-11 (position rounding); 9.5e-10 rad; 6.7e-6 |
| Gaussian beam: waist field, 1/e radius, Rayleigh amplitude A/sqrt2, Gouy phase, width law w(z), tilted frame | exact to 1e-12 |
| Gaussian beam power through transverse planes (Simpson), z = 0..3 z_R | 2.2e-11 |
| Gaussian Helmholtz residual (paraxial, nonzero by nature), w0 = 2 / 3 / 20 lambda | 2.9e-4 / 5.2e-5 / 6.8e-6 (finite-difference floor; falls with w0) |
| Spectrum: line intensities add to amplitude^2; weights are power fractions | 1e-14 |
| Polarization basis: orthonormal, right-handed with d, e1 = projection of pol_ref (3000 random directions + parallel fallback) | 6e-15 |
| Scene closest hit over sphere / plane / mesh, 30000 rays (12238 hit) vs per-object reference | 0 mismatches; insertion order irrelevant |
| 4-thread vs serial scene queries | 0 differences |
| Solids: open mesh with a dielectric material rejected; inside-out closed mesh rejected; open mirror sheet accepted | pass |
| Loader: 54 distinct error cases, each with the right message and line number; CRLF, comments, blank lines | 0 wrong |
| Loader: values, defaults, normalization (direction, plane normal, Jones), transform order (translate-then-scale != scale-then-translate), mirror transform keeps outward winding | pass |
| demo.scene end to end: all 5 objects visible, aimed rays hit lens / prism face / wall at analytic distances, 4 vs 1 thread byte-identical | pass; image in build/scene_check.ppm |
| Whole suite under AddressSanitizer + UBSan with leak detection | clean |

Findings from the tests (kept visible on purpose):
- **Phase floor of the representation.** A double at coordinate x has spacing ~eps |x|, so any phase k s carries an error ~ k eps |x|.
  At 500 nm and |x| ~ 0.4 m that is ~1e-9 rad (observed 1.0e-9; the same test with the source at the origin and probes within 1 mm gives 7e-12).
  Not a code error, but it matters: keep the scene origin near the region of interest, and do not expect better phase than this from
  absolute world coordinates. Transport code should compute path differences from local offsets where it can.
- **A real bug caught:** the Gaussian beam computed rho^2 as |p|^2 - z^2, losing ~11 digits for points close to the axis far from the
  waist (2e-11 width error). Fixed by computing the transverse part directly.
- The box generator wound its faces inward (found in 4b); the loader's validation now also refuses such meshes as solids.

Limitations: media are not nested (inside = the object's material, outside = ambient); overlapping solids are not detected; Cauchy and
Sellmeier have constant extinction (no absorption dispersion); the OBJ loader triangulates polygons as fans (convex faces only);
object and light velocities are stored and parsed but nothing uses them until Phase 8; coherence length and Jones vectors are stored and
validated but only interpreted in Phases 5 and 7; the scene is rebuilt, not updated, when geometry changes; one camera per scene.
