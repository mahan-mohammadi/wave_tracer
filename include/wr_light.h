/* wr_light.h - light source descriptions and exact source fields (Phase 4c).
 *
 * A light is a description plus the analytic scalar field it radiates in
 * vacuum. Propagation through a scene is the transport layer's job.
 *
 * Spectrum: a set of spectral lines (vacuum wavelength, relative power weight).
 * Lines have different frequencies and are therefore mutually incoherent:
 * intensities of different lines ADD (roadmap rule 6), while the field within
 * one line is a complex amplitude. Line i carries power fraction w_i / sum(w),
 * so its amplitude is amplitude * sqrt(w_i / sum(w)). Continuous spectra come
 * with the spectral phase (Phase 6).
 *
 * Amplitude (field units, intensity = |E|^2):
 *   POINT      : field amplitude at 1 m distance; intensity amplitude^2 / r^2
 *   PLANE_WAVE : constant amplitude everywhere
 *   GAUSSIAN   : on-axis amplitude at the waist (waist centre is `pos`)
 * Convention exp(-i w t): phase grows along propagation: exp(i (k s + phase)).
 *
 * Fields stored now and interpreted by later phases (parsed and validated,
 * but not used by any field function here):
 *   coherence_length (m): INFINITY = fully coherent line, 0 = incoherent, in between = partial (Phase 5)
 *   jones / pol_ref     : Jones vector (unit) in the transverse basis built by wr_light_pol_basis (Phase 7)
 *   velocity (m/s)      : source motion, Doppler (Phase 8)
 */
#ifndef WR_LIGHT_H
#define WR_LIGHT_H

#include "wr_vec3.h"
#include "wr_complex.h"
#include "wr_material.h"

#define WR_MAX_LINES 16

typedef enum { WR_LIGHT_POINT, WR_LIGHT_PLANE_WAVE, WR_LIGHT_GAUSSIAN } wr_light_type;

typedef struct {
    int    n;
    double lambda[WR_MAX_LINES];    /* vacuum wavelengths, m */
    double weight[WR_MAX_LINES];    /* relative power weights, > 0 */
} wr_spectrum;

typedef struct {
    char          name[WR_NAME_LEN];
    wr_light_type type;
    wr_vec3       pos;              /* point: location; plane wave: phase reference point; Gaussian: waist centre */
    wr_vec3       dir;              /* unit propagation direction (plane wave, Gaussian) */
    double        waist;            /* Gaussian 1/e field radius w0 at the waist, m */
    double        amplitude;
    double        phase;            /* initial phase, rad */
    wr_spectrum   spectrum;
    double        coherence_length;
    wr_cplx       jones[2];
    wr_vec3       pol_ref;
    wr_vec3       velocity;
} wr_light;

/* Constructors give one spectral line, phase 0, fully coherent, Jones (1,0) along pol_ref (default +y), at rest. */
wr_light wr_light_point(const char *name, wr_vec3 pos, double lambda0, double amplitude);
wr_light wr_light_plane_wave(const char *name, wr_vec3 pos, wr_vec3 dir, double lambda0, double amplitude);
wr_light wr_light_gaussian(const char *name, wr_vec3 pos, wr_vec3 dir, double waist, double lambda0, double amplitude);

int    wr_light_validate(const wr_light *l, char *err, size_t errlen);
double wr_light_line_amplitude(const wr_light *l, int line);

/* Scalar field of spectral line `line` at point p, in vacuum. A point source returns NaN at its own position. */
wr_cplx wr_light_field(const wr_light *l, wr_vec3 p, int line);

/* Transverse polarization basis for propagation direction d (unit): e1 = pol_ref projected perpendicular to d
 * (any perpendicular if pol_ref is parallel to d), e2 = d x e1, so (e1, e2, d) is right-handed. */
void wr_light_pol_basis(const wr_light *l, wr_vec3 d, wr_vec3 *e1, wr_vec3 *e2);

#endif /* WR_LIGHT_H */
