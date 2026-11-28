#include <math.h>
#include <stdio.h>
#include <string.h>
#include "wr_light.h"
#include "wr_optics.h"

static wr_light base(const char *name, wr_light_type type, wr_vec3 pos, double lambda0, double amp)
{
    wr_light l;
    memset(&l, 0, sizeof l);
    snprintf(l.name, sizeof l.name, "%s", name ? name : "");
    l.type = type;
    l.pos = pos;
    l.dir = wr_v3(0, 0, 1);
    l.amplitude = amp;
    l.spectrum.n = 1;
    l.spectrum.lambda[0] = lambda0;
    l.spectrum.weight[0] = 1.0;
    l.coherence_length = INFINITY;
    l.jones[0] = wr_c(1, 0);
    l.jones[1] = wr_c(0, 0);
    l.pol_ref = wr_v3(0, 1, 0);
    return l;
}

wr_light wr_light_point(const char *name, wr_vec3 pos, double lambda0, double amplitude)
{
    return base(name, WR_LIGHT_POINT, pos, lambda0, amplitude);
}

wr_light wr_light_plane_wave(const char *name, wr_vec3 pos, wr_vec3 dir, double lambda0, double amplitude)
{
    wr_light l = base(name, WR_LIGHT_PLANE_WAVE, pos, lambda0, amplitude);
    l.dir = wr_v3_normalize(dir);
    return l;
}

wr_light wr_light_gaussian(const char *name, wr_vec3 pos, wr_vec3 dir, double waist, double lambda0, double amplitude)
{
    wr_light l = base(name, WR_LIGHT_GAUSSIAN, pos, lambda0, amplitude);
    l.dir = wr_v3_normalize(dir);
    l.waist = waist;
    return l;
}

static int bad(char *err, size_t n, const char *name, const char *msg)
{
    if (err && n) snprintf(err, n, "light '%s': %s", name, msg);
    return 0;
}

int wr_light_validate(const wr_light *l, char *err, size_t errlen)
{
    if (l->name[0] == '\0') return bad(err, errlen, "?", "empty name");
    if (!(l->amplitude > 0.0) || !isfinite(l->amplitude)) return bad(err, errlen, l->name, "amplitude must be finite and > 0");
    if (!isfinite(l->phase)) return bad(err, errlen, l->name, "phase must be finite");
    if (l->spectrum.n < 1 || l->spectrum.n > WR_MAX_LINES) return bad(err, errlen, l->name, "spectrum needs 1..16 lines");
    for (int i = 0; i < l->spectrum.n; i++) {
        if (!(l->spectrum.lambda[i] > 0.0) || !isfinite(l->spectrum.lambda[i])) return bad(err, errlen, l->name, "wavelength must be finite and > 0");
        if (!(l->spectrum.weight[i] > 0.0) || !isfinite(l->spectrum.weight[i])) return bad(err, errlen, l->name, "spectral weight must be finite and > 0");
    }
    if (l->type != WR_LIGHT_POINT && fabs(wr_v3_len(l->dir) - 1.0) > 1e-9) return bad(err, errlen, l->name, "direction must be a unit vector");
    if (l->type == WR_LIGHT_GAUSSIAN && (!(l->waist > 0.0) || !isfinite(l->waist))) return bad(err, errlen, l->name, "waist must be finite and > 0");
    if (!(l->coherence_length >= 0.0)) return bad(err, errlen, l->name, "coherence length must be >= 0 (or inf)");
    double j2 = wr_cabs2(l->jones[0]) + wr_cabs2(l->jones[1]);
    if (!(fabs(j2 - 1.0) < 1e-9)) return bad(err, errlen, l->name, "Jones vector must have unit norm");
    if (!(wr_v3_len2(l->pol_ref) > 0.0)) return bad(err, errlen, l->name, "pol_ref must be non-zero");
    if (!isfinite(wr_v3_len2(l->velocity))) return bad(err, errlen, l->name, "velocity must be finite");
    return 1;
}

double wr_light_line_amplitude(const wr_light *l, int line)
{
    double s = 0.0;
    for (int i = 0; i < l->spectrum.n; i++) s += l->spectrum.weight[i];
    return l->amplitude * sqrt(l->spectrum.weight[line] / s);
}

wr_cplx wr_light_field(const wr_light *l, wr_vec3 p, int line)
{
    double lambda = l->spectrum.lambda[line];
    double k = wr_wavenumber(lambda);
    double A = wr_light_line_amplitude(l, line);
    wr_vec3 rel = wr_v3_sub(p, l->pos);
    switch (l->type) {
    case WR_LIGHT_POINT: {
        double r = wr_v3_len(rel);
        if (r == 0.0) return wr_c(NAN, NAN);
        return wr_cpolar(A / r, wr_wrap_phase(k * r + l->phase));
    }
    case WR_LIGHT_PLANE_WAVE:
        return wr_cpolar(A, wr_wrap_phase(k * wr_v3_dot(rel, l->dir) + l->phase));
    case WR_LIGHT_GAUSSIAN: {
        double z = wr_v3_dot(rel, l->dir);
        wr_vec3 perp = wr_v3_madd(rel, l->dir, -z);                       /* transverse part, no |rel|^2 - z^2 cancellation */
        double rho2 = wr_v3_len2(perp);
        double zr = WR_PI * l->waist * l->waist / lambda;                 /* Rayleigh range */
        double w2 = l->waist * l->waist * (1.0 + (z / zr) * (z / zr));    /* w(z)^2 */
        double inv_R = z / (z * z + zr * zr);                             /* 1/R(z), finite at the waist */
        double amp = A * l->waist / sqrt(w2) * exp(-rho2 / w2);
        double ph = k * z + 0.5 * k * rho2 * inv_R - atan(z / zr) + l->phase;   /* Gouy phase -atan(z/zr) */
        return wr_cpolar(amp, wr_wrap_phase(ph));
    }
    }
    return wr_c(NAN, NAN);
}

void wr_light_pol_basis(const wr_light *l, wr_vec3 d, wr_vec3 *e1, wr_vec3 *e2)
{
    wr_vec3 e = wr_v3_madd(l->pol_ref, d, -wr_v3_dot(l->pol_ref, d));
    double len = wr_v3_len(e);
    if (len > 1e-12 * wr_v3_len(l->pol_ref)) {
        e = wr_v3_scale(e, 1.0 / len);
    } else {
        wr_vec3 t, b;
        wr_v3_basis(d, &t, &b);
        e = t;
    }
    *e1 = e;
    *e2 = wr_v3_cross(d, e);
}
