/* wr_vec3.h - 3D vectors (Phase 4a). Pure geometry: no optics in here.
 *
 * Plain struct of three doubles, value semantics, all inline.
 * Directions passed to reflect/refract must be unit length (caller's job;
 * the tests check the unit-length results).
 */
#ifndef WR_VEC3_H
#define WR_VEC3_H

#include <math.h>

typedef struct { double x, y, z; } wr_vec3;

static inline wr_vec3 wr_v3(double x, double y, double z) { wr_vec3 v = {x, y, z}; return v; }
static inline wr_vec3 wr_v3_zero(void) { return wr_v3(0.0, 0.0, 0.0); }

static inline wr_vec3 wr_v3_add(wr_vec3 a, wr_vec3 b) { return wr_v3(a.x + b.x, a.y + b.y, a.z + b.z); }
static inline wr_vec3 wr_v3_sub(wr_vec3 a, wr_vec3 b) { return wr_v3(a.x - b.x, a.y - b.y, a.z - b.z); }
static inline wr_vec3 wr_v3_neg(wr_vec3 a) { return wr_v3(-a.x, -a.y, -a.z); }
static inline wr_vec3 wr_v3_scale(wr_vec3 a, double s) { return wr_v3(a.x * s, a.y * s, a.z * s); }
static inline wr_vec3 wr_v3_mul(wr_vec3 a, wr_vec3 b) { return wr_v3(a.x * b.x, a.y * b.y, a.z * b.z); }
/* a + b * s */
static inline wr_vec3 wr_v3_madd(wr_vec3 a, wr_vec3 b, double s)
{
    return wr_v3(a.x + b.x * s, a.y + b.y * s, a.z + b.z * s);
}
static inline wr_vec3 wr_v3_lerp(wr_vec3 a, wr_vec3 b, double t)
{
    return wr_v3(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, a.z + (b.z - a.z) * t);
}

static inline double wr_v3_dot(wr_vec3 a, wr_vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline wr_vec3 wr_v3_cross(wr_vec3 a, wr_vec3 b)
{
    return wr_v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
static inline double wr_v3_len2(wr_vec3 a) { return wr_v3_dot(a, a); }
static inline double wr_v3_len(wr_vec3 a) { return sqrt(wr_v3_dot(a, a)); }
static inline double wr_v3_dist(wr_vec3 a, wr_vec3 b) { return wr_v3_len(wr_v3_sub(a, b)); }

/* Unit vector. The zero vector maps to the zero vector (no NaNs). */
static inline wr_vec3 wr_v3_normalize(wr_vec3 a)
{
    double l = wr_v3_len(a);
    return l > 0.0 ? wr_v3_scale(a, 1.0 / l) : wr_v3_zero();
}

static inline wr_vec3 wr_v3_min(wr_vec3 a, wr_vec3 b)
{
    return wr_v3(fmin(a.x, b.x), fmin(a.y, b.y), fmin(a.z, b.z));
}
static inline wr_vec3 wr_v3_max(wr_vec3 a, wr_vec3 b)
{
    return wr_v3(fmax(a.x, b.x), fmax(a.y, b.y), fmax(a.z, b.z));
}
/* Component by axis index 0,1,2. */
static inline double wr_v3_get(wr_vec3 a, int axis) { return axis == 0 ? a.x : (axis == 1 ? a.y : a.z); }

/* Mirror direction: d - 2 (d.n) n.  d and n unit. */
static inline wr_vec3 wr_v3_reflect(wr_vec3 d, wr_vec3 n)
{
    return wr_v3_madd(d, n, -2.0 * wr_v3_dot(d, n));
}

/* Refracted direction by Snell's law (vector form).
 *   d   : unit incident direction
 *   n   : unit normal on the side the ray comes from, i.e. dot(d, n) <= 0
 *   eta : n_incident / n_transmitted (real indices)
 * Returns 1 and writes *out (unit) on transmission, 0 on total internal
 * reflection (*out untouched). This is only the direction; amplitudes and
 * phase belong to the optics layer, not to geometry. */
static inline int wr_v3_refract(wr_vec3 d, wr_vec3 n, double eta, wr_vec3 *out)
{
    double cosi = -wr_v3_dot(d, n);
    double k = 1.0 - eta * eta * (1.0 - cosi * cosi);
    if (k < 0.0) return 0;
    *out = wr_v3_madd(wr_v3_scale(d, eta), n, eta * cosi - sqrt(k));
    return 1;
}

/* Orthonormal basis (t, b, n) around unit n, branch-free and singularity-free
 * (Duff et al., "Building an Orthonormal Basis, Revisited", JCGT 2017). */
static inline void wr_v3_basis(wr_vec3 n, wr_vec3 *t, wr_vec3 *b)
{
    double s = copysign(1.0, n.z);
    double a = -1.0 / (s + n.z);
    double c = n.x * n.y * a;
    *t = wr_v3(1.0 + s * n.x * n.x * a, s * c, -s * n.x);
    *b = wr_v3(c, s + n.y * n.y * a, -n.y);
}

#endif /* WR_VEC3_H */
