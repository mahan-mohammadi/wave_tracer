/* wr_aabb.h - axis-aligned bounding box with slab ray test (Phase 4a). */
#ifndef WR_AABB_H
#define WR_AABB_H

#include <math.h>
#include "wr_vec3.h"
#include "wr_ray.h"

typedef struct { wr_vec3 lo, hi; } wr_aabb;

/* Empty box: contains nothing; the identity for union/expand. */
static inline wr_aabb wr_aabb_empty(void)
{
    wr_aabb b = {wr_v3(INFINITY, INFINITY, INFINITY), wr_v3(-INFINITY, -INFINITY, -INFINITY)};
    return b;
}
static inline wr_aabb wr_aabb_expand(wr_aabb b, wr_vec3 p)
{
    b.lo = wr_v3_min(b.lo, p);
    b.hi = wr_v3_max(b.hi, p);
    return b;
}
static inline wr_aabb wr_aabb_union(wr_aabb a, wr_aabb b)
{
    a.lo = wr_v3_min(a.lo, b.lo);
    a.hi = wr_v3_max(a.hi, b.hi);
    return a;
}
static inline wr_vec3 wr_aabb_centroid(const wr_aabb *b) { return wr_v3_scale(wr_v3_add(b->lo, b->hi), 0.5); }
static inline double wr_aabb_surface_area(const wr_aabb *b)
{
    wr_vec3 e = wr_v3_sub(b->hi, b->lo);
    return 2.0 * (e.x * e.y + e.y * e.z + e.z * e.x);
}
/* Axis (0,1,2) with the largest extent. */
static inline int wr_aabb_longest_axis(const wr_aabb *b)
{
    wr_vec3 e = wr_v3_sub(b->hi, b->lo);
    return (e.x >= e.y && e.x >= e.z) ? 0 : (e.y >= e.z ? 1 : 2);
}

/* Slab test. Clips [ray->tmin, ray->tmax] against the box; returns 1 and the
 * clipped interval if it is non-empty (touching counts). A ray parallel to a
 * slab is handled explicitly (no 0 * inf NaNs). No robust-rounding padding
 * of tfar: for BVH use, expand the boxes by a few ulps instead. */
static inline int wr_aabb_hit(const wr_aabb *b, const wr_ray *r, double *tnear, double *tfar)
{
    double t0 = r->tmin, t1 = r->tmax;
    for (int a = 0; a < 3; a++) {
        double o = wr_v3_get(r->o, a), d = wr_v3_get(r->d, a);
        double lo = wr_v3_get(b->lo, a), hi = wr_v3_get(b->hi, a);
        if (d == 0.0) {
            if (o < lo || o > hi) return 0;
            continue;
        }
        double inv = 1.0 / d;
        double ta = (lo - o) * inv, tb = (hi - o) * inv;
        if (ta > tb) { double s = ta; ta = tb; tb = s; }
        if (ta > t0) t0 = ta;
        if (tb < t1) t1 = tb;
        if (t0 > t1) return 0;
    }
    if (tnear) *tnear = t0;
    if (tfar) *tfar = t1;
    return 1;
}

#endif /* WR_AABB_H */
