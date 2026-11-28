/* wr_ray.h - ray and hit record (Phase 4a).
 *
 * Ray: p(t) = o + t d, valid for t in [tmin, tmax]. d need NOT be unit: the
 * intersection routines work for any non-zero d and report t in units of
 * |d| (this makes object-space instancing with scaled transforms trivial,
 * because t is the same in both spaces). wr_ray_make normalizes d.
 *
 * tmin default 1e-9 (1 nm, project units are metres): it only exists to stop
 * a ray leaving a surface from re-hitting that same surface through rounding.
 * It is scale dependent; a scene whose features approach 1 nm needs another
 * value. Spawning with an explicit offset is preferable in the transport code.
 *
 * Hit: n is the GEOMETRIC unit normal pointing to the outside of the
 * surface (for triangles: the side from which the vertices look
 * counter-clockwise). front = 1 when the ray arrives from the outside
 * (dot(d, n) < 0). The optics layer uses front to decide which side of an
 * interface is vacuum and which is glass, so the geometry must never flip n
 * to face the ray; use wr_hit_facing_normal for that.
 */
#ifndef WR_RAY_H
#define WR_RAY_H

#include <math.h>
#include "wr_vec3.h"

#define WR_RAY_TMIN 1e-9

typedef struct {
    wr_vec3 o, d;
    double tmin, tmax;
} wr_ray;

typedef struct {
    double  t;
    wr_vec3 p;        /* hit position */
    wr_vec3 n;        /* outward geometric normal, unit */
    double  u, v;     /* triangle: barycentric weights of v1 and v2; else 0 */
    int     front;    /* 1: ray hit the outside, 0: the inside */
    int     prim;     /* triangle index within a mesh, filled by containers; -1 = unset */
    int     obj;      /* scene object index, filled by wr_scene_intersect; -1 = unset */
} wr_hit;

static inline wr_ray wr_ray_make(wr_vec3 o, wr_vec3 d)
{
    wr_ray r;
    r.o = o;
    r.d = wr_v3_normalize(d);
    r.tmin = WR_RAY_TMIN;
    r.tmax = INFINITY;
    return r;
}

static inline wr_vec3 wr_ray_at(const wr_ray *r, double t) { return wr_v3_madd(r->o, r->d, t); }

/* Normal flipped to face against the ray: dot(d, result) <= 0. */
static inline wr_vec3 wr_hit_facing_normal(const wr_hit *h) { return h->front ? h->n : wr_v3_neg(h->n); }

#endif /* WR_RAY_H */
