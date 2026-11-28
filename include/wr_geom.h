/* wr_geom.h - primitives and ray intersection (Phase 4a).
 *
 * Geometry only. Nothing here knows about wavelength, phase, materials or
 * the transport method (roadmap rule: keep geometry independent of wave
 * transport). Each intersect function returns 1 and fills *hit with the
 * NEAREST intersection whose t lies in [ray->tmin, ray->tmax], else 0 and
 * *hit is untouched. hit->prim is not touched either (containers set it).
 *
 * Numerical notes (all in double):
 *  - Sphere uses the cancellation-free discriminant a*(r^2 - |v|^2) with
 *    v = the part of (o - c) perpendicular to d, and the Press stable pair
 *    of roots, so tiny spheres far from the origin stay accurate.
 *  - Triangle is Moller-Trumbore, double sided, edges inclusive. It is NOT
 *    watertight: a ray that passes exactly through a shared edge can in rare
 *    rounding cases miss both neighbours. Acceptable for a wave renderer that
 *    integrates over many rays; a watertight test can replace it later.
 */
#ifndef WR_GEOM_H
#define WR_GEOM_H

#include "wr_vec3.h"
#include "wr_ray.h"
#include "wr_aabb.h"

/* Sphere; the inside is the ball. */
typedef struct { wr_vec3 center; double radius; } wr_sphere;

/* Infinite plane through `point` with unit normal `normal`. The half-space
 * the normal points away from is the "inside" (front = 1 when hit from the
 * normal side). */
typedef struct { wr_vec3 point, normal; } wr_plane;

/* Triangle; counter-clockwise (v0, v1, v2) seen from outside gives the
 * outward normal normalize(cross(v1 - v0, v2 - v0)). */
typedef struct { wr_vec3 v0, v1, v2; } wr_triangle;

wr_sphere   wr_sphere_make(wr_vec3 center, double radius);
wr_plane    wr_plane_make(wr_vec3 point, wr_vec3 normal);   /* normalizes normal */
wr_triangle wr_triangle_make(wr_vec3 v0, wr_vec3 v1, wr_vec3 v2);

int wr_sphere_intersect(const wr_sphere *s, const wr_ray *ray, wr_hit *hit);
int wr_plane_intersect(const wr_plane *p, const wr_ray *ray, wr_hit *hit);
int wr_triangle_intersect(const wr_triangle *tri, const wr_ray *ray, wr_hit *hit);

wr_aabb wr_sphere_bounds(const wr_sphere *s);
wr_aabb wr_triangle_bounds(const wr_triangle *tri);
/* Unit geometric normal and area of a triangle (zero normal if degenerate). */
wr_vec3 wr_triangle_normal(const wr_triangle *tri);
double  wr_triangle_area(const wr_triangle *tri);

#endif /* WR_GEOM_H */
