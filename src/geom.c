#include <math.h>
#include "wr_geom.h"

wr_sphere wr_sphere_make(wr_vec3 center, double radius)
{
    wr_sphere s = {center, radius};
    return s;
}

wr_plane wr_plane_make(wr_vec3 point, wr_vec3 normal)
{
    wr_plane p = {point, wr_v3_normalize(normal)};
    return p;
}

wr_triangle wr_triangle_make(wr_vec3 v0, wr_vec3 v1, wr_vec3 v2)
{
    wr_triangle t = {v0, v1, v2};
    return t;
}

int wr_sphere_intersect(const wr_sphere *s, const wr_ray *ray, wr_hit *hit)
{
    wr_vec3 oc = wr_v3_sub(ray->o, s->center);
    double a = wr_v3_dot(ray->d, ray->d);
    if (!(a > 0.0)) return 0;
    double hb = wr_v3_dot(oc, ray->d);                    /* half of the linear coefficient */
    wr_vec3 v = wr_v3_madd(oc, ray->d, -hb / a);          /* oc perpendicular to d */
    double disc = a * (s->radius * s->radius - wr_v3_len2(v));
    if (disc < 0.0) return 0;
    double sd = sqrt(disc);
    double q = -(hb + copysign(sd, hb));
    double t0, t1;
    if (q == 0.0) {
        t0 = t1 = 0.0;                                    /* hb = 0 and tangent at the origin */
    } else {
        double cc = wr_v3_len2(oc) - s->radius * s->radius;
        t0 = q / a;
        t1 = cc / q;
        if (t0 > t1) { double tmp = t0; t0 = t1; t1 = tmp; }
    }
    double t;
    if (t0 >= ray->tmin && t0 <= ray->tmax) t = t0;
    else if (t1 >= ray->tmin && t1 <= ray->tmax) t = t1;
    else return 0;

    wr_vec3 rel = wr_v3_sub(wr_ray_at(ray, t), s->center);
    wr_vec3 n = wr_v3_normalize(rel);
    hit->t = t;
    hit->p = wr_v3_madd(s->center, n, s->radius);         /* re-projected onto the surface */
    hit->n = n;
    hit->u = hit->v = 0.0;
    hit->front = wr_v3_dot(ray->d, n) < 0.0;
    return 1;
}

int wr_plane_intersect(const wr_plane *pl, const wr_ray *ray, wr_hit *hit)
{
    double denom = wr_v3_dot(ray->d, pl->normal);
    if (fabs(denom) <= 1e-15 * wr_v3_len(ray->d)) return 0;   /* parallel */
    double t = wr_v3_dot(wr_v3_sub(pl->point, ray->o), pl->normal) / denom;
    if (!(t >= ray->tmin && t <= ray->tmax)) return 0;
    hit->t = t;
    hit->p = wr_ray_at(ray, t);
    hit->n = pl->normal;
    hit->u = hit->v = 0.0;
    hit->front = denom < 0.0;
    return 1;
}

int wr_triangle_intersect(const wr_triangle *tri, const wr_ray *ray, wr_hit *hit)
{
    wr_vec3 e1 = wr_v3_sub(tri->v1, tri->v0);
    wr_vec3 e2 = wr_v3_sub(tri->v2, tri->v0);
    wr_vec3 pv = wr_v3_cross(ray->d, e2);
    double det = wr_v3_dot(e1, pv);
    /* |det| <= |e1||e2||d| always; skip near-parallel or degenerate cases. */
    if (det * det <= 1e-30 * wr_v3_len2(e1) * wr_v3_len2(e2) * wr_v3_len2(ray->d)) return 0;
    double inv = 1.0 / det;
    wr_vec3 tv = wr_v3_sub(ray->o, tri->v0);
    double u = wr_v3_dot(tv, pv) * inv;
    if (u < 0.0 || u > 1.0) return 0;
    wr_vec3 qv = wr_v3_cross(tv, e1);
    double v = wr_v3_dot(ray->d, qv) * inv;
    if (v < 0.0 || u + v > 1.0) return 0;
    double t = wr_v3_dot(e2, qv) * inv;
    if (!(t >= ray->tmin && t <= ray->tmax)) return 0;

    wr_vec3 n = wr_v3_normalize(wr_v3_cross(e1, e2));
    hit->t = t;
    hit->p = wr_v3_add(tri->v0, wr_v3_add(wr_v3_scale(e1, u), wr_v3_scale(e2, v)));
    hit->n = n;
    hit->u = u;
    hit->v = v;
    hit->front = wr_v3_dot(ray->d, n) < 0.0;
    return 1;
}

wr_aabb wr_sphere_bounds(const wr_sphere *s)
{
    wr_vec3 r = wr_v3(s->radius, s->radius, s->radius);
    wr_aabb b = {wr_v3_sub(s->center, r), wr_v3_add(s->center, r)};
    return b;
}

wr_aabb wr_triangle_bounds(const wr_triangle *tri)
{
    wr_aabb b = wr_aabb_empty();
    b = wr_aabb_expand(b, tri->v0);
    b = wr_aabb_expand(b, tri->v1);
    b = wr_aabb_expand(b, tri->v2);
    return b;
}

wr_vec3 wr_triangle_normal(const wr_triangle *tri)
{
    return wr_v3_normalize(wr_v3_cross(wr_v3_sub(tri->v1, tri->v0), wr_v3_sub(tri->v2, tri->v0)));
}

double wr_triangle_area(const wr_triangle *tri)
{
    return 0.5 * wr_v3_len(wr_v3_cross(wr_v3_sub(tri->v1, tri->v0), wr_v3_sub(tri->v2, tri->v0)));
}
