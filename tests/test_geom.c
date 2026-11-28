#include "wr_test.h"
#include "wr_rng.h"
#include "wr_geom.h"

static wr_vec3 rv(wr_rng *r, double lo, double hi)
{
    return wr_v3(wr_rng_range(r, lo, hi), wr_rng_range(r, lo, hi), wr_rng_range(r, lo, hi));
}
static wr_ray ray_raw(wr_vec3 o, wr_vec3 d)   /* keeps d as given (not normalized) */
{
    wr_ray r = {o, d, WR_RAY_TMIN, INFINITY};
    return r;
}
static void near3(wr_vec3 a, wr_vec3 b, double tol, int line)
{
    wr_checks++;
    if (!(wr_v3_dist(a, b) <= tol)) {
        wr_failures++;
        fprintf(stderr, "FAIL %s:%d: (%g,%g,%g) vs (%g,%g,%g)\n", __FILE__, line, a.x, a.y, a.z, b.x, b.y, b.z);
    }
}
#define NEAR3(a, b, tol) near3((a), (b), (tol), __LINE__)

static void test_sphere(void)
{
    wr_sphere s = wr_sphere_make(wr_v3(0, 0, 0), 1.0);
    wr_hit h;
    wr_ray r = wr_ray_make(wr_v3(0, 0, -5), wr_v3(0, 0, 1));
    CHECK(wr_sphere_intersect(&s, &r, &h));
    CHECK_NEAR(h.t, 4, 1e-15); NEAR3(h.p, wr_v3(0, 0, -1), 1e-15); NEAR3(h.n, wr_v3(0, 0, -1), 1e-15);
    CHECK(h.front == 1);

    r = wr_ray_make(wr_v3(0, 0, 0), wr_v3(0, 0, 1));            /* from inside: far wall, back side */
    CHECK(wr_sphere_intersect(&s, &r, &h));
    CHECK_NEAR(h.t, 1, 1e-15); NEAR3(h.n, wr_v3(0, 0, 1), 1e-15); CHECK(h.front == 0);
    NEAR3(wr_hit_facing_normal(&h), wr_v3(0, 0, -1), 1e-15);

    r = wr_ray_make(wr_v3(1, 0, -5), wr_v3(0, 0, 1));           /* tangent */
    CHECK(wr_sphere_intersect(&s, &r, &h)); CHECK_NEAR(h.t, 5, 1e-15);
    r = wr_ray_make(wr_v3(1.0001, 0, -5), wr_v3(0, 0, 1));      /* just misses */
    CHECK(!wr_sphere_intersect(&s, &r, &h));
    r = wr_ray_make(wr_v3(0, 0, -5), wr_v3(0, 0, -1));          /* sphere behind */
    CHECK(!wr_sphere_intersect(&s, &r, &h));

    r = wr_ray_make(wr_v3(0, 0, -5), wr_v3(0, 0, 1));
    r.tmax = 3.5; CHECK(!wr_sphere_intersect(&s, &r, &h));      /* clipped before the sphere */
    r.tmax = INFINITY; r.tmin = 4.5;                            /* near root rejected: back wall */
    CHECK(wr_sphere_intersect(&s, &r, &h));
    CHECK_NEAR(h.t, 6, 1e-15); NEAR3(h.p, wr_v3(0, 0, 1), 1e-15); CHECK(h.front == 0);

    r = ray_raw(wr_v3(0, 0, -5), wr_v3(0, 0, 2));               /* non-unit direction: t in units of |d| */
    CHECK(wr_sphere_intersect(&s, &r, &h));
    CHECK_NEAR(h.t, 2, 1e-15); NEAR3(h.p, wr_v3(0, 0, -1), 1e-15);

    /* ray leaving a surface point must not re-hit that point (tmin) */
    r = wr_ray_make(wr_v3(0, 0, 1), wr_v3(0, 0, 1));
    CHECK(!wr_sphere_intersect(&s, &r, &h));

    /* tiny sphere far away: a naive discriminant loses all digits here */
    wr_sphere tiny = wr_sphere_make(wr_v3(0, 0, 1e6), 1e-3);
    r = wr_ray_make(wr_v3(0, 0, 0), wr_v3(0, 0, 1));
    CHECK(wr_sphere_intersect(&tiny, &r, &h));
    CHECK_NEAR(h.t, 1e6 - 1e-3, 1e-8);
    r = wr_ray_make(wr_v3(0.9e-3, 0, 0), wr_v3(0, 0, 1));
    CHECK(wr_sphere_intersect(&tiny, &r, &h));
    CHECK_NEAR(h.t, 1e6 - sqrt(1e-6 - 0.81e-6), 1e-8);          /* chord offset 0.9 mm */
    r = wr_ray_make(wr_v3(1.001e-3, 0, 0), wr_v3(0, 0, 1));
    CHECK(!wr_sphere_intersect(&tiny, &r, &h));

    /* randomized: hit lies on the sphere, on the ray, with consistent normal and front flag;
       and the ray aimed at a point inside always hits */
    wr_rng rng = {12345};
    double w_on = 0, w_ray = 0, w_nrm = 0;
    int n_hit = 0;
    for (int i = 0; i < 20000; i++) {
        wr_sphere sp = wr_sphere_make(rv(&rng, -2, 2), wr_rng_range(&rng, 0.1, 2));
        wr_vec3 target = wr_v3_madd(sp.center, wr_v3_normalize(rv(&rng, -1, 1)), 0.99 * sp.radius * wr_rng_uniform(&rng));
        wr_vec3 o = rv(&rng, -10, 10);
        if (wr_v3_dist(o, sp.center) < 2 * sp.radius) continue;
        wr_ray rr = wr_ray_make(o, wr_v3_sub(target, o));
        CHECK(wr_sphere_intersect(&sp, &rr, &h));
        n_hit++;
        w_on = fmax(w_on, fabs(wr_v3_dist(h.p, sp.center) - sp.radius) / sp.radius);
        w_ray = fmax(w_ray, wr_v3_dist(h.p, wr_ray_at(&rr, h.t)));
        w_nrm = fmax(w_nrm, wr_v3_dist(h.n, wr_v3_scale(wr_v3_sub(h.p, sp.center), 1.0 / sp.radius)));
        CHECK(h.front == 1);                                    /* outside origin: always the near wall */
        CHECK(wr_v3_dot(rr.d, wr_v3_sub(target, o)) > 0);
    }
    CHECK_NEAR(w_on, 0, 1e-14); CHECK_NEAR(w_ray, 0, 1e-13); CHECK_NEAR(w_nrm, 0, 1e-14);
    printf("sphere: %d random hits, on-surface %.1e (rel), on-ray %.1e, normal %.1e\n", n_hit, w_on, w_ray, w_nrm);
}

static void test_plane(void)
{
    wr_plane p = wr_plane_make(wr_v3(0, 0, 2), wr_v3(0, 0, 5));  /* normal is normalized by make */
    wr_hit h;
    wr_ray r = wr_ray_make(wr_v3(0, 0, 5), wr_v3(0, 0, -1));
    CHECK(wr_plane_intersect(&p, &r, &h));
    CHECK_NEAR(h.t, 3, 1e-15); NEAR3(h.p, wr_v3(0, 0, 2), 1e-15); NEAR3(h.n, wr_v3(0, 0, 1), 1e-15); CHECK(h.front == 1);
    r = wr_ray_make(wr_v3(0, 0, -1), wr_v3(0, 0, 1));
    CHECK(wr_plane_intersect(&p, &r, &h)); CHECK_NEAR(h.t, 3, 1e-15); CHECK(h.front == 0);
    r = ray_raw(wr_v3(0, 0, 5), wr_v3(1, 0, -1));
    CHECK(wr_plane_intersect(&p, &r, &h)); CHECK_NEAR(h.t, 3, 1e-15); NEAR3(h.p, wr_v3(3, 0, 2), 1e-15);
    r = wr_ray_make(wr_v3(0, 0, 5), wr_v3(1, 0, 0));
    CHECK(!wr_plane_intersect(&p, &r, &h));                      /* parallel */
    r = wr_ray_make(wr_v3(0, 0, 5), wr_v3(0, 0, 1));
    CHECK(!wr_plane_intersect(&p, &r, &h));                      /* plane behind */
    r = wr_ray_make(wr_v3(0, 0, 5), wr_v3(0, 0, -1)); r.tmax = 2.9;
    CHECK(!wr_plane_intersect(&p, &r, &h));
}

static void test_triangle(void)
{
    wr_triangle tri = wr_triangle_make(wr_v3(0, 0, 0), wr_v3(1, 0, 0), wr_v3(0, 1, 0));
    wr_hit h;
    NEAR3(wr_triangle_normal(&tri), wr_v3(0, 0, 1), 1e-16);
    CHECK_NEAR(wr_triangle_area(&tri), 0.5, 1e-16);
    wr_ray r = wr_ray_make(wr_v3(0.25, 0.25, 1), wr_v3(0, 0, -1));
    CHECK(wr_triangle_intersect(&tri, &r, &h));
    CHECK_NEAR(h.t, 1, 1e-15); CHECK_NEAR(h.u, 0.25, 1e-15); CHECK_NEAR(h.v, 0.25, 1e-15);
    CHECK(h.front == 1); NEAR3(h.n, wr_v3(0, 0, 1), 1e-16);
    r = wr_ray_make(wr_v3(0.25, 0.25, -1), wr_v3(0, 0, 1));       /* from behind: hit, front = 0, n unchanged */
    CHECK(wr_triangle_intersect(&tri, &r, &h));
    CHECK(h.front == 0); NEAR3(h.n, wr_v3(0, 0, 1), 1e-16);
    r = wr_ray_make(wr_v3(0.8, 0.8, 1), wr_v3(0, 0, -1));
    CHECK(!wr_triangle_intersect(&tri, &r, &h));
    r = wr_ray_make(wr_v3(0.5, 0.5, 1), wr_v3(0, 0, -1));         /* on the hypotenuse: inclusive */
    CHECK(wr_triangle_intersect(&tri, &r, &h));
    r = wr_ray_make(wr_v3(-1e-9, 0.25, 1), wr_v3(0, 0, -1));      /* just outside an edge */
    CHECK(!wr_triangle_intersect(&tri, &r, &h));
    r = wr_ray_make(wr_v3(-1, 0.25, 1), wr_v3(1, 0, 0));          /* parallel to the plane */
    CHECK(!wr_triangle_intersect(&tri, &r, &h));
    wr_triangle deg = wr_triangle_make(wr_v3(0, 0, 0), wr_v3(1, 0, 0), wr_v3(2, 0, 0));
    r = wr_ray_make(wr_v3(0.5, 0, 1), wr_v3(0, 0, -1));
    CHECK(!wr_triangle_intersect(&deg, &r, &h));                  /* degenerate: miss, no NaN */
    wr_triangle cw = wr_triangle_make(wr_v3(0, 0, 0), wr_v3(0, 1, 0), wr_v3(1, 0, 0));   /* reversed winding */
    NEAR3(wr_triangle_normal(&cw), wr_v3(0, 0, -1), 1e-16);

    /* randomized: ray aimed at a barycentric point must report those barycentrics;
       ray aimed outside (margin) must miss */
    wr_rng rng = {777};
    double w_bary = 0, w_pos = 0;
    int n_in = 0, n_out = 0;
    for (int i = 0; i < 20000; i++) {
        wr_triangle t = wr_triangle_make(rv(&rng, -3, 3), rv(&rng, -3, 3), rv(&rng, -3, 3));
        if (wr_triangle_area(&t) < 0.2) continue;
        wr_vec3 nrm = wr_triangle_normal(&t);
        wr_vec3 e1 = wr_v3_sub(t.v1, t.v0), e2 = wr_v3_sub(t.v2, t.v0);
        wr_vec3 o = rv(&rng, -8, 8);
        double u, v;
        int inside = (i & 1);
        if (inside) {                                             /* strictly inside, margin 1% */
            u = wr_rng_range(&rng, 0.01, 0.98); v = wr_rng_range(&rng, 0.01, 0.99 - u);
        } else {                                                  /* clearly outside */
            u = wr_rng_range(&rng, -0.5, 1.5); v = wr_rng_range(&rng, -0.5, 1.5);
            if (u >= -0.02 && v >= -0.02 && u + v <= 1.02) continue;
        }
        wr_vec3 target = wr_v3_add(t.v0, wr_v3_add(wr_v3_scale(e1, u), wr_v3_scale(e2, v)));
        wr_vec3 d = wr_v3_normalize(wr_v3_sub(target, o));
        if (fabs(wr_v3_dot(d, nrm)) < 0.05) continue;             /* skip grazing rays */
        if (fabs(wr_v3_dot(wr_v3_sub(o, t.v0), nrm)) < 0.05) continue;
        wr_ray rr = wr_ray_make(o, d);
        int hit = wr_triangle_intersect(&t, &rr, &h);
        if (inside) {
            CHECK(hit);
            if (hit) {
                n_in++;
                w_bary = fmax(w_bary, fmax(fabs(h.u - u), fabs(h.v - v)));
                w_pos = fmax(w_pos, wr_v3_dist(h.p, target));
                w_pos = fmax(w_pos, wr_v3_dist(wr_ray_at(&rr, h.t), target));
                CHECK(h.front == (wr_v3_dot(d, nrm) < 0));
            }
        } else {
            CHECK(!hit);
            n_out++;
        }
    }
    CHECK(n_in > 2000 && n_out > 2000);
    CHECK_NEAR(w_bary, 0, 1e-9); CHECK_NEAR(w_pos, 0, 1e-12);
    printf("triangle: %d inside / %d outside rays, barycentric error %.1e, position error %.1e\n",
           n_in, n_out, w_bary, w_pos);
}

static void test_aabb(void)
{
    wr_aabb b = {wr_v3(-1, -1, -1), wr_v3(1, 1, 1)};
    double tn, tf;
    wr_ray r = wr_ray_make(wr_v3(-5, 0, 0), wr_v3(1, 0, 0));
    CHECK(wr_aabb_hit(&b, &r, &tn, &tf)); CHECK_NEAR(tn, 4, 1e-15); CHECK_NEAR(tf, 6, 1e-15);
    r = wr_ray_make(wr_v3(0, 0, 0), wr_v3(0, 1, 0));               /* origin inside: tnear = tmin */
    CHECK(wr_aabb_hit(&b, &r, &tn, &tf)); CHECK_NEAR(tn, WR_RAY_TMIN, 0); CHECK_NEAR(tf, 1, 1e-15);
    r = wr_ray_make(wr_v3(-5, 2, 0), wr_v3(1, 0, 0));              /* parallel to slabs, outside */
    CHECK(!wr_aabb_hit(&b, &r, &tn, &tf));
    r = wr_ray_make(wr_v3(-5, 1, 0), wr_v3(1, 0, 0));              /* grazing a face: touching counts */
    CHECK(wr_aabb_hit(&b, &r, &tn, &tf));
    r = wr_ray_make(wr_v3(5, 0, 0), wr_v3(1, 0, 0));               /* box behind */
    CHECK(!wr_aabb_hit(&b, &r, &tn, &tf));
    r = wr_ray_make(wr_v3(-5, 0, 0), wr_v3(1, 0, 0)); r.tmax = 3;
    CHECK(!wr_aabb_hit(&b, &r, &tn, &tf));
    r = wr_ray_make(wr_v3(-3, -3, 0), wr_v3(1, 1, 0));             /* diagonal through the corner region */
    CHECK(wr_aabb_hit(&b, &r, &tn, &tf)); CHECK_NEAR(tn, 2 * sqrt(2.0), 1e-14);

    wr_aabb e = wr_aabb_empty();
    wr_aabb u = wr_aabb_union(e, b);
    CHECK(u.lo.x == -1 && u.hi.z == 1);
    r = wr_ray_make(wr_v3(-5, 0, 0), wr_v3(1, 0, 0));
    CHECK(!wr_aabb_hit(&e, &r, &tn, &tf));                          /* empty box never hit */
    CHECK_NEAR(wr_aabb_surface_area(&b), 24, 1e-15);
    CHECK(wr_aabb_longest_axis(&b) == 0);
    wr_aabb slab = {wr_v3(0, 0, 0), wr_v3(4, 1, 2)};
    CHECK(wr_aabb_longest_axis(&slab) == 0);
    wr_aabb tall = {wr_v3(0, 0, 0), wr_v3(1, 1, 5)};
    CHECK(wr_aabb_longest_axis(&tall) == 2);

    /* bounds must be conservative: every primitive hit also hits its box, inside [tnear, tfar] */
    wr_rng rng = {4242};
    int n = 0;
    for (int i = 0; i < 100000; i++) {
        wr_hit h;
        wr_ray rr = wr_ray_make(rv(&rng, -8, 8), rv(&rng, -1, 1));
        wr_sphere sp = wr_sphere_make(rv(&rng, -3, 3), wr_rng_range(&rng, 0.2, 2));
        wr_aabb sb = wr_sphere_bounds(&sp);
        if (wr_sphere_intersect(&sp, &rr, &h)) {
            n++;
            CHECK(wr_aabb_hit(&sb, &rr, &tn, &tf));
            CHECK(h.t >= tn - 1e-12 && h.t <= tf + 1e-12);
        }
        wr_triangle t = wr_triangle_make(rv(&rng, -3, 3), rv(&rng, -3, 3), rv(&rng, -3, 3));
        wr_aabb tb = wr_triangle_bounds(&t);
        if (wr_triangle_intersect(&t, &rr, &h)) {
            n++;
            CHECK(wr_aabb_hit(&tb, &rr, &tn, &tf));
            CHECK(h.t >= tn - 1e-12 && h.t <= tf + 1e-12);
        }
    }
    CHECK(n > 1000);          /* ~2% of random rays hit; 100000 samples give ~1900 hits */
    printf("aabb: %d primitive hits all inside their boxes' slab intervals\n", n);
}

int main(void)
{
    test_sphere();
    test_plane();
    test_triangle();
    test_aabb();
    TEST_REPORT("geom");
}
