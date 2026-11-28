#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include "wr_test.h"
#include "wr_rng.h"
#include "wr_scene.h"
#include "wr_optics.h"

static wr_vec3 rv(wr_rng *r, double lo, double hi)
{
    return wr_v3(wr_rng_range(r, lo, hi), wr_rng_range(r, lo, hi), wr_rng_range(r, lo, hi));
}
static int has(const char *hay, const char *needle) { return strstr(hay, needle) != NULL; }

/* sphere (glass) at the origin, absorbing plane behind it, box mesh to the side */
static int build(wr_scene *s, int sphere_first)
{
    char e[160];
    if (!wr_scene_init(s)) return 0;
    wr_material glass = wr_material_constant("glass", 1.5, 0.0), wall = wr_material_absorber("wall");
    int g = wr_scene_add_material(s, &glass, e, sizeof e), w = wr_scene_add_material(s, &wall, e, sizeof e);
    wr_mesh box;
    wr_mesh_make_box(wr_v3(4, -1, -1), wr_v3(6, 1, 1), &box);
    if (sphere_first) {
        wr_scene_add_sphere(s, "ball", g, wr_v3_zero(), 1.0, e, sizeof e);
        wr_scene_add_plane(s, "wall", w, wr_v3(0, 0, -3), wr_v3(0, 0, 2), e, sizeof e);
        wr_scene_add_mesh(s, "block", g, &box, e, sizeof e);
    } else {
        wr_scene_add_mesh(s, "block", g, &box, e, sizeof e);
        wr_scene_add_plane(s, "wall", w, wr_v3(0, 0, -3), wr_v3(0, 0, 2), e, sizeof e);
        wr_scene_add_sphere(s, "ball", g, wr_v3_zero(), 1.0, e, sizeof e);
    }
    return 1;
}

int main(void)
{
    char err[200] = "";
    wr_scene s;

    CHECK(wr_scene_init(&s));
    CHECK(s.nmat == 1 && strcmp(s.materials[0].name, "vacuum") == 0 && s.ambient == 0 && s.nobj == 0 && s.time == 0.0);
    wr_aabb eb = wr_scene_bounds(&s);
    CHECK(eb.lo.x > eb.hi.x);                                          /* empty scene: empty box */
    CHECK(wr_scene_validate(&s, err, sizeof err));
    wr_hit h;
    wr_ray ray0 = wr_ray_make(wr_v3_zero(), wr_v3(0, 0, 1));
    CHECK(!wr_scene_intersect(&s, &ray0, &h) && !wr_scene_occluded(&s, &ray0));

    /* materials */
    wr_material glass = wr_material_constant("glass", 1.5, 0.0);
    int g = wr_scene_add_material(&s, &glass, err, sizeof err);
    CHECK(g == 1 && wr_scene_find_material(&s, "glass") == 1 && wr_scene_find_material(&s, "nope") == -1);
    CHECK(wr_scene_add_material(&s, &glass, err, sizeof err) == -1 && has(err, "duplicate"));
    wr_material badm = wr_material_constant("bad", -1.0, 0.0);
    CHECK(wr_scene_add_material(&s, &badm, err, sizeof err) == -1 && has(err, "bad"));
    wr_material mir = wr_material_mirror("mir");
    int mi = wr_scene_add_material(&s, &mir, err, sizeof err);
    CHECK(mi == 2);
    for (int i = 0; i < 20; i++) {                                     /* growth of the array */
        char nm[16]; snprintf(nm, sizeof nm, "m%d", i);
        wr_material m = wr_material_constant(nm, 1.1, 0.0);
        CHECK(wr_scene_add_material(&s, &m, err, sizeof err) == 3 + i);
    }
    CHECK(s.nmat == 23 && strcmp(s.materials[22].name, "m19") == 0 && s.materials[1].n0 == 1.5);

    /* objects: rejected inputs leave the scene unchanged */
    CHECK(wr_scene_add_sphere(&s, "a", 99, wr_v3_zero(), 1, err, sizeof err) == -1 && has(err, "material"));
    CHECK(wr_scene_add_sphere(&s, "a", g, wr_v3_zero(), 0, err, sizeof err) == -1 && has(err, "radius"));
    CHECK(wr_scene_add_sphere(&s, "", g, wr_v3_zero(), 1, err, sizeof err) == -1);
    CHECK(wr_scene_add_plane(&s, "p", g, wr_v3_zero(), wr_v3_zero(), err, sizeof err) == -1 && has(err, "normal"));
    wr_mesh none = wr_mesh_empty();
    CHECK(wr_scene_add_mesh(&s, "m", g, &none, err, sizeof err) == -1 && has(err, "no triangles"));
    CHECK(s.nobj == 0);
    CHECK(wr_scene_add_sphere(&s, "a", g, wr_v3_zero(), 1, err, sizeof err) == 0);
    CHECK(wr_scene_add_sphere(&s, "a", g, wr_v3_zero(), 2, err, sizeof err) == -1 && has(err, "duplicate"));
    CHECK(wr_scene_find_object(&s, "a") == 0 && wr_scene_find_object(&s, "b") == -1);
    wr_mesh box;
    wr_mesh_make_box(wr_v3(0, 0, 0), wr_v3(1, 1, 1), &box);
    CHECK(wr_scene_add_mesh(&s, "a", g, &box, err, sizeof err) == -1 && box.nt == 12);       /* failure: caller keeps the mesh */
    CHECK(wr_scene_add_mesh(&s, "box", g, &box, err, sizeof err) == 1 && box.nt == 0 && box.v == NULL);   /* success: moved */
    CHECK(s.objects[1].mesh != NULL && s.objects[1].bvh.nnodes > 0);

    /* ambient */
    CHECK(wr_scene_set_ambient(&s, mi, err, sizeof err) == -1 && has(err, "dielectric") && s.ambient == 0);
    CHECK(wr_scene_set_ambient(&s, 77, err, sizeof err) == -1);
    wr_material air = wr_material_constant("air", 1.000277, 0.0);
    int ai = wr_scene_add_material(&s, &air, err, sizeof err);
    CHECK(wr_scene_set_ambient(&s, ai, err, sizeof err) == 0 && s.ambient == ai);
    wr_scene_free(&s);
    CHECK(s.nobj == 0 && s.materials == NULL);
    wr_scene_free(&s);                                                  /* double free is harmless */

    /* intersection semantics */
    CHECK(build(&s, 1));
    wr_ray r = wr_ray_make(wr_v3(0, 0, 10), wr_v3(0, 0, -1));
    CHECK(wr_scene_intersect(&s, &r, &h));
    CHECK(h.obj == wr_scene_find_object(&s, "ball") && h.prim == -1 && h.front == 1);
    CHECK_NEAR(h.t, 9.0, 1e-14);
    CHECK(wr_scene_medium_before(&s, &h) == s.ambient && wr_scene_medium_after(&s, &h) == 1);       /* entering glass */
    r = wr_ray_make(wr_v3_zero(), wr_v3(0, 0, 1));                      /* from inside the sphere */
    CHECK(wr_scene_intersect(&s, &r, &h) && h.obj == 0 && h.front == 0 && h.t > 0.999999 && h.t < 1.000001);
    CHECK(wr_scene_medium_before(&s, &h) == 1 && wr_scene_medium_after(&s, &h) == s.ambient);       /* leaving glass */
    r = wr_ray_make(wr_v3(2, 0, 10), wr_v3(0, 0, -1));                  /* misses the sphere, hits the absorbing wall */
    CHECK(wr_scene_intersect(&s, &r, &h) && h.obj == wr_scene_find_object(&s, "wall") && h.prim == -1);
    CHECK_NEAR(h.t, 13.0, 1e-14);
    wr_ray r2 = wr_ray_make(wr_v3(5, 0, 10), wr_v3(0, 0, -1));
    CHECK(wr_scene_intersect(&s, &r2, &h) && h.obj == wr_scene_find_object(&s, "block") && h.prim >= 0 && h.prim < 12);
    CHECK(h.front == 1 && h.n.z == 1);
    r2.tmax = 8.0;                                                      /* clipped before the block (t = 9) */
    CHECK(!wr_scene_intersect(&s, &r2, &h) && !wr_scene_occluded(&s, &r2));
    r2.tmax = INFINITY;
    CHECK(wr_scene_occluded(&s, &r2));
    wr_ray r3 = wr_ray_make(wr_v3(0, 0, 10), wr_v3(0, 0, 1));           /* looking away from everything */
    CHECK(!wr_scene_intersect(&s, &r3, &h) && !wr_scene_occluded(&s, &r3));
    wr_aabb b = wr_scene_bounds(&s);
    CHECK(b.lo.x == -1 && b.hi.x == 6 && b.lo.z == -1 && b.hi.z == 1);   /* plane ignored */

    /* closest hit over all objects == reference computed per object; independent of insertion order */
    wr_scene s2;
    CHECK(build(&s2, 0));
    wr_rng rng = {31415};
    int mismatch = 0, nhit = 0, order_diff = 0;
    for (int i = 0; i < 30000; i++) {
        wr_ray q = wr_ray_make(rv(&rng, -8, 8), rv(&rng, -1, 1));
        if (i % 3 == 0) q.tmax = wr_rng_range(&rng, 2, 12);
        wr_hit a, o2;
        int ha = wr_scene_intersect(&s, &q, &a);
        /* reference: ask every object separately */
        int found = 0; double bt = INFINITY; int bobj = -1;
        for (int k = 0; k < s.nobj; k++) {
            wr_hit hk; int ok = 0;
            const wr_object *ob = &s.objects[k];
            if (ob->kind == WR_OBJ_SPHERE) ok = wr_sphere_intersect(&ob->sphere, &q, &hk);
            else if (ob->kind == WR_OBJ_PLANE) ok = wr_plane_intersect(&ob->plane, &q, &hk);
            else ok = wr_mesh_intersect_brute(ob->mesh, &q, &hk);
            if (ok && hk.t < bt) { bt = hk.t; bobj = k; found = 1; }
        }
        if (ha != found || (ha && (a.t != bt))) mismatch++;
        if (ha) {
            nhit++;
            if (a.obj != bobj) mismatch++;       /* objects are disjoint here, so exact-t ties between objects do not occur */
        }
        if (wr_scene_occluded(&s, &q) != found) mismatch++;
        int h2 = wr_scene_intersect(&s2, &q, &o2);
        if (h2 != ha || (ha && o2.t != a.t)) order_diff++;
        if (ha && h2 && strcmp(s.objects[a.obj].name, s2.objects[o2.obj].name) != 0) order_diff++;
    }
    CHECK(mismatch == 0 && order_diff == 0 && nhit > 5000);
    printf("scene: closest hit over 3 object kinds, 30000 rays (%d hit): mismatches vs per-object reference %d, insertion-order differences %d\n",
           nhit, mismatch, order_diff);

    /* threads */
    {
        int N = 30000;
        wr_ray *rr = malloc(N * sizeof *rr);
        wr_hit *hs = malloc(N * sizeof *hs), *hp = malloc(N * sizeof *hp);
        unsigned char *ks = malloc(N), *kp = malloc(N);
        for (int i = 0; i < N; i++) rr[i] = wr_ray_make(rv(&rng, -8, 8), rv(&rng, -1, 1));
        for (int i = 0; i < N; i++) ks[i] = (unsigned char)wr_scene_intersect(&s, &rr[i], &hs[i]);
        omp_set_num_threads(4);
        #pragma omp parallel for schedule(dynamic, 128)
        for (int i = 0; i < N; i++) kp[i] = (unsigned char)wr_scene_intersect(&s, &rr[i], &hp[i]);
        int diff = 0;
        for (int i = 0; i < N; i++) if (ks[i] != kp[i] || (ks[i] && (hs[i].t != hp[i].t || hs[i].obj != hp[i].obj || hs[i].prim != hp[i].prim))) diff++;
        CHECK(diff == 0);
        printf("scene: 4-thread vs serial scene queries, %d rays: %d differences\n", N, diff);
        free(rr); free(hs); free(hp); free(ks); free(kp);
    }
    wr_scene_free(&s2);
    wr_scene_free(&s);

    /* validation of solids */
    wr_scene v;
    CHECK(wr_scene_init(&v));
    wr_material bk = wr_material_constant("bk7", 1.5168, 0.0), mr = wr_material_mirror("mr");
    int bk_i = wr_scene_add_material(&v, &bk, err, sizeof err), mr_i = wr_scene_add_material(&v, &mr, err, sizeof err);
    wr_mesh tri = wr_mesh_empty();
    tri.v = malloc(3 * sizeof *tri.v); tri.idx = malloc(3 * sizeof *tri.idx);
    tri.v[0] = wr_v3(0, 0, 0); tri.v[1] = wr_v3(1, 0, 0); tri.v[2] = wr_v3(0, 1, 0);
    tri.idx[0] = 0; tri.idx[1] = 1; tri.idx[2] = 2; tri.nv = 3; tri.nt = 1;
    CHECK(wr_scene_add_mesh(&v, "sheet", mr_i, &tri, err, sizeof err) == 0);            /* an open mirror sheet is fine */
    CHECK(wr_scene_validate(&v, err, sizeof err));
    wr_mesh tri2 = wr_mesh_empty();
    tri2.v = malloc(3 * sizeof *tri2.v); tri2.idx = malloc(3 * sizeof *tri2.idx);
    memcpy(tri2.v, (wr_vec3[3]){wr_v3(0, 0, 0), wr_v3(1, 0, 0), wr_v3(0, 1, 0)}, 3 * sizeof *tri2.v);
    tri2.idx[0] = 0; tri2.idx[1] = 1; tri2.idx[2] = 2; tri2.nv = 3; tri2.nt = 1;
    CHECK(wr_scene_add_mesh(&v, "open_glass", bk_i, &tri2, err, sizeof err) == 1);
    CHECK(!wr_scene_validate(&v, err, sizeof err) && has(err, "open_glass") && has(err, "closed"));
    wr_scene_free(&v);

    CHECK(wr_scene_init(&v));
    bk_i = wr_scene_add_material(&v, &bk, err, sizeof err);
    wr_mesh inward;
    wr_mesh_make_box(wr_v3(0, 0, 0), wr_v3(1, 1, 1), &inward);
    for (size_t i = 0; i < inward.nt; i++) { uint32_t t = inward.idx[3 * i + 1]; inward.idx[3 * i + 1] = inward.idx[3 * i + 2]; inward.idx[3 * i + 2] = t; }
    CHECK(wr_mesh_is_closed_oriented(&inward) && wr_mesh_volume(&inward) < 0);
    CHECK(wr_scene_add_mesh(&v, "inside_out", bk_i, &inward, err, sizeof err) == 0);
    CHECK(!wr_scene_validate(&v, err, sizeof err) && has(err, "inside_out") && has(err, "wound inward"));
    wr_scene_free(&v);

    /* closed-mesh predicate */
    wr_mesh ok;
    wr_mesh_make_box(wr_v3(0, 0, 0), wr_v3(1, 1, 1), &ok);
    CHECK(wr_mesh_is_closed_oriented(&ok));
    ok.nt = 11;                                                         /* drop one triangle: open surface */
    CHECK(!wr_mesh_is_closed_oriented(&ok));
    ok.nt = 12;
    ok.idx[4] = ok.idx[1];                                              /* corrupt: duplicate directed edge / orientation */
    CHECK(!wr_mesh_is_closed_oriented(&ok));
    wr_mesh_free(&ok);
    wr_mesh nomesh = wr_mesh_empty();
    CHECK(!wr_mesh_is_closed_oriented(&nomesh));

    /* lights and camera */
    CHECK(wr_scene_init(&v));
    wr_light l1 = wr_light_point("lamp", wr_v3(0, 1, 0), 550e-9, 1.0);
    CHECK(wr_scene_add_light(&v, &l1, err, sizeof err) == 0 && wr_scene_find_light(&v, "lamp") == 0);
    CHECK(wr_scene_add_light(&v, &l1, err, sizeof err) == -1 && has(err, "duplicate"));
    wr_light l2 = l1; snprintf(l2.name, sizeof l2.name, "bad"); l2.amplitude = -1;
    CHECK(wr_scene_add_light(&v, &l2, err, sizeof err) == -1 && has(err, "amplitude") && v.nlight == 1);
    for (int i = 0; i < 20; i++) {
        wr_light li = l1; snprintf(li.name, sizeof li.name, "l%d", i);
        CHECK(wr_scene_add_light(&v, &li, err, sizeof err) == 1 + i);
    }
    CHECK(!v.has_camera);
    wr_camera cam;
    wr_camera_pinhole(&cam, wr_v3(0, 0, 5), wr_v3_zero(), wr_v3(0, 1, 0), 1.0, 64, 48);
    wr_scene_set_camera(&v, &cam);
    CHECK(v.has_camera && v.camera.width == 64);
    CHECK(wr_scene_validate(&v, err, sizeof err));
    wr_scene_free(&v);

    TEST_REPORT("scene");
}
