#include <stdlib.h>
#include <string.h>
#include "wr_test.h"
#include "wr_mesh.h"
#include "wr_optics.h"

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : (x > y);
}

/* Closed, consistently oriented 2-manifold: every directed edge appears exactly once and its reverse exactly once. */
static int is_closed_oriented(const wr_mesh *m)
{
    size_t ne = 3 * m->nt;
    uint64_t *k = malloc(ne * sizeof *k);
    for (size_t t = 0; t < m->nt; t++)
        for (int e = 0; e < 3; e++)
            k[3 * t + e] = (uint64_t)m->idx[3 * t + e] * m->nv + m->idx[3 * t + (e + 1) % 3];
    qsort(k, ne, sizeof *k, cmp_u64);
    int ok = 1;
    for (size_t i = 1; i < ne; i++) if (k[i] == k[i - 1]) ok = 0;            /* directed edge used twice */
    for (size_t t = 0; t < m->nt && ok; t++)
        for (int e = 0; e < 3; e++) {
            uint64_t rev = (uint64_t)m->idx[3 * t + (e + 1) % 3] * m->nv + m->idx[3 * t + e];
            if (!bsearch(&rev, k, ne, sizeof *k, cmp_u64)) { ok = 0; break; }
        }
    free(k);
    return ok;
}

static int all_outward(const wr_mesh *m, wr_vec3 c)
{
    for (size_t i = 0; i < m->nt; i++) {
        wr_triangle t = wr_mesh_triangle(m, i);
        wr_vec3 cen = wr_v3_scale(wr_v3_add(wr_v3_add(t.v0, t.v1), t.v2), 1.0 / 3.0);
        if (!(wr_v3_dot(wr_triangle_normal(&t), wr_v3_sub(cen, c)) > 0.0)) return 0;
    }
    return 1;
}

static const char *CUBE_QUADS =
    "# unit cube, quads\n"
    "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv 0 0 1\nv 1 0 1\nv 1 1 1\nv 0 1 1\n"
    "f 1 4 3 2\nf 5 6 7 8\nf 1 2 6 5\nf 2 3 7 6\nf 3 4 8 7\nf 4 1 5 8\n";

int main(void)
{
    char err[128] = "";
    wr_mesh m;

    /* OBJ: quads, fan triangulation, volume/area/bounds */
    CHECK(wr_mesh_parse_obj(CUBE_QUADS, &m, err, sizeof err));
    CHECK(m.nv == 8 && m.nt == 12);
    CHECK_NEAR(wr_mesh_volume(&m), 1.0, 1e-15);
    CHECK_NEAR(wr_mesh_area(&m), 6.0, 1e-15);
    wr_aabb b = wr_mesh_bounds(&m);
    CHECK(b.lo.x == 0 && b.hi.z == 1);
    CHECK(is_closed_oriented(&m));
    wr_mesh_free(&m);
    CHECK(m.v == NULL && m.nt == 0);

    /* token forms, comments, CRLF, tabs, ignored statements, negative indices */
    const char *fancy =
        "# comment\r\nmtllib x.mtl\r\no thing\r\ng grp\r\nusemtl red\r\ns off\r\n\r\n"
        "v 0 0 0\r\nv\t1 0 0\r\nv 0 1 0\r\nv 1 1 0 1.0\r\n"      /* optional w ignored */
        "vt 0 0\r\nvt 1 0\r\nvt 0 1\r\nvn 0 0 1\r\n"
        "f 1/1/1 2/2/1 3/3/1\r\n"                                  /* v/vt/vn */
        "f 2//1 4//1 3//1\r\n"                                     /* v//vn */
        "f -4 -3 -2\r\n"                                           /* relative */
        "f 1/1 2/2 3/3 \r\n";                                      /* v/vt, trailing space */
    CHECK(wr_mesh_parse_obj(fancy, &m, err, sizeof err));
    CHECK(m.nv == 4 && m.nt == 4);
    CHECK(m.idx[0] == 0 && m.idx[1] == 1 && m.idx[2] == 2);
    CHECK(m.idx[3] == 1 && m.idx[4] == 3 && m.idx[5] == 2);
    CHECK(m.idx[6] == 0 && m.idx[7] == 1 && m.idx[8] == 2);       /* -4,-3,-2 with 4 vertices defined */
    wr_mesh_free(&m);

    /* negative indices are relative to the vertices defined SO FAR */
    CHECK(wr_mesh_parse_obj("v 0 0 0\nv 1 0 0\nv 0 1 0\nf -3 -2 -1\nv 5 5 5\nf -1 -2 -3\n", &m, err, sizeof err));
    CHECK(m.nt == 2 && m.idx[0] == 0 && m.idx[2] == 2 && m.idx[3] == 3 && m.idx[5] == 1);
    wr_mesh_free(&m);

    /* regular pentagon, fan triangulation: area = (5/2) R^2 sin(72 deg), no file-order assumption */
    {
        char buf[512] = "";
        for (int i = 0; i < 5; i++) {
            char l[64];
            snprintf(l, sizeof l, "v %.17g %.17g 0\n", cos(WR_TWO_PI * i / 5), sin(WR_TWO_PI * i / 5));
            strcat(buf, l);
        }
        strcat(buf, "f 1 2 3 4 5\n");
        CHECK(wr_mesh_parse_obj(buf, &m, err, sizeof err));
        CHECK(m.nt == 3);
        CHECK_NEAR(wr_mesh_area(&m), 2.5 * sin(WR_TWO_PI / 5), 1e-14);
        wr_mesh_free(&m);
    }

    /* forward reference (face before its vertices) is legal and resolved at the end */
    CHECK(wr_mesh_parse_obj("f 1 2 3\nv 0 0 0\nv 1 0 0\nv 0 1 0\n", &m, err, sizeof err));
    CHECK(m.nt == 1);
    wr_mesh_free(&m);

    /* empty input is a valid empty mesh */
    CHECK(wr_mesh_parse_obj("", &m, err, sizeof err) && m.nt == 0 && m.nv == 0);
    CHECK(wr_mesh_parse_obj("# nothing\n", &m, err, sizeof err) && m.nt == 0);

    /* malformed input: reported, no leak, mesh left empty */
    static const char *bad[] = {
        "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 0 1 2\n",           /* index 0 */
        "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 4\n",           /* beyond vertex count */
        "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2\n",             /* two vertices */
        "v 0 0 0\nv 1 0 0\nv 0 1 0\nf -1 -2 -4\n",        /* relative beyond start */
        "v 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n",             /* two coordinates */
        "v 0 0 nan\nv 1 0 0\nv 0 1 0\nf 1 2 3\n",         /* non-finite */
        "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 x\n"};          /* garbage token */
    for (size_t i = 0; i < sizeof bad / sizeof *bad; i++) {
        err[0] = '\0';
        int ok = wr_mesh_parse_obj(bad[i], &m, err, sizeof err);
        CHECK(!ok && m.v == NULL && m.idx == NULL && m.nt == 0 && err[0] != '\0');
    }

    /* file round trip */
    FILE *f = fopen("build/test_cube.obj", "wb");
    CHECK(f != NULL);
    if (f) {
        fputs(CUBE_QUADS, f); fclose(f);
        CHECK(wr_mesh_load_obj("build/test_cube.obj", &m, err, sizeof err) && m.nt == 12);
        wr_mesh_free(&m);
    }
    CHECK(!wr_mesh_load_obj("build/does_not_exist.obj", &m, err, sizeof err) && m.nt == 0);

    /* generators */
    CHECK(wr_mesh_make_box(wr_v3(-1, -2, -3), wr_v3(1, 2, 3), &m));
    CHECK_NEAR(wr_mesh_volume(&m), 2 * 4 * 6, 1e-13);
    CHECK_NEAR(wr_mesh_area(&m), 2 * (2 * 4 + 2 * 6 + 4 * 6), 1e-13);
    CHECK(is_closed_oriented(&m) && all_outward(&m, wr_v3_zero()));
    wr_mesh_free(&m);

    wr_vec3 c = wr_v3(0.3, -0.2, 0.7);
    double R = 1.7;
    CHECK(wr_mesh_make_uv_sphere(c, R, 64, 32, &m));
    CHECK(m.nv == 2 + 31 * 64 && m.nt == 2 * 64 + 2 * 30 * 64);
    CHECK(is_closed_oriented(&m) && all_outward(&m, c));
    double worst = 0;
    for (size_t i = 0; i < m.nv; i++) worst = fmax(worst, fabs(wr_v3_dist(m.v[i], c) - R));
    CHECK_NEAR(worst, 0, 1e-14);                                 /* every vertex on the sphere */
    double Vex = 4.0 / 3.0 * WR_PI * R * R * R, Aex = 4.0 * WR_PI * R * R;
    double vr = wr_mesh_volume(&m) / Vex, ar = wr_mesh_area(&m) / Aex;
    CHECK(vr > 0.99 && vr < 1.0);                               /* inscribed: smaller than the sphere, never larger */
    CHECK(ar > 0.99 && ar < 1.0);
    printf("mesh: UV sphere 64x32 volume/exact %.5f, area/exact %.5f (inscribed polyhedron)\n", vr, ar);

    /* transforms: scale, rotation, mirror (winding must flip), singular */
    double V0 = wr_mesh_volume(&m);
    wr_xform S = wr_xform_scale(wr_v3(2, 3, 4));
    CHECK(wr_mesh_transform(&m, &S));
    CHECK_NEAR(wr_mesh_volume(&m) / V0, 24.0, 1e-12);
    wr_xform Sinv; CHECK(wr_xform_inverse(&S, &Sinv)); CHECK(wr_mesh_transform(&m, &Sinv));
    wr_xform Rot = wr_xform_rotate(wr_v3(1, 2, 3), 0.8);
    CHECK(wr_mesh_transform(&m, &Rot));
    CHECK_NEAR(wr_mesh_volume(&m) / V0, 1.0, 1e-12);
    wr_vec3 c2 = wr_xform_point(&Rot, c);
    CHECK(all_outward(&m, c2));
    wr_xform Mir = wr_xform_scale(wr_v3(-1, 1, 1));
    CHECK(wr_mesh_transform(&m, &Mir));
    CHECK_NEAR(wr_mesh_volume(&m) / V0, 1.0, 1e-12);            /* would be -1 without the winding flip */
    CHECK(all_outward(&m, wr_xform_point(&Mir, c2)) && is_closed_oriented(&m));
    wr_vec3 before = m.v[5];
    wr_xform flat = wr_xform_scale(wr_v3(1, 0, 1));
    CHECK(!wr_mesh_transform(&m, &flat));
    CHECK(m.v[5].x == before.x && m.v[5].y == before.y);        /* untouched on failure */
    wr_mesh_free(&m);

    CHECK(!wr_mesh_make_uv_sphere(c, 1.0, 2, 8, &m));
    CHECK(!wr_mesh_make_uv_sphere(c, 1.0, 8, 1, &m));
    CHECK(!wr_mesh_make_uv_sphere(c, -1.0, 8, 8, &m));
    CHECK(wr_mesh_make_uv_sphere(c, 1.0, 3, 2, &m) && m.nt == 6);   /* smallest legal: a triangular bipyramid */
    CHECK(is_closed_oriented(&m) && all_outward(&m, c));
    wr_mesh_free(&m);

    TEST_REPORT("mesh");
}
