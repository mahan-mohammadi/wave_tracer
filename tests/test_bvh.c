#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include "wr_test.h"
#include "wr_rng.h"
#include "wr_bvh.h"
#include "wr_optics.h"

static wr_vec3 rv(wr_rng *r, double lo, double hi)
{
    return wr_v3(wr_rng_range(r, lo, hi), wr_rng_range(r, lo, hi), wr_rng_range(r, lo, hi));
}

static int alloc_soup(wr_mesh *m, size_t nt)
{
    m->v = malloc(3 * nt * sizeof *m->v);
    m->idx = malloc(3 * nt * sizeof *m->idx);
    m->nv = 3 * nt; m->nt = nt;
    for (size_t i = 0; i < 3 * nt; i++) m->idx[i] = (uint32_t)i;
    return m->v && m->idx;
}

/* Compare BVH against brute force. Returns number of rays that hit. */
static size_t compare(const wr_mesh *m, const wr_bvh *bvh, const wr_ray *rays, size_t n, int *mismatch)
{
    size_t hits = 0;
    for (size_t i = 0; i < n; i++) {
        wr_hit a, b;
        int ha = wr_bvh_intersect(bvh, &rays[i], &a);
        int hb = wr_mesh_intersect_brute(m, &rays[i], &b);
        int oc = wr_bvh_occluded(bvh, &rays[i]);
        if (ha != hb || oc != hb) { (*mismatch)++; continue; }
        if (ha) {
            hits++;
            /* the same routine decides every hit, so t must be bit-identical; prim may differ only on exact ties */
            if (a.t != b.t) (*mismatch)++;
            else if (a.prim != b.prim) {
                wr_triangle t1 = wr_mesh_triangle(m, (size_t)a.prim), t2 = wr_mesh_triangle(m, (size_t)b.prim);
                wr_hit h1, h2;
                if (!(wr_triangle_intersect(&t1, &rays[i], &h1) && wr_triangle_intersect(&t2, &rays[i], &h2) && h1.t == h2.t))
                    (*mismatch)++;
            }
            if (a.prim < 0 || (size_t)a.prim >= m->nt) (*mismatch)++;
        }
    }
    return hits;
}

int main(void)
{
    wr_rng rng = {0xC0FFEEULL};
    wr_mesh m;
    wr_bvh bvh;

    /* empty mesh: valid BVH, never hit */
    m = wr_mesh_empty();
    CHECK(wr_bvh_build(&bvh, &m));
    wr_ray r0 = wr_ray_make(wr_v3(0, 0, 5), wr_v3(0, 0, -1));
    wr_hit h;
    CHECK(!wr_bvh_intersect(&bvh, &r0, &h) && !wr_bvh_occluded(&bvh, &r0));
    wr_bvh_free(&bvh);

    /* random triangle soup vs brute force: random, tmax-clipped and axis-aligned rays */
    CHECK(alloc_soup(&m, 3000));
    for (size_t i = 0; i < m.nt; i++) {
        wr_vec3 c = rv(&rng, -5, 5);
        for (int k = 0; k < 3; k++) m.v[3 * i + k] = wr_v3_add(c, rv(&rng, -0.6, 0.6));
    }
    CHECK(wr_bvh_build(&bvh, &m));
    size_t N = 40000;
    wr_ray *rays = malloc(N * sizeof *rays);
    for (size_t i = 0; i < N; i++) {
        wr_vec3 o = rv(&rng, -8, 8);
        wr_vec3 d = rv(&rng, -1, 1);
        if (i % 5 == 1) { d = wr_v3_zero(); (&d.x)[i % 3] = (i & 8) ? 1.0 : -1.0; }       /* axis aligned: exact zeros in d */
        if (i % 5 == 2) { d.x = 0.0; }                                                    /* parallel to one slab */
        if (wr_v3_len2(d) < 1e-6) d = wr_v3(0, 0, 1);
        rays[i] = wr_ray_make(o, d);
        if (i % 4 == 3) rays[i].tmax = wr_rng_range(&rng, 0.5, 6.0);                      /* clipped rays */
        if (i % 7 == 6) rays[i].tmin = wr_rng_range(&rng, 0.0, 3.0);
    }
    int mismatch = 0;
    size_t hits = compare(&m, &bvh, rays, N, &mismatch);
    CHECK(mismatch == 0);
    CHECK(hits > 3000);
    printf("bvh: random soup, %zu triangles, %zu rays (%zu hit): mismatches vs brute force = %d\n", m.nt, N, hits, mismatch);
    printf("bvh: nodes %zu, leaves %zu, max depth %d, max leaf %u, SAH cost %.2f triangle tests/ray (brute force: %zu)\n",
           bvh.nnodes, bvh.nleaves, bvh.max_depth, bvh.max_leaf, bvh.sah_cost, m.nt);
    CHECK(bvh.max_leaf <= 4 && bvh.max_depth < 64);
    CHECK(bvh.nnodes <= 2 * m.nt);

    /* every triangle appears in exactly one leaf */
    {
        unsigned char *seen = calloc(m.nt, 1);
        size_t covered = 0, dup = 0;
        for (size_t i = 0; i < bvh.nnodes; i++)
            if (bvh.nodes[i].n)
                for (uint32_t k = 0; k < bvh.nodes[i].n; k++) {
                    uint32_t t = bvh.order[bvh.nodes[i].a + k];
                    if (seen[t]) dup++; else { seen[t] = 1; covered++; }
                }
        CHECK(covered == m.nt && dup == 0);
        free(seen);
    }

    /* concurrent queries (OpenMP) must give bit-identical results to a serial run */
    {
        wr_hit *ser = malloc(N * sizeof *ser), *par = malloc(N * sizeof *par);
        unsigned char *hs = malloc(N), *hp = malloc(N);
        for (size_t i = 0; i < N; i++) hs[i] = (unsigned char)wr_bvh_intersect(&bvh, &rays[i], &ser[i]);
        omp_set_num_threads(4);
        #pragma omp parallel for schedule(dynamic, 256)
        for (long i = 0; i < (long)N; i++) hp[i] = (unsigned char)wr_bvh_intersect(&bvh, &rays[i], &par[i]);
        int diff = 0;
        for (size_t i = 0; i < N; i++) {
            if (hs[i] != hp[i]) diff++;
            else if (hs[i] && (ser[i].t != par[i].t || ser[i].prim != par[i].prim)) diff++;
        }
        CHECK(diff == 0);
        printf("bvh: 4-thread vs serial traversal, %zu rays: %d differences\n", N, diff);
        free(ser); free(par); free(hs); free(hp);
    }
    free(rays);
    wr_bvh_free(&bvh);
    wr_mesh_free(&m);

    /* coincident centroids: 500 identical triangles. No spatial split exists; the median fallback must keep
       leaves small and the tree shallow. */
    CHECK(alloc_soup(&m, 500));
    for (size_t i = 0; i < m.nt; i++) {
        m.v[3 * i] = wr_v3(0, 0, 0); m.v[3 * i + 1] = wr_v3(1, 0, 0); m.v[3 * i + 2] = wr_v3(0, 1, 0);
    }
    CHECK(wr_bvh_build(&bvh, &m));
    CHECK(bvh.max_leaf <= 4 && bvh.max_depth <= 12);
    wr_ray rd = wr_ray_make(wr_v3(0.2, 0.2, 1), wr_v3(0, 0, -1));
    CHECK(wr_bvh_intersect(&bvh, &rd, &h) && fabs(h.t - 1.0) < 1e-15 && wr_bvh_occluded(&bvh, &rd));
    printf("bvh: 500 identical triangles -> %zu nodes, depth %d, max leaf %u\n", bvh.nnodes, bvh.max_depth, bvh.max_leaf);
    wr_bvh_free(&bvh);
    wr_mesh_free(&m);

    /* single triangle */
    CHECK(alloc_soup(&m, 1));
    m.v[0] = wr_v3(0, 0, 0); m.v[1] = wr_v3(1, 0, 0); m.v[2] = wr_v3(0, 1, 0);
    CHECK(wr_bvh_build(&bvh, &m) && bvh.nnodes == 1 && bvh.nleaves == 1);
    CHECK(wr_bvh_intersect(&bvh, &rd, &h) && h.prim == 0);
    wr_ray miss = wr_ray_make(wr_v3(2, 2, 1), wr_v3(0, 0, -1));
    CHECK(!wr_bvh_intersect(&bvh, &miss, &h));
    wr_bvh_free(&bvh);
    wr_mesh_free(&m);

    /* coplanar sheet (z = 0): flat boxes, rays through shared edges and vertices. The BVH must agree with
       brute force exactly, including where the non-watertight triangle test itself lets a ray leak. */
    {
        int G = 40;
        CHECK(alloc_soup(&m, (size_t)G * G * 2));
        size_t t = 0;
        for (int i = 0; i < G; i++)
            for (int j = 0; j < G; j++) {
                wr_vec3 p00 = wr_v3(i * 0.1, j * 0.1, 0), p10 = wr_v3((i + 1) * 0.1, j * 0.1, 0);
                wr_vec3 p01 = wr_v3(i * 0.1, (j + 1) * 0.1, 0), p11 = wr_v3((i + 1) * 0.1, (j + 1) * 0.1, 0);
                m.v[3 * t] = p00; m.v[3 * t + 1] = p10; m.v[3 * t + 2] = p11; t++;
                m.v[3 * t] = p00; m.v[3 * t + 1] = p11; m.v[3 * t + 2] = p01; t++;
            }
        CHECK(wr_bvh_build(&bvh, &m));
        size_t NR = (size_t)(G + 1) * (G + 1) + 20000;
        wr_ray *rr = malloc(NR * sizeof *rr);
        size_t n = 0, vertex_rays = 0;
        for (int i = 0; i <= G; i++)
            for (int j = 0; j <= G; j++) {                          /* exactly through grid vertices */
                rr[n++] = wr_ray_make(wr_v3(i * 0.1, j * 0.1, 1.0), wr_v3(0, 0, -1)); vertex_rays++;
            }
        while (n < NR) {                                            /* random interior points, mixed directions */
            wr_vec3 tgt = wr_v3(wr_rng_range(&rng, 0.001, 3.999), wr_rng_range(&rng, 0.001, 3.999), 0);
            wr_vec3 o = wr_v3(tgt.x + wr_rng_range(&rng, -1, 1), tgt.y + wr_rng_range(&rng, -1, 1), wr_rng_range(&rng, 0.2, 2));
            if (n & 1) o.z = -o.z;
            rr[n++] = wr_ray_make(o, wr_v3_sub(tgt, o));
        }
        int mm = 0;
        size_t hh = compare(&m, &bvh, rr, n, &mm);
        CHECK(mm == 0);
        size_t vhits = 0, rhits = 0;
        for (size_t i = 0; i < n; i++) {
            wr_hit hx;
            int hit = wr_bvh_intersect(&bvh, &rr[i], &hx);
            if (i < vertex_rays) vhits += hit; else rhits += hit;
        }
        CHECK(rhits == n - vertex_rays);                            /* every interior ray (>= 1 mm from the border) hits */
        printf("bvh: coplanar sheet, %zu rays: mismatches vs brute force %d; interior rays hit %zu/%zu; rays through exact grid vertices hit %zu/%zu (rest leak through the non-watertight test)\n",
               n, mm, rhits, n - vertex_rays, vhits, vertex_rays);
        (void)hh;
        free(rr);
        wr_bvh_free(&bvh);
        wr_mesh_free(&m);
    }

    /* closed sphere mesh: rays aimed inside the silhouette always hit the outside, with an outward normal */
    {
        wr_vec3 c = wr_v3(0.5, -0.5, 0.2);
        double R = 2.0;
        CHECK(wr_mesh_make_uv_sphere(c, R, 64, 32, &m));
        CHECK(wr_bvh_build(&bvh, &m));
        int bad = 0, n = 0;
        double rmin = 1e9, rmax = 0;
        for (int i = 0; i < 20000; i++) {
            wr_vec3 dir = wr_v3_normalize(rv(&rng, -1, 1));
            wr_vec3 o = wr_v3_madd(c, dir, 10.0);
            wr_vec3 off = rv(&rng, -1, 1);
            wr_vec3 tgt = wr_v3_madd(c, off, 0.5 * R);                         /* well inside the silhouette */
            wr_ray r = wr_ray_make(o, wr_v3_sub(tgt, o));
            if (wr_v3_len(off) > 1.0) continue;
            n++;
            if (!wr_bvh_intersect(&bvh, &r, &h) || !h.front) { bad++; continue; }
            if (wr_v3_dot(h.n, wr_v3_sub(h.p, c)) <= 0.0) bad++;
            double rad = wr_v3_dist(h.p, c);
            rmin = fmin(rmin, rad); rmax = fmax(rmax, rad);
        }
        CHECK(bad == 0);
        CHECK(rmax <= R * (1 + 1e-12));                                         /* inscribed mesh never pokes out */
        CHECK(rmin > 0.99 * R * cos(WR_PI / 64));                               /* and stays within the polyhedron depth */
        printf("bvh: sphere mesh, %d rays, hit radius in [%.5f, %.5f] x R, bad = %d\n", n, rmin / R, rmax / R, bad);

        /* timing (informational only, not a pass/fail criterion) */
        size_t NT = 200000;
        wr_ray *tr = malloc(NT * sizeof *tr);
        for (size_t i = 0; i < NT; i++) {
            wr_vec3 o = wr_v3(wr_rng_range(&rng, -1.5, 1.5), wr_rng_range(&rng, -1.5, 1.5), 10);
            tr[i] = wr_ray_make(o, wr_v3(0, 0, -1));
        }
        double t0 = omp_get_wtime();
        size_t th = 0;
        for (size_t i = 0; i < NT; i++) th += wr_bvh_intersect(&bvh, &tr[i], &h);
        double tb = omp_get_wtime() - t0;
        size_t NB = 2000;
        t0 = omp_get_wtime();
        for (size_t i = 0; i < NB; i++) th += wr_mesh_intersect_brute(&m, &tr[i], &h);
        double tbr = omp_get_wtime() - t0;
        printf("bvh: timing (1 thread, %zu triangles): BVH %.2f Mray/s, brute force %.4f Mray/s, speedup %.0fx (informational)\n",
               m.nt, NT / tb / 1e6, NB / tbr / 1e6, (NT / tb) / (NB / tbr));
        (void)th;
        free(tr);
        wr_bvh_free(&bvh);
        wr_mesh_free(&m);
    }

    TEST_REPORT("bvh");
}
