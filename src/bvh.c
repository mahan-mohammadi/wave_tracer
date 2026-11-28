#include <math.h>
#include <stdlib.h>
#include "wr_bvh.h"

#define NBINS      16
#define MAX_LEAF   4
#define TRAV_COST  0.5
#define FORCE_MEDIAN_DEPTH 64   /* beyond this depth split by median so depth stays < 64 + log2(n) < 128 */
#define STACK_SIZE 128

typedef struct {
    wr_bvh      *bvh;
    const wr_aabb *tb;       /* per-triangle bounds */
    const wr_vec3 *tc;       /* per-triangle centroids */
    wr_bvh_node *nodes;
    size_t       nn;
    double       pad;
} build_ctx;

static double box_area(wr_aabb b)
{
    return (b.hi.x < b.lo.x) ? 0.0 : wr_aabb_surface_area(&b);
}

static wr_aabb padded(wr_aabb b, double pad)
{
    b.lo = wr_v3_sub(b.lo, wr_v3(pad, pad, pad));
    b.hi = wr_v3_add(b.hi, wr_v3(pad, pad, pad));
    return b;
}

/* Quickselect on order[lo..hi) by centroid coordinate; puts the nth element in place. */
static void select_nth(const build_ctx *c, uint32_t *a, long lo, long hi, long nth, int axis)
{
    while (hi - lo > 1) {
        double pivot = wr_v3_get(c->tc[a[(lo + hi) / 2]], axis);
        long i = lo, j = hi - 1;
        while (i <= j) {
            while (wr_v3_get(c->tc[a[i]], axis) < pivot) i++;
            while (wr_v3_get(c->tc[a[j]], axis) > pivot) j--;
            if (i <= j) { uint32_t t = a[i]; a[i] = a[j]; a[j] = t; i++; j--; }
        }
        if (nth <= j) hi = j + 1;
        else if (nth >= i) lo = i;
        else break;
    }
}

static uint32_t build_rec(build_ctx *c, uint32_t s, uint32_t e, int depth)
{
    uint32_t self = (uint32_t)c->nn++;
    uint32_t *ord = c->bvh->order;
    uint32_t n = e - s;
    wr_aabb box = wr_aabb_empty(), cbox = wr_aabb_empty();
    for (uint32_t i = s; i < e; i++) {
        box = wr_aabb_union(box, c->tb[ord[i]]);
        cbox = wr_aabb_expand(cbox, c->tc[ord[i]]);
    }
    wr_bvh_node *node = &c->nodes[self];
    node->box = padded(box, c->pad);
    if (depth > c->bvh->max_depth) c->bvh->max_depth = depth;

    int axis = wr_aabb_longest_axis(&cbox);
    double clo = wr_v3_get(cbox.lo, axis), ext = wr_v3_get(cbox.hi, axis) - clo;
    uint32_t mid = 0;                       /* split position, 0 = make a leaf */

    if (n > 1 && ext > 0.0 && depth < FORCE_MEDIAN_DEPTH) {
        wr_aabb bb[NBINS];
        uint32_t cnt[NBINS] = {0};
        for (int k = 0; k < NBINS; k++) bb[k] = wr_aabb_empty();
        for (uint32_t i = s; i < e; i++) {
            int k = (int)(NBINS * (wr_v3_get(c->tc[ord[i]], axis) - clo) / ext);
            if (k > NBINS - 1) k = NBINS - 1;
            if (k < 0) k = 0;
            cnt[k]++;
            bb[k] = wr_aabb_union(bb[k], c->tb[ord[i]]);
        }
        double left_cost[NBINS - 1];
        wr_aabb acc = wr_aabb_empty();
        uint32_t nacc = 0;
        for (int k = 0; k < NBINS - 1; k++) {
            acc = wr_aabb_union(acc, bb[k]); nacc += cnt[k];
            left_cost[k] = box_area(acc) * nacc;
        }
        acc = wr_aabb_empty(); nacc = 0;
        double best = INFINITY; int bestk = -1;
        for (int k = NBINS - 1; k >= 1; k--) {
            acc = wr_aabb_union(acc, bb[k]); nacc += cnt[k];
            double cost = left_cost[k - 1] + box_area(acc) * nacc;      /* split between bin k-1 and k */
            uint32_t nleft = n - nacc;
            if (nleft > 0 && nacc > 0 && cost < best) { best = cost; bestk = k; }
        }
        double sa = box_area(box);
        double split_cost = TRAV_COST + (sa > 0.0 ? best / sa : (double)n);
        if (bestk >= 0 && !(n <= MAX_LEAF && split_cost >= (double)n)) {
            uint32_t i = s, j = e;                                       /* partition: bin < bestk goes left */
            while (i < j) {
                int k = (int)(NBINS * (wr_v3_get(c->tc[ord[i]], axis) - clo) / ext);
                if (k > NBINS - 1) k = NBINS - 1;
                if (k < bestk) i++;
                else { j--; uint32_t t = ord[i]; ord[i] = ord[j]; ord[j] = t; }
            }
            if (i > s && i < e) mid = i;
        }
    }
    if (mid == 0 && n > MAX_LEAF) {          /* no useful SAH split (coincident centroids, depth guard): halve by median */
        uint32_t m = s + n / 2;
        if (ext > 0.0) select_nth(c, ord, (long)s, (long)e, (long)m, axis);
        mid = m;
    }
    if (mid == 0) {                           /* leaf */
        node->a = s; node->n = n; node->axis = 0;
        c->bvh->nleaves++;
        if (n > c->bvh->max_leaf) c->bvh->max_leaf = n;
        return self;
    }
    node->n = 0; node->axis = (uint32_t)axis;
    build_rec(c, s, mid, depth + 1);
    uint32_t right = build_rec(c, mid, e, depth + 1);
    c->nodes[self].a = right;                 /* nodes may not move: array is preallocated */
    return self;
}

int wr_bvh_build(wr_bvh *bvh, const wr_mesh *mesh)
{
    bvh->mesh = mesh; bvh->nodes = NULL; bvh->nnodes = 0; bvh->order = NULL;
    bvh->nleaves = 0; bvh->max_depth = 0; bvh->max_leaf = 0; bvh->sah_cost = 0.0;
    size_t nt = mesh->nt;
    if (nt == 0) return 1;
    if (nt > 0x7FFFFFFFu) return 0;
    bvh->order = malloc(nt * sizeof *bvh->order);
    bvh->nodes = malloc(2 * nt * sizeof *bvh->nodes);
    wr_aabb *tb = malloc(nt * sizeof *tb);
    wr_vec3 *tc = malloc(nt * sizeof *tc);
    if (!bvh->order || !bvh->nodes || !tb || !tc) { free(tb); free(tc); wr_bvh_free(bvh); return 0; }
    wr_aabb all = wr_aabb_empty();
    for (size_t i = 0; i < nt; i++) {
        wr_triangle t = wr_mesh_triangle(mesh, i);
        tb[i] = wr_triangle_bounds(&t);
        tc[i] = wr_aabb_centroid(&tb[i]);
        all = wr_aabb_union(all, tb[i]);
        bvh->order[i] = (uint32_t)i;
    }
    double scale = fmax(wr_v3_len(all.lo), wr_v3_len(all.hi));
    build_ctx c = {bvh, tb, tc, bvh->nodes, 0, 1e-12 * scale};
    build_rec(&c, 0, (uint32_t)nt, 0);
    bvh->nnodes = c.nn;
    /* SAH cost estimate: sum over nodes of area/root_area * (trav or leaf count) */
    double root = box_area(bvh->nodes[0].box), cost = 0.0;
    for (size_t i = 0; i < bvh->nnodes; i++) {
        double w = root > 0.0 ? box_area(bvh->nodes[i].box) / root : 1.0;
        cost += w * (bvh->nodes[i].n ? (double)bvh->nodes[i].n : TRAV_COST);
    }
    bvh->sah_cost = cost;
    free(tb); free(tc);
    return 1;
}

void wr_bvh_free(wr_bvh *bvh)
{
    free(bvh->nodes); free(bvh->order);
    bvh->nodes = NULL; bvh->order = NULL; bvh->nnodes = 0; bvh->nleaves = 0;
}

/* Slab test with precomputed 1/d. zero[a] marks d[a] == 0 (parallel to that slab). */
static inline int box_hit(const wr_aabb *b, const wr_ray *r, wr_vec3 inv, const int zero[3], double tmax)
{
    double t0 = r->tmin, t1 = tmax;
    for (int a = 0; a < 3; a++) {
        double o = wr_v3_get(r->o, a), lo = wr_v3_get(b->lo, a), hi = wr_v3_get(b->hi, a);
        if (zero[a]) {
            if (o < lo || o > hi) return 0;
            continue;
        }
        double iv = wr_v3_get(inv, a);
        double ta = (lo - o) * iv, tb = (hi - o) * iv;
        if (ta > tb) { double s = ta; ta = tb; tb = s; }
        if (ta > t0) t0 = ta;
        if (tb < t1) t1 = tb;
        if (t0 > t1) return 0;
    }
    return 1;
}

static int traverse(const wr_bvh *bvh, const wr_ray *ray, wr_hit *hit, int any)
{
    if (bvh->nnodes == 0) return 0;
    wr_ray r = *ray;
    int zero[3] = {r.d.x == 0.0, r.d.y == 0.0, r.d.z == 0.0};
    wr_vec3 inv = wr_v3(zero[0] ? 0.0 : 1.0 / r.d.x, zero[1] ? 0.0 : 1.0 / r.d.y, zero[2] ? 0.0 : 1.0 / r.d.z);
    uint32_t stack[STACK_SIZE];
    int sp = 0;
    uint32_t cur = 0;
    int found = 0;
    wr_hit best;
    for (;;) {
        const wr_bvh_node *nd = &bvh->nodes[cur];
        if (box_hit(&nd->box, &r, inv, zero, r.tmax)) {
            if (nd->n) {
                for (uint32_t i = 0; i < nd->n; i++) {
                    uint32_t tri = bvh->order[nd->a + i];
                    wr_triangle t = wr_mesh_triangle(bvh->mesh, tri);
                    wr_hit h;
                    if (wr_triangle_intersect(&t, &r, &h)) {
                        if (any) return 1;
                        best = h; best.prim = (int)tri;
                        found = 1;
                        r.tmax = h.t;                                   /* only closer hits from now on */
                    }
                }
            } else {
                uint32_t left = cur + 1, right = nd->a;
                if (wr_v3_get(r.d, (int)nd->axis) < 0.0) { uint32_t t = left; left = right; right = t; }
                stack[sp++] = right;                                    /* far child later */
                cur = left;                                             /* near child now */
                continue;
            }
        }
        if (sp == 0) break;
        cur = stack[--sp];
    }
    if (found && hit) *hit = best;
    return found;
}

int wr_bvh_intersect(const wr_bvh *bvh, const wr_ray *ray, wr_hit *hit) { return traverse(bvh, ray, hit, 0); }
int wr_bvh_occluded(const wr_bvh *bvh, const wr_ray *ray) { return traverse(bvh, ray, NULL, 1); }

int wr_mesh_intersect_brute(const wr_mesh *mesh, const wr_ray *ray, wr_hit *hit)
{
    wr_ray r = *ray;
    int found = 0;
    for (size_t i = 0; i < mesh->nt; i++) {
        wr_triangle t = wr_mesh_triangle(mesh, i);
        wr_hit h;
        if (wr_triangle_intersect(&t, &r, &h)) { *hit = h; hit->prim = (int)i; found = 1; r.tmax = h.t; }
    }
    return found;
}
