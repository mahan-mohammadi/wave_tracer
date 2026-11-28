#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wr_mesh.h"
#include "wr_optics.h"

wr_mesh wr_mesh_empty(void)
{
    wr_mesh m = {NULL, 0, NULL, 0};
    return m;
}

void wr_mesh_free(wr_mesh *m)
{
    free(m->v);
    free(m->idx);
    *m = wr_mesh_empty();
}

wr_triangle wr_mesh_triangle(const wr_mesh *m, size_t i)
{
    return wr_triangle_make(m->v[m->idx[3 * i]], m->v[m->idx[3 * i + 1]], m->v[m->idx[3 * i + 2]]);
}

wr_aabb wr_mesh_bounds(const wr_mesh *m)
{
    wr_aabb b = wr_aabb_empty();
    for (size_t i = 0; i < m->nv; i++) b = wr_aabb_expand(b, m->v[i]);
    return b;
}

double wr_mesh_area(const wr_mesh *m)
{
    double a = 0.0;
    for (size_t i = 0; i < m->nt; i++) { wr_triangle t = wr_mesh_triangle(m, i); a += wr_triangle_area(&t); }
    return a;
}

double wr_mesh_volume(const wr_mesh *m)
{
    double s = 0.0;
    for (size_t i = 0; i < m->nt; i++) {
        wr_triangle t = wr_mesh_triangle(m, i);
        s += wr_v3_dot(t.v0, wr_v3_cross(t.v1, t.v2));
    }
    return s / 6.0;
}

int wr_mesh_transform(wr_mesh *m, const wr_xform *x)
{
    wr_mat3 inv;
    if (!wr_mat3_inverse(&x->m, &inv)) return 0;
    for (size_t i = 0; i < m->nv; i++) m->v[i] = wr_xform_point(x, m->v[i]);
    if (wr_mat3_det(&x->m) < 0.0)
        for (size_t i = 0; i < m->nt; i++) { uint32_t t = m->idx[3 * i + 1]; m->idx[3 * i + 1] = m->idx[3 * i + 2]; m->idx[3 * i + 2] = t; }
    return 1;
}

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : (x > y);
}

int wr_mesh_is_closed_oriented(const wr_mesh *m)
{
    if (m->nt == 0) return 0;
    size_t ne = 3 * m->nt;
    uint64_t *k = malloc(ne * sizeof *k);
    if (!k) return 0;
    for (size_t t = 0; t < m->nt; t++)
        for (int e = 0; e < 3; e++)
            k[3 * t + e] = (uint64_t)m->idx[3 * t + e] * m->nv + m->idx[3 * t + (e + 1) % 3];
    qsort(k, ne, sizeof *k, cmp_u64);
    int ok = 1;
    for (size_t i = 1; i < ne; i++) if (k[i] == k[i - 1]) ok = 0;
    for (size_t t = 0; t < m->nt && ok; t++)
        for (int e = 0; e < 3; e++) {
            uint64_t rev = (uint64_t)m->idx[3 * t + (e + 1) % 3] * m->nv + m->idx[3 * t + e];
            if (!bsearch(&rev, k, ne, sizeof *k, cmp_u64)) { ok = 0; break; }
        }
    free(k);
    return ok;
}

static int fail(char *err, size_t errlen, const char *fmt, long a)
{
    if (err && errlen) snprintf(err, errlen, fmt, a);
    return 0;
}

/* ---- OBJ parsing ---- */

typedef struct {
    wr_vec3  *v; size_t nv, capv;
    uint32_t *idx; size_t ni, capi;
} obj_builder;

static int push_v(obj_builder *b, wr_vec3 p)
{
    if (b->nv == b->capv) {
        size_t c = b->capv ? 2 * b->capv : 64;
        wr_vec3 *n = realloc(b->v, c * sizeof *n);
        if (!n) return 0;
        b->v = n; b->capv = c;
    }
    b->v[b->nv++] = p;
    return 1;
}
static int push_i(obj_builder *b, uint32_t i)
{
    if (b->ni == b->capi) {
        size_t c = b->capi ? 2 * b->capi : 192;
        uint32_t *n = realloc(b->idx, c * sizeof *n);
        if (!n) return 0;
        b->idx = n; b->capi = c;
    }
    b->idx[b->ni++] = i;
    return 1;
}

int wr_mesh_parse_obj(const char *text, wr_mesh *out, char *err, size_t errlen)
{
    *out = wr_mesh_empty();
    obj_builder b = {NULL, 0, 0, NULL, 0, 0};
    long lineno = 0;
    const char *p = text;
    int ok = 1;
    while (ok && *p) {
        const char *eol = strchr(p, '\n');
        size_t len = eol ? (size_t)(eol - p) : strlen(p);
        lineno++;
        char *line = malloc(len + 1);
        if (!line) { ok = fail(err, errlen, "out of memory%ld", 0); break; }
        memcpy(line, p, len); line[len] = '\0';
        p = eol ? eol + 1 : p + len;

        char *s = line;
        while (*s == ' ' || *s == '\t' || *s == '\r') s++;
        if (s[0] == 'v' && (s[1] == ' ' || s[1] == '\t')) {
            char *e;
            double x = strtod(s + 1, &e); char *e1 = e;
            double y = strtod(e, &e); char *e2 = e;
            double z = strtod(e, &e);
            if (e1 == s + 1 || e2 == e1 || e == e2 || !isfinite(x) || !isfinite(y) || !isfinite(z)) {
                ok = fail(err, errlen, "bad vertex on line %ld", lineno);
            } else if (!push_v(&b, wr_v3(x, y, z))) {
                ok = fail(err, errlen, "out of memory%ld", 0);
            }
        } else if (s[0] == 'f' && (s[1] == ' ' || s[1] == '\t')) {
            long first = 0, prev = 0;
            int count = 0;
            char *q = s + 1;
            for (;;) {
                while (*q == ' ' || *q == '\t' || *q == '\r') q++;
                if (!*q) break;
                char *e;
                long k = strtol(q, &e, 10);
                if (e == q || k == 0) { ok = fail(err, errlen, "bad face index on line %ld", lineno); break; }
                if (*e && *e != '/' && *e != ' ' && *e != '\t' && *e != '\r') { ok = fail(err, errlen, "bad face token on line %ld", lineno); break; }
                while (*e && *e != ' ' && *e != '\t' && *e != '\r') e++;   /* skip /vt/vn */
                q = e;
                long vi = k > 0 ? k - 1 : (long)b.nv + k;                 /* negative: relative to vertices so far */
                if (vi < 0 || (size_t)vi > 0xFFFFFFF0u) { ok = fail(err, errlen, "face index out of range on line %ld", lineno); break; }
                if (count == 0) first = vi;
                else if (count >= 2) {                                    /* fan: (first, prev, current) */
                    if (!push_i(&b, (uint32_t)first) || !push_i(&b, (uint32_t)prev) || !push_i(&b, (uint32_t)vi)) {
                        ok = fail(err, errlen, "out of memory%ld", 0); break;
                    }
                }
                prev = vi;
                count++;
            }
            if (ok && count < 3) ok = fail(err, errlen, "face with fewer than 3 vertices on line %ld", lineno);
        }
        free(line);
    }
    if (ok) {                                                              /* forward references resolved: validate now */
        for (size_t i = 0; i < b.ni; i++)
            if (b.idx[i] >= b.nv) { ok = fail(err, errlen, "face index beyond vertex count (triangle %ld)", (long)(i / 3)); break; }
    }
    if (!ok) { free(b.v); free(b.idx); return 0; }
    out->v = b.v; out->nv = b.nv; out->idx = b.idx; out->nt = b.ni / 3;
    return 1;
}

int wr_mesh_load_obj(const char *path, wr_mesh *out, char *err, size_t errlen)
{
    *out = wr_mesh_empty();
    FILE *f = fopen(path, "rb");
    if (!f) { if (err && errlen) snprintf(err, errlen, "cannot open %s", path); return 0; }
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *buf = n >= 0 ? malloc((size_t)n + 1) : NULL;
    if (!buf) { fclose(f); if (err && errlen) snprintf(err, errlen, "cannot read %s", path); return 0; }
    size_t got = fread(buf, 1, (size_t)n, f);
    fclose(f);
    buf[got] = '\0';
    int ok = wr_mesh_parse_obj(buf, out, err, errlen);
    free(buf);
    return ok;
}

/* ---- generators ---- */

static int alloc_mesh(wr_mesh *m, size_t nv, size_t nt)
{
    m->v = malloc(nv * sizeof *m->v);
    m->idx = malloc(3 * nt * sizeof *m->idx);
    m->nv = nv; m->nt = nt;
    if (!m->v || !m->idx) { wr_mesh_free(m); return 0; }
    return 1;
}

int wr_mesh_make_box(wr_vec3 lo, wr_vec3 hi, wr_mesh *out)
{
    *out = wr_mesh_empty();
    if (!alloc_mesh(out, 8, 12)) return 0;
    for (int i = 0; i < 8; i++)
        out->v[i] = wr_v3((i & 1) ? hi.x : lo.x, (i & 2) ? hi.y : lo.y, (i & 4) ? hi.z : lo.z);
    /* vertex index = x + 2y + 4z bits; faces listed as quads, clockwise seen from outside */
    static const uint32_t q[6][4] = {
        {0, 2, 6, 4},   /* -x */  {1, 5, 7, 3},   /* +x */
        {0, 4, 5, 1},   /* -y */  {2, 3, 7, 6},   /* +y */
        {0, 1, 3, 2},   /* -z */  {4, 6, 7, 5}};  /* +z */
    for (int f = 0; f < 6; f++) {
        uint32_t *t = out->idx + 6 * f;
        t[0] = q[f][0]; t[1] = q[f][2]; t[2] = q[f][1];          /* quads below are listed CW from outside: reverse */
        t[3] = q[f][0]; t[4] = q[f][3]; t[5] = q[f][2];
    }
    return 1;
}

int wr_mesh_make_uv_sphere(wr_vec3 c, double r, int nseg, int nstack, wr_mesh *out)
{
    *out = wr_mesh_empty();
    if (nseg < 3 || nstack < 2 || !(r > 0.0)) return 0;
    size_t nv = 2 + (size_t)(nstack - 1) * nseg;
    size_t nt = 2 * (size_t)nseg + 2 * (size_t)(nstack - 2) * nseg;
    if (!alloc_mesh(out, nv, nt)) return 0;
    out->v[0] = wr_v3_add(c, wr_v3(0, r, 0));                       /* north pole, theta = 0 */
    out->v[1] = wr_v3_add(c, wr_v3(0, -r, 0));                      /* south pole, theta = pi */
    for (int i = 1; i < nstack; i++) {
        double th = WR_PI * i / nstack;
        for (int j = 0; j < nseg; j++) {
            double ph = WR_TWO_PI * j / nseg;
            out->v[2 + (size_t)(i - 1) * nseg + j] =
                wr_v3(c.x + r * sin(th) * cos(ph), c.y + r * cos(th), c.z + r * sin(th) * sin(ph));
        }
    }
#define RING(i, j) ((uint32_t)(2 + (size_t)((i) - 1) * nseg + ((j) % nseg)))
    uint32_t *t = out->idx;
    for (int j = 0; j < nseg; j++) { *t++ = 0; *t++ = RING(1, j + 1); *t++ = RING(1, j); }
    for (int i = 1; i < nstack - 1; i++)
        for (int j = 0; j < nseg; j++) {
            uint32_t a = RING(i, j), b = RING(i, j + 1), cc = RING(i + 1, j), d = RING(i + 1, j + 1);
            *t++ = a; *t++ = b; *t++ = cc;
            *t++ = b; *t++ = d; *t++ = cc;
        }
    for (int j = 0; j < nseg; j++) { *t++ = RING(nstack - 1, j); *t++ = RING(nstack - 1, j + 1); *t++ = 1; }
#undef RING
    return 1;
}
