#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wr_sceneio.h"
#include "wr_optics.h"

#define MAXTOK 80
#define MAX_OPS 16

typedef struct { char *tok[MAXTOK]; int n; } line_t;
typedef struct { int kind; double v[4]; } mesh_op;      /* 0 translate, 1 scale, 2 rotate (v[3] = degrees) */

typedef struct {
    char *err; size_t errlen; int lineno;
} ctx_t;

static int fail(ctx_t *c, const char *fmt, ...)
{
    if (c->err && c->errlen) {
        int n = snprintf(c->err, c->errlen, "line %d: ", c->lineno);
        va_list ap;
        va_start(ap, fmt);
        if (n > 0 && (size_t)n < c->errlen) vsnprintf(c->err + n, c->errlen - n, fmt, ap);
        va_end(ap);
    }
    return 0;
}

static int parse_num(const char *s, double *out)
{
    char *e;
    double v = strtod(s, &e);
    if (e == s || *e != '\0' || isnan(v)) return 0;
    *out = v;
    return 1;
}

/* read `count` numbers following the keyword at *i; advances *i past them */
static int take(ctx_t *c, const line_t *L, int *i, int count, double *out)
{
    const char *kw = L->tok[*i];
    if (*i + count >= L->n) return fail(c, "'%s' needs %d number(s)", kw, count);
    for (int k = 0; k < count; k++)
        if (!parse_num(L->tok[*i + 1 + k], &out[k])) return fail(c, "'%s': bad number '%s'", kw, L->tok[*i + 1 + k]);
    *i += count + 1;
    return 1;
}

static int take_vec(ctx_t *c, const line_t *L, int *i, wr_vec3 *v)
{
    double d[3];
    if (!take(c, L, i, 3, d)) return 0;
    *v = wr_v3(d[0], d[1], d[2]);
    return 1;
}

static int valid_name(ctx_t *c, const char *name)
{
    if (strlen(name) >= WR_NAME_LEN) return fail(c, "name '%s' is too long (max %d)", name, WR_NAME_LEN - 1);
    return 1;
}

static int need_material(ctx_t *c, const wr_scene *s, const char *name)
{
    int m = wr_scene_find_material(s, name);
    if (m < 0) fail(c, "unknown material '%s'", name);
    return m;
}

/* ---- directives ---- */

static int do_material(ctx_t *c, wr_scene *s, const line_t *L)
{
    if (L->n < 3) return fail(c, "material needs a name and a model");
    const char *name = L->tok[1], *model = L->tok[2];
    if (!valid_name(c, name)) return 0;
    wr_material m;
    int i = 3;
    double d[6];
    if (strcmp(model, "constant") == 0) {
        if (!take(c, L, &(int){2}, 1, d)) return 0;
        i = 4; m = wr_material_constant(name, d[0], 0.0);
    } else if (strcmp(model, "cauchy") == 0) {
        if (!take(c, L, &(int){2}, 3, d)) return 0;
        i = 6; m = wr_material_cauchy(name, d[0], d[1], d[2]);
    } else if (strcmp(model, "sellmeier") == 0) {
        if (!take(c, L, &(int){2}, 6, d)) return 0;
        i = 9; m = wr_material_sellmeier(name, d, d + 3);
    } else if (strcmp(model, "ideal_mirror") == 0) {
        m = wr_material_mirror(name);
    } else if (strcmp(model, "absorber") == 0) {
        m = wr_material_absorber(name);
    } else {
        return fail(c, "unknown material model '%s'", model);
    }
    while (i < L->n) {
        const char *kw = L->tok[i];
        if (strcmp(kw, "extinction") == 0) { if (!take(c, L, &i, 1, d)) return 0; m.kappa = d[0]; }
        else if (strcmp(kw, "range") == 0) { if (!take(c, L, &i, 2, d)) return 0; m.lambda_min = d[0]; m.lambda_max = d[1]; }
        else return fail(c, "unexpected '%s' in material", kw);
    }
    char e[160];
    if (wr_scene_add_material(s, &m, e, sizeof e) < 0) return fail(c, "%s", e);
    return 1;
}

static int do_sphere_plane(ctx_t *c, wr_scene *s, const line_t *L, int is_sphere)
{
    if (L->n < 2) return fail(c, "%s needs a name", L->tok[0]);
    const char *name = L->tok[1];
    if (!valid_name(c, name)) return 0;
    int mat = -1, have_a = 0, have_b = 0;
    wr_vec3 a = wr_v3_zero(), nrm = wr_v3_zero(), vel = wr_v3_zero();
    double r = 0, d[1];
    for (int i = 2; i < L->n;) {
        const char *kw = L->tok[i];
        if (strcmp(kw, "material") == 0) {
            if (i + 1 >= L->n) return fail(c, "'material' needs a name");
            mat = need_material(c, s, L->tok[i + 1]);
            if (mat < 0) return 0;
            i += 2;
        }
        else if (is_sphere && strcmp(kw, "center") == 0) { if (!take_vec(c, L, &i, &a)) return 0; have_a = 1; }
        else if (is_sphere && strcmp(kw, "radius") == 0) { if (!take(c, L, &i, 1, d)) return 0; r = d[0]; have_b = 1; }
        else if (!is_sphere && strcmp(kw, "point") == 0) { if (!take_vec(c, L, &i, &a)) return 0; have_a = 1; }
        else if (!is_sphere && strcmp(kw, "normal") == 0) { if (!take_vec(c, L, &i, &nrm)) return 0; have_b = 1; }
        else if (strcmp(kw, "velocity") == 0) { if (!take_vec(c, L, &i, &vel)) return 0; }
        else return fail(c, "unexpected '%s' in %s", kw, L->tok[0]);
    }
    if (mat < 0) return fail(c, "%s '%s' needs 'material'", L->tok[0], name);
    if (!have_a || !have_b) return fail(c, "%s '%s' needs %s", L->tok[0], name, is_sphere ? "center and radius" : "point and normal");
    char e[160];
    int id = is_sphere ? wr_scene_add_sphere(s, name, mat, a, r, e, sizeof e) : wr_scene_add_plane(s, name, mat, a, nrm, e, sizeof e);
    if (id < 0) return fail(c, "%s", e);
    s->objects[id].velocity = vel;
    return 1;
}

static int do_mesh(ctx_t *c, wr_scene *s, const line_t *L, const char *base_dir)
{
    if (L->n < 2) return fail(c, "mesh needs a name");
    const char *name = L->tok[1];
    if (!valid_name(c, name)) return 0;
    int mat = -1, have_src = 0;
    wr_mesh mesh = wr_mesh_empty();
    mesh_op ops[MAX_OPS];
    int nops = 0;
    wr_vec3 vel = wr_v3_zero();
    double d[6];
    char e[200];
    for (int i = 2; i < L->n;) {
        const char *kw = L->tok[i];
        if (strcmp(kw, "material") == 0) {
            if (i + 1 >= L->n) { wr_mesh_free(&mesh); return fail(c, "'material' needs a name"); }
            mat = need_material(c, s, L->tok[i + 1]);
            if (mat < 0) { wr_mesh_free(&mesh); return 0; }
            i += 2;
        } else if (strcmp(kw, "obj") == 0 || strcmp(kw, "box") == 0 || strcmp(kw, "uvsphere") == 0) {
            if (have_src) { wr_mesh_free(&mesh); return fail(c, "mesh '%s' has more than one source", name); }
            have_src = 1;
            int ok;
            if (strcmp(kw, "obj") == 0) {
                if (i + 1 >= L->n) return fail(c, "'obj' needs a path");
                char path[1024];
                const char *rel = L->tok[i + 1];
                if (base_dir && base_dir[0] && rel[0] != '/') snprintf(path, sizeof path, "%s/%s", base_dir, rel);
                else snprintf(path, sizeof path, "%s", rel);
                ok = wr_mesh_load_obj(path, &mesh, e, sizeof e);
                if (!ok) return fail(c, "obj '%s': %s", rel, e);
                i += 2;
            } else if (strcmp(kw, "box") == 0) {
                if (!take(c, L, &i, 6, d)) return 0;
                if (!(d[3] > d[0] && d[4] > d[1] && d[5] > d[2])) return fail(c, "box needs lo < hi on every axis");
                if (!wr_mesh_make_box(wr_v3(d[0], d[1], d[2]), wr_v3(d[3], d[4], d[5]), &mesh)) return fail(c, "box generation failed");
            } else {
                if (!take(c, L, &i, 6, d)) return 0;
                if (d[4] != floor(d[4]) || d[5] != floor(d[5]) || d[4] > 4096 || d[5] > 4096)
                    return fail(c, "uvsphere segments and stacks must be integers <= 4096");
                if (!wr_mesh_make_uv_sphere(wr_v3(d[0], d[1], d[2]), d[3], (int)d[4], (int)d[5], &mesh))
                    return fail(c, "uvsphere needs radius > 0, segments >= 3, stacks >= 2");
            }
        } else if (strcmp(kw, "translate") == 0 || strcmp(kw, "scale") == 0 || strcmp(kw, "rotate") == 0) {
            int kind = kw[0] == 't' ? 0 : (kw[1] == 'c' ? 1 : 2);
            int cnt = kind == 2 ? 4 : 3;
            if (nops >= MAX_OPS) { wr_mesh_free(&mesh); return fail(c, "too many transforms (max %d)", MAX_OPS); }
            if (!take(c, L, &i, cnt, d)) { wr_mesh_free(&mesh); return 0; }
            ops[nops].kind = kind;
            for (int k = 0; k < 4; k++) ops[nops].v[k] = k < cnt ? d[k] : 0.0;
            nops++;
        } else if (strcmp(kw, "velocity") == 0) {
            if (!take_vec(c, L, &i, &vel)) { wr_mesh_free(&mesh); return 0; }
        } else {
            wr_mesh_free(&mesh);
            return fail(c, "unexpected '%s' in mesh", kw);
        }
    }
    if (mat < 0) { wr_mesh_free(&mesh); return fail(c, "mesh '%s' needs 'material'", name); }
    if (!have_src) return fail(c, "mesh '%s' needs a source (obj, box or uvsphere)", name);
    for (int k = 0; k < nops; k++) {
        wr_xform x;
        if (ops[k].kind == 0) x = wr_xform_translate(wr_v3(ops[k].v[0], ops[k].v[1], ops[k].v[2]));
        else if (ops[k].kind == 1) x = wr_xform_scale(wr_v3(ops[k].v[0], ops[k].v[1], ops[k].v[2]));
        else {
            if (wr_v3_len2(wr_v3(ops[k].v[0], ops[k].v[1], ops[k].v[2])) == 0.0) { wr_mesh_free(&mesh); return fail(c, "rotate axis must be non-zero"); }
            x = wr_xform_rotate(wr_v3(ops[k].v[0], ops[k].v[1], ops[k].v[2]), ops[k].v[3] * WR_PI / 180.0);
        }
        if (!wr_mesh_transform(&mesh, &x)) { wr_mesh_free(&mesh); return fail(c, "singular transform on mesh '%s'", name); }
    }
    int id = wr_scene_add_mesh(s, name, mat, &mesh, e, sizeof e);
    if (id < 0) { wr_mesh_free(&mesh); return fail(c, "%s", e); }
    s->objects[id].velocity = vel;
    return 1;
}

static int do_light(ctx_t *c, wr_scene *s, const line_t *L)
{
    if (L->n < 3) return fail(c, "light needs a type and a name");
    const char *type = L->tok[1], *name = L->tok[2];
    wr_light_type lt;
    if (strcmp(type, "point") == 0) lt = WR_LIGHT_POINT;
    else if (strcmp(type, "plane_wave") == 0) lt = WR_LIGHT_PLANE_WAVE;
    else if (strcmp(type, "gaussian") == 0) lt = WR_LIGHT_GAUSSIAN;
    else return fail(c, "unknown light type '%s'", type);
    if (!valid_name(c, name)) return 0;

    wr_vec3 pos = wr_v3_zero(), dir = wr_v3_zero(), pol_ref = wr_v3(0, 1, 0), vel = wr_v3_zero();
    int have_dir = 0, have_waist = 0, have_spec = 0;
    double waist = 0, amp = 1, phase = 0, coh = INFINITY, d[4];
    wr_cplx jones[2] = {wr_c(1, 0), wr_c(0, 0)};
    wr_spectrum sp;
    memset(&sp, 0, sizeof sp);
    for (int i = 3; i < L->n;) {
        const char *kw = L->tok[i];
        if (strcmp(kw, "position") == 0) { if (!take_vec(c, L, &i, &pos)) return 0; }
        else if (strcmp(kw, "direction") == 0) { if (!take_vec(c, L, &i, &dir)) return 0; have_dir = 1; }
        else if (strcmp(kw, "waist") == 0) { if (!take(c, L, &i, 1, d)) return 0; waist = d[0]; have_waist = 1; }
        else if (strcmp(kw, "wavelength") == 0) {
            if (!take(c, L, &i, 1, d)) return 0;
            sp.n = 1; sp.lambda[0] = d[0]; sp.weight[0] = 1.0; have_spec = 1;
        } else if (strcmp(kw, "lines") == 0) {
            double nn;
            if (i + 1 >= L->n || !parse_num(L->tok[i + 1], &nn) || nn != floor(nn) || nn < 1 || nn > WR_MAX_LINES)
                return fail(c, "'lines' needs a count between 1 and %d", WR_MAX_LINES);
            int cnt = (int)nn;
            if (i + 1 + 2 * cnt >= L->n) return fail(c, "'lines %d' needs %d wavelength/weight pairs", cnt, cnt);
            for (int k = 0; k < cnt; k++) {
                if (!parse_num(L->tok[i + 2 + 2 * k], &sp.lambda[k]) || !parse_num(L->tok[i + 3 + 2 * k], &sp.weight[k]))
                    return fail(c, "'lines': bad number in pair %d", k + 1);
            }
            sp.n = cnt; have_spec = 1;
            i += 2 + 2 * cnt;
        }
        else if (strcmp(kw, "amplitude") == 0) { if (!take(c, L, &i, 1, d)) return 0; amp = d[0]; }
        else if (strcmp(kw, "phase") == 0) { if (!take(c, L, &i, 1, d)) return 0; phase = d[0]; }
        else if (strcmp(kw, "coherence_length") == 0) { if (!take(c, L, &i, 1, d)) return 0; coh = d[0]; }
        else if (strcmp(kw, "jones") == 0) {
            if (!take(c, L, &i, 4, d)) return 0;
            double nrm = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + d[3] * d[3]);
            if (!(nrm > 0.0) || !isfinite(nrm)) return fail(c, "jones vector must be non-zero and finite");
            jones[0] = wr_c(d[0] / nrm, d[1] / nrm); jones[1] = wr_c(d[2] / nrm, d[3] / nrm);
        }
        else if (strcmp(kw, "pol_ref") == 0) { if (!take_vec(c, L, &i, &pol_ref)) return 0; }
        else if (strcmp(kw, "velocity") == 0) { if (!take_vec(c, L, &i, &vel)) return 0; }
        else return fail(c, "unexpected '%s' in light", kw);
    }
    if (!have_spec) return fail(c, "light '%s' needs 'wavelength' or 'lines'", name);
    if (lt != WR_LIGHT_POINT) {
        if (!have_dir) return fail(c, "light '%s' needs 'direction'", name);
        if (!(wr_v3_len2(dir) > 0.0) || !isfinite(wr_v3_len2(dir))) return fail(c, "light '%s': direction must be non-zero", name);
    }
    if (lt == WR_LIGHT_GAUSSIAN && !have_waist) return fail(c, "light '%s' needs 'waist'", name);
    wr_light l = lt == WR_LIGHT_POINT ? wr_light_point(name, pos, sp.lambda[0], amp)
               : lt == WR_LIGHT_PLANE_WAVE ? wr_light_plane_wave(name, pos, dir, sp.lambda[0], amp)
               : wr_light_gaussian(name, pos, dir, waist, sp.lambda[0], amp);
    l.spectrum = sp;
    l.phase = phase; l.coherence_length = coh; l.jones[0] = jones[0]; l.jones[1] = jones[1];
    l.pol_ref = pol_ref; l.velocity = vel;
    char e[160];
    if (wr_scene_add_light(s, &l, e, sizeof e) < 0) return fail(c, "%s", e);
    return 1;
}

static int do_camera(ctx_t *c, wr_scene *s, const line_t *L)
{
    if (L->n < 2) return fail(c, "camera needs a type");
    const char *type = L->tok[1];
    int kind = strcmp(type, "pinhole") == 0 ? 0 : strcmp(type, "orthographic") == 0 ? 1 : strcmp(type, "thin_lens") == 0 ? 2 : -1;
    if (kind < 0) return fail(c, "unknown camera type '%s'", type);
    wr_vec3 pos = wr_v3_zero(), look = wr_v3_zero(), up = wr_v3(0, 1, 0);
    int have_pos = 0, have_look = 0, have_size = 0, have_fov = 0, have_vh = 0, have_ap = 0, have_foc = 0;
    double vfov = 0, vh = 0, ap = 0, foc = 0, d[2];
    int W = 0, H = 0;
    for (int i = 2; i < L->n;) {
        const char *kw = L->tok[i];
        if (strcmp(kw, "position") == 0) { if (!take_vec(c, L, &i, &pos)) return 0; have_pos = 1; }
        else if (strcmp(kw, "look_at") == 0) { if (!take_vec(c, L, &i, &look)) return 0; have_look = 1; }
        else if (strcmp(kw, "up") == 0) { if (!take_vec(c, L, &i, &up)) return 0; }
        else if (strcmp(kw, "size") == 0) {
            if (!take(c, L, &i, 2, d)) return 0;
            if (d[0] != floor(d[0]) || d[1] != floor(d[1]) || d[0] < 1 || d[1] < 1 || d[0] > 65536 || d[1] > 65536) return fail(c, "size needs two integers in 1..65536");
            W = (int)d[0]; H = (int)d[1]; have_size = 1;
        }
        else if (strcmp(kw, "vfov_deg") == 0) { if (!take(c, L, &i, 1, d)) return 0; vfov = d[0] * WR_PI / 180.0; have_fov = 1; }
        else if (strcmp(kw, "view_height") == 0) { if (!take(c, L, &i, 1, d)) return 0; vh = d[0]; have_vh = 1; }
        else if (strcmp(kw, "aperture") == 0) { if (!take(c, L, &i, 1, d)) return 0; ap = d[0]; have_ap = 1; }
        else if (strcmp(kw, "focus") == 0) { if (!take(c, L, &i, 1, d)) return 0; foc = d[0]; have_foc = 1; }
        else return fail(c, "unexpected '%s' in camera", kw);
    }
    if (!have_pos || !have_look || !have_size) return fail(c, "camera needs position, look_at and size");
    if (kind != 1 && !have_fov) return fail(c, "%s camera needs vfov_deg", type);
    if (kind == 1 && !have_vh) return fail(c, "orthographic camera needs view_height");
    if (kind == 2 && !(have_ap && have_foc)) return fail(c, "thin_lens camera needs aperture and focus");
    wr_camera cam;
    int ok = kind == 0 ? wr_camera_pinhole(&cam, pos, look, up, vfov, W, H)
           : kind == 1 ? wr_camera_orthographic(&cam, pos, look, up, vh, W, H)
           : wr_camera_thin_lens(&cam, pos, look, up, vfov, ap, foc, W, H);
    if (!ok) return fail(c, "invalid camera (check up vs. view direction, field of view in (0, 180) deg, positive sizes)");
    wr_scene_set_camera(s, &cam);
    return 1;
}

int wr_scene_parse(const char *text, const char *base_dir, wr_scene *scene, char *err, size_t errlen)
{
    if (!wr_scene_init(scene)) { if (err && errlen) snprintf(err, errlen, "out of memory"); return 0; }
    ctx_t c = {err, errlen, 0};
    if (err && errlen) err[0] = '\0';
    int ok = 1;
    const char *p = text;
    while (ok && *p) {
        const char *eol = strchr(p, '\n');
        size_t len = eol ? (size_t)(eol - p) : strlen(p);
        c.lineno++;
        char *buf = malloc(len + 1);
        if (!buf) { ok = fail(&c, "out of memory"); break; }
        memcpy(buf, p, len); buf[len] = '\0';
        p = eol ? eol + 1 : p + len;
        char *hash = strchr(buf, '#');
        if (hash) *hash = '\0';
        line_t L; L.n = 0;
        for (char *t = strtok(buf, " \t\r"); t; t = strtok(NULL, " \t\r")) {
            if (L.n == MAXTOK) { ok = fail(&c, "too many tokens on one line (max %d)", MAXTOK); break; }
            L.tok[L.n++] = t;
        }
        if (ok && L.n > 0) {
            const char *d = L.tok[0];
            double v;
            if (strcmp(d, "time") == 0) {
                if (L.n != 2 || !parse_num(L.tok[1], &v) || !isfinite(v)) ok = fail(&c, "time needs one finite number");
                else scene->time = v;
            } else if (strcmp(d, "ambient") == 0) {
                if (L.n != 2) ok = fail(&c, "ambient needs one material name");
                else {
                    int m = need_material(&c, scene, L.tok[1]);
                    char e[160];
                    if (m < 0) ok = 0;
                    else if (wr_scene_set_ambient(scene, m, e, sizeof e) < 0) ok = fail(&c, "%s", e);
                }
            }
            else if (strcmp(d, "material") == 0) ok = do_material(&c, scene, &L);
            else if (strcmp(d, "sphere") == 0) ok = do_sphere_plane(&c, scene, &L, 1);
            else if (strcmp(d, "plane") == 0) ok = do_sphere_plane(&c, scene, &L, 0);
            else if (strcmp(d, "mesh") == 0) ok = do_mesh(&c, scene, &L, base_dir);
            else if (strcmp(d, "light") == 0) ok = do_light(&c, scene, &L);
            else if (strcmp(d, "camera") == 0) ok = do_camera(&c, scene, &L);
            else ok = fail(&c, "unknown directive '%s'", d);
        }
        free(buf);
    }
    if (ok) {
        char e[200];
        if (!wr_scene_validate(scene, e, sizeof e)) {
            if (err && errlen) snprintf(err, errlen, "scene: %s", e);
            ok = 0;
        }
    }
    if (!ok) wr_scene_free(scene);
    return ok;
}

int wr_scene_load(const char *path, wr_scene *scene, char *err, size_t errlen)
{
    memset(scene, 0, sizeof *scene);
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
    char dir[1024];
    snprintf(dir, sizeof dir, "%s", path);
    char *slash = strrchr(dir, '/');
    if (slash) *slash = '\0'; else strcpy(dir, ".");
    int ok = wr_scene_parse(buf, dir, scene, err, errlen);
    free(buf);
    return ok;
}
