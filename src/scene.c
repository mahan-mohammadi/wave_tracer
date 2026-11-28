#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wr_scene.h"

static int fail(char *err, size_t n, const char *fmt, const char *a)
{
    if (err && n) snprintf(err, n, fmt, a ? a : "");
    return -1;
}

int wr_scene_init(wr_scene *s)
{
    memset(s, 0, sizeof *s);
    wr_material vac = wr_material_constant("vacuum", 1.0, 0.0);
    s->materials = malloc(4 * sizeof *s->materials);
    if (!s->materials) return 0;
    s->capmat = 4;
    s->materials[0] = vac;
    s->nmat = 1;
    s->ambient = 0;
    return 1;
}

void wr_scene_free(wr_scene *s)
{
    for (int i = 0; i < s->nobj; i++)
        if (s->objects[i].kind == WR_OBJ_MESH) {
            wr_bvh_free(&s->objects[i].bvh);
            wr_mesh_free(s->objects[i].mesh);
            free(s->objects[i].mesh);
        }
    free(s->materials); free(s->objects); free(s->lights);
    memset(s, 0, sizeof *s);
}

/* Ensure room for one more element; returns 0 on failure. */
static int grow(void **arr, int n, int *cap, size_t size)
{
    if (n < *cap) return 1;
    int nc = *cap ? 2 * *cap : 8;
    void *p = realloc(*arr, (size_t)nc * size);
    if (!p) return 0;
    *arr = p; *cap = nc;
    return 1;
}

int wr_scene_find_material(const wr_scene *s, const char *name)
{
    for (int i = 0; i < s->nmat; i++) if (strcmp(s->materials[i].name, name) == 0) return i;
    return -1;
}
int wr_scene_find_object(const wr_scene *s, const char *name)
{
    for (int i = 0; i < s->nobj; i++) if (strcmp(s->objects[i].name, name) == 0) return i;
    return -1;
}
int wr_scene_find_light(const wr_scene *s, const char *name)
{
    for (int i = 0; i < s->nlight; i++) if (strcmp(s->lights[i].name, name) == 0) return i;
    return -1;
}

int wr_scene_add_material(wr_scene *s, const wr_material *m, char *err, size_t errlen)
{
    if (!wr_material_validate(m, err, errlen)) return -1;
    if (wr_scene_find_material(s, m->name) >= 0) return fail(err, errlen, "duplicate material name '%s'", m->name);
    if (!grow((void **)&s->materials, s->nmat, &s->capmat, sizeof *s->materials)) return fail(err, errlen, "out of memory%s", NULL);
    s->materials[s->nmat] = *m;
    return s->nmat++;
}

static int check_object(const wr_scene *s, const char *name, int material, char *err, size_t errlen)
{
    if (!name || !name[0]) return fail(err, errlen, "object needs a name%s", NULL);
    if (strlen(name) >= WR_NAME_LEN) return fail(err, errlen, "object name too long: '%s'", name);
    if (wr_scene_find_object(s, name) >= 0) return fail(err, errlen, "duplicate object name '%s'", name);
    if (material < 0 || material >= s->nmat) return fail(err, errlen, "object '%s': unknown material index", name);
    return 0;
}

static wr_object *push_object(wr_scene *s, const char *name, wr_obj_kind kind, int material, char *err, size_t errlen)
{
    if (!grow((void **)&s->objects, s->nobj, &s->capobj, sizeof *s->objects)) { fail(err, errlen, "out of memory%s", NULL); return NULL; }
    wr_object *o = &s->objects[s->nobj];
    memset(o, 0, sizeof *o);
    snprintf(o->name, sizeof o->name, "%s", name);
    o->kind = kind;
    o->material = material;
    return o;
}

int wr_scene_add_sphere(wr_scene *s, const char *name, int material, wr_vec3 c, double r, char *err, size_t errlen)
{
    if (check_object(s, name, material, err, errlen) < 0) return -1;
    if (!(r > 0.0) || !isfinite(r) || !isfinite(wr_v3_len2(c))) return fail(err, errlen, "sphere '%s': radius must be > 0 and center finite", name);
    wr_object *o = push_object(s, name, WR_OBJ_SPHERE, material, err, errlen);
    if (!o) return -1;
    o->sphere = wr_sphere_make(c, r);
    return s->nobj++;
}

int wr_scene_add_plane(wr_scene *s, const char *name, int material, wr_vec3 p, wr_vec3 n, char *err, size_t errlen)
{
    if (check_object(s, name, material, err, errlen) < 0) return -1;
    if (!(wr_v3_len2(n) > 0.0) || !isfinite(wr_v3_len2(n)) || !isfinite(wr_v3_len2(p))) return fail(err, errlen, "plane '%s': needs a finite point and a non-zero normal", name);
    wr_object *o = push_object(s, name, WR_OBJ_PLANE, material, err, errlen);
    if (!o) return -1;
    o->plane = wr_plane_make(p, n);
    return s->nobj++;
}

int wr_scene_add_mesh(wr_scene *s, const char *name, int material, wr_mesh *mesh, char *err, size_t errlen)
{
    if (check_object(s, name, material, err, errlen) < 0) return -1;
    if (mesh->nt == 0) return fail(err, errlen, "mesh '%s' has no triangles", name);
    wr_object *o = push_object(s, name, WR_OBJ_MESH, material, err, errlen);
    if (!o) return -1;
    o->mesh = malloc(sizeof *o->mesh);
    if (!o->mesh) return fail(err, errlen, "out of memory%s", NULL);
    *o->mesh = *mesh;
    if (!wr_bvh_build(&o->bvh, o->mesh)) { free(o->mesh); o->mesh = NULL; return fail(err, errlen, "BVH build failed for '%s'", name); }
    *mesh = wr_mesh_empty();                                  /* ownership moved */
    return s->nobj++;
}

int wr_scene_add_light(wr_scene *s, const wr_light *l, char *err, size_t errlen)
{
    if (!wr_light_validate(l, err, errlen)) return -1;
    if (wr_scene_find_light(s, l->name) >= 0) return fail(err, errlen, "duplicate light name '%s'", l->name);
    if (!grow((void **)&s->lights, s->nlight, &s->caplight, sizeof *s->lights)) return fail(err, errlen, "out of memory%s", NULL);
    s->lights[s->nlight] = *l;
    return s->nlight++;
}

int wr_scene_set_ambient(wr_scene *s, int material, char *err, size_t errlen)
{
    if (material < 0 || material >= s->nmat) return fail(err, errlen, "ambient: unknown material index%s", NULL);
    if (s->materials[material].surface != WR_SURF_DIELECTRIC) return fail(err, errlen, "ambient medium '%s' must be a dielectric", s->materials[material].name);
    s->ambient = material;
    return 0;
}

void wr_scene_set_camera(wr_scene *s, const wr_camera *c) { s->camera = *c; s->has_camera = 1; }

static int closest(const wr_object *o, const wr_ray *r, wr_hit *h)
{
    switch (o->kind) {
    case WR_OBJ_SPHERE: { int ok = wr_sphere_intersect(&o->sphere, r, h); if (ok) h->prim = -1; return ok; }
    case WR_OBJ_PLANE:  { int ok = wr_plane_intersect(&o->plane, r, h);   if (ok) h->prim = -1; return ok; }
    case WR_OBJ_MESH:   return wr_bvh_intersect(&o->bvh, r, h);
    }
    return 0;
}

int wr_scene_intersect(const wr_scene *s, const wr_ray *ray, wr_hit *hit)
{
    wr_ray r = *ray;
    int found = 0;
    for (int i = 0; i < s->nobj; i++) {
        wr_hit h;
        if (closest(&s->objects[i], &r, &h)) {
            *hit = h; hit->obj = i; found = 1;
            r.tmax = h.t;
        }
    }
    return found;
}

int wr_scene_occluded(const wr_scene *s, const wr_ray *ray)
{
    for (int i = 0; i < s->nobj; i++) {
        const wr_object *o = &s->objects[i];
        wr_hit h;
        switch (o->kind) {
        case WR_OBJ_SPHERE: if (wr_sphere_intersect(&o->sphere, ray, &h)) return 1; break;
        case WR_OBJ_PLANE:  if (wr_plane_intersect(&o->plane, ray, &h)) return 1; break;
        case WR_OBJ_MESH:   if (wr_bvh_occluded(&o->bvh, ray)) return 1; break;
        }
    }
    return 0;
}

int wr_scene_medium_before(const wr_scene *s, const wr_hit *hit)
{
    return hit->front ? s->ambient : s->objects[hit->obj].material;
}
int wr_scene_medium_after(const wr_scene *s, const wr_hit *hit)
{
    return hit->front ? s->objects[hit->obj].material : s->ambient;
}

wr_aabb wr_scene_bounds(const wr_scene *s)
{
    wr_aabb b = wr_aabb_empty();
    for (int i = 0; i < s->nobj; i++) {
        const wr_object *o = &s->objects[i];
        if (o->kind == WR_OBJ_SPHERE) { wr_aabb sb = wr_sphere_bounds(&o->sphere); b = wr_aabb_union(b, sb); }
        else if (o->kind == WR_OBJ_MESH) { wr_aabb mb = wr_mesh_bounds(o->mesh); b = wr_aabb_union(b, mb); }
    }
    return b;
}

static int bad2(char *err, size_t n, const char *fmt, const char *a, const char *b)
{
    if (err && n) snprintf(err, n, fmt, a, b);
    return 0;
}

int wr_scene_validate(const wr_scene *s, char *err, size_t errlen)
{
    if (s->ambient < 0 || s->ambient >= s->nmat) return bad2(err, errlen, "ambient material index invalid%s%s", "", "");
    if (s->materials[s->ambient].surface != WR_SURF_DIELECTRIC)
        return bad2(err, errlen, "ambient medium '%s' must be a dielectric%s", s->materials[s->ambient].name, "");
    for (int i = 0; i < s->nmat; i++) if (!wr_material_validate(&s->materials[i], err, errlen)) return 0;
    for (int i = 0; i < s->nlight; i++) if (!wr_light_validate(&s->lights[i], err, errlen)) return 0;
    for (int i = 0; i < s->nobj; i++) {
        const wr_object *o = &s->objects[i];
        if (o->material < 0 || o->material >= s->nmat) return bad2(err, errlen, "object '%s': unknown material%s", o->name, "");
        const wr_material *m = &s->materials[o->material];
        if (o->kind == WR_OBJ_MESH && m->surface == WR_SURF_DIELECTRIC && o->material != s->ambient) {
            if (!wr_mesh_is_closed_oriented(o->mesh))
                return bad2(err, errlen, "mesh '%s' is a solid of '%s' but is not a closed, consistently oriented surface", o->name, m->name);
            if (!(wr_mesh_volume(o->mesh) > 0.0))
                return bad2(err, errlen, "mesh '%s' (solid of '%s') has non-positive volume: triangles wound inward", o->name, m->name);
        }
    }
    return 1;
}
