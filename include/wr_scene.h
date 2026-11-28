/* wr_scene.h - scene container (Phase 4c).
 *
 * Owns materials, objects (sphere / plane / triangle mesh with BVH), lights
 * and one optional camera, plus the ambient medium and the simulation time.
 * Geometry queries only; no light transport.
 *
 * Media model: every object has one material = the medium INSIDE it. Outside
 * is the scene's ambient material. Nested or touching media are NOT supported
 * (a ray leaving an object always enters the ambient medium); wr_scene_validate
 * cannot detect overlapping objects, so scenes must keep solids disjoint.
 * wr_scene_medium_before/after encode exactly this rule, so the transport code
 * has one place to change when nested media arrive.
 *
 * A mesh with a dielectric material (other than the ambient) is a solid and
 * must be closed, consistently oriented and outward wound (positive volume);
 * wr_scene_validate checks it. Mirror/absorber meshes may be open sheets.
 * Plane + dielectric means a half-space (below the normal).
 *
 * Object velocity is stored for the time-dependent phase (Phase 8) and ignored
 * by every query here: objects do not move yet.
 *
 * Add functions return the new index (>= 0) or -1 with a message in err
 * (err may be NULL). wr_scene_add_mesh takes ownership of *mesh on success
 * only (the mesh is moved to the heap; the caller's struct is reset to empty).
 */
#ifndef WR_SCENE_H
#define WR_SCENE_H

#include "wr_material.h"
#include "wr_light.h"
#include "wr_camera.h"
#include "wr_mesh.h"
#include "wr_bvh.h"

typedef enum { WR_OBJ_SPHERE, WR_OBJ_PLANE, WR_OBJ_MESH } wr_obj_kind;

typedef struct {
    char        name[WR_NAME_LEN];
    wr_obj_kind kind;
    int         material;       /* index into scene materials: medium inside */
    wr_vec3     velocity;       /* m/s, stored only */
    wr_sphere   sphere;
    wr_plane    plane;
    wr_mesh    *mesh;           /* heap, owned (stable address: the BVH points at it) */
    wr_bvh      bvh;
} wr_object;

typedef struct {
    wr_material *materials; int nmat, capmat;
    wr_object   *objects;   int nobj, capobj;
    wr_light    *lights;    int nlight, caplight;
    wr_camera    camera;    int has_camera;
    int          ambient;   /* material index of the surrounding medium */
    double       time;      /* simulation time, s */
} wr_scene;

/* Creates the built-in material 0 "vacuum" (n = 1) and makes it the ambient medium. */
int  wr_scene_init(wr_scene *s);
void wr_scene_free(wr_scene *s);

int wr_scene_add_material(wr_scene *s, const wr_material *m, char *err, size_t errlen);
int wr_scene_add_sphere(wr_scene *s, const char *name, int material, wr_vec3 center, double radius, char *err, size_t errlen);
int wr_scene_add_plane(wr_scene *s, const char *name, int material, wr_vec3 point, wr_vec3 normal, char *err, size_t errlen);
int wr_scene_add_mesh(wr_scene *s, const char *name, int material, wr_mesh *mesh, char *err, size_t errlen);
int wr_scene_add_light(wr_scene *s, const wr_light *l, char *err, size_t errlen);
int wr_scene_set_ambient(wr_scene *s, int material, char *err, size_t errlen);
void wr_scene_set_camera(wr_scene *s, const wr_camera *c);

int wr_scene_find_material(const wr_scene *s, const char *name);   /* -1 if absent */
int wr_scene_find_object(const wr_scene *s, const char *name);
int wr_scene_find_light(const wr_scene *s, const char *name);

/* Nearest hit over all objects with t in [ray->tmin, ray->tmax]. Fills hit->obj, and hit->prim for meshes
 * (-1 otherwise). Thread-safe (read only). */
int wr_scene_intersect(const wr_scene *s, const wr_ray *ray, wr_hit *hit);
int wr_scene_occluded(const wr_scene *s, const wr_ray *ray);

/* Medium the ray is in just before / just after the hit: entering (front = 1) goes ambient -> object,
 * exiting (front = 0) goes object -> ambient. Returns material indices. */
int wr_scene_medium_before(const wr_scene *s, const wr_hit *hit);
int wr_scene_medium_after(const wr_scene *s, const wr_hit *hit);

/* Bounding box of all finite objects (spheres, meshes); planes are infinite and ignored.
 * Empty box (lo > hi) if there are none. */
wr_aabb wr_scene_bounds(const wr_scene *s);

/* 1 if consistent, else 0 with the first problem in err. */
int wr_scene_validate(const wr_scene *s, char *err, size_t errlen);

#endif /* WR_SCENE_H */
