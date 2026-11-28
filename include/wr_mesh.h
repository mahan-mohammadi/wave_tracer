/* wr_mesh.h - indexed triangle meshes and an OBJ loader (Phase 4b).
 *
 * Geometry only. Triangles are counter-clockwise seen from outside (same
 * convention as wr_triangle), so a closed mesh has positive signed volume.
 *
 * OBJ support: v, f only (v, v/vt, v/vt/vn, v//vn tokens; negative indices;
 * polygons are fan-triangulated, which is correct for convex faces only).
 * vt, vn, o, g, s, usemtl, mtllib and unknown lines are ignored (materials
 * are assigned by the scene, not by the file). Line continuation with '\' is
 * not supported. Any index that is 0 or out of range is an error.
 */
#ifndef WR_MESH_H
#define WR_MESH_H

#include <stddef.h>
#include <stdint.h>
#include "wr_vec3.h"
#include "wr_xform.h"
#include "wr_geom.h"

typedef struct {
    wr_vec3  *v;      /* nv vertices */
    size_t    nv;
    uint32_t *idx;    /* 3 * nt vertex indices */
    size_t    nt;
} wr_mesh;

wr_mesh     wr_mesh_empty(void);
void        wr_mesh_free(wr_mesh *m);                /* safe on empty; resets to empty */
wr_triangle wr_mesh_triangle(const wr_mesh *m, size_t i);
wr_aabb     wr_mesh_bounds(const wr_mesh *m);
double      wr_mesh_area(const wr_mesh *m);
/* Signed volume sum(v0 . (v1 x v2)) / 6; positive for a closed, outward-wound mesh. */
double      wr_mesh_volume(const wr_mesh *m);
/* 1 if the mesh is a closed, consistently oriented 2-manifold surface: every directed edge occurs exactly
 * once and so does its reverse. (Says nothing about inside/outside: check wr_mesh_volume > 0 for that.) */
int         wr_mesh_is_closed_oriented(const wr_mesh *m);

/* In-place affine transform. A transform with negative determinant (mirror)
 * reverses the triangle winding so normals stay outward. Returns 0 (mesh
 * untouched) if the linear part is singular. */
int         wr_mesh_transform(wr_mesh *m, const wr_xform *x);

/* All loaders/generators return 1 on success, 0 on failure (*out is then
 * empty; err, if non-NULL, receives a message). */
int wr_mesh_parse_obj(const char *text, wr_mesh *out, char *err, size_t errlen);
int wr_mesh_load_obj(const char *path, wr_mesh *out, char *err, size_t errlen);
int wr_mesh_make_box(wr_vec3 lo, wr_vec3 hi, wr_mesh *out);                 /* 8 vertices, 12 triangles */
/* UV sphere, y up; nseg >= 3 around, nstack >= 2 from pole to pole. Inscribed: every
 * vertex lies on the true sphere, so the mesh is slightly smaller than the sphere. */
int wr_mesh_make_uv_sphere(wr_vec3 center, double radius, int nseg, int nstack, wr_mesh *out);

#endif /* WR_MESH_H */
