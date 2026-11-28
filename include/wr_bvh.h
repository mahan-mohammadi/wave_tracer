/* wr_bvh.h - bounding volume hierarchy over a triangle mesh (Phase 4b).
 *
 * Built with a binned surface-area heuristic (16 bins, leaves of at most 4
 * triangles). Build is serial; queries are read-only, so any number of
 * threads may trace concurrently (OpenMP over pixels is safe).
 *
 * Boxes are padded by 1e-12 * (largest coordinate magnitude of the mesh) so
 * that rounding in the slab test cannot drop a hit that the exact triangle
 * test would find. Closest-hit results are identical to brute force: the
 * same wr_triangle_intersect decides every hit, the BVH only skips triangles.
 *
 * The BVH keeps a pointer to the mesh; the mesh must outlive it and must not
 * be modified (vertices moved) without rebuilding.
 */
#ifndef WR_BVH_H
#define WR_BVH_H

#include <stdint.h>
#include "wr_mesh.h"

typedef struct {
    wr_aabb  box;
    uint32_t a;      /* leaf: first entry in order[]; interior: index of right child (left = self + 1) */
    uint32_t n;      /* leaf: triangle count (> 0); interior: 0 */
    uint32_t axis;   /* interior: split axis, used for near-child-first traversal */
} wr_bvh_node;

typedef struct {
    const wr_mesh *mesh;
    wr_bvh_node   *nodes;
    size_t         nnodes;
    uint32_t      *order;       /* triangle ids, leaves index into this */
    /* statistics */
    size_t         nleaves;
    int            max_depth;
    uint32_t       max_leaf;
    double         sah_cost;    /* expected cost per ray, in triangle tests, trav cost 0.5 */
} wr_bvh;

/* Returns 1 on success. An empty mesh gives a valid BVH that is never hit. */
int  wr_bvh_build(wr_bvh *bvh, const wr_mesh *mesh);
void wr_bvh_free(wr_bvh *bvh);

/* Nearest hit with t in [ray->tmin, ray->tmax]; fills hit (hit->prim = triangle index). */
int  wr_bvh_intersect(const wr_bvh *bvh, const wr_ray *ray, wr_hit *hit);
/* True if anything is hit in [tmin, tmax]; cheaper, no hit data. */
int  wr_bvh_occluded(const wr_bvh *bvh, const wr_ray *ray);
/* Reference implementation: tests every triangle. Used by the test suite. */
int  wr_mesh_intersect_brute(const wr_mesh *mesh, const wr_ray *ray, wr_hit *hit);

#endif /* WR_BVH_H */
