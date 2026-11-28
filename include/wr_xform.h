/* wr_xform.h - 3x3 matrices and affine transforms (Phase 4a).
 *
 * wr_mat3 is row-major. wr_xform is the affine map  p' = M p + t.
 * Points get the translation, vectors do not, normals use the inverse
 * transpose (so they stay perpendicular to the surface under non-uniform
 * scale). Perspective is not needed here; cameras are handled separately.
 */
#ifndef WR_XFORM_H
#define WR_XFORM_H

#include "wr_vec3.h"

typedef struct { double m[3][3]; } wr_mat3;
typedef struct { wr_mat3 m; wr_vec3 t; } wr_xform;

wr_mat3 wr_mat3_identity(void);
wr_mat3 wr_mat3_mul(const wr_mat3 *a, const wr_mat3 *b);   /* a * b */
wr_mat3 wr_mat3_transpose(const wr_mat3 *a);
double  wr_mat3_det(const wr_mat3 *a);
/* Returns 1 on success, 0 if the matrix is (numerically) singular, in which
 * case *out is left untouched. "Singular" means |det| <= 1e-14 * maxabs^3. */
int     wr_mat3_inverse(const wr_mat3 *a, wr_mat3 *out);
wr_vec3 wr_mat3_apply(const wr_mat3 *a, wr_vec3 v);
/* Rotation about a (not necessarily unit) axis by angle in radians,
 * right-handed (Rodrigues). Zero axis gives identity. */
wr_mat3 wr_mat3_rotation(wr_vec3 axis, double angle);

wr_xform wr_xform_identity(void);
wr_xform wr_xform_translate(wr_vec3 t);
wr_xform wr_xform_scale(wr_vec3 s);
wr_xform wr_xform_rotate(wr_vec3 axis, double angle);
/* Apply b first, then a. */
wr_xform wr_xform_compose(const wr_xform *a, const wr_xform *b);
/* Returns 1 on success, 0 if the linear part is singular. */
int      wr_xform_inverse(const wr_xform *x, wr_xform *out);

wr_vec3 wr_xform_point(const wr_xform *x, wr_vec3 p);
wr_vec3 wr_xform_vector(const wr_xform *x, wr_vec3 v);
/* Transform a normal. Pass the INVERSE of the transform that moves the
 * surface; result is the normalized (M^-1)^T n. */
wr_vec3 wr_xform_normal(const wr_xform *inverse, wr_vec3 n);

#endif /* WR_XFORM_H */
