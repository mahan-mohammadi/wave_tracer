#include <math.h>
#include "wr_xform.h"

wr_mat3 wr_mat3_identity(void)
{
    wr_mat3 r = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    return r;
}

wr_mat3 wr_mat3_mul(const wr_mat3 *a, const wr_mat3 *b)
{
    wr_mat3 r;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            r.m[i][j] = a->m[i][0] * b->m[0][j] + a->m[i][1] * b->m[1][j] + a->m[i][2] * b->m[2][j];
    return r;
}

wr_mat3 wr_mat3_transpose(const wr_mat3 *a)
{
    wr_mat3 r;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) r.m[i][j] = a->m[j][i];
    return r;
}

double wr_mat3_det(const wr_mat3 *a)
{
    const double (*m)[3] = a->m;
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
         + m[0][1] * (m[1][2] * m[2][0] - m[1][0] * m[2][2])
         + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

int wr_mat3_inverse(const wr_mat3 *a, wr_mat3 *out)
{
    const double (*m)[3] = a->m;
    double c00 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
    double c01 = m[1][2] * m[2][0] - m[1][0] * m[2][2];
    double c02 = m[1][0] * m[2][1] - m[1][1] * m[2][0];
    double det = m[0][0] * c00 + m[0][1] * c01 + m[0][2] * c02;
    double mx = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) mx = fmax(mx, fabs(m[i][j]));
    if (!(fabs(det) > 1e-14 * mx * mx * mx)) return 0;   /* also catches NaN */
    double id = 1.0 / det;
    wr_mat3 r;
    r.m[0][0] = c00 * id;
    r.m[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * id;
    r.m[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * id;
    r.m[1][0] = c01 * id;
    r.m[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * id;
    r.m[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * id;
    r.m[2][0] = c02 * id;
    r.m[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * id;
    r.m[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * id;
    *out = r;
    return 1;
}

wr_vec3 wr_mat3_apply(const wr_mat3 *a, wr_vec3 v)
{
    return wr_v3(a->m[0][0] * v.x + a->m[0][1] * v.y + a->m[0][2] * v.z,
                 a->m[1][0] * v.x + a->m[1][1] * v.y + a->m[1][2] * v.z,
                 a->m[2][0] * v.x + a->m[2][1] * v.y + a->m[2][2] * v.z);
}

wr_mat3 wr_mat3_rotation(wr_vec3 axis, double angle)
{
    wr_vec3 k = wr_v3_normalize(axis);
    double c = cos(angle), s = sin(angle), w = 1.0 - c;
    wr_mat3 r = {{
        {c + w * k.x * k.x,         w * k.x * k.y - s * k.z,  w * k.x * k.z + s * k.y},
        {w * k.y * k.x + s * k.z,   c + w * k.y * k.y,        w * k.y * k.z - s * k.x},
        {w * k.z * k.x - s * k.y,   w * k.z * k.y + s * k.x,  c + w * k.z * k.z}}};
    if (k.x == 0.0 && k.y == 0.0 && k.z == 0.0) return wr_mat3_identity();
    return r;
}

wr_xform wr_xform_identity(void)
{
    wr_xform x = {wr_mat3_identity(), wr_v3_zero()};
    return x;
}

wr_xform wr_xform_translate(wr_vec3 t)
{
    wr_xform x = {wr_mat3_identity(), t};
    return x;
}

wr_xform wr_xform_scale(wr_vec3 s)
{
    wr_xform x = {{{{s.x, 0, 0}, {0, s.y, 0}, {0, 0, s.z}}}, wr_v3_zero()};
    return x;
}

wr_xform wr_xform_rotate(wr_vec3 axis, double angle)
{
    wr_xform x = {wr_mat3_rotation(axis, angle), wr_v3_zero()};
    return x;
}

wr_xform wr_xform_compose(const wr_xform *a, const wr_xform *b)
{
    wr_xform r;
    r.m = wr_mat3_mul(&a->m, &b->m);
    r.t = wr_v3_add(wr_mat3_apply(&a->m, b->t), a->t);
    return r;
}

int wr_xform_inverse(const wr_xform *x, wr_xform *out)
{
    wr_mat3 mi;
    if (!wr_mat3_inverse(&x->m, &mi)) return 0;
    out->m = mi;
    out->t = wr_v3_neg(wr_mat3_apply(&mi, x->t));
    return 1;
}

wr_vec3 wr_xform_point(const wr_xform *x, wr_vec3 p)
{
    return wr_v3_add(wr_mat3_apply(&x->m, p), x->t);
}

wr_vec3 wr_xform_vector(const wr_xform *x, wr_vec3 v) { return wr_mat3_apply(&x->m, v); }

wr_vec3 wr_xform_normal(const wr_xform *inverse, wr_vec3 n)
{
    wr_mat3 it = wr_mat3_transpose(&inverse->m);
    return wr_v3_normalize(wr_mat3_apply(&it, n));
}
