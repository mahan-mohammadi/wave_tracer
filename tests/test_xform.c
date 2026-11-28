#include "wr_test.h"
#include "wr_rng.h"
#include "wr_xform.h"
#include "wr_optics.h"

static double mat_dist(const wr_mat3 *a, const wr_mat3 *b)
{
    double d = 0;
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) d = fmax(d, fabs(a->m[i][j] - b->m[i][j]));
    return d;
}
static wr_vec3 rv(wr_rng *r, double lo, double hi)
{
    return wr_v3(wr_rng_range(r, lo, hi), wr_rng_range(r, lo, hi), wr_rng_range(r, lo, hi));
}

int main(void)
{
    wr_mat3 I = wr_mat3_identity();

    /* 90 deg about z maps x -> y, y -> -x (right-handed) */
    wr_mat3 rz = wr_mat3_rotation(wr_v3(0, 0, 1), (WR_PI / 2));
    wr_vec3 v = wr_mat3_apply(&rz, wr_v3(1, 0, 0));
    CHECK_NEAR(v.x, 0, 1e-16); CHECK_NEAR(v.y, 1, 1e-16); CHECK_NEAR(v.z, 0, 0);
    v = wr_mat3_apply(&rz, wr_v3(0, 1, 0));
    CHECK_NEAR(v.x, -1, 1e-16); CHECK_NEAR(v.y, 0, 1e-16);
    wr_mat3 r0 = wr_mat3_rotation(wr_v3_zero(), 1.0);          /* degenerate axis */
    CHECK_NEAR(mat_dist(&r0, &I), 0, 0);

    wr_rng rng = {0x9E3779B97F4A7C15ULL};
    double w_orth = 0, w_det = 0, w_comp = 0, w_inv = 0, w_rt = 0, w_nrm = 0, w_cmp = 0;
    for (int i = 0; i < 5000; i++) {
        wr_vec3 axis = rv(&rng, -1, 1);
        if (wr_v3_len(axis) < 0.05) continue;
        double a1 = wr_rng_range(&rng, -7, 7), a2 = wr_rng_range(&rng, -7, 7);
        wr_mat3 R = wr_mat3_rotation(axis, a1);
        wr_mat3 Rt = wr_mat3_transpose(&R);
        wr_mat3 RRt = wr_mat3_mul(&R, &Rt);
        w_orth = fmax(w_orth, mat_dist(&RRt, &I));
        w_det = fmax(w_det, fabs(wr_mat3_det(&R) - 1.0));
        /* same-axis angles add */
        wr_mat3 Ra = wr_mat3_rotation(axis, a1), Rb = wr_mat3_rotation(axis, a2), Rs = wr_mat3_rotation(axis, a1 + a2);
        wr_mat3 Rab = wr_mat3_mul(&Ra, &Rb);
        w_comp = fmax(w_comp, mat_dist(&Rab, &Rs));
        /* axis is fixed */
        wr_vec3 k = wr_v3_normalize(axis);
        w_comp = fmax(w_comp, wr_v3_dist(wr_mat3_apply(&R, k), k));

        /* random affine: rotate, anisotropic scale, translate */
        wr_xform T = wr_xform_translate(rv(&rng, -5, 5));
        wr_xform Rx = wr_xform_rotate(axis, a1);
        wr_xform S = wr_xform_scale(wr_v3(wr_rng_range(&rng, 0.3, 3), wr_rng_range(&rng, 0.3, 3), wr_rng_range(&rng, 0.3, 3)));
        wr_xform RS = wr_xform_compose(&Rx, &S);
        wr_xform X = wr_xform_compose(&T, &RS);
        wr_xform Xi;
        CHECK(wr_xform_inverse(&X, &Xi));
        wr_vec3 p = rv(&rng, -3, 3);
        w_rt = fmax(w_rt, wr_v3_dist(wr_xform_point(&Xi, wr_xform_point(&X, p)), p));
        w_cmp = fmax(w_cmp, wr_v3_dist(wr_xform_point(&X, p),
                         wr_xform_point(&T, wr_xform_point(&Rx, wr_xform_point(&S, p)))));  /* compose order */
        wr_mat3 MMi = wr_mat3_mul(&X.m, &Xi.m);
        w_inv = fmax(w_inv, mat_dist(&MMi, &I));

        /* vectors ignore translation; points do not */
        wr_vec3 dv = wr_v3_sub(wr_xform_point(&X, wr_v3_add(p, wr_v3(1, 2, 3))), wr_xform_point(&X, p));
        w_cmp = fmax(w_cmp, wr_v3_dist(dv, wr_xform_vector(&X, wr_v3(1, 2, 3))));

        /* normals stay perpendicular to the transformed surface under non-uniform scale */
        wr_vec3 n = wr_v3_normalize(rv(&rng, -1, 1));
        wr_vec3 t1, t2;
        wr_v3_basis(n, &t1, &t2);
        wr_vec3 nn = wr_xform_normal(&Xi, n);
        w_nrm = fmax(w_nrm, fabs(wr_v3_len(nn) - 1));
        w_nrm = fmax(w_nrm, fabs(wr_v3_dot(nn, wr_xform_vector(&X, t1))));
        w_nrm = fmax(w_nrm, fabs(wr_v3_dot(nn, wr_xform_vector(&X, t2))));
    }
    CHECK_NEAR(w_orth, 0, 5e-15);
    CHECK_NEAR(w_det, 0, 5e-15);
    CHECK_NEAR(w_comp, 0, 1e-14);
    CHECK_NEAR(w_inv, 0, 1e-14);
    CHECK_NEAR(w_rt, 0, 1e-13);
    CHECK_NEAR(w_cmp, 0, 1e-13);
    CHECK_NEAR(w_nrm, 0, 1e-14);
    printf("xform: orthonormality %.1e, det-1 %.1e, angle-add %.1e, M*Minv %.1e, roundtrip %.1e, normals %.1e\n",
           w_orth, w_det, w_comp, w_inv, w_rt, w_nrm);

    /* singular transforms are reported, not inverted */
    wr_xform flat = wr_xform_scale(wr_v3(1, 0, 1)), out;
    CHECK(!wr_xform_inverse(&flat, &out));
    wr_mat3 rank1 = {{{1, 2, 3}, {2, 4, 6}, {3, 6, 9}}};
    wr_mat3 mo;
    CHECK(!wr_mat3_inverse(&rank1, &mo));
    /* known inverse */
    wr_mat3 A = {{{2, 0, 0}, {0, 4, 0}, {0, 0, 8}}}, Ai;
    CHECK(wr_mat3_inverse(&A, &Ai));
    CHECK_NEAR(Ai.m[0][0], 0.5, 0); CHECK_NEAR(Ai.m[1][1], 0.25, 0); CHECK_NEAR(Ai.m[2][2], 0.125, 0);

    TEST_REPORT("xform");
}
