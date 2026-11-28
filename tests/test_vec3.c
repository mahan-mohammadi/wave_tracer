#include "wr_test.h"
#include "wr_rng.h"
#include "wr_vec3.h"

static wr_vec3 rand_unit(wr_rng *r)
{
    for (;;) {
        wr_vec3 v = wr_v3(wr_rng_range(r, -1, 1), wr_rng_range(r, -1, 1), wr_rng_range(r, -1, 1));
        double l = wr_v3_len(v);
        if (l > 0.1 && l <= 1.0) return wr_v3_scale(v, 1.0 / l);
    }
}

int main(void)
{
    wr_vec3 a = wr_v3(1, 2, 3), b = wr_v3(-4, 5, 0.5);
    CHECK_NEAR(wr_v3_dot(a, b), -4 + 10 + 1.5, 1e-15);
    wr_vec3 c = wr_v3_cross(wr_v3(1, 0, 0), wr_v3(0, 1, 0));
    CHECK_NEAR(c.x, 0, 0); CHECK_NEAR(c.y, 0, 0); CHECK_NEAR(c.z, 1, 0);
    CHECK_NEAR(wr_v3_len(wr_v3(3, 4, 12)), 13, 1e-15);
    CHECK_NEAR(wr_v3_dist(wr_v3(1, 1, 1), wr_v3(2, 3, 3)), 3, 1e-15);
    wr_vec3 z = wr_v3_normalize(wr_v3_zero());
    CHECK(z.x == 0 && z.y == 0 && z.z == 0);               /* no NaN */
    wr_vec3 l = wr_v3_lerp(a, b, 0.25);
    CHECK_NEAR(l.x, 0.75 * 1 + 0.25 * -4, 1e-15);

    wr_rng rng = {88172645463325252ULL};
    double worst_cross = 0, worst_basis = 0, worst_refl = 0, worst_snell = 0, worst_unit = 0;
    int tir = 0, trans = 0;
    for (int i = 0; i < 20000; i++) {
        wr_vec3 p = rand_unit(&rng), q = rand_unit(&rng);
        wr_vec3 cr = wr_v3_cross(p, q);
        worst_cross = fmax(worst_cross, fmax(fabs(wr_v3_dot(cr, p)), fabs(wr_v3_dot(cr, q))));
        /* Lagrange identity |p x q|^2 = |p|^2|q|^2 - (p.q)^2 */
        double lag = wr_v3_len2(cr) - (1.0 - wr_v3_dot(p, q) * wr_v3_dot(p, q));
        worst_cross = fmax(worst_cross, fabs(lag));

        wr_vec3 t, bb;
        wr_v3_basis(p, &t, &bb);
        worst_basis = fmax(worst_basis, fabs(wr_v3_len(t) - 1));
        worst_basis = fmax(worst_basis, fabs(wr_v3_len(bb) - 1));
        worst_basis = fmax(worst_basis, fabs(wr_v3_dot(t, bb)));
        worst_basis = fmax(worst_basis, fabs(wr_v3_dot(t, p)));
        worst_basis = fmax(worst_basis, fabs(wr_v3_dot(bb, p)));
        worst_basis = fmax(worst_basis, fabs(wr_v3_dot(wr_v3_cross(t, bb), p) - 1.0)); /* right-handed */

        /* reflection: unit length, angle preserved, tangential part kept */
        wr_vec3 n = q;
        wr_vec3 d = p;
        if (wr_v3_dot(d, n) > 0) n = wr_v3_neg(n);          /* n faces the incoming ray */
        wr_vec3 r = wr_v3_reflect(d, n);
        worst_refl = fmax(worst_refl, fabs(wr_v3_len(r) - 1));
        worst_refl = fmax(worst_refl, fabs(wr_v3_dot(r, n) + wr_v3_dot(d, n)));

        /* refraction vs Snell's law: eta sin(theta_i) = sin(theta_t) */
        double eta = wr_rng_range(&rng, 0.4, 2.5);
        wr_vec3 tdir;
        double cosi = -wr_v3_dot(d, n);
        double sini = sqrt(fmax(0.0, 1 - cosi * cosi));
        if (wr_v3_refract(d, n, eta, &tdir)) {
            trans++;
            double cost = -wr_v3_dot(tdir, n);
            double sint = sqrt(fmax(0.0, 1 - cost * cost));
            worst_snell = fmax(worst_snell, fabs(eta * sini - sint));
            worst_unit = fmax(worst_unit, fabs(wr_v3_len(tdir) - 1));
            CHECK(cost >= 0);                                 /* continues through the surface */
            /* coplanarity: d, n, tdir lie in one plane */
            worst_snell = fmax(worst_snell, fabs(wr_v3_dot(wr_v3_cross(d, n), tdir)));
        } else {
            tir++;
            CHECK(eta * sini > 1.0 - 1e-12);                  /* TIR only past the critical angle */
        }
    }
    CHECK(trans > 1000 && tir > 1000);                        /* both branches exercised */
    CHECK_NEAR(worst_cross, 0, 5e-15);
    CHECK_NEAR(worst_basis, 0, 5e-15);
    CHECK_NEAR(worst_refl, 0, 5e-15);
    CHECK_NEAR(worst_snell, 0, 1e-14);
    CHECK_NEAR(worst_unit, 0, 5e-15);
    printf("vec3: cross %.1e, basis %.1e, reflect %.1e, Snell %.1e (transmitted %d, TIR %d)\n",
           worst_cross, worst_basis, worst_refl, worst_snell, trans, tir);

    /* exact cases: normal incidence passes straight; 45 deg into glass n=1.5 */
    wr_vec3 out;
    CHECK(wr_v3_refract(wr_v3(0, 0, -1), wr_v3(0, 0, 1), 1.0 / 1.5, &out));
    CHECK_NEAR(out.x, 0, 1e-16); CHECK_NEAR(out.z, -1, 1e-15);
    double s = sqrt(0.5);
    CHECK(wr_v3_refract(wr_v3(s, 0, -s), wr_v3(0, 0, 1), 1.0 / 1.5, &out));
    CHECK_NEAR(out.x, s / 1.5, 1e-15);                         /* sin(theta_t) = sin(45 deg)/1.5 */
    /* glass to air at 45 deg: critical angle is 41.8 deg -> TIR */
    CHECK(!wr_v3_refract(wr_v3(s, 0, -s), wr_v3(0, 0, 1), 1.5, &out));
    /* just inside/outside the critical angle */
    double crit = asin(1.0 / 1.5);
    CHECK(wr_v3_refract(wr_v3(sin(crit - 1e-6), 0, -cos(crit - 1e-6)), wr_v3(0, 0, 1), 1.5, &out));
    CHECK(!wr_v3_refract(wr_v3(sin(crit + 1e-6), 0, -cos(crit + 1e-6)), wr_v3(0, 0, 1), 1.5, &out));

    TEST_REPORT("vec3");
}
