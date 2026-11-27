/* Validation: Newton's rings (roadmap Phase 3, scene 5).
 *
 * Scene: a plano-convex lens (radius of curvature R, glass n = 1.5) rests on a
 * flat glass plate; light of wavelength lambda hits it at normal incidence from
 * inside the lens. The gap between lens and plate has thickness
 *     h(r) = R - sqrt(R^2 - r^2)      (exact sphere, not the parabola)
 * and is filled with a medium of index n_g. The reflected field at each pixel is
 * the coherent sum of all partial waves of the layer (wr_layer_series): no ring
 * formula appears anywhere in the code; rings, their radii and the dark centre
 * emerge from phase accumulation and the sign (pi shift) of the Fresnel
 * coefficients.
 *
 * Theory (reflection, n_lens = n_plate): the two interfaces have r12 = -r23, so
 * the reflected intensity is monotone in cos(delta), delta = 2 k n_g h:
 *     dark   at 2 n_g h = m lambda      -> r_m^2 = m lambda R / n_g          (m = 0: dark centre)
 *     bright at 2 n_g h = (m+1/2) lambda -> r^2 = (m + 1/2) lambda R / n_g
 *     I_max = (2 r12 / (1 + r12^2))^2  (all partial waves), not 4 r12^2 (two-beam)
 * (paraxial radii; the sphere correction is a relative r^2/(4R^2) ~ 2e-6 here.)
 * Limits: normal incidence, scalar, lossless, thin-element (local-film)
 * approximation: valid while the gap slope r/R << 1 (here <= 3.2e-3) and
 * transverse spreading across the gap is negligible (gap <= 5 um).
 * Measured along the x axis by parabolic refinement of the intensity extrema.
 * Writes build/newton_rings.ppm (P6, sqrt tone map) for viewing.
 */
#include <stdio.h>
#include <stdlib.h>
#include "wr_test.h"
#include "wr_field.h"
#include "wr_film.h"
#include "wr_optics.h"

#define N 512
static const double DX = 12.5e-6, RC = 1.0, LAM = 500e-9, NG_LENS = 1.5;

/* returns max |R+T-1| over the image through *e_err */
static void build(wr_field *f, double ng, double *e_err)
{
    double worst = 0.0;
    #pragma omp parallel for schedule(static) reduction(max : worst)
    for (long iy = 0; iy < N; iy++)
        for (long ix = 0; ix < N; ix++) {
            double x = wr_field_x(f, (size_t)ix), y = wr_field_y(f, (size_t)iy);
            double r2 = x * x + y * y;
            double h = r2 / (RC + sqrt(RC * RC - r2));           /* exact sphere, no cancellation */
            wr_rt rt = wr_layer_series(NG_LENS, ng, NG_LENS, h, LAM, 30);
            *wr_field_at(f, (size_t)ix, (size_t)iy) = rt.r;
            worst = fmax(worst, fabs(wr_layer_reflectance(rt) + wr_layer_transmittance(rt, NG_LENS, NG_LENS) - 1.0));
        }
    *e_err = worst;
}

/* position of the intensity extremum near x0 on row `row` (positive x side) */
static double extremum(const wr_field *f, size_t row, double x0, double half, int want_max)
{
    long lo = (long)floor((x0 - half) / f->dx + (double)(f->nx / 2));
    long hi = (long)ceil((x0 + half) / f->dx + (double)(f->nx / 2));
    long best = lo;
    double bv = want_max ? -1.0 : 1e30;
    for (long i = lo; i <= hi; i++) {
        double v = wr_cabs2(*wr_field_at((wr_field *)f, (size_t)i, row));
        if (want_max ? v > bv : v < bv) { bv = v; best = i; }
    }
    double a = wr_cabs2(*wr_field_at((wr_field *)f, (size_t)best - 1, row));
    double c = wr_cabs2(*wr_field_at((wr_field *)f, (size_t)best + 1, row));
    double den = a - 2.0 * bv + c;
    return wr_field_x(f, (size_t)best) + ((den != 0.0) ? 0.5 * (a - c) / den : 0.0) * f->dx;
}

static double run_case(double ng, const char *label, int write_ppm)
{
    wr_field *f = wr_field_create(N, N, DX, DX);
    double e_err;
    build(f, ng, &e_err);
    const size_t row = N / 2;       /* y = +dx/2: 6 um off the axis, negligible */
    const double y0 = wr_field_y(f, row);
    double imax = 0.0, imin = 1e30;
    for (size_t i = 0; i < f->nx * f->ny; i++) { double v = wr_cabs2(f->data[i]); imax = fmax(imax, v); imin = fmin(imin, v); }

    printf("\n%s: n_gap = %.2f, R = %.1f m, lambda = %.0f nm, grid %dx%d at %.1f um\n", label, ng, RC, LAM * 1e9, N, N, DX * 1e6);
    printf("  R+T-1 max error %.1e; I_max = %.5f (all partial waves: %.5f, two-beam 4r^2: %.5f)\n", e_err, imax,
           pow(2.0 * wr_fresnel_r_normal(NG_LENS, ng) / (1.0 + pow(wr_fresnel_r_normal(NG_LENS, ng), 2.0)), 2.0),
           4.0 * pow(wr_fresnel_r_normal(NG_LENS, ng), 2.0));
    CHECK_NEAR(e_err, 0.0, 1e-13);
    double r12 = wr_fresnel_r_normal(NG_LENS, ng);
    CHECK_NEAR(imax, pow(2.0 * r12 / (1.0 + r12 * r12), 2.0), 5e-5);
    CHECK(imin < 1e-4 * imax);                        /* dark centre (pi shift), no ring formula used */
    double centre = wr_cabs2(*wr_field_at(f, N / 2, row)) ;
    CHECK(centre < 1e-4 * imax);

    printf("  ring   expected r (um)  measured r (um)   error (um)\n");
    double worst = 0.0;
    const int mmax = (ng < 1.2) ? 14 : 16;
    for (int m = 1; m <= mmax; m++) {
        for (int bright = 0; bright <= 1; bright++) {
            double rexp = sqrt((m - (bright ? 0.5 : 0.0)) * LAM * RC / ng);        /* bright: (m-1/2), m >= 1 */
            double sp = 0.5 * LAM * RC / (ng * rexp);                                  /* local ring spacing / 1 */
            double rmeas = extremum(f, row, rexp, 0.3 * sp, bright);
            double err = rmeas - rexp;
            worst = fmax(worst, fabs(err));
            if (m == 1 || m == 5 || m == 10 || m == mmax)
                printf("  %s%2d   %12.3f   %12.3f   %+9.3f\n", bright ? "B" : "D", m, rexp * 1e6, rmeas * 1e6, err * 1e6);
        }
    }
    (void)y0;
    printf("  worst radius error over all measured rings: %.3f um\n", worst * 1e6);

    if (write_ppm) {
        FILE *fp = fopen("build/newton_rings.ppm", "wb");
        if (fp) {
            fprintf(fp, "P6\n%d %d\n255\n", N, N);
            for (size_t i = 0; i < f->nx * f->ny; i++) {
                unsigned char g = (unsigned char)(255.0 * sqrt(wr_cabs2(f->data[i]) / imax) + 0.5);
                unsigned char px[3] = {g, g, g};
                fwrite(px, 1, 3, fp);
            }
            fclose(fp);
        }
    }
    wr_field_destroy(f);
    return worst;
}

int main(void)
{
    double w_air = run_case(1.0, "air gap", 1);
    double w_wat = run_case(1.33, "water-filled gap", 0);
    CHECK(w_air < 0.5e-6);       /* < 0.5 um on rings ~ 90-700 um apart */
    CHECK(w_wat < 0.5e-6);
    TEST_REPORT("newton");
}
