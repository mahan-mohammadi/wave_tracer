#include "wr_test.h"
#include "wr_rng.h"
#include "wr_camera.h"
#include "wr_geom.h"
#include "wr_optics.h"

static wr_vec3 rv(wr_rng *r, double lo, double hi)
{
    return wr_v3(wr_rng_range(r, lo, hi), wr_rng_range(r, lo, hi), wr_rng_range(r, lo, hi));
}

int main(void)
{
    wr_camera cam;
    double vfov = WR_PI / 3.0;     /* 60 degrees */

    /* frame and orientation: forward -z, up +y  =>  right +x, top of the image is +y */
    CHECK(wr_camera_pinhole(&cam, wr_v3(0, 0, 0), wr_v3(0, 0, -1), wr_v3(0, 1, 0), vfov, 200, 100));
    CHECK_NEAR(cam.right.x, 1, 1e-16); CHECK_NEAR(cam.up.y, 1, 1e-16);
    wr_ray r = wr_camera_ray(&cam, 100, 50, 0, 0, 0, 0);               /* exact image centre */
    CHECK_NEAR(r.d.x, 0, 1e-16); CHECK_NEAR(r.d.y, 0, 1e-16); CHECK_NEAR(r.d.z, -1, 1e-16);
    r = wr_camera_ray(&cam, 150, 25, 0.5, 0.5, 0, 0);                  /* right of centre, above centre */
    CHECK(r.d.x > 0 && r.d.y > 0);
    r = wr_camera_ray(&cam, 20, 80, 0.5, 0.5, 0, 0);                   /* left, below */
    CHECK(r.d.x < 0 && r.d.y < 0);

    /* field of view: top edge makes exactly vfov/2 with the axis; right edge atan(tan(vfov/2) * aspect) */
    r = wr_camera_ray(&cam, 100, 0, 0, 0, 0, 0);
    CHECK_NEAR(acos(wr_v3_dot(r.d, cam.forward)), vfov / 2, 1e-14);
    r = wr_camera_ray(&cam, 199, 50, 1.0, 0, 0, 0);
    CHECK_NEAR(atan2(wr_v3_dot(r.d, cam.right), wr_v3_dot(r.d, cam.forward)), atan(tan(vfov / 2) * 2.0), 1e-14);
    /* the pixel grid is uniform on the image plane: neighbouring directions differ by constant tan-space steps */
    {
        wr_ray a = wr_camera_ray(&cam, 10, 50, 0, 0, 0, 0), b = wr_camera_ray(&cam, 11, 50, 0, 0, 0, 0);
        wr_ray c = wr_camera_ray(&cam, 150, 50, 0, 0, 0, 0), d = wr_camera_ray(&cam, 151, 50, 0, 0, 0, 0);
        double s1 = wr_v3_dot(b.d, cam.right) / wr_v3_dot(b.d, cam.forward) - wr_v3_dot(a.d, cam.right) / wr_v3_dot(a.d, cam.forward);
        double s2 = wr_v3_dot(d.d, cam.right) / wr_v3_dot(d.d, cam.forward) - wr_v3_dot(c.d, cam.right) / wr_v3_dot(c.d, cam.forward);
        CHECK_NEAR(s1, s2, 1e-15);
        CHECK_NEAR(s1, 2.0 * cam.tan_half_h / 200.0, 1e-15);
    }

    /* random poses: orthonormal right-handed frame (right x up = -forward), ray -> point -> project round trip */
    wr_rng rng = {2024};
    double w_frame = 0, w_rt = 0;
    int n_rt = 0;
    for (int i = 0; i < 5000; i++) {
        wr_vec3 pos = rv(&rng, -5, 5), look = rv(&rng, -5, 5), up = rv(&rng, -1, 1);
        int W = 16 + (int)(wr_rng_uniform(&rng) * 300), H = 16 + (int)(wr_rng_uniform(&rng) * 300);
        double fov = wr_rng_range(&rng, 0.2, 2.5);
        if (!wr_camera_pinhole(&cam, pos, look, up, fov, W, H)) continue;
        w_frame = fmax(w_frame, fabs(wr_v3_len(cam.right) - 1));
        w_frame = fmax(w_frame, fabs(wr_v3_len(cam.up) - 1));
        w_frame = fmax(w_frame, fabs(wr_v3_dot(cam.right, cam.up)));
        w_frame = fmax(w_frame, fabs(wr_v3_dot(cam.right, cam.forward)));
        w_frame = fmax(w_frame, fabs(wr_v3_dot(cam.up, cam.forward)));
        w_frame = fmax(w_frame, wr_v3_dist(wr_v3_cross(cam.right, cam.up), wr_v3_neg(cam.forward)));
        int px = (int)(wr_rng_uniform(&rng) * W), py = (int)(wr_rng_uniform(&rng) * H);
        double jx = wr_rng_uniform(&rng), jy = wr_rng_uniform(&rng);
        wr_ray ray = wr_camera_ray(&cam, px, py, jx, jy, 0, 0);
        w_frame = fmax(w_frame, fabs(wr_v3_len(ray.d) - 1));
        wr_vec3 p = wr_ray_at(&ray, wr_rng_range(&rng, 0.1, 50));
        double qx, qy;
        CHECK(wr_camera_project(&cam, p, &qx, &qy));
        w_rt = fmax(w_rt, fmax(fabs(qx - (px + jx)), fabs(qy - (py + jy))));
        n_rt++;
        wr_vec3 behind = wr_v3_madd(pos, cam.forward, -1.0);
        CHECK(!wr_camera_project(&cam, behind, &qx, &qy));
    }
    CHECK(n_rt > 4000);
    CHECK_NEAR(w_frame, 0, 5e-15);
    CHECK_NEAR(w_rt, 0, 1e-9);
    printf("camera: %d random poses, frame orthonormality %.1e, ray->point->pixel round trip %.1e px\n", n_rt, w_frame, w_rt);

    /* orthographic: parallel rays, origins on a regular grid covering the window */
    CHECK(wr_camera_orthographic(&cam, wr_v3(1, 2, 3), wr_v3(1, 2, 0), wr_v3(0, 1, 0), 4.0, 200, 100));
    wr_ray a = wr_camera_ray(&cam, 0, 0, 0, 0, 0, 0), b = wr_camera_ray(&cam, 199, 99, 1, 1, 0, 0), c = wr_camera_ray(&cam, 100, 50, 0, 0, 0, 0);
    CHECK_NEAR(wr_v3_dist(a.d, b.d), 0, 1e-16);
    CHECK_NEAR(wr_v3_dist(c.o, cam.pos), 0, 1e-16);                    /* centre ray starts at the position */
    wr_vec3 span = wr_v3_sub(b.o, a.o);
    CHECK_NEAR(wr_v3_dot(span, cam.right), 8.0, 1e-14);                /* window width = 4 * aspect(2) */
    CHECK_NEAR(wr_v3_dot(span, cam.up), -4.0, 1e-14);                  /* top-left to bottom-right goes down */
    wr_ray a2 = wr_camera_ray(&cam, 1, 0, 0, 0, 0, 0);
    CHECK_NEAR(wr_v3_dist(a2.o, a.o), 4.0 / 100.0, 1e-15);            /* pixel pitch = view_height / height */

    /* thin lens: every lens sample sees the same point on the focus plane (sharp there), origins stay inside
       the aperture, and off-focus planes show the geometric-optics circle of confusion */
    double F = 3.0, A = 0.05;
    CHECK(wr_camera_thin_lens(&cam, wr_v3(0.2, -0.1, 0.4), wr_v3(1, 1, -2), wr_v3(0, 1, 0), 0.8, A, F, 128, 96));
    wr_plane focal = wr_plane_make(wr_v3_madd(cam.pos, cam.forward, F), cam.forward);
    double w_focus = 0, w_lens = 0, w_coc = 0;
    for (int i = 0; i < 3000; i++) {
        int px = (int)(wr_rng_uniform(&rng) * 128), py = (int)(wr_rng_uniform(&rng) * 96);
        double jx = wr_rng_uniform(&rng), jy = wr_rng_uniform(&rng);
        wr_camera pin = cam; pin.type = WR_CAM_PINHOLE;
        wr_ray ref = wr_camera_ray(&pin, px, py, jx, jy, 0, 0);
        wr_hit hr;
        wr_plane_intersect(&focal, &ref, &hr);
        wr_vec3 pf = hr.p;
        for (int k = 0; k < 8; k++) {
            wr_ray lr = wr_camera_ray(&cam, px, py, jx, jy, wr_rng_uniform(&rng), wr_rng_uniform(&rng));
            wr_hit hl;
            lr.tmin = 0.0;
            CHECK(wr_plane_intersect(&focal, &lr, &hl));
            w_focus = fmax(w_focus, wr_v3_dist(hl.p, pf));
            wr_vec3 off = wr_v3_sub(lr.o, cam.pos);
            w_lens = fmax(w_lens, fabs(wr_v3_dot(off, cam.forward)));                 /* lens plane is perpendicular to forward */
            CHECK(wr_v3_len(off) <= A * (1 + 1e-14));
            /* circle of confusion on the plane at depth z: radius = A |1 - z/F| */
            double z = 2.0 * F;
            wr_plane pz = wr_plane_make(wr_v3_madd(cam.pos, cam.forward, z), cam.forward);
            wr_hit hz;
            lr.tmin = 0.0;
            wr_plane_intersect(&pz, &lr, &hz);
            wr_vec3 chief;                                                           /* where the pinhole ray crosses that plane */
            wr_hit hc; wr_plane_intersect(&pz, &ref, &hc); chief = hc.p;
            w_coc = fmax(w_coc, wr_v3_dist(hz.p, chief));
        }
    }
    CHECK_NEAR(w_focus, 0, 1e-13);
    CHECK_NEAR(w_lens, 0, 1e-16);
    CHECK(w_coc <= A * (1 + 1e-12) && w_coc > 0.9 * A);                               /* blur disc radius A at z = 2F */
    printf("camera: thin lens, focus-plane spread %.1e m, max blur radius at 2F = %.5f m (aperture %.3f)\n", w_focus, w_coc, A);

    /* aperture 0 is exactly a pinhole */
    CHECK(wr_camera_thin_lens(&cam, wr_v3_zero(), wr_v3(0, 0, -1), wr_v3(0, 1, 0), 1.0, 0.0, 2.0, 64, 64));
    wr_ray t1 = wr_camera_ray(&cam, 5, 9, 0.3, 0.7, 0.123, 0.456);
    wr_camera pin2 = cam; pin2.type = WR_CAM_PINHOLE;
    wr_ray t2 = wr_camera_ray(&pin2, 5, 9, 0.3, 0.7, 0, 0);
    CHECK(wr_v3_dist(t1.d, t2.d) == 0.0 && wr_v3_dist(t1.o, t2.o) == 0.0);

    /* concentric disc map: stays in the unit disc, area uniform (mean r^2 = 1/2, quadrant balance) */
    {
        int n = 256, q[4] = {0, 0, 0, 0};
        double sum = 0, worst = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double x, y;
                wr_sample_disc_concentric((i + 0.5) / n, (j + 0.5) / n, &x, &y);
                double rr = x * x + y * y;
                worst = fmax(worst, rr);
                sum += rr;
                q[(x >= 0) + 2 * (y >= 0)]++;
            }
        CHECK(worst <= 1.0 + 1e-14);
        CHECK_NEAR(sum / ((double)n * n), 0.5, 2e-4);
        for (int k = 0; k < 4; k++) CHECK_NEAR((double)q[k] / ((double)n * n), 0.25, 1e-3);
        double x, y; wr_sample_disc_concentric(0.5, 0.5, &x, &y);
        CHECK(x == 0.0 && y == 0.0);
    }

    /* invalid construction is rejected, not silently accepted */
    wr_camera bad;
    CHECK(!wr_camera_pinhole(&bad, wr_v3_zero(), wr_v3(0, 1, 0), wr_v3(0, 1, 0), 1.0, 10, 10));   /* up parallel to view */
    CHECK(!wr_camera_pinhole(&bad, wr_v3_zero(), wr_v3_zero(), wr_v3(0, 1, 0), 1.0, 10, 10));     /* look_at == pos */
    CHECK(!wr_camera_pinhole(&bad, wr_v3_zero(), wr_v3(0, 0, -1), wr_v3_zero(), 1.0, 10, 10));    /* zero up */
    CHECK(!wr_camera_pinhole(&bad, wr_v3_zero(), wr_v3(0, 0, -1), wr_v3(0, 1, 0), 0.0, 10, 10));
    CHECK(!wr_camera_pinhole(&bad, wr_v3_zero(), wr_v3(0, 0, -1), wr_v3(0, 1, 0), WR_PI, 10, 10));
    CHECK(!wr_camera_pinhole(&bad, wr_v3_zero(), wr_v3(0, 0, -1), wr_v3(0, 1, 0), 1.0, 0, 10));
    CHECK(!wr_camera_orthographic(&bad, wr_v3_zero(), wr_v3(0, 0, -1), wr_v3(0, 1, 0), 0.0, 10, 10));
    CHECK(!wr_camera_thin_lens(&bad, wr_v3_zero(), wr_v3(0, 0, -1), wr_v3(0, 1, 0), 1.0, -0.1, 1.0, 10, 10));
    CHECK(!wr_camera_thin_lens(&bad, wr_v3_zero(), wr_v3(0, 0, -1), wr_v3(0, 1, 0), 1.0, 0.1, 0.0, 10, 10));
    CHECK(!wr_camera_project(&cam, wr_v3(0, 0, -1), &a.d.x, &a.d.y));                            /* thin lens: no projection */

    TEST_REPORT("camera");
}
