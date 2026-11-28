/* End-to-end geometry check: pinhole camera -> analytic sphere and BVH sphere mesh -> normal-shaded PPM.
 * Not an optics test: it only shows that camera, BVH, primitives and OpenMP work together. */
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wr_test.h"
#include "wr_camera.h"
#include "wr_bvh.h"
#include "wr_optics.h"

#define W 256
#define H 256

typedef struct { unsigned char rgb[W * H * 3]; unsigned char hit[W * H]; } image;

static void render(image *img, const wr_camera *cam, const wr_sphere *sph, const wr_bvh *bvh)
{
    #pragma omp parallel for schedule(dynamic, 4)
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            wr_ray ray = wr_camera_ray(cam, x, y, 0.5, 0.5, 0.0, 0.0);
            wr_hit h;
            int ok = sph ? wr_sphere_intersect(sph, &ray, &h) : wr_bvh_intersect(bvh, &ray, &h);
            unsigned char *p = &img->rgb[3 * (y * W + x)];
            img->hit[y * W + x] = (unsigned char)ok;
            if (ok) {
                double s = -wr_v3_dot(ray.d, wr_hit_facing_normal(&h));       /* headlight shading */
                unsigned char v = (unsigned char)(255.0 * fmax(0.0, fmin(1.0, s)) + 0.5);
                p[0] = v; p[1] = v; p[2] = v;
            } else { p[0] = 20; p[1] = 24; p[2] = 40; }
        }
}

int main(void)
{
    wr_camera cam;
    double D = 5.0, R = 1.0, vfov = 40.0 * WR_PI / 180.0;
    CHECK(wr_camera_pinhole(&cam, wr_v3(0, 0, D), wr_v3_zero(), wr_v3(0, 1, 0), vfov, W, H));
    wr_sphere sph = wr_sphere_make(wr_v3_zero(), R);
    wr_mesh m;
    CHECK(wr_mesh_make_uv_sphere(wr_v3_zero(), R, 64, 32, &m));
    wr_bvh bvh;
    CHECK(wr_bvh_build(&bvh, &m));

    image *a = malloc(sizeof *a), *b = malloc(sizeof *b), *a1 = malloc(sizeof *a1), *b1 = malloc(sizeof *b1);
    omp_set_num_threads(4);
    render(a, &cam, &sph, NULL);
    render(b, &cam, NULL, &bvh);
    omp_set_num_threads(1);
    render(a1, &cam, &sph, NULL);
    render(b1, &cam, NULL, &bvh);
    CHECK(memcmp(a->rgb, a1->rgb, sizeof a->rgb) == 0);                      /* 4 threads == 1 thread, byte for byte */
    CHECK(memcmp(b->rgb, b1->rgb, sizeof b->rgb) == 0);

    /* analytic silhouette: the cone of half-angle asin(R/D) meets the image plane in a circle of radius
       f * tan(alpha), f = (H/2) / tan(vfov/2) pixels. Exact for a sphere on the optical axis. */
    double f = (H / 2.0) / tan(vfov / 2.0);
    double rpx = f * tan(asin(R / D));
    double area_exact = WR_PI * rpx * rpx;
    long na = 0, nb = 0, ndiff = 0;
    for (int i = 0; i < W * H; i++) { na += a->hit[i]; nb += b->hit[i]; ndiff += a->hit[i] != b->hit[i]; }
    double ea = fabs(na - area_exact) / area_exact, eb = fabs(nb - area_exact) / area_exact;
    CHECK(ea < 0.005);                                                       /* pixel-centre sampling of a circle of ~72 px radius */
    CHECK(nb <= na);                                                         /* inscribed mesh is never larger than the sphere */
    CHECK(eb < 0.015);
    printf("render: silhouette radius %.2f px, exact area %.0f px; analytic %ld (%.2f%%), mesh+BVH %ld (%.2f%%), pixels that differ %ld\n",
           rpx, area_exact, na, 100 * ea, nb, 100 * eb, ndiff);

    /* centre pixel looks straight at the sphere: headlight shading is 1 there, and every lit pixel is in [0,1] */
    CHECK(a->rgb[3 * ((H / 2) * W + W / 2)] >= 254 && b->rgb[3 * ((H / 2) * W + W / 2)] >= 254);

    /* side-by-side PPM (P6): analytic | mesh */
    FILE *fp = fopen("build/geometry_check.ppm", "wb");
    CHECK(fp != NULL);
    if (fp) {
        fprintf(fp, "P6\n%d %d\n255\n", 2 * W, H);
        for (int y = 0; y < H; y++) {
            fwrite(&a->rgb[3 * y * W], 1, 3 * W, fp);
            fwrite(&b->rgb[3 * y * W], 1, 3 * W, fp);
        }
        fclose(fp);
    }
    free(a); free(b); free(a1); free(b1);
    wr_bvh_free(&bvh);
    wr_mesh_free(&m);
    TEST_REPORT("render");
}
