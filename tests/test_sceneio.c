#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include "wr_test.h"
#include "wr_sceneio.h"
#include "wr_optics.h"

static int has(const char *hay, const char *needle) { return strstr(hay, needle) != NULL; }

static const char *FULL =
    "# every directive once\r\n"
    "time 1.5e-9\r\n"
    "material air constant 1.000277 range 3e-7 1e-6   # trailing comment\r\n"
    "ambient air\n"
    "\n"
    "material bk7  sellmeier 1.03961212 0.231792344 1.01046945 0.00600069867 0.0200179144 103.560653\n"
    "material crown cauchy 1.5046 0.0042 0 range 4e-7 7e-7\n"
    "material tint constant 1.5 extinction 2e-6\n"
    "material mir ideal_mirror\n"
    "material blk absorber\n"
    "sphere ball material bk7 center 1 2 3 radius 0.5 velocity 10 0 0\n"
    "plane  floor material mir point 0 -1 0 normal 0 3 4\n"
    "mesh   cube material bk7 box 0 0 0 1 2 3 scale 2 3 4 translate 1 1 1\n"
    "mesh   mirrored material bk7 box 0 0 0 1 1 1 scale -1 1 1 velocity 0 0 -2\n"
    "mesh   rot material tint box 0 0 0 1 2 3 rotate 0 0 1 90\n"
    "mesh   ts material bk7 box 0 0 0 1 1 1 translate 1 0 0 scale 2 1 1\n"
    "mesh   st material bk7 box 0 0 0 1 1 1 scale 2 1 1 translate 1 0 0\n"
    "mesh   orb material bk7 uvsphere 0 0 0 1 64 32\n"
    "light point p1 position 1 2 3 wavelength 500e-9\n"
    "light plane_wave pw direction 0 0 2 lines 3 450e-9 1 550e-9 2 650e-9 1 amplitude 3 phase 0.5 coherence_length 0.1 velocity 0 0 5\n"
    "light gaussian gb position 0 0 0 direction 1 0 0 waist 2e-3 wavelength 633e-9 jones 1 0 0 1 pol_ref 0 0 1\n"
    "camera pinhole position 0 0 5 look_at 0 0 0 vfov_deg 60 size 200 100\n";

static wr_vec3 lo_of(const wr_scene *s, const char *name) { wr_aabb b = wr_mesh_bounds(s->objects[wr_scene_find_object(s, name)].mesh); return b.lo; }
static wr_vec3 hi_of(const wr_scene *s, const char *name) { wr_aabb b = wr_mesh_bounds(s->objects[wr_scene_find_object(s, name)].mesh); return b.hi; }

typedef struct { const char *text; const char *msg; int line; } errcase;

int main(void)
{
    char err[300] = "";
    wr_scene s;

    CHECK(wr_scene_parse(FULL, NULL, &s, err, sizeof err));
    if (err[0]) printf("unexpected: %s\n", err);
    CHECK(s.nmat == 7 && s.nobj == 8 && s.nlight == 3 && s.has_camera);
    CHECK(s.time == 1.5e-9);
    CHECK(strcmp(s.materials[s.ambient].name, "air") == 0);

    /* materials */
    const wr_material *air = &s.materials[wr_scene_find_material(&s, "air")];
    CHECK(air->n0 == 1.000277 && air->lambda_min == 3e-7 && air->lambda_max == 1e-6);
    int ok;
    wr_cplx n587 = wr_material_index(&s.materials[wr_scene_find_material(&s, "bk7")], 587.56e-9, &ok);
    CHECK(ok); CHECK_NEAR(n587.re, 1.5168, 2e-5);                          /* the file reproduces the catalogue index */
    const wr_material *cr = &s.materials[wr_scene_find_material(&s, "crown")];
    CHECK(cr->model == WR_IOR_CAUCHY && cr->cauchy[1] == 0.0042 && cr->lambda_min == 4e-7 && cr->lambda_max == 7e-7);
    CHECK(s.materials[wr_scene_find_material(&s, "tint")].kappa == 2e-6);
    CHECK(s.materials[wr_scene_find_material(&s, "mir")].surface == WR_SURF_IDEAL_MIRROR);
    CHECK(s.materials[wr_scene_find_material(&s, "blk")].surface == WR_SURF_ABSORBER);

    /* analytic objects */
    const wr_object *ball = &s.objects[wr_scene_find_object(&s, "ball")];
    CHECK(ball->kind == WR_OBJ_SPHERE && ball->sphere.radius == 0.5 && ball->sphere.center.z == 3 && ball->velocity.x == 10);
    const wr_object *fl = &s.objects[wr_scene_find_object(&s, "floor")];
    CHECK_NEAR(fl->plane.normal.y, 0.6, 1e-15); CHECK_NEAR(fl->plane.normal.z, 0.8, 1e-15);   /* (0,3,4) normalized */

    /* mesh sources and transform order */
    wr_vec3 lo = lo_of(&s, "cube"), hi = hi_of(&s, "cube");
    CHECK(lo.x == 1 && lo.y == 1 && lo.z == 1 && hi.x == 3 && hi.y == 7 && hi.z == 13);
    CHECK_NEAR(wr_mesh_volume(s.objects[wr_scene_find_object(&s, "cube")].mesh), 144.0, 1e-12);
    CHECK_NEAR(wr_mesh_volume(s.objects[wr_scene_find_object(&s, "mirrored")].mesh), 1.0, 1e-14);        /* mirror kept outward */
    CHECK(s.objects[wr_scene_find_object(&s, "mirrored")].velocity.z == -2);
    lo = lo_of(&s, "rot"); hi = hi_of(&s, "rot");                          /* 90 deg about z: x in [-2,0], y in [0,1] */
    CHECK_NEAR(lo.x, -2, 1e-15); CHECK_NEAR(hi.x, 0, 1e-15); CHECK_NEAR(lo.y, 0, 1e-15); CHECK_NEAR(hi.y, 1, 1e-15);
    CHECK_NEAR(wr_mesh_volume(s.objects[wr_scene_find_object(&s, "rot")].mesh), 6.0, 1e-13);
    CHECK(lo_of(&s, "ts").x == 2 && hi_of(&s, "ts").x == 4);               /* translate then scale */
    CHECK(lo_of(&s, "st").x == 1 && hi_of(&s, "st").x == 3);               /* scale then translate: order matters */
    double vr = wr_mesh_volume(s.objects[wr_scene_find_object(&s, "orb")].mesh) / (4.0 / 3.0 * WR_PI);
    CHECK(vr > 0.99 && vr < 1.0);

    /* lights */
    const wr_light *p1 = &s.lights[0], *pw = &s.lights[1], *gb = &s.lights[2];
    CHECK(p1->type == WR_LIGHT_POINT && p1->amplitude == 1.0 && p1->spectrum.n == 1 && p1->spectrum.lambda[0] == 500e-9);
    CHECK(isinf(p1->coherence_length) && p1->pos.z == 3);
    CHECK(pw->type == WR_LIGHT_PLANE_WAVE && pw->spectrum.n == 3 && pw->spectrum.lambda[1] == 550e-9 && pw->spectrum.weight[1] == 2);
    CHECK_NEAR(pw->dir.z, 1.0, 1e-15);                                      /* (0,0,2) normalized */
    CHECK(pw->amplitude == 3 && pw->phase == 0.5 && pw->coherence_length == 0.1 && pw->velocity.z == 5);
    CHECK(gb->type == WR_LIGHT_GAUSSIAN && gb->waist == 2e-3 && gb->pol_ref.z == 1);
    CHECK_NEAR(gb->jones[0].re, sqrt(0.5), 1e-15); CHECK_NEAR(gb->jones[1].im, sqrt(0.5), 1e-15);   /* "1 0 0 1" -> (1, i)/sqrt2 */
    CHECK(wr_light_validate(gb, err, sizeof err));

    /* camera */
    CHECK(s.camera.type == WR_CAM_PINHOLE && s.camera.width == 200 && s.camera.height == 100);
    CHECK_NEAR(s.camera.tan_half_v, tan(WR_PI / 6), 1e-15);                 /* 60 deg field of view, converted from degrees */
    wr_scene_free(&s);

    /* the other camera types */
    CHECK(wr_scene_parse("camera orthographic position 0 0 5 look_at 0 0 0 view_height 2 size 40 20\n", NULL, &s, err, sizeof err));
    CHECK(s.camera.type == WR_CAM_ORTHOGRAPHIC && s.camera.half_h_ortho == 1.0 && s.camera.half_w_ortho == 2.0);
    wr_scene_free(&s);
    CHECK(wr_scene_parse("camera thin_lens position 0 0 5 look_at 0 0 0 up 0 0 1 vfov_deg 40 aperture 0.01 focus 5 size 8 8\n"
                         "light point x wavelength 5e-7\n", NULL, &s, err, sizeof err) == 0);   /* up parallel to view direction */
    CHECK(has(err, "line 1") && has(err, "invalid camera"));
    CHECK(wr_scene_parse("camera thin_lens position 0 0 5 look_at 0 0 0 vfov_deg 40 aperture 0.01 focus 5 size 8 8\n", NULL, &s, err, sizeof err));
    CHECK(s.camera.type == WR_CAM_THIN_LENS && s.camera.aperture == 0.01 && s.camera.focus == 5);
    wr_scene_free(&s);

    /* empty / comment-only text is a valid empty scene */
    CHECK(wr_scene_parse("", NULL, &s, err, sizeof err) && s.nobj == 0 && s.nmat == 1);
    wr_scene_free(&s);
    CHECK(wr_scene_parse("# nothing\n   \n\t\n", NULL, &s, err, sizeof err));
    wr_scene_free(&s);

    /* open and inside-out solids written to files, referenced relative to a base directory */
    FILE *f = fopen("build/open_tri.obj", "wb");
    CHECK(f != NULL);
    if (f) { fputs("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n", f); fclose(f); }
    f = fopen("build/inward.obj", "wb");
    if (f) {
        fputs("v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\nf 1 2 3\nf 1 4 2\nf 1 3 4\nf 2 4 3\n", f);   /* tetrahedron wound inward */
        fclose(f);
    }
    CHECK(!wr_scene_parse("material glass constant 1.5\nmesh t material glass obj open_tri.obj\n", "build", &s, err, sizeof err));
    CHECK(has(err, "scene:") && has(err, "closed") && s.materials == NULL);
    CHECK(!wr_scene_parse("material glass constant 1.5\nmesh t material glass obj inward.obj\n", "build", &s, err, sizeof err));
    CHECK(has(err, "wound inward"));
    CHECK(wr_scene_parse("material m ideal_mirror\nmesh t material m obj open_tri.obj\n", "build", &s, err, sizeof err));   /* open mirror sheet is fine */
    wr_scene_free(&s);

    /* every error path: message, and the line number it is reported on */
    static const errcase cases[] = {
        {"bogus 1 2\n", "unknown directive", 1},
        {"material a constant\n", "needs 1 number", 1},
        {"material a constant 1.5x\n", "bad number", 1},
        {"material a constant nan\n", "bad number", 1},
        {"material a constant -1\n", "index must be", 1},
        {"material a foo 1\n", "unknown material model", 1},
        {"material vacuum constant 1\n", "duplicate material", 1},
        {"material a constant 1.5 range 5e-7\n", "needs 2 number", 1},
        {"material a constant 1.5 extinction -1\n", "extinction", 1},
        {"material a constant 1.5 wibble\n", "unexpected 'wibble'", 1},
        {"material\n", "needs a name and a model", 1},
        {"material a sellmeier 1 2 3\n", "needs 6 number", 1},
        {"sphere s material nope center 0 0 0 radius 1\n", "unknown material 'nope'", 1},
        {"sphere s center 0 0 0 radius 1\n", "needs 'material'", 1},
        {"sphere s material vacuum center 0 0 0\n", "needs center and radius", 1},
        {"sphere s material vacuum center 0 0 radius 1\n", "bad number 'radius'", 1},
        {"sphere s material vacuum center 0 0 0 radius -1\n", "radius must be", 1},
        {"sphere s material vacuum center 0 0 0 radius 1\nsphere s material vacuum center 1 0 0 radius 1\n", "duplicate object", 2},
        {"plane p material vacuum point 0 0 0 normal 0 0 0\n", "non-zero normal", 1},
        {"plane p material vacuum point 0 0 0\n", "needs point and normal", 1},
        {"mesh m material vacuum box 0 0 0 1 1\n", "needs 6 number", 1},
        {"mesh m material vacuum box 1 1 1 0 0 0\n", "lo < hi", 1},
        {"mesh m material vacuum\n", "needs a source", 1},
        {"mesh m material vacuum obj /nonexistent/x.obj\n", "cannot open", 1},
        {"mesh m material vacuum box 0 0 0 1 1 1 box 0 0 0 1 1 1\n", "more than one source", 1},
        {"mesh m material vacuum box 0 0 0 1 1 1 rotate 0 0 0 30\n", "axis must be non-zero", 1},
        {"mesh m material vacuum box 0 0 0 1 1 1 scale 1 0 1\n", "singular", 1},
        {"mesh m material vacuum uvsphere 0 0 0 1 2.5 8\n", "integers", 1},
        {"mesh m material vacuum uvsphere 0 0 0 1 2 8\n", "segments >= 3", 1},
        {"mesh m material vacuum box 0 0 0 1 1 1 shear 1\n", "unexpected 'shear'", 1},
        {"light spot l wavelength 5e-7\n", "unknown light type", 1},
        {"light point l\n", "needs 'wavelength' or 'lines'", 1},
        {"light plane_wave l wavelength 5e-7\n", "needs 'direction'", 1},
        {"light plane_wave l wavelength 5e-7 direction 0 0 0\n", "non-zero", 1},
        {"light gaussian l wavelength 5e-7 direction 0 0 1\n", "needs 'waist'", 1},
        {"light point l wavelength -5e-7\n", "wavelength", 1},
        {"light point l wavelength 5e-7 jones 0 0 0 0\n", "non-zero", 1},
        {"light point l lines 2 4e-7 1\n", "wavelength/weight pairs", 1},
        {"light point l lines 17 4e-7 1\n", "between 1 and 16", 1},
        {"light point l wavelength 5e-7 amplitude 0\n", "amplitude", 1},
        {"light point l wavelength 5e-7\nlight point l wavelength 6e-7\n", "duplicate light", 2},
        {"camera fisheye position 0 0 0\n", "unknown camera type", 1},
        {"camera pinhole position 0 0 0 look_at 0 0 -1 size 10 10\n", "needs vfov_deg", 1},
        {"camera pinhole position 0 0 0 look_at 0 0 -1 vfov_deg 60 size 10.5 10\n", "integers", 1},
        {"camera pinhole position 0 0 0 look_at 0 0 -1 vfov_deg 180 size 10 10\n", "invalid camera", 1},
        {"camera pinhole position 0 0 0 look_at 0 0 -1 vfov_deg 60\n", "position, look_at and size", 1},
        {"camera thin_lens position 0 0 0 look_at 0 0 -1 vfov_deg 60 size 10 10 aperture 0.1\n", "aperture and focus", 1},
        {"camera orthographic position 0 0 0 look_at 0 0 -1 size 10 10\n", "needs view_height", 1},
        {"ambient nope\n", "unknown material", 1},
        {"material m ideal_mirror\nambient m\n", "must be a dielectric", 2},
        {"time\n", "time needs", 1},
        {"time inf\n", "time needs", 1},
        {"# a\n\n  # b\nmaterial a constant 1.5\n\n\nsphere s material a center 0 0 0 radius x\n", "bad number 'x'", 7},
        {"material a constant 1.5\r\nsphere s material a center 0 0 0 radius -2\r\n", "radius must be", 2},
    };
    int nbad = 0;
    for (size_t i = 0; i < sizeof cases / sizeof *cases; i++) {
        char want[16];
        snprintf(want, sizeof want, "line %d:", cases[i].line);
        err[0] = '\0';
        int rc = wr_scene_parse(cases[i].text, NULL, &s, err, sizeof err);
        int good = rc == 0 && strncmp(err, want, strlen(want)) == 0 && has(err, cases[i].msg) && s.materials == NULL && s.nobj == 0;
        if (!good) { nbad++; printf("  error case %zu failed: rc=%d err='%s' (wanted '%s' at line %d)\n", i, rc, err, cases[i].msg, cases[i].line); }
        CHECK(good);
        if (rc) wr_scene_free(&s);
    }
    printf("sceneio: %zu error cases reported with the right message and line (%d wrong)\n", sizeof cases / sizeof *cases, nbad);

    /* the demo scene file, end to end: load, query, render */
    CHECK(wr_scene_load("scenes/demo.scene", &s, err, sizeof err));
    if (err[0]) printf("demo: %s\n", err);
    CHECK(s.nmat == 6 && s.nobj == 5 && s.nlight == 2 && s.has_camera);
    CHECK(strcmp(s.materials[s.ambient].name, "air") == 0);
    CHECK(s.lights[0].spectrum.n == 3 && s.lights[0].amplitude == 2 && s.lights[1].coherence_length == 0.3);
    CHECK(!wr_scene_load("scenes/no_such_file.scene", &(wr_scene){0}, err, sizeof err) && has(err, "cannot open"));

    int W = s.camera.width, H = s.camera.height;
    unsigned char *img = calloc((size_t)W * H, 3), *img1 = calloc((size_t)W * H, 3);
    int *ids = malloc((size_t)W * H * sizeof *ids);
    static const unsigned char pal[8][3] = {{230, 90, 80}, {90, 200, 120}, {90, 140, 240}, {230, 200, 90}, {120, 120, 130}, {200, 100, 200}, {80, 200, 200}, {240, 240, 240}};
    for (int pass = 0; pass < 2; pass++) {
        omp_set_num_threads(pass == 0 ? 4 : 1);
        unsigned char *out = pass == 0 ? img : img1;
        #pragma omp parallel for schedule(dynamic, 4)
        for (int y = 0; y < H; y++)
            for (int x = 0; x < W; x++) {
                wr_ray ray = wr_camera_ray(&s.camera, x, y, 0.5, 0.5, 0, 0);
                wr_hit h;
                unsigned char *px = &out[3 * (y * W + x)];
                if (wr_scene_intersect(&s, &ray, &h)) {
                    double shade = 0.35 + 0.65 * fmax(0.0, -wr_v3_dot(ray.d, wr_hit_facing_normal(&h)));
                    for (int c = 0; c < 3; c++) px[c] = (unsigned char)(pal[h.obj % 8][c] * shade);
                    if (pass == 0) ids[y * W + x] = h.obj;
                } else { px[0] = 15; px[1] = 18; px[2] = 30; if (pass == 0) ids[y * W + x] = -1; }
            }
    }
    CHECK(memcmp(img, img1, (size_t)W * H * 3) == 0);                      /* 4 threads == 1 thread */
    int cnt[8] = {0}, miss = 0;
    for (int i = 0; i < W * H; i++) { if (ids[i] < 0) miss++; else cnt[ids[i]]++; }
    printf("sceneio: demo.scene render %dx%d, visible pixels per object:", W, H);
    for (int i = 0; i < s.nobj; i++) { printf(" %s=%d", s.objects[i].name, cnt[i]); CHECK(cnt[i] > 0); }
    printf(" (background %d)\n", miss);

    /* aimed rays: unobstructed paths with analytic answers */
    wr_vec3 cp = s.camera.pos;
    int lens = wr_scene_find_object(&s, "lens"), prism = wr_scene_find_object(&s, "prism"), back = wr_scene_find_object(&s, "back");
    wr_hit h;
    wr_vec3 c = s.objects[lens].sphere.center;
    wr_ray r = wr_ray_make(cp, wr_v3_sub(c, cp));
    CHECK(wr_scene_intersect(&s, &r, &h) && h.obj == lens && h.front == 1);
    CHECK_NEAR(h.t, wr_v3_dist(cp, c) - 0.02, 1e-12);                       /* distance to the centre minus the radius */
    CHECK(wr_scene_medium_before(&s, &h) == s.ambient && strcmp(s.materials[wr_scene_medium_after(&s, &h)].name, "bk7") == 0);
    wr_vec3 tgt = wr_v3(0, 0.0346410 / 3.0, 0.02);                          /* centroid of the prism's +z end face */
    r = wr_ray_make(cp, wr_v3_sub(tgt, cp));
    CHECK(wr_scene_intersect(&s, &r, &h) && h.obj == prism && h.prim == 1 && h.front == 1);
    CHECK_NEAR(h.n.z, 1.0, 1e-15); CHECK_NEAR(h.t, wr_v3_dist(cp, tgt), 1e-9);
    r = wr_ray_make(cp, wr_v3_sub(wr_v3(0.0, 0.1, -0.3), cp));
    CHECK(wr_scene_intersect(&s, &r, &h) && h.obj == back && h.front == 1);
    CHECK_NEAR(h.t, wr_v3_dist(cp, wr_v3(0.0, 0.1, -0.3)), 1e-12);
    wr_aabb bb = wr_scene_bounds(&s);
    CHECK(bb.lo.x <= -0.02 && bb.hi.x >= 0.12 && bb.hi.y >= 0.0346 && bb.lo.z <= -0.05);

    FILE *pf = fopen("build/scene_check.ppm", "wb");
    CHECK(pf != NULL);
    if (pf) { fprintf(pf, "P6\n%d %d\n255\n", W, H); fwrite(img, 3, (size_t)W * H, pf); fclose(pf); }
    free(img); free(img1); free(ids);
    wr_scene_free(&s);

    TEST_REPORT("sceneio");
}
