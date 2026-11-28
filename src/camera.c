#include <math.h>
#include "wr_camera.h"
#include "wr_optics.h"

static int make_frame(wr_camera *c, wr_vec3 pos, wr_vec3 look_at, wr_vec3 up_hint, int w, int h)
{
    if (w <= 0 || h <= 0) return 0;
    wr_vec3 f = wr_v3_normalize(wr_v3_sub(look_at, pos));
    if (wr_v3_len2(f) == 0.0) return 0;
    wr_vec3 r = wr_v3_cross(f, up_hint);
    double rl = wr_v3_len(r);
    if (!(rl > 1e-12 * wr_v3_len(up_hint))) return 0;           /* up_hint parallel to the view direction or zero */
    r = wr_v3_scale(r, 1.0 / rl);
    c->pos = pos; c->forward = f; c->right = r; c->up = wr_v3_cross(r, f);
    c->width = w; c->height = h;
    c->aperture = 0.0; c->focus = 1.0;
    c->tan_half_v = c->tan_half_h = c->half_h_ortho = c->half_w_ortho = 0.0;
    return 1;
}

int wr_camera_pinhole(wr_camera *c, wr_vec3 pos, wr_vec3 look_at, wr_vec3 up_hint,
                      double vfov, int width, int height)
{
    if (!(vfov > 0.0 && vfov < WR_PI)) return 0;
    if (!make_frame(c, pos, look_at, up_hint, width, height)) return 0;
    c->type = WR_CAM_PINHOLE;
    c->tan_half_v = tan(0.5 * vfov);
    c->tan_half_h = c->tan_half_v * (double)width / (double)height;
    return 1;
}

int wr_camera_orthographic(wr_camera *c, wr_vec3 pos, wr_vec3 look_at, wr_vec3 up_hint,
                           double view_height, int width, int height)
{
    if (!(view_height > 0.0) || !isfinite(view_height)) return 0;
    if (!make_frame(c, pos, look_at, up_hint, width, height)) return 0;
    c->type = WR_CAM_ORTHOGRAPHIC;
    c->half_h_ortho = 0.5 * view_height;
    c->half_w_ortho = c->half_h_ortho * (double)width / (double)height;
    return 1;
}

int wr_camera_thin_lens(wr_camera *c, wr_vec3 pos, wr_vec3 look_at, wr_vec3 up_hint,
                        double vfov, double aperture, double focus, int width, int height)
{
    if (!(aperture >= 0.0) || !(focus > 0.0) || !isfinite(focus)) return 0;
    if (!wr_camera_pinhole(c, pos, look_at, up_hint, vfov, width, height)) return 0;
    c->type = WR_CAM_THIN_LENS;
    c->aperture = aperture;
    c->focus = focus;
    return 1;
}

void wr_sample_disc_concentric(double u, double v, double *x, double *y)
{
    double a = 2.0 * u - 1.0, b = 2.0 * v - 1.0;
    if (a == 0.0 && b == 0.0) { *x = *y = 0.0; return; }
    double r, th;
    if (fabs(a) > fabs(b)) { r = a; th = (WR_PI / 4.0) * (b / a); }
    else { r = b; th = WR_PI / 2.0 - (WR_PI / 4.0) * (a / b); }
    *x = r * cos(th);
    *y = r * sin(th);
}

wr_ray wr_camera_ray(const wr_camera *c, int px, int py, double jx, double jy, double lu, double lv)
{
    double sx = 2.0 * ((px + jx) / (double)c->width) - 1.0;      /* +1 = right edge  */
    double sy = 1.0 - 2.0 * ((py + jy) / (double)c->height);     /* +1 = top edge    */
    if (c->type == WR_CAM_ORTHOGRAPHIC) {
        wr_vec3 o = wr_v3_madd(wr_v3_madd(c->pos, c->right, sx * c->half_w_ortho), c->up, sy * c->half_h_ortho);
        return wr_ray_make(o, c->forward);
    }
    /* pinhole direction, scaled so its component along forward is exactly 1 */
    wr_vec3 dir = wr_v3_madd(wr_v3_madd(c->forward, c->right, sx * c->tan_half_h), c->up, sy * c->tan_half_v);
    if (c->type == WR_CAM_PINHOLE || c->aperture == 0.0) return wr_ray_make(c->pos, dir);
    wr_vec3 focus_pt = wr_v3_madd(c->pos, dir, c->focus);        /* on the plane at distance `focus` along forward */
    double dx, dy;
    wr_sample_disc_concentric(lu, lv, &dx, &dy);
    wr_vec3 o = wr_v3_madd(wr_v3_madd(c->pos, c->right, dx * c->aperture), c->up, dy * c->aperture);
    return wr_ray_make(o, wr_v3_sub(focus_pt, o));
}

int wr_camera_project(const wr_camera *c, wr_vec3 p, double *px, double *py)
{
    if (c->type != WR_CAM_PINHOLE) return 0;
    wr_vec3 d = wr_v3_sub(p, c->pos);
    double z = wr_v3_dot(d, c->forward);
    if (!(z > 0.0)) return 0;
    double sx = wr_v3_dot(d, c->right) / (z * c->tan_half_h);
    double sy = wr_v3_dot(d, c->up) / (z * c->tan_half_v);
    *px = 0.5 * (sx + 1.0) * c->width;
    *py = 0.5 * (1.0 - sy) * c->height;
    return 1;
}
