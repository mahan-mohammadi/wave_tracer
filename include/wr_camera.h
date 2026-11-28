/* wr_camera.h - geometric camera models (Phase 4c, camera part).
 *
 * Generates rays for pixels. This is the GEOMETRIC camera only: where rays
 * start and which way they go. The wave-optical sensor model (aperture
 * diffraction, point-spread function, spectral response) comes with the
 * transport phase and must not be faked by this layer.
 *
 * Image conventions: pixel (px, py) has its origin at the top-left corner,
 * px grows to the right, py grows downward. Sub-pixel offsets (jx, jy) are in
 * [0,1); 0.5, 0.5 is the pixel centre. Camera space: forward is the viewing
 * direction, right = normalize(forward x up_hint), up = right x forward
 * (right-handed; forward = -z, up_hint = +y gives right = +x).
 *
 * Models:
 *   PINHOLE      : all rays start at the position. vfov is the full vertical field of view.
 *   ORTHOGRAPHIC : parallel rays; `height` is the full vertical extent of the view window in metres.
 *   THIN_LENS    : pinhole rays refocused through a circular lens of radius `aperture` so that the plane at
 *                  distance `focus` (measured along forward) is sharp. aperture = 0 reduces to PINHOLE.
 *                  Geometric-optics model: defocus blur is the ray-optics circle of confusion, no diffraction.
 */
#ifndef WR_CAMERA_H
#define WR_CAMERA_H

#include "wr_vec3.h"
#include "wr_ray.h"

typedef enum { WR_CAM_PINHOLE, WR_CAM_ORTHOGRAPHIC, WR_CAM_THIN_LENS } wr_cam_type;

typedef struct {
    wr_cam_type type;
    wr_vec3 pos, forward, right, up;   /* orthonormal frame */
    int     width, height;             /* image size in pixels */
    double  tan_half_v, tan_half_h;    /* pinhole / thin lens: tan(vfov/2) and tan(hfov/2) (square pixels) */
    double  half_h_ortho, half_w_ortho;/* orthographic: half extents of the view window */
    double  aperture, focus;           /* thin lens: lens radius and focus distance */
} wr_camera;

/* Return 1 on success, 0 on invalid input (degenerate frame, vfov outside (0, pi), non-positive size/height/focus). */
int wr_camera_pinhole(wr_camera *c, wr_vec3 pos, wr_vec3 look_at, wr_vec3 up_hint,
                      double vfov_rad, int width, int height);
int wr_camera_orthographic(wr_camera *c, wr_vec3 pos, wr_vec3 look_at, wr_vec3 up_hint,
                           double view_height, int width, int height);
int wr_camera_thin_lens(wr_camera *c, wr_vec3 pos, wr_vec3 look_at, wr_vec3 up_hint,
                        double vfov_rad, double aperture_radius, double focus_distance,
                        int width, int height);

/* Ray through pixel (px, py) at sub-pixel offset (jx, jy). (lu, lv) in [0,1)^2 select the lens point
 * (ignored unless THIN_LENS; they are mapped to the disc area-uniformly by the concentric map). The
 * returned ray has unit direction and default tmin/tmax. */
wr_ray wr_camera_ray(const wr_camera *c, int px, int py, double jx, double jy, double lu, double lv);

/* Pinhole only: project a world point to continuous pixel coordinates (centre of pixel (i,j) is
 * (i + 0.5, j + 0.5)). Returns 0 if the point is behind the camera (or the camera is not pinhole). */
int wr_camera_project(const wr_camera *c, wr_vec3 p, double *px, double *py);

/* Shirley-Chiu concentric map from the unit square to the unit disc (area preserving). */
void wr_sample_disc_concentric(double u, double v, double *x, double *y);

#endif /* WR_CAMERA_H */
