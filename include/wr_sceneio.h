/* wr_sceneio.h - text scene format (Phase 4c).
 *
 * One directive per line; '#' starts a comment; tokens are whitespace separated; numbers are decimal
 * (500e-9, inf allowed; nan rejected). Units: metres, seconds, radians, EXCEPT where a keyword says deg.
 * Names are at most 47 bytes, no spaces. Anything may refer only to things defined on earlier lines.
 *
 *   time T
 *   ambient MATERIAL
 *   material NAME constant N | cauchy A B C | sellmeier B1 B2 B3 C1 C2 C3 | ideal_mirror | absorber
 *            [extinction K] [range LMIN LMAX]            (Cauchy B in um^2, C in um^4; Sellmeier C in um^2; range in m)
 *   sphere NAME material M center X Y Z radius R [velocity X Y Z]
 *   plane  NAME material M point X Y Z normal X Y Z [velocity X Y Z]
 *   mesh   NAME material M  SOURCE  [OP]...  [velocity X Y Z]
 *            SOURCE: obj PATH | box LX LY LZ HX HY HZ | uvsphere CX CY CZ R SEG STACK
 *            OP (applied in the order written): translate X Y Z | scale X Y Z | rotate AX AY AZ DEGREES
 *            obj paths are relative to the scene file's directory (or the base directory given to the parser)
 *   light point|plane_wave|gaussian NAME  [position X Y Z] [direction X Y Z] [waist W]
 *            (wavelength L | lines N L1 W1 ... LN WN) [amplitude A] [phase P] [coherence_length L]
 *            [jones RE1 IM1 RE2 IM2] [pol_ref X Y Z] [velocity X Y Z]
 *            direction required for plane_wave/gaussian, waist for gaussian, wavelength or lines for all;
 *            amplitude defaults to 1; jones is normalized to unit length on load
 *   camera pinhole|orthographic|thin_lens  position X Y Z  look_at X Y Z  [up X Y Z]  size W H
 *            pinhole/thin_lens: vfov_deg V;  orthographic: view_height H;  thin_lens also: aperture A  focus F
 *
 * Every failure returns 0 with "line N: message" in err; the finished scene is validated with wr_scene_validate.
 */
#ifndef WR_SCENEIO_H
#define WR_SCENEIO_H

#include "wr_scene.h"

/* base_dir may be NULL (obj paths then used as written). On success the scene is initialized and must be
 * released with wr_scene_free; on failure it is released already and left zeroed. */
int wr_scene_parse(const char *text, const char *base_dir, wr_scene *scene, char *err, size_t errlen);
int wr_scene_load(const char *path, wr_scene *scene, char *err, size_t errlen);

#endif /* WR_SCENEIO_H */
