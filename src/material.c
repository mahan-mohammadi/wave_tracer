#include <math.h>
#include <stdio.h>
#include <string.h>
#include "wr_material.h"
#include "wr_optics.h"

static wr_material base(const char *name)
{
    wr_material m;
    memset(&m, 0, sizeof m);
    snprintf(m.name, sizeof m.name, "%s", name ? name : "");
    m.surface = WR_SURF_DIELECTRIC;
    m.model = WR_IOR_CONSTANT;
    m.n0 = 1.0;
    m.lambda_min = 0.0;
    m.lambda_max = INFINITY;
    return m;
}

wr_material wr_material_constant(const char *name, double n, double kappa)
{
    wr_material m = base(name);
    m.n0 = n; m.kappa = kappa;
    return m;
}

wr_material wr_material_cauchy(const char *name, double a, double b, double c)
{
    wr_material m = base(name);
    m.model = WR_IOR_CAUCHY;
    m.cauchy[0] = a; m.cauchy[1] = b; m.cauchy[2] = c;
    m.lambda_min = 0.4e-6; m.lambda_max = 0.8e-6;
    return m;
}

wr_material wr_material_sellmeier(const char *name, const double b[3], const double c[3])
{
    wr_material m = base(name);
    m.model = WR_IOR_SELLMEIER;
    for (int i = 0; i < 3; i++) { m.sell_b[i] = b[i]; m.sell_c[i] = c[i]; }
    m.lambda_min = 0.3e-6; m.lambda_max = 2.5e-6;
    return m;
}

wr_material wr_material_mirror(const char *name)
{
    wr_material m = base(name);
    m.surface = WR_SURF_IDEAL_MIRROR;
    return m;
}

wr_material wr_material_absorber(const char *name)
{
    wr_material m = base(name);
    m.surface = WR_SURF_ABSORBER;
    return m;
}

wr_cplx wr_material_index(const wr_material *m, double lambda0, int *in_range)
{
    int ok = lambda0 > 0.0 && lambda0 >= m->lambda_min && lambda0 <= m->lambda_max;
    double n = 1.0, kappa = 0.0;
    if (m->surface == WR_SURF_DIELECTRIC) {
        kappa = m->kappa;
        double l = lambda0 * 1e6, l2 = l * l;                      /* micrometres */
        switch (m->model) {
        case WR_IOR_CONSTANT: n = m->n0; break;
        case WR_IOR_CAUCHY:   n = m->cauchy[0] + m->cauchy[1] / l2 + m->cauchy[2] / (l2 * l2); break;
        case WR_IOR_SELLMEIER: {
            double n2 = 1.0;
            for (int i = 0; i < 3; i++) n2 += m->sell_b[i] * l2 / (l2 - m->sell_c[i]);
            if (n2 > 0.0 && isfinite(n2)) n = sqrt(n2); else { n = NAN; ok = 0; }
            break;
        }
        }
        if (!(n > 0.0)) ok = 0;
    }
    if (in_range) *in_range = ok;
    return wr_c(n, kappa);
}

double wr_material_absorption_coeff(const wr_material *m, double lambda0)
{
    return m->surface == WR_SURF_DIELECTRIC ? 2.0 * WR_TWO_PI * m->kappa / lambda0 : 0.0;   /* 4 pi kappa / lambda0 */
}

static int bad(char *err, size_t n, const char *msg, const char *name)
{
    if (err && n) snprintf(err, n, "material '%s': %s", name, msg);
    return 0;
}

int wr_material_validate(const wr_material *m, char *err, size_t errlen)
{
    if (m->name[0] == '\0') return bad(err, errlen, "empty name", "?");
    if (!(m->kappa >= 0.0) || !isfinite(m->kappa)) return bad(err, errlen, "extinction must be finite and >= 0", m->name);
    if (!(m->lambda_min >= 0.0) || !(m->lambda_max > m->lambda_min)) return bad(err, errlen, "invalid wavelength range", m->name);
    if (m->surface != WR_SURF_DIELECTRIC) return 1;
    switch (m->model) {
    case WR_IOR_CONSTANT:
        if (!(m->n0 > 0.0) || !isfinite(m->n0)) return bad(err, errlen, "index must be finite and > 0", m->name);
        break;
    case WR_IOR_CAUCHY:
        for (int i = 0; i < 3; i++) if (!isfinite(m->cauchy[i])) return bad(err, errlen, "non-finite Cauchy coefficient", m->name);
        break;
    case WR_IOR_SELLMEIER:
        for (int i = 0; i < 3; i++)
            if (!isfinite(m->sell_b[i]) || !isfinite(m->sell_c[i]) || m->sell_c[i] < 0.0)
                return bad(err, errlen, "Sellmeier coefficients must be finite and C_i >= 0", m->name);
        break;
    }
    return 1;
}
