#include "wr_film.h"
#include "wr_optics.h"

double wr_fresnel_r_normal(double n_i, double n_t) { return (n_i - n_t) / (n_i + n_t); }
double wr_fresnel_t_normal(double n_i, double n_t) { return 2.0 * n_i / (n_i + n_t); }

/* one-way phase factor exp(i k0 n2 d); phase taken from the OPL so that
 * thick films keep full precision. e = half^2 = exp(i delta). */
static wr_cplx one_way(double n2, double d, double lambda0)
{
    return wr_cexp_i(wr_phase_from_opl(n2 * d, lambda0));
}

wr_rt wr_layer_series(double n1, double n2, double n3, double d, double lambda0, int nterms)
{
    double r12 = wr_fresnel_r_normal(n1, n2), r21 = -r12, r23 = wr_fresnel_r_normal(n2, n3);
    double t12 = wr_fresnel_t_normal(n1, n2), t21 = wr_fresnel_t_normal(n2, n1), t23 = wr_fresnel_t_normal(n2, n3);
    wr_cplx half = one_way(n2, d, lambda0), e = wr_cmul(half, half);
    wr_cplx q = wr_cscale(e, r21 * r23);           /* ratio between successive partial waves */
    wr_cplx term = wr_c(1.0, 0.0), sum = wr_c(0.0, 0.0);
    for (int j = 0; j < nterms; j++) { sum = wr_cadd(sum, term); term = wr_cmul(term, q); }
    wr_rt out;
    out.r = wr_cadd(wr_c(r12, 0.0), wr_cmul(wr_cscale(e, t12 * t21 * r23), sum));
    out.t = wr_cmul(wr_cscale(half, t12 * t23), sum);
    return out;
}

wr_rt wr_layer_airy(double n1, double n2, double n3, double d, double lambda0)
{
    double r12 = wr_fresnel_r_normal(n1, n2), r23 = wr_fresnel_r_normal(n2, n3);
    double t12 = wr_fresnel_t_normal(n1, n2), t23 = wr_fresnel_t_normal(n2, n3);
    wr_cplx half = one_way(n2, d, lambda0), e = wr_cmul(half, half);
    wr_cplx den = wr_cadd(wr_c(1.0, 0.0), wr_cscale(e, r12 * r23));
    wr_rt out;
    out.r = wr_cdiv(wr_cadd(wr_c(r12, 0.0), wr_cscale(e, r23)), den);
    out.t = wr_cdiv(wr_cscale(half, t12 * t23), den);
    return out;
}

double wr_layer_reflectance(wr_rt x) { return wr_cabs2(x.r); }
double wr_layer_transmittance(wr_rt x, double n1, double n3) { return (n3 / n1) * wr_cabs2(x.t); }
