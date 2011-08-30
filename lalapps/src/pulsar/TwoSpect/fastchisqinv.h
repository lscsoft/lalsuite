
#ifndef __FASTCHISQINV_H__
#define __FASTCHISQINV_H__

#include <lal/LALStdlib.h>

struct cheb_series_struct {
   double * c;   /* coefficients                */
   int order;    /* order of expansion          */
   double a;     /* lower interval point        */
   double b;     /* upper interval point        */
   int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;

REAL8 cdf_chisq_Pinv(REAL8 P, REAL8 nu);
REAL8 cdf_chisq_Qinv(REAL8 Q, REAL8 nu);
REAL8 cdf_gamma_Pinv(REAL8 P, REAL8 a, REAL8 b);
REAL8 cdf_gamma_Qinv(REAL8 Q, REAL8 a, REAL8 b);
REAL8 cdf_ugaussian_Pinv(REAL8 P);
REAL8 cdf_ugaussian_Qinv(REAL8 Q);
REAL8 small(REAL8 q);
REAL8 intermediate(REAL8 r);
REAL8 tail(REAL8 r);
REAL8 rat_eval(const REAL8 a[], const size_t na, const REAL8 b[], const size_t nb, const REAL8 x);
REAL8 cdf_gamma_P(REAL8 x, REAL8 a, REAL8 b);
REAL8 cdf_gamma_Q(REAL8 x, REAL8 a, REAL8 b);
REAL8 ran_gamma_pdf(REAL8 x, REAL8 a, REAL8 b);
REAL8 sf_gamma_inc_P(REAL8 a, REAL8 x);
REAL8 sf_gamma_inc_Q(REAL8 a, REAL8 x);
REAL8 gamma_inc_P_series(REAL8 a, REAL8 x);
REAL8 gamma_inc_Q_series(REAL8 a, REAL8 x);
REAL8 gamma_inc_D(REAL8 a, REAL8 x);
REAL8 twospect_sf_gammastar(REAL8 x);
REAL8 twospect_cheb_eval(const cheb_series * cs, REAL8 x);
REAL8 gammastar_ser(REAL8 x);
REAL8 sf_exprel_n_CF(REAL8 N, REAL8 x);
REAL8 gamma_inc_Q_asymp_unif(REAL8 a, REAL8 x);
REAL8 gamma_inc_Q_CF(REAL8 a, REAL8 x);
REAL8 gamma_inc_F_CF(REAL8 a, REAL8 x);
REAL8 gamma_inc_Q_large_x(REAL8 a, REAL8 x);



#endif
