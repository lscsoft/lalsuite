/*
 *  Copyright (C) 2011 Evan Goetz
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

//Based on GSL functions to determine chi-squared inversions
//Some functions based from Matab 2012a functions, but optimized for TwoSpect analysis

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>
#include <lal/LALConstants.h>

#include "TwoSpectSpecFunc.h"

/**
 * Compute the Dirichlet kernel for large N values
 * \param [in] delta The delta value as the arguement
 * \return Complex valued Dirichlet kernel
 */
COMPLEX16 DirichletKernelLargeN(const REAL8 delta)
{
   if (fabs(delta)<1.0e-6) return crect(1.0, 0.0);

   COMPLEX16 val = crect(0.0, LAL_TWOPI*delta);
   return (cexp(val)-1.0)/val;
}

/**
 * Compute the Dirichlet kernel for large N values and the Hann window
 * \param [in] delta The delta value as the arguement
 * \return Complex valued Dirichlet kernel
 */
COMPLEX16 DirichletKernelLargeNHann(const REAL8 delta)
{
   //return DirichletKernelLargeN(delta) - 0.5*DirichletKernelLargeN(delta+1.0) - 0.5*DirichletKernelLargeN(delta-1.0);
   if (fabs(delta)<1.0e-6) return crect(0.5, 0.0);
   else if (fabs(delta*delta - 1.0)<1.0e-6) return crect(-0.25, 0.0);
   else return crect(0.0, 1.0)*(cpolar(1.0, LAL_TWOPI*delta) - 1.0)/(2.0*LAL_TWOPI*delta*(delta*delta - 1.0));
}

/**
 * Compute the argument of the ratio of Dirichlet kernels for large N values and the Hann window
 * \param [out] ratio   The complex ratio
 * \param [in]  delta0  The delta0 value (denominator)
 * \param [in]  delta1  The delta1 value (numerator)
 * \param [in]  scaling The real-valued scaling
 * \return Error code: 0 = no error, 1 = divergent value
 */
INT4 DirichletKernalLargeNHannRatio(COMPLEX8 *ratio, const REAL4 delta0, const REAL4 delta1, const REAL4 scaling)
{
   //if (fabsf((REAL4)delta1)<(REAL4)1.0e-6 || fabsf((REAL4)(delta1*delta1-1.0))<(REAL4)1.0e-6 || fabsf((REAL4)(delta0-roundf(delta0)))<(REAL4)1.0e-6) return 1;
   //else *ratio = scaling*(cpolarf(1.0, (REAL4)(LAL_TWOPI*delta1))-1.0)/(cpolarf(1.0, (REAL4)(LAL_TWOPI*delta0))-1.0);
   *ratio = crectf(0.0, 0.0);
   if (fabsf(delta1)<(REAL4)1.0e-6) {
      if (fabsf(delta0)<(REAL4)1.0e-6) *ratio = crectf(1.0, 0.0);
      else if (fabsf((REAL4)(delta0*delta0 - 1.0))<(REAL4)1.0e-6) *ratio = crectf(-2.0, 0.0);
      else if (fabsf(delta0 - roundf(delta0))<(REAL4)1.0e-6) return 1;
      else *ratio = 0.5/(crectf(0.0, 1.0)*(cpolarf(1.0, LAL_TWOPI*delta0) - 1.0)/(2.0*LAL_TWOPI*delta0*(delta0*delta0 - 1.0)));
   } else if (fabsf((REAL4)(delta1*delta1 - 1.0))<(REAL4)1.0e-6) {
      if (fabsf(delta0)<(REAL4)1.0e-6) *ratio = crectf(-0.5, 0.0);
      else if (fabsf((REAL4)(delta0*delta0 - 1.0))<(REAL4)1.0e-6) *ratio = crectf(1.0, 0.0);
      else if (fabsf(delta0 - roundf(delta0))<(REAL4)1.0e-6) return 1;
      else *ratio = -0.25/(crectf(0.0, 1.0)*(cpolarf(1.0, LAL_TWOPI*delta0) - 1.0)/(2.0*LAL_TWOPI*delta0*(delta0*delta0 - 1.0)));
   } else if (fabsf(delta0 - roundf(delta0))<(REAL4)1.0e-6) return 1;
   else *ratio = scaling*sinf((REAL4)LAL_PI*delta1)/sinf((REAL4)LAL_PI*delta0)*cpolarf(1.0, (REAL4)(LAL_PI*(delta1 - delta0)));
   return 0;
}

/* Chebyshev coefficients for Gamma*(3/4(t+1)+1/2), -1<t<1
 */
static REAL8 gstar_a_data[30] = {
   2.16786447866463034423060819465,
   -0.05533249018745584258035832802,
   0.01800392431460719960888319748,
   -0.00580919269468937714480019814,
   0.00186523689488400339978881560,
   -0.00059746524113955531852595159,
   0.00019125169907783353925426722,
   -0.00006124996546944685735909697,
   0.00001963889633130842586440945,
   -6.3067741254637180272515795142e-06,
   2.0288698405861392526872789863e-06,
   -6.5384896660838465981983750582e-07,
   2.1108698058908865476480734911e-07,
   -6.8260714912274941677892994580e-08,
   2.2108560875880560555583978510e-08,
   -7.1710331930255456643627187187e-09,
   2.3290892983985406754602564745e-09,
   -7.5740371598505586754890405359e-10,
   2.4658267222594334398525312084e-10,
   -8.0362243171659883803428749516e-11,
   2.6215616826341594653521346229e-11,
   -8.5596155025948750540420068109e-12,
   2.7970831499487963614315315444e-12,
   -9.1471771211886202805502562414e-13,
   2.9934720198063397094916415927e-13,
   -9.8026575909753445931073620469e-14,
   3.2116773667767153777571410671e-14,
   -1.0518035333878147029650507254e-14,
   3.4144405720185253938994854173e-15,
   -1.0115153943081187052322643819e-15
};
static cheb_series gstar_a_cs = {
   gstar_a_data,
   29,
   -1, 1,
   17
};
/* Chebyshev coefficients for
 * x^2(Gamma*(x) - 1 - 1/(12x)), x = 4(t+1)+2, -1 < t < 1
 */
static REAL8 gstar_b_data[] = {
   0.0057502277273114339831606096782,
   0.0004496689534965685038254147807,
   -0.0001672763153188717308905047405,
   0.0000615137014913154794776670946,
   -0.0000223726551711525016380862195,
   8.0507405356647954540694800545e-06,
   -2.8671077107583395569766746448e-06,
   1.0106727053742747568362254106e-06,
   -3.5265558477595061262310873482e-07,
   1.2179216046419401193247254591e-07,
   -4.1619640180795366971160162267e-08,
   1.4066283500795206892487241294e-08,
   -4.6982570380537099016106141654e-09,
   1.5491248664620612686423108936e-09,
   -5.0340936319394885789686867772e-10,
   1.6084448673736032249959475006e-10,
   -5.0349733196835456497619787559e-11,
   1.5357154939762136997591808461e-11,
   -4.5233809655775649997667176224e-12,
   1.2664429179254447281068538964e-12,
   -3.2648287937449326771785041692e-13,
   7.1528272726086133795579071407e-14,
   -9.4831735252566034505739531258e-15,
   -2.3124001991413207293120906691e-15,
   2.8406613277170391482590129474e-15,
   -1.7245370321618816421281770927e-15,
   8.6507923128671112154695006592e-16,
   -3.9506563665427555895391869919e-16,
   1.6779342132074761078792361165e-16,
   -6.0483153034414765129837716260e-17
};
static cheb_series gstar_b_cs = {
   gstar_b_data,
   29,
   -1, 1,
   18
};


REAL8 twospect_small(REAL8 q)
{
   const REAL8 a[8] = { 3.387132872796366608, 133.14166789178437745,
      1971.5909503065514427, 13731.693765509461125,
      45921.953931549871457, 67265.770927008700853,
      33430.575583588128105, 2509.0809287301226727
   };

   const REAL8 b[8] = { 1.0, 42.313330701600911252,
      687.1870074920579083, 5394.1960214247511077,
      21213.794301586595867, 39307.89580009271061,
      28729.085735721942674, 5226.495278852854561
   };

   REAL8 r = 0.180625 - q * q;

   REAL8 x = q * rat_eval(a, 8, b, 8, r);

   return x;
}
REAL8 twospect_intermediate(REAL8 r)
{
   const REAL8 a[] = { 1.42343711074968357734, 4.6303378461565452959,
      5.7694972214606914055, 3.64784832476320460504,
      1.27045825245236838258, 0.24178072517745061177,
      0.0227238449892691845833, 7.7454501427834140764e-4
   };

   const REAL8 b[] = { 1.0, 2.05319162663775882187,
      1.6763848301838038494, 0.68976733498510000455,
      0.14810397642748007459, 0.0151986665636164571966,
      5.475938084995344946e-4, 1.05075007164441684324e-9
   };

   REAL8 x = rat_eval(a, 8, b, 8, (r - 1.6));

   return x;
}
REAL8 twospect_tail(REAL8 r)
{
   const REAL8 a[] = { 6.6579046435011037772, 5.4637849111641143699,
      1.7848265399172913358, 0.29656057182850489123,
      0.026532189526576123093, 0.0012426609473880784386,
      2.71155556874348757815e-5, 2.01033439929228813265e-7
   };

   const REAL8 b[] = { 1.0, 0.59983220655588793769,
      0.13692988092273580531, 0.0148753612908506148525,
      7.868691311456132591e-4, 1.8463183175100546818e-5,
      1.4215117583164458887e-7, 2.04426310338993978564e-15
   };

   REAL8 x = rat_eval(a, 8, b, 8, (r - 5.0));

   return x;
}
REAL8 rat_eval(const REAL8 a[], const size_t na, const REAL8 b[], const size_t nb, const REAL8 x)
{
   size_t i, j;
   REAL8 u, v, r;

   u = a[na - 1];

   for (i = na - 1; i > 0; i--) u = x * u + a[i - 1];

   v = b[nb - 1];

   for (j = nb - 1; j > 0; j--) v = x * v + b[j - 1];

   r = u / v;

   return r;
}

REAL8 ran_gamma_pdf(REAL8 x, REAL8 a, REAL8 b)
{
   if (x < 0.0) return 0.0;
   else if (x == 0.0) {
      if (a == 1.0) return 1.0/b ;
      else return 0.0;
   } else if (a == 1.0) return (exp(-x/b)/b);
   else return (exp((a - 1.0) * log(x/b) - x/b - lgamma(a))/b);
}
REAL8 matlab_gamma_inc(REAL8 x, REAL8 a, INT4 upper)
{
   const REAL8 amax = 1048576.0;
   const REAL8 amaxthird = amax-1.0/3.0;
   REAL8 xint = x;
   REAL8 aint = a;
   if (aint>amax) {
      xint = fmax(amaxthird+sqrt(amax/a)*(xint-(aint-1.0/3.0)), 0.0);
      aint = amax;
   }
   if (aint==0.0) {
      if (upper==0) return 1.0;
      else return 0.0;
   } else if (xint==0.0) {
      if (upper==0) return 0.0;
      else return 1.0;
   } else if (xint<aint+1.0) {
      REAL8 ap = aint;
      REAL8 del = 1.0;
      REAL8 sum = del;
      while (fabs(del)>=100.0*epsval(fabs(sum))) {
         ap += 1.0;
         del = xint*del/ap;
         sum += del;
      }
      REAL8 b = sum * exp(-xint + aint*log(xint) - lgamma(aint+1.0));
      if (xint>0.0 && b>1.0) b = 1.0;
      if (upper==0) return b;
      else return 1.0-b;
   } else {
      REAL8 a0 = 1.0;
      REAL8 a1 = xint;
      REAL8 b0 = 0.0;
      REAL8 b1 = a0;
      REAL8 fac = 1.0/a1;
      INT8 n = 1;
      REAL8 g = b1 * fac;
      REAL8 gold = b0;
      while (fabs(g-gold) >= 100*epsval(fabs(g))) {
         gold = g;
         REAL8 ana = n - a;
         a0 = (a1 + a0 * ana) * fac;
         b0 = (b1 + b0 * ana) * fac;
         REAL8 anf = n*fac;
         a1 = xint * a0 + anf * a1;
         b1 = xint * b0 + anf * b1;
         fac = 1.0/a1;
         g = b1 * fac;
         n++;
      }
      REAL8 b = exp(-xint + aint*log(xint) - lgamma(aint)) * g;
      if (upper==0) return 1.0-b;
      else return b;
   }

}
INT4 sf_gamma_inc_P(REAL8 *out, REAL8 a, REAL8 x)
{
   XLAL_CHECK( a > 0.0 && x >= 0.0, XLAL_EINVAL, "Invalid input of zero (a = %f), or less than zero (x = %f)\n", a, x );

   if(x == 0.0) *out = 0.0;
   else if(x < 20.0 || x < 0.5*a) {
      REAL8 val;
      XLAL_CHECK( gamma_inc_P_series(&val, a, x) == XLAL_SUCCESS, XLAL_EFUNC );
      *out = val;
   } else if(a > 1.0e+06 && (x-a)*(x-a) < a) {
      /* Crossover region. Note that Q and P are
       * roughly the same order of magnitude here,
       * so the subtraction is stable.
       */
      *out = 1.0 - gamma_inc_Q_asymp_unif(a, x);
   } else if(a <= x) {
      /* Q <~ P in this area, so the
       * subtractions are stable.
       */
      REAL8 Q;
      if (a > 0.2*x) {
         XLAL_CHECK( gamma_inc_Q_CF(&Q, a, x) == XLAL_SUCCESS, XLAL_EFUNC );
      } else {
         Q = gamma_inc_Q_large_x(a, x);
         XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
      }
      *out = 1.0 - Q;
   } else {
      if ((x-a)*(x-a) < a) {
         /* This condition is meant to insure
          * that Q is not very close to 1,
          * so the subtraction is stable.
          */
         REAL8 Q;
         XLAL_CHECK( gamma_inc_Q_CF(&Q, a, x) == XLAL_SUCCESS, XLAL_EFUNC );
         *out = 1.0 - Q;
      } else {
         REAL8 val;
         XLAL_CHECK( gamma_inc_P_series(&val, a, x) == XLAL_SUCCESS, XLAL_EFUNC );
         *out = val;
      }
   }
   return XLAL_SUCCESS;
}
INT4 sf_gamma_inc_Q(REAL8 *out, REAL8 a, REAL8 x)
{
   XLAL_CHECK( a >= 0.0 && x >= 0.0, XLAL_EINVAL, "Invalid input of less than zero (a = %f), or less than zero (x = %f)\n", a, x );
   if (x == 0.0) *out = 1.0;
   else if (a == 0.0) *out = 0.0;
   else if(x <= 0.5*a) {
      /* If the series is quick, do that. It is
       * robust and simple.
       */
      REAL8 P;
      XLAL_CHECK( gamma_inc_P_series(&P, a, x) == XLAL_SUCCESS, XLAL_EFUNC );
      *out = 1.0 - P;
   } else if(a >= 1.0e+06 && (x-a)*(x-a) < a) {
      /* Then try the difficult asymptotic regime.
       * This is the only way to do this region.
       */
      *out = gamma_inc_Q_asymp_unif(a, x);
   } else if(a < 0.2 && x < 5.0) {
      /* Cancellations at small a must be handled
       * analytically; x should not be too big
       * either since the series terms grow
       * with x and log(x).
       */
      REAL8 val;
      XLAL_CHECK( gamma_inc_Q_series(&val, a, x) == XLAL_SUCCESS, XLAL_EFUNC );
      *out = val;
   } else if(a <= x) {
      if(x <= 1.0e+06) {
         /* Continued fraction is excellent for x >~ a.
          * We do not let x be too large when x > a since
          * it is somewhat pointless to try this there;
          * the function is rapidly decreasing for
          * x large and x > a, and it will just
          * underflow in that region anyway. We
          * catch that case in the standard
          * large-x method.
          */
         REAL8 val;
         XLAL_CHECK( gamma_inc_Q_CF(&val, a, x) == XLAL_SUCCESS, XLAL_EFUNC );
         *out = val;
      } else {
         REAL8 val = gamma_inc_Q_large_x(a, x);
         XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
         *out = val;
      }
   } else {
      if(x > a - sqrt(a)) {
         /* Continued fraction again. The convergence
          * is a little slower here, but that is fine.
          * We have to trade that off against the slow
          * convergence of the series, which is the
          * only other option.
          */
         REAL8 val;
         XLAL_CHECK( gamma_inc_Q_CF(&val, a, x) == XLAL_SUCCESS, XLAL_EFUNC );
         *out = val;
      } else {
         REAL8 P;
         XLAL_CHECK( gamma_inc_P_series(&P, a, x) == XLAL_SUCCESS, XLAL_EFUNC );
         *out = 1.0 - P;
      }
   }
   return XLAL_SUCCESS;
}
INT4 gamma_inc_P_series(REAL8 *out, REAL8 a, REAL8 x)
{

   INT4 nmax = 10000;

   REAL8 D;
   XLAL_CHECK( gamma_inc_D(&D, a, x) == XLAL_SUCCESS, XLAL_EFUNC );

   /* Approximating the terms of the series using Stirling's
    approximation gives t_n = (x/a)^n * exp(-n(n+1)/(2a)), so the
    convergence condition is n^2 / (2a) + (1-(x/a) + (1/2a)) n >>
    -log(GSL_DBL_EPS) if we want t_n < O(1e-16) t_0. The condition
    below detects cases where the minimum value of n is > 5000 */

   if (x > 0.995 * a && a > 1.0e5) { /* Difficult case: try continued fraction */
      REAL8 cf_res = sf_exprel_n_CF(a, x);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
      *out = D * cf_res;
      return XLAL_SUCCESS;
   }

   /* Series would require excessive number of terms */

   XLAL_CHECK( x <= (a + nmax), XLAL_EMAXITER, "gamma_inc_P_series x>>a exceeds range\n" );

   /* Normal case: sum the series */
   REAL8 sum  = 1.0;
   REAL8 term = 1.0;
   REAL8 remainderval;
   INT4 n;

   /* Handle lower part of the series where t_n is increasing, |x| > a+n */

   INT4 nlow = (x > a) ? (x - a): 0;

   for (n=1; n < nlow; n++) {
      term *= x/(a+n);
      sum  += term;
   }

   /* Handle upper part of the series where t_n is decreasing, |x| < a+n */

   for (/* n = previous n */ ; n<nmax; n++)  {
      term *= x/(a+n);
      sum  += term;
      if (fabs(term/sum) < LAL_REAL4_EPS) break;
   }

   /*  Estimate remainder of series ~ t_(n+1)/(1-x/(a+n+1)) */
   REAL8 tnp1 = (x/(a+n)) * term;
   remainderval =  tnp1 / (1.0 - x/(a + n + 1.0));

   REAL8 val = D * sum;

   XLAL_CHECK( !(n == nmax && fabs(remainderval/sum) > sqrt(LAL_REAL4_EPS)), XLAL_EMAXITER, "gamma_inc_P_series_float failed to converge\n" );

   *out = val;
   return XLAL_SUCCESS;

}
INT4 gamma_inc_Q_series(REAL8 *out, REAL8 a, REAL8 x)
{
   REAL8 term1;  /* 1 - x^a/Gamma(a+1) */
   REAL8 sum;    /* 1 + (a+1)/(a+2)(-x)/2! + (a+1)/(a+3)(-x)^2/3! + ... */
   REAL8 term2;  /* a temporary variable used at the end */

   {
      /* Evaluate series for 1 - x^a/Gamma(a+1), small a
       */
      const REAL8 pg21 = -2.404113806319188570799476;  /* PolyGamma[2,1] */
      const REAL8 lnx  = log(x);
      const REAL8 el   = M_EULER+lnx;
      const REAL8 c1 = -el;
      const REAL8 c2 = M_PI*M_PI/12.0 - 0.5*el*el;
      const REAL8 c3 = el*(M_PI*M_PI/12.0 - el*el/6.0) + pg21/6.0;
      const REAL8 c4 = -0.04166666666666666667
      * (-1.758243446661483480 + lnx)
      * (-0.764428657272716373 + lnx)
      * ( 0.723980571623507657 + lnx)
      * ( 4.107554191916823640 + lnx);
      const REAL8 c5 = -0.0083333333333333333
      * (-2.06563396085715900 + lnx)
      * (-1.28459889470864700 + lnx)
      * (-0.27583535756454143 + lnx)
      * ( 1.33677371336239618 + lnx)
      * ( 5.17537282427561550 + lnx);
      const REAL8 c6 = -0.0013888888888888889
      * (-2.30814336454783200 + lnx)
      * (-1.65846557706987300 + lnx)
      * (-0.88768082560020400 + lnx)
      * ( 0.17043847751371778 + lnx)
      * ( 1.92135970115863890 + lnx)
      * ( 6.22578557795474900 + lnx);
      const REAL8 c7 = -0.00019841269841269841
      * (-2.5078657901291800 + lnx)
      * (-1.9478900888958200 + lnx)
      * (-1.3194837322612730 + lnx)
      * (-0.5281322700249279 + lnx)
      * ( 0.5913834939078759 + lnx)
      * ( 2.4876819633378140 + lnx)
      * ( 7.2648160783762400 + lnx);
      const REAL8 c8 = -0.00002480158730158730
      * (-2.677341544966400 + lnx)
      * (-2.182810448271700 + lnx)
      * (-1.649350342277400 + lnx)
      * (-1.014099048290790 + lnx)
      * (-0.191366955370652 + lnx)
      * ( 0.995403817918724 + lnx)
      * ( 3.041323283529310 + lnx)
      * ( 8.295966556941250 + lnx);
      const REAL8 c9 = -2.75573192239859e-6
      * (-2.8243487670469080 + lnx)
      * (-2.3798494322701120 + lnx)
      * (-1.9143674728689960 + lnx)
      * (-1.3814529102920370 + lnx)
      * (-0.7294312810261694 + lnx)
      * ( 0.1299079285269565 + lnx)
      * ( 1.3873333251885240 + lnx)
      * ( 3.5857258865210760 + lnx)
      * ( 9.3214237073814600 + lnx);
      const REAL8 c10 = -2.75573192239859e-7
      * (-2.9540329644556910 + lnx)
      * (-2.5491366926991850 + lnx)
      * (-2.1348279229279880 + lnx)
      * (-1.6741881076349450 + lnx)
      * (-1.1325949616098420 + lnx)
      * (-0.4590034650618494 + lnx)
      * ( 0.4399352987435699 + lnx)
      * ( 1.7702236517651670 + lnx)
      * ( 4.1231539047474080 + lnx)
      * ( 10.342627908148680 + lnx);

      term1 = a*(c1+a*(c2+a*(c3+a*(c4+a*(c5+a*(c6+a*(c7+a*(c8+a*(c9+a*c10)))))))));
   }

   {
      /* Evaluate the sum.
       */
      const INT4 nmax = 5000;
      REAL8 t = 1.0;
      INT4 n;
      sum = 1.0;

      for (n=1; n<nmax; n++) {
         t *= -x/(n+1.0);
         sum += (a+1.0)/(a+n+1.0)*t;
         if(fabs(t/sum) < LAL_REAL4_EPS) break;
      }

      XLAL_CHECK( n != nmax, XLAL_EMAXITER, "Maximum iterations reached\n" );
   }

   term2 = (1.0 - term1) * a/(a+1.0) * x * sum;
   *out = (term1 + term2);
   return XLAL_SUCCESS;

}
REAL8 twospect_cheb_eval(const cheb_series * cs, REAL8 x)
{

   INT4 j;
   REAL8 d  = 0.0;
   REAL8 dd = 0.0;

   REAL8 y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
   REAL8 y2 = 2.0 * y;

   for (j = cs->order; j>=1; j--) {
      REAL8 temp = d;
      d = y2*d - dd + cs->c[j];
      dd = temp;
   }

   {
      d = y*d - dd + 0.5 * cs->c[0];
   }

   return d;

}
INT4 gamma_inc_D(REAL8 *out, REAL8 a, REAL8 x)
{

   if (a < 10.0) {
      REAL8 lnr = a * log(x) - x - lgamma(a+1.0);
      *out = exp(lnr);
   } else {
      REAL8 gstar;
      REAL8 ln_term;
      REAL8 term1;
      if (x < 0.5*a) {
         REAL8 u = x/a;
         REAL8 ln_u = log(u);
         ln_term = ln_u - u + 1.0;
      } else {
         REAL8 mu = (x-a)/a;
         //ln_term = gsl_sf_log_1plusx_mx(mu);  /* log(1+mu) - mu */
         ln_term = log1p(mu) - mu;  /* log(1+mu) - mu */
      }
      XLAL_CHECK( twospect_sf_gammastar(&gstar, a) == XLAL_SUCCESS, XLAL_EFUNC );
      term1 = exp(a*ln_term)/sqrt(2.0*LAL_PI*a);
      *out = term1/gstar;
   }
   return XLAL_SUCCESS;

}
INT4 twospect_sf_gammastar(REAL8 *out, REAL8 x)
{
   XLAL_CHECK( x > 0.0, XLAL_EINVAL, "Invalid input of zero or less: %f\n", x );

   if (x < 0.5) {
      REAL8 lg = lgamma(x);
      REAL8 lx = log(x);
      REAL8 c  = 0.5*(LAL_LN2+M_LNPI);
      REAL8 lnr_val = lg - (x-0.5)*lx + x - c;
      *out = exp(lnr_val);
   } else if(x < 2.0) {
      REAL8 t = 4.0/3.0*(x-0.5) - 1.0;
      REAL8 val = twospect_cheb_eval(&gstar_a_cs, t);
      *out = val;
   } else if(x < 10.0) {
      REAL8 t = 0.25*(x-2.0) - 1.0;
      REAL8 c = twospect_cheb_eval(&gstar_b_cs, t);
      *out = c/(x*x) + 1.0 + 1.0/(12.0*x);
   } else if(x < 1.0/1.2207031250000000e-04) {
      REAL8 val = gammastar_ser(x);
      *out = val;
   } else if(x < 1.0/LAL_REAL8_EPS) {
      /* Use Stirling formula for Gamma(x).
       */
      REAL8 xi = 1.0/x;
      *out = 1.0 + xi/12.0*(1.0 + xi/24.0*(1.0 - xi*(139.0/180.0 + 571.0/8640.0*xi)));
   } else *out = 1.0;

   return XLAL_SUCCESS;
}
REAL8 gammastar_ser(REAL8 x)
{
   /* Use the Stirling series for the correction to Log(Gamma(x)),
    * which is better behaved and easier to compute than the
    * regular Stirling series for Gamma(x).
    */
   const REAL8 y = 1.0/(x*x);
   const REAL8 c0 =  1.0/12.0;
   const REAL8 c1 = -1.0/360.0;
   const REAL8 c2 =  1.0/1260.0;
   const REAL8 c3 = -1.0/1680.0;
   const REAL8 c4 =  1.0/1188.0;
   const REAL8 c5 = -691.0/360360.0;
   const REAL8 c6 =  1.0/156.0;
   const REAL8 c7 = -3617.0/122400.0;
   const REAL8 ser = c0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*c7))))));
   return exp(ser/x);
}
REAL8 sf_exprel_n_CF(REAL8 N, REAL8 x)
{
   const REAL8 RECUR_BIG = sqrt(LAL_REAL8_MAX);
   INT4 maxiter = 5000;
   INT4 n = 1;
   REAL8 Anm2 = 1.0;
   REAL8 Bnm2 = 0.0;
   REAL8 Anm1 = 0.0;
   REAL8 Bnm1 = 1.0;
   REAL8 a1 = 1.0;
   REAL8 b1 = 1.0;
   REAL8 a2 = -x;
   REAL8 b2 = N+1;
   REAL8 an, bn;

   REAL8 fn;

   REAL8 An = b1*Anm1 + a1*Anm2;   /* A1 */
   REAL8 Bn = b1*Bnm1 + a1*Bnm2;   /* B1 */

   /* One explicit step, before we get to the main pattern. */
   n++;
   Anm2 = Anm1;
   Bnm2 = Bnm1;
   Anm1 = An;
   Bnm1 = Bn;
   An = b2*Anm1 + a2*Anm2;   /* A2 */
   Bn = b2*Bnm1 + a2*Bnm2;   /* B2 */

   fn = An/Bn;

   while(n < maxiter) {
      REAL8 old_fn;
      REAL8 del;
      n++;
      Anm2 = Anm1;
      Bnm2 = Bnm1;
      Anm1 = An;
      Bnm1 = Bn;
      an = ( GSL_IS_ODD(n) ? ((n-1)/2)*x : -(N+(n/2)-1)*x );
      bn = N + n - 1;
      An = bn*Anm1 + an*Anm2;
      Bn = bn*Bnm1 + an*Bnm2;

      if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
         An /= RECUR_BIG;
         Bn /= RECUR_BIG;
         Anm1 /= RECUR_BIG;
         Bnm1 /= RECUR_BIG;
         Anm2 /= RECUR_BIG;
         Bnm2 /= RECUR_BIG;
      }

      old_fn = fn;
      fn = An/Bn;
      del = old_fn/fn;

      if (fabs(del - 1.0) < 2.0*LAL_REAL4_EPS) break;
   }

   XLAL_CHECK_REAL8( n < maxiter, XLAL_EMAXITER, "Reached maximum number of iterations (5000)\n" );

   return fn;

}
REAL8 gamma_inc_Q_asymp_unif(REAL8 a, REAL8 x)
{

   REAL8 rta = sqrt(a);
   REAL8 eps = (x-a)/a;

   REAL8 ln_term = gsl_sf_log_1plusx_mx(eps);  /* log(1+eps) - eps */
   REAL8 eta  = GSL_SIGN(eps) * sqrt(-2.0*ln_term);

   REAL8 R;
   REAL8 c0, c1;

   REAL8 erfcval = erfc(eta*rta/LAL_SQRT2);

   if(fabs(eps) < 7.4009597974140505e-04) {
      c0 = -1.0/3.0 + eps*(1.0/12.0 - eps*(23.0/540.0 - eps*(353.0/12960.0 - eps*589.0/30240.0)));
      c1 = -1.0/540.0 - eps/288.0;
   } else {
      REAL8 rt_term = sqrt(-2.0 * ln_term/(eps*eps));
      REAL8 lam = x/a;
      c0 = (1.0 - 1.0/rt_term)/eps;
      c1 = -(eta*eta*eta * (lam*lam + 10.0*lam + 1.0) - 12.0 * eps*eps*eps) / (12.0 * eta*eta*eta*eps*eps*eps);
   }

   R = exp(-0.5*a*eta*eta)/(LAL_SQRT2*M_SQRTPI*rta) * (c0 + c1/a);

   return (0.5 * erfcval + R);

}
INT4 gamma_inc_Q_CF(REAL8 *out, REAL8 a, REAL8 x)
{
   REAL8 D, F;
   XLAL_CHECK( gamma_inc_D(&D, a, x) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( gamma_inc_F_CF(&F, a, x) == XLAL_SUCCESS, XLAL_EFUNC );

   *out = (D * (a/x) * F);

   return XLAL_SUCCESS;

}
INT4 gamma_inc_F_CF(REAL8 *out, REAL8 a, REAL8 x)
{
   INT4 nmax  =  5000;
   const REAL8 smallval =  LAL_REAL8_EPS*LAL_REAL8_EPS*LAL_REAL8_EPS;

   REAL8 hn = 1.0;           /* convergent */
   REAL8 Cn = 1.0 / smallval;
   REAL8 Dn = 1.0;
   INT4 n;

   /* n == 1 has a_1, b_1, b_0 independent of a,x,
    so that has been done by hand                */
   for ( n = 2 ; n < nmax ; n++ ) {
      REAL8 an;
      REAL8 delta;

      if (GSL_IS_ODD(n)) an = 0.5*(n-1)/x;
      else an = (0.5*n-a)/x;

      Dn = 1.0 + an * Dn;
      if ( fabs(Dn) < smallval ) Dn = smallval;
      Cn = 1.0 + an/Cn;
      if ( fabs(Cn) < smallval ) Cn = smallval;
      Dn = 1.0 / Dn;
      delta = Cn * Dn;
      hn *= delta;
      if (fabs(delta-1.0) < LAL_REAL4_EPS) break;
   }

   XLAL_CHECK( n < nmax, XLAL_EMAXITER );

   *out = hn;

   return XLAL_SUCCESS;

}
REAL8 gamma_inc_Q_large_x(REAL8 a, REAL8 x)
{
   const INT4 nmax = 100000;

   REAL8 D;
   XLAL_CHECK_REAL8( gamma_inc_D(&D, a, x) == XLAL_SUCCESS, XLAL_EFUNC );

   REAL8 sum  = 1.0;
   REAL8 term = 1.0;
   REAL8 last = 1.0;
   INT4 n;
   for(n=1; n<nmax; n++) {
      term *= (a-n)/x;
      if(fabs(term/last) > 1.0) break;
      if(fabs(term/sum) < LAL_REAL4_EPS) break;
      sum  += term;
      last = term;
   }

   XLAL_CHECK_REAL8( n < nmax, XLAL_EMAXITER );

   return (D * (a/x) * sum);

}
//Matlab's eps function for REAL8, but written in C
REAL8 epsval(REAL8 val)
{
   //Same as matlab
   REAL8 absval = fabs(val);
   int exponentval = 0;
   frexp(absval, &exponentval);
   exponentval -= LAL_REAL8_MANT;
   return ldexp(1.0, exponentval);
} /* epsval() */

//Matlab's eps function for REAL4, but written in C
REAL4 epsval_float(REAL4 val)
{
   //Same as matlab
   REAL4 absval = fabsf(val);
   int exponentval = 0;
   frexpf(absval, &exponentval);
   exponentval -= LAL_REAL4_MANT;
   return ldexpf(1.0, exponentval);
} /* epsval_float() */
//Matlab's sumseries function
void sumseries(REAL8 *computedprob, REAL8 P, REAL8 C, REAL8 E, INT8 counter, REAL8 x, REAL8 dof, REAL8 halfdelta, REAL8 err, INT4 countdown)
{

   //Exit with error if halfdelta = 0.0
   XLAL_CHECK_VOID( halfdelta != 0.0, XLAL_EINVAL );

   REAL8 Pint = P, Cint = C, Eint = E;
   INT8 counterint = counter;
   INT8 j = 0;
   if (countdown!=0) {
      if (counterint>=0) j = 1;
      if (j==1) {
         Pint *= (counterint+1.0)/halfdelta;
         Cint += E;
      } else counterint = -1;
   }

   while (counterint!=-1) {
      REAL8 pplus = Pint*Cint;
      *(computedprob) += pplus;

      if (pplus > *(computedprob)*err) j = 1;
      else j = 0;
      if (countdown!=0 && counterint<0) j = 0;
      if (j==0) return;

      if (countdown!=0) {
         counterint--;
         Pint *= (counterint+1.0)/halfdelta;
         Eint *= (0.5*dof + counterint+1.0)/(x*0.5);
         Cint += Eint;
      } else {
         counterint++;
         Pint *= halfdelta/counterint;
         Eint *= (0.5*x)/(0.5*dof+counterint-1.0);
         Cint -= Eint;
      }
   }

} /* sumseries() */


//Evan's sumseries function based on matlab's sumseries() version above, but faster
void sumseries_eg(REAL8 *computedprob, REAL8 P, REAL8 C, REAL8 E, INT8 counter, REAL8 x, REAL8 dof, REAL8 halfdelta, REAL8 err, INT4 countdown)
{

   //If halfdelta = 0.0, then exit with error
   XLAL_CHECK_VOID( halfdelta != 0.0, XLAL_EINVAL );

   REAL8 Pint = P, Cint = C, Eint = E;
   INT8 counterint = counter;
   REAL8 oneoverhalfdelta = 1.0/halfdelta;   //pre-compute
   REAL8 halfdof = 0.5*dof;                  //pre-compute
   REAL8 halfx = 0.5*x;                      //pre-compute

   if (countdown!=0) {
      if (counterint>=0) {
         Pint *= (counterint+1.0)*oneoverhalfdelta;
         Cint += E;
      } else counterint = -1;
   }

   if (counterint==-1) return;
   else if (countdown!=0) {
      REAL8 oneoverhalfx = 1.0/halfx;
      while (counterint!=-1) {
         REAL8 pplus = Pint*Cint;
         *(computedprob) += pplus;

         if (pplus<=*(computedprob)*err || counterint<0) return;

         counterint--;
         Pint *= (counterint+1)*oneoverhalfdelta;
         Eint *= (halfdof + counterint+1)*oneoverhalfx;
         Cint += Eint;
      }
   } else {
      while (counterint!=-1) {
         REAL8 pplus = Pint*Cint;
         *(computedprob) += pplus;

         if (pplus<=*(computedprob)*err) return;

         counterint++;
         Pint *= halfdelta/counterint;
         Eint *= halfx/(halfdof+counterint-1);
         Cint -= Eint;
      }
   }

} /* sumseries_eg() */

//Matlab's binodeviance, a "hidden" function
REAL8 binodeviance(REAL8 x, REAL8 np)
{
   if (fabs(x-np)<0.1*(x+np)) {
      REAL8 s = (x-np)*(x-np)/(x+np);
      REAL8 v = (x-np)/(x+np);
      REAL8 ej = 2.0*x*v;
      REAL8 s1 = 0.0;
      INT4 jj = 0;
      INT4 ok = 1;
      while (ok==1) {
         ej *= v*v;
         jj++;
         s1 = s + ej/(2.0*jj+1.0);
         if (s1!=s) s = s1;
         else ok = 0;
      }
      return s;
   } else {
      return x*log(x/np)+np-x;
   }
} /* binodeviance() */
