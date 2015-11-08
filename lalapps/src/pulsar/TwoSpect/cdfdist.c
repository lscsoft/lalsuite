/*
 *  Copyright (C) 2014 Evan Goetz
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
#include <lal/LALConstants.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include "cdfdist.h"
#include "TwoSpectSpecFunc.h"

REAL8 cdf_chisq_Pinv(REAL8 P, REAL8 nu)
{
   REAL8 val = cdf_gamma_Pinv(P, 0.5*nu, 2.0);
   XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
   return val;
}
REAL8 cdf_chisq_Qinv(REAL8 Q, REAL8 nu)
{
   REAL8 val;
   XLAL_CHECK_REAL8( cdf_gamma_Qinv(&val, Q, 0.5*nu, 2.0) == XLAL_SUCCESS, XLAL_EFUNC );
   return val;
}
REAL8 cdf_gamma_Pinv(REAL8 P, REAL8 a, REAL8 b)
{
   REAL8 x;

   XLAL_CHECK_REAL8( P < 1.0, XLAL_EFPOVRFLW, "Input P of 1.0 or larger returns infinity\n" );
   if (P == 0.0) return 0.0;

   /* Consider, small, large and intermediate cases separately.  The
    boundaries at 0.05 and 0.95 have not been optimised, but seem ok
    for an initial approximation.

    BJG: These approximations aren't really valid, the relevant
    criterion is P*gamma(a+1) < 1. Need to rework these routines and
    use a single bisection style solver for all the inverse
    functions.
    */

   if (P < 0.05) x = exp((lgamma(a) + log(P))/a);
   else if (P > 0.95) x = -log1p(-P) + lgamma(a);
   else {
      REAL8 xg = cdf_ugaussian_Pinv(P);
      XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
      x = (xg < -0.5*sqrt(a)) ? a : sqrt(a) * xg + a;
   }

   /* Use Lagrange's interpolation for E(x)/phi(x0) to work backwards
    to an improved value of x (Abramowitz & Stegun, 3.6.6)

    where E(x)=P-integ(phi(u),u,x0,x) and phi(u) is the pdf.
    */

   REAL8 lambda, dP, phi;
   UINT4 n = 0;

   INT4 keepgoing = 1;
   while (keepgoing == 1) {
      REAL8 val = cdf_gamma_P(x, a, 1.0);
      XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
      dP = P - val;
      phi = ran_gamma_pdf(x, a, 1.0);

      if (dP == 0.0 || n++ > 32) {
         XLAL_CHECK_REAL8( fabs(dP) <= sqrt(LAL_REAL4_EPS) * P, XLAL_EFPINEXCT, "Inverse failed to converge\n" );
         return b * x;
      }

      lambda = dP / fmax (2.0 * fabs (dP / x), phi);


      REAL8 step0 = lambda;
      REAL8 step1 = -((a - 1.0) / x - 1.0) * lambda * lambda / 4.0;

      REAL8 step = step0;
      if (fabs (step1) < 0.5 * fabs (step0)) step += step1;

      if (x + step > 0) x += step;
      else x *= 0.5;

      if (fabs (step0) > 1e-6 * x || fabs(step0 * phi) > 1e-6 * P) keepgoing = 1;
      else keepgoing = 0;
   }

   XLAL_CHECK_REAL8( fabs(dP) <= sqrt(LAL_REAL4_EPS) * P, XLAL_EFPINEXCT, "Inverse failed to converge\n" );
   return b * x;

}
INT4 cdf_gamma_Qinv(REAL8 *out, REAL8 Q, REAL8 a, REAL8 b)
{
   REAL8 x;

   if (Q == 1.0) {
      *out = 0.0;
      return XLAL_SUCCESS;
   }
   XLAL_CHECK( Q > 0.0 && Q < 1.0, XLAL_EFPOVRFLW, "Input P of 0.0 returns infinity\n" );

   /* Consider, small, large and intermediate cases separately.  The
    boundaries at 0.05 and 0.95 have not been optimised, but seem ok
    for an initial approximation. */

   if (Q < 0.05) x = -log(Q) + lgamma(a);
   else if (Q > 0.95) x = exp((lgamma(a) + log1p(-Q)) / a);
   else {
      REAL8 xg;
      XLAL_CHECK( cdf_ugaussian_Qinv(&xg, Q) == XLAL_SUCCESS, XLAL_EFUNC );
      x = (xg < -0.5*sqrt (a)) ? a : sqrt (a) * xg + a;
   }

   /* Use Lagrange's interpolation for E(x)/phi(x0) to work backwards
    to an improved value of x (Abramowitz & Stegun, 3.6.6)

    where E(x)=P-integ(phi(u),u,x0,x) and phi(u) is the pdf.
    */

   REAL8 lambda, dQ, phi;
   UINT4 n = 0;

   BOOLEAN keepgoing = 1;
   while (keepgoing == 1) {
      REAL8 val;
      XLAL_CHECK( cdf_gamma_Q(&val, x, a, 1.0) == XLAL_SUCCESS, XLAL_EFUNC );
      dQ = Q - val;
      phi = ran_gamma_pdf(x, a, 1.0);

      if (dQ == 0.0 || n++ > 32) {
         *out = b * x;
         return XLAL_SUCCESS;
      }

      lambda = -dQ / fmax (2 * fabs (dQ / x), phi);

      REAL8 step0 = lambda;
      REAL8 step1 = -((a - 1) / x - 1) * lambda * lambda / 4.0;

      REAL8 step = step0;
      if (fabs (step1) < 0.5 * fabs (step0)) step += step1;

      if (x + step > 0) x += step;
      else x /= 2.0;

      if (fabs (step0) > 1e-6 * x) keepgoing = 1;
      else keepgoing = 0;
   }

   *out = b * x;
   return XLAL_SUCCESS;
}
REAL8 cdf_ugaussian_Pinv(REAL8 P)
{
   REAL8 r, x, pp;

   REAL8 dP = P - 0.5;

   XLAL_CHECK_REAL8( P != 1.0, XLAL_EFPOVRFLW );
   XLAL_CHECK_REAL8( P != 0.0, XLAL_EFPOVRFLW );

   if (fabs(dP) <= 0.425) {
      x = twospect_small(dP);
      return x;
   }

   pp = (P < 0.5) ? P : 1.0 - P;

   r = sqrt(-log(pp));

   if (r <= 5.0) x = twospect_intermediate(r);
   else x = twospect_tail(r);

   if (P < 0.5) return -x;
   else return x;

}
INT4 cdf_ugaussian_Qinv(REAL8 *out, REAL8 Q)
{
   REAL8 r, x, pp;

   REAL8 dQ = Q - 0.5;

   XLAL_CHECK( Q != 1.0, XLAL_EFPOVRFLW );
   XLAL_CHECK( Q != 0.0, XLAL_EFPOVRFLW );

   if (fabs(dQ) <= 0.425) {
      x = twospect_small(dQ);
      *out = -x;
      return XLAL_SUCCESS;
   }

   pp = (Q < 0.5) ? Q : 1.0 - Q;

   r = sqrt(-log(pp));

   if (r <= 5.0) x = twospect_intermediate(r);
   else x = twospect_tail(r);

   if (Q < 0.5) *out = x;
   else *out = -x;
   return XLAL_SUCCESS;
}
REAL8 cdf_gamma_P(REAL8 x, REAL8 a, REAL8 b)
{
   REAL8 P;
   REAL8 y = x / b;

   if (x <= 0.0) return 0.0;

   if (y > a) {
      REAL8 val;
      XLAL_CHECK_REAL8( sf_gamma_inc_Q(&val, a, y) == XLAL_SUCCESS, XLAL_EFUNC );
      P = 1.0 - val;
   } else {
      XLAL_CHECK_REAL8( sf_gamma_inc_P(&P, a, y) == XLAL_SUCCESS, XLAL_EFUNC );
   }

   return P;
}
REAL8 cdf_gamma_P_usingmatlab(REAL8 x, REAL8 a, REAL8 b)
{

   REAL8 P;
   REAL8 y = x / b;

   if (x <= 0.0) return 0.0;

   if (y > a) {
      REAL8 val = matlab_gamma_inc(y, a, 1);
      P = 1.0 - val;
   } else {
      P = matlab_gamma_inc(y, a, 0);
   }

   return P;
}
INT4 cdf_gamma_Q(REAL8 *out, REAL8 x, REAL8 a, REAL8 b)
{
   REAL8 Q;
   REAL8 y = x / b;

   if (x <= 0.0) return 1.0;

   if (y < a) {
      REAL8 val;
      XLAL_CHECK( sf_gamma_inc_P(&val, a, y) == XLAL_SUCCESS, XLAL_EFUNC );
      Q = 1.0 - val;
   } else {
      XLAL_CHECK( sf_gamma_inc_Q(&Q, a, y) == XLAL_SUCCESS, XLAL_EFUNC );
   }

   *out = Q;
   return XLAL_SUCCESS;
}
REAL8 cdf_gamma_Q_usingmatlab(REAL8 x, REAL8 a, REAL8 b)
{

   REAL8 Q;
   REAL8 y = x / b;

   if (x <= 0.0) return 1.0;

   if (y < a) {
      REAL8 val = matlab_gamma_inc(y, a, 0);
      Q = 1.0 - val;
   } else {
      Q = matlab_gamma_inc(y, a, 1);
   }

   return Q;
}

/**
 * Compute the CDF P value at value x of a chi squared distribution with nu degrees of freedom
 * Rougly REAL4 precision
 * \param [in] x  CDF value at value x
 * \param [in] nu Number of degrees of freedom
 * \return CDF value
 */
REAL8 twospect_cdf_chisq_P(REAL8 x, REAL8 nu)
{
   REAL8 val = cdf_gamma_P(x, 0.5*nu, 2.0);
   XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
   return val;
} /* twospect_cdf_chisq_P() */


/**
 * Compute the CDF P value at value x of a chi squared distrubution with nu degrees of freedom using the Matlab-based function
 * \param [in] x  CDF value at value x
 * \param [in] nu Number of degrees of freedom
 * \return CDF value
 */
REAL8 matlab_cdf_chisq_P(REAL8 x, REAL8 nu)
{
   REAL8 val = cdf_gamma_P_usingmatlab(x, 0.5*nu, 2.0);
   XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
   return val;
} /* matlab_cdf_chisq_P() */


/**
 * Matlab's version of the non-central chi-squared CDF with nu degrees of freedom and non-centrality delta at value x
 * \param [in] x     Value at which to compute the CDF
 * \param [in] dof   Number of degrees of freedom
 * \param [in] delta Non-centrality parameter
 * \return CDF value
 */
REAL8 ncx2cdf(REAL8 x, REAL8 dof, REAL8 delta)
{

   REAL8 prob = 0.0;

   //Fail for bad inputs or return 0 if x<=0
   XLAL_CHECK_REAL8( dof >= 0.0 && delta >= 0.0, XLAL_EINVAL );
   if (x<=0.0) {
      return prob;
   }

   REAL8 err = LAL_REAL8_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT8 counter = (INT8)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = gsl_cdf_chisq_P(x, dof+2.0*counter);
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
   counter--;
   if (counter<0) return fmin(prob, 1.0);

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );

   //This part computes small probabilities
   INT4 fromzero = 0;
   if (prob==0.0) fromzero = 1;
   if (fromzero==1) {
      counter = 0;
      REAL8 pk = gsl_ran_poisson_pdf(0, halfdelta)*gsl_cdf_chisq_P(x, dof);
      REAL8 dp = 0.0;
      INT4 ok = 0;
      if ((REAL8)counter<halfdelta) ok = 1;
      while (ok==1) {
         counter++;
         P = gsl_ran_poisson_pdf(counter, halfdelta);
         C = gsl_cdf_chisq_P(x, dof+2.0*counter);
         dp = P*C;
         pk += dp;
         if (!(ok==1 && (REAL8)counter<halfdelta && dp>=err*pk)) ok = 0;
      }
      prob = pk;
   }

   return fmin(prob, 1.0);

} /* ncx2cdf() */

//Matlab's non-central chi square CDF up to REAL4 precision
REAL4 ncx2cdf_float(REAL4 x, REAL4 dof, REAL4 delta)
{

   REAL8 prob = 0.0;

   //Fail for bad inputs or return 0 if x<=0
   XLAL_CHECK_REAL4( dof >= 0.0 && delta >= 0.0, XLAL_EINVAL );
   if (x<=0.0) {
      return (REAL4)prob;
   }

   REAL8 err = (REAL8)LAL_REAL4_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT8 counter = (INT8)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = twospect_cdf_chisq_P((REAL8)x, (REAL8)(dof+2.0*counter));
   XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );
   counter--;
   if (counter<0) return (REAL4)fmin(prob, 1.0);

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );

   //This part computes small probabilities
   INT4 fromzero = 0;
   if (prob==0.0) fromzero = 1;
   if (fromzero==1) {
      counter = 0;
      REAL8 pk = gsl_ran_poisson_pdf(0, halfdelta)*twospect_cdf_chisq_P(x, dof);
      REAL8 dp = 0.0;
      INT4 ok = 0;
      if ((REAL8)counter<halfdelta) ok = 1;
      while (ok==1) {
         counter++;
         P = gsl_ran_poisson_pdf(counter, halfdelta);
         C = twospect_cdf_chisq_P(x, dof+2.0*counter);
         XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );
         dp = P*C;
         pk += dp;
         if (!(ok==1 && (REAL8)counter<halfdelta && dp>=err*pk)) ok = 0;
      }
      prob = pk;
   }

   return (REAL4)fmin(prob, 1.0);

} /* ncx2cdf_float() */

//Matlab's non-central chi-square tries to compute very small probabilities. We don't normally need this,
//so this function leaves out the last part to compute small probabilities.
REAL8 ncx2cdf_withouttinyprob(REAL8 x, REAL8 dof, REAL8 delta)
{

   REAL8 prob = 0.0;

   //Fail for bad inputs or return 0 if x<=0
   XLAL_CHECK_REAL8( dof >= 0.0 && delta >= 0.0, XLAL_EINVAL );
   if (x<=0.0) return prob;

   REAL8 err = LAL_REAL8_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT8 counter = (INT8)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = gsl_cdf_chisq_P(x, dof+2.0*counter);
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
   counter--;
   if (counter<0) return fmin(prob, 1.0);

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );

   return fmin(prob, 1.0);

} /* ncx2cdf_withouttinyprob() */

//Without small probabilities up to REAL4 precision
REAL4 ncx2cdf_float_withouttinyprob(REAL4 x, REAL4 dof, REAL4 delta)
{

   REAL8 prob = 0.0;

   //Fail for bad inputs or return 0 if x<=0
   XLAL_CHECK_REAL4( dof >= 0.0 && delta >= 0.0, XLAL_EINVAL );
   if (x<=0.0) return (REAL4)prob;

   REAL8 err = (REAL8)LAL_REAL4_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT8 counter = (INT8)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = twospect_cdf_chisq_P((REAL8)x, (REAL8)(dof+2.0*counter));
   XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );
   counter--;
   if (counter<0) return (REAL4)fmin(prob, 1.0);

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );

   return (REAL4)fmin(prob, 1.0);

} /* ncx2cdf_float_withouttinyprob() */


//This is ncx2cdf function like in Matlab, but using the Matlab version of the central chi square calculation instead of the GSL version
REAL8 ncx2cdf_withouttinyprob_withmatlabchi2cdf(REAL8 x, REAL8 dof, REAL8 delta)
{

   REAL8 prob = 0.0;

   //Fail for bad inputs or return 0 if x<=0
   XLAL_CHECK_REAL8( dof >= 0.0 && delta >= 0.0, XLAL_EINVAL );
   if (x<=0.0) return prob;

   REAL8 err = LAL_REAL8_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT8 counter = (INT8)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = matlab_cdf_chisq_P(x, dof+2.0*counter);  //Matlab chi square cdf calculation
   XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
   counter--;
   if (counter<0) return fmin(prob, 1.0);

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );

   return fmin(prob, 1.0);

} /* ncx2cdf_withouttinyprob_withmatlabchi2cdf() */

//This is ncx2cdf function like in Matlab, but using the Matlab version of the central chi square calculation instead of the GSL version; up to REAL4 precision
REAL4 ncx2cdf_float_withouttinyprob_withmatlabchi2cdf(REAL4 x, REAL4 dof, REAL4 delta)
{

   REAL8 prob = 0.0;

   //Fail for bad inputs or return 0 if x<=0
   XLAL_CHECK_REAL4( dof >= 0.0 && delta >= 0.0, XLAL_EINVAL );
   if (x<=0.0) return (REAL4)prob;

   REAL8 err = (REAL8)LAL_REAL4_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT8 counter = (INT8)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = matlab_cdf_chisq_P((REAL8)x, (REAL8)(dof+2.0*counter));  //Matlab chi2cdf
   XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );
   counter--;
   if (counter<0) return (REAL4)fmin(prob, 1.0);

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );

   return (REAL4)fmin(prob, 1.0);

} /* ncx2cdf_float_withouttinyprob_withmatlabchi2cdf() */


//Like Matlabs ncx2pdf
REAL8 ncx2pdf(REAL8 x, REAL8 dof, REAL8 delta)
{

   REAL8 dofint = 0.5*dof-1.0;
   REAL8 x1 = sqrt(x);
   REAL8 delta1 = sqrt(delta);
   REAL8 logreal8min = -708.3964185322641;

   REAL8 ul = 0.0;
   if (dofint<=-0.5) ul = -0.5*(delta+x) + 0.5*x1*delta1/(dofint+1.0) + dofint*(log(x)-LAL_LN2) - LAL_LN2 - lgamma(dofint+1.0);
   else ul = -0.5*(delta1-x1)*(delta1-x1) + dofint*(log(x)-LAL_LN2) - LAL_LN2 - lgamma(dofint+1.0) + (dofint+0.5)*log((dofint+0.5)/(x1*delta1+dofint+0.5));

   if (ul<logreal8min) return 0.0;

   //Scaled Bessel function?
   gsl_sf_result sbes = {0,0};
   INT4 status = gsl_sf_bessel_Inu_scaled_e(dofint, delta1*x1, &sbes);
   if (status==GSL_SUCCESS && sbes.val>0.0) return exp(-LAL_LN2 - 0.5*(x1-delta1)*(x1-delta1) + dofint*log(x1/delta1))*sbes.val;

   //Bessel function without scaling?
   gsl_sf_result bes;
   status = gsl_sf_bessel_Inu_e(dofint, delta1*x1, &bes);
   if (status==GSL_SUCCESS && bes.val>0.0) return exp(-LAL_LN2 - 0.5*(x+delta) + dofint*log(x1/delta1))*bes.val;

   //Okay, now recursion
   REAL8 lnsr2pi = log(sqrt(LAL_TWOPI));
   REAL8 dx = delta*x*0.25;
   INT8 K = GSL_MAX_INT(0, (INT8)floor(0.5*(sqrt(dofint*dofint+4.0*dx) - dofint)));
   REAL8 lntK = 0.0;
   if (K==0) lntK = -lnsr2pi - 0.5*(delta+log(dofint)) - (lgamma(dofint+1)-0.5*log(LAL_TWOPI*dofint)+dofint*log(dofint)-dofint) - binodeviance(dofint, 0.5*x);
   else lntK = -2.0*lnsr2pi - 0.5*(log(K) + log(dofint+K)) - (lgamma(K+1)-0.5*log(LAL_TWOPI*K)+K*log(K)-K) - (lgamma(dofint+K+1)-0.5*log(LAL_TWOPI*(dofint+K))+(dofint+K)*log(dofint+K)-(dofint+K)) - binodeviance(K, 0.5*delta) - binodeviance(dofint+K, 0.5*x);
   REAL8 sumK = 1.0;
   INT4 keep = 0;
   if (K>0) keep = 1;
   REAL8 term = 1.0;
   REAL8 k = K;
   while (keep==1) {
      term *= (dofint+k)*k/dx;
      sumK += term;
      if (k<=0 || term<=epsval(sumK) || keep!=1) keep = 0;
      k--;
   }
   keep = 1;
   term = 1.0;
   k = K+1;
   while (keep==1) {
      term /= (dofint+k)*k/dx;
      sumK += term;
      if (term<=epsval(sumK) || keep!=1) keep = 0;
      k++;
   }
   return 0.5*exp(lntK + log(sumK));

} /* ncx2pdf() */

/**
 * Matlab's ncx2inv function
 * \param [in] p     CDF P value from which to compute the inversion
 * \param [in] dof   Number of degrees of freedom
 * \param [in] delta Non-centrality parameter
 * \return The x value that corresponds to the P value
 */
REAL8 ncx2inv(REAL8 p, REAL8 dof, REAL8 delta)
{

   //Fail if bad input
   XLAL_CHECK_REAL8( delta >= 0.0, XLAL_EINVAL );

   if (delta==0.0) return gsl_cdf_chisq_Pinv(p, dof);

   REAL8 pk = p;
   INT4 count_limit = 100;
   INT4 count = 0;
   REAL8 crit = sqrt(LAL_REAL8_EPS);
   REAL8 mn = dof + delta;
   REAL8 variance = 2.0*(dof + 2.0*delta);
   REAL8 temp = log(variance + mn*mn);
   REAL8 mu = 2.0*log(mn) - 0.5*temp;
   REAL8 sigma = -2.0*log(mn) + temp;
   REAL8 xk = exp(norminv(pk, mu, sigma));
   REAL8 h = 0.0;
   REAL8 F = ncx2cdf(xk, dof, delta);
   XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
   while (count < count_limit) {
      count++;
      REAL8 f = ncx2pdf(xk, dof, delta);
      XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
      h = (F-pk)/f;
      REAL8 xnew = fmax(0.2*xk, fmin(5.0*xk, xk-h));
      REAL8 newF = ncx2cdf(xnew, dof, delta);
      XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
      INT4 worse = 0;
      while (worse==0) {
         if (!(fabs(newF-pk)>fabs(F-pk)*(1.0+crit) && fabs(xk-xnew)>crit*xk)) worse = 1;
         else {
            xnew = 0.5*(xnew + xk);
            newF = ncx2cdf(xnew, dof, delta);
            XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
         }
      }
      h = xk-xnew;
      if (!(fabs(h)>crit*fabs(xk) && fabs(h)>crit)) return xk;
      xk = xnew;
      F = newF;
   }

   fprintf(stderr, "%s: Warning! ncx2inv(%g, %g, %g) failed to converge!\n", __func__, p, dof, delta);
   return xk;

} /* ncx2inv() */


//Matlab's ncx2inv() function to REAL4 precision
REAL4 ncx2inv_float(REAL8 p, REAL8 dof, REAL8 delta)
{

   //Fail if bad input
   XLAL_CHECK_REAL4( delta >= 0.0, XLAL_EINVAL );

   if (delta==0.0) return (REAL4)gsl_cdf_chisq_Pinv(p, dof);

   REAL8 pk = p;
   INT4 count_limit = 100;
   INT4 count = 0;
   REAL8 crit = sqrt(LAL_REAL4_EPS);
   REAL8 mn = dof + delta;
   REAL8 variance = 2.0*(dof + 2.0*delta);
   REAL8 temp = log(variance + mn*mn);
   REAL8 mu = 2.0*log(mn) - 0.5*temp;
   REAL8 sigma = -2.0*log(mn) + temp;
   REAL8 xk = exp(norminv(pk, mu, sigma));
   REAL8 h = 0.0;
   REAL8 F = ncx2cdf_float_withouttinyprob_withmatlabchi2cdf(xk, dof, delta);
   XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );
   while (count < count_limit) {
      count++;
      REAL8 f = ncx2pdf(xk, dof, delta);
      XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );
      h = (F-pk)/f;
      REAL8 xnew = fmax(0.2*xk, fmin(5.0*xk, xk-h));
      REAL8 newF = ncx2cdf_float_withouttinyprob_withmatlabchi2cdf(xnew, dof, delta);
      XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );
      INT4 worse = 0;
      while (worse==0) {
         if (!(fabs(newF-pk)>fabs(F-pk)*(1.0+crit) && fabs(xk-xnew)>crit*xk)) worse = 1;
         else {
            xnew = 0.5*(xnew + xk);
            newF = ncx2cdf_float_withouttinyprob_withmatlabchi2cdf(xnew, dof, delta);
            XLAL_CHECK_REAL4( xlalErrno == 0, XLAL_EFUNC );
         }
      }
      h = xk-xnew;
      if (!(fabs(h)>crit*fabs(xk) && fabs(h)>crit)) return xk;
      xk = xnew;
      F = newF;
   }

   fprintf(stderr, "%s: Warning! ncx2inv_float() failed to converge!\n", __func__);
   return xk;

} /* ncx2inv_float() */

//Matlab's norminv function
REAL8 norminv(REAL8 p, REAL8 mu, REAL8 sigma)
{
   return mu - sigma*gsl_cdf_ugaussian_Qinv(p);
} /* norminv() */


//For the normal distribution, what is the SNR of a given value
REAL8 unitGaussianSNR(REAL8 value, REAL8 dof)
{
   REAL8 snr = (value - dof) / sqrt(2.0*dof);
   return snr;
} /* unitGaussianSNR() */


