/*
*  Copyright (C) 2010, 2011, 2014 Evan Goetz
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

// Program based on Robert Davies C algorithm from his 1980 paper
// "The distribution of a linear combination of chi^2 random variables"
// Journal of the Royal Statistical Society. Series C (Applied Statistics)
// Vol. 29, No. 3, 1980, pp. 323-33


#include <math.h>
#include <stdio.h>

#include <lal/LALConstants.h>
#include <lal/Sort.h>

#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_machine.h>

#include "cdfwchisq.h"
#include "vectormath.h"


//Exp function to avoid underflows
REAL8 exp1(REAL8 x)
{
   if (x<-700.0) return 0.0;
   else return exp(x);
} /* exp1() */

//Special functions
REAL8 twospect_log_1plusx(REAL8 x)
{
   return log1p(x);
} /* twospect_log_1plusx() */
REAL8 twospect_log_1plusx_mx(REAL8 x)
{
   return log1p(x)-x;
} /* twospect_log_1plusx_mx() */


//find order of absolute values of weights
void order(qfvars *vars)
{

   INT4 ascend = 1;     //To sort descending, set ascend to zero
   XLAL_CHECK_VOID( XLALHeapIndex(vars->sorting->data, vars->weights->data, vars->weights->length, sizeof(REAL8), &ascend, compar) == XLAL_SUCCESS, XLAL_EFUNC );

   vars->arrayNotSorted = 0; //Signify that we have done the sorting

} /* order() */

//Comparison routine for sorting algorithm (NOTE: Ascending order p=1, descending p=0)
int compar(void *p, const void *a, const void *b)
{
   REAL8 x = *((const REAL8 *)a);
   REAL8 y = *((const REAL8 *)b);
   int ascend = *(int *)p;

   if (ascend) {
      if (x < y) return -1;
      if (x > y) return 1;
      return 0;
   }

   if (x > y) return -1;
   if (x < y) return 1;
   return 0;

} /* compar() */

//find bound on tail probability using mgf, cutoff point returned to *cx
REAL8 errbound(qfvars *vars, REAL8 u, REAL8* cx)
{

   REAL8 sum1, x, y, xconst;
   INT4 ii;

   (vars->count)++;           //Increase counter

   xconst = u * vars->sigsq;  //xconst = u * sigma**2 + sum{ }
   sum1 = u * xconst;         //sum1 = u**2 * sigma**2 + sum{ } this is almost the equation after eq 9 in Davies 1973
                              //without the factor of 1/2 (applied at the end of this function)
   u *= 2.0;
   for (ii=vars->weights->length-1; ii>=0; ii--) {
      x = u * vars->weights->data[ii];       //x=2*u*lambda_j
      y = 1.0 - x;                           //y=1-2*u*lambda_j
      xconst += vars->weights->data[ii] * (vars->noncentrality->data[ii] / y + vars->dofs->data[ii]) / y;
      sum1 += vars->noncentrality->data[ii] * (x*x/(y*y)) + vars->dofs->data[ii] * (x*x / y + gsl_sf_log_1plusx_mx(-x));
   }
   *cx = xconst;

   return exp1(-0.5 * sum1);

}/* errbound() */
//Very similar to errbound() but this has no non-centrality parameters and assumes all dofs are 2.0
REAL8 errbound_twospect(qfvars *vars, REAL8 u, REAL8* cx)
{

   REAL8 sum1, x, y, xconst;
   INT4 ii;

   (vars->count)++;           //Increase counter

   xconst = u * vars->sigsq;
   sum1 = u * xconst;
   u *= 2.0;
   for (ii=vars->weights->length-1; ii>=0; ii--) {
      x = u * vars->weights->data[ii];
      y = 1.0 - x;
      xconst += 2.0*vars->weights->data[ii]/y;
      sum1 += 2 * (x*x / y + (log1p(-x)+x));
   }
   *cx = xconst;

   return exp1(-0.5 * sum1);

} /* errbound_twospect() */

//find cutoff so that p(qf > cutoff) < accx  if (upn > 0, p(qf < cutoff) < accx otherwise
REAL8 cutoff(qfvars *vars, REAL8 accx, REAL8* upn)
{

   REAL8 u1, u2, u, rb, xconst, c1, c2;

   u2 = *upn;
   u1 = 0.0;
   c1 = vars->wnmean;
   if (u2>0.0) rb = 2.0*vars->wnmax;
   else rb = 2.0*vars->wnmin;

   u = u2/(1.0 + u2*rb);
   while (errbound(vars, u, &c2)>accx) {
      u1 = u2;
      c1 = c2;
      u2 *= 2.0;
      u = u2/(1.0 + u2*rb);
   }

   u = (c1-vars->wnmean)/(c2-vars->wnmean);
   while (u<0.9) {
      u = 0.5*(u1 + u2);
      if (errbound(vars, u/(1.0+u*rb), &xconst)>accx) {
         u1 = u;
         c1 = xconst;
      } else {
         u2 = u;
         c2 = xconst;
      }
      u = (c1-vars->wnmean)/(c2-vars->wnmean);
   }
   *upn = u2;

   return c2;

} /* cutoff() */
//Same as cutoff() but uses *_twospect() functions
REAL8 cutoff_twospect(qfvars *vars, REAL8 accx, REAL8* upn)
{

   REAL8 u1, u2, u, rb, xconst, c1, c2;

   u2 = *upn;
   u1 = 0.0;
   c1 = vars->wnmean;
   if (u2>0.0) rb = 2.0*vars->wnmax;
   else rb = 2.0*vars->wnmin;

   u = u2/(1.0 + u2*rb);
   while (errbound_twospect(vars, u, &c2)>accx) {
      u1 = u2;
      c1 = c2;
      u2 *= 2.0;
      u = u2/(1.0 + u2*rb);
   }

   u = (c1-vars->wnmean)/(c2-vars->wnmean);
   while (u<0.9) {
      u = 0.5*(u1 + u2);
      if (errbound_twospect(vars, u/(1.0+u*rb), &xconst)>accx) {
         u1 = u;
         c1 = xconst;
      } else {
         u2 = u;
         c2 = xconst;
      }
      u = (c1-vars->wnmean)/(c2-vars->wnmean);
   }
   *upn = u2;

   return c2;

} /* cutoff_twospect() */

//bound integration error due to truncation at u
//Eq. 6, 7, 8 of Davies 1980
REAL8 truncation(qfvars *vars, REAL8 u, REAL8 tausq)
{
   REAL8 sum1, sum2, prod1, prod2, prod3, x, y, err1, err2;
   INT4 ii, s;

   //counter(vars);
   (vars->count)++;           //Increase counter

   sum1  = 0.0;   //Calculating N(u) = exp(-2u**2 sum_j(lambda_j**2 delta_j**2/(1+4u**2 lambda_j**2)))
   prod2 = 0.0;   //Calculating product (i)
   prod3 = 0.0;   //Calculating product (ii)
   s = 0;   //Sum of degrees of freedom

   sum2 = (vars->sigsq + tausq) * u*u;
   prod1 = 2.0 * sum2;
   u *= 2.0;      //This produces the factor of 4 in front of the products (U*lambda_j)**2 (i and ii) in Davies 1980

   for (ii=0; ii<(INT4)vars->weights->length; ii++ ) {
      x = (u * vars->weights->data[ii])*(u * vars->weights->data[ii]);  //(2*U*lambda_j)**2
      sum1 += vars->noncentrality->data[ii] * x / (1.0 + x);   //Sum after eq 4 in Davies 1980
      if (x > 1.0) {
         prod2 += vars->dofs->data[ii] * log(x);      //Logarithim of product (ii) produces sum of logorithms
         prod3 += vars->dofs->data[ii] * gsl_sf_log_1plusx(x);    //Logarithim of product (i) produces sum of logorithms
         s += vars->dofs->data[ii];    //sum of degrees of freedom
      }
      else prod1 += vars->dofs->data[ii] * gsl_sf_log_1plusx(x);
   } /* for ii < vars->weights->length */

   sum1 *= 0.5;      //Remove the extra prefactor of 2 before taking the exponential
   prod2 += prod1;
   prod3 += prod1;
   x = exp1(-sum1 - 0.25*prod2)*LAL_1_PI;    //Now remove logarithm by computing exponential (eq 6)
   y = exp1(-sum1 - 0.25*prod3)*LAL_1_PI;    //Now remove logarithm by computing exponential (eq 8)

   if (s==0) err1 = 1.0;
   else err1 = 2.0*x/s;

   if (prod3>1.0) err2 = 2.5*y;  //eq 8
   else err2 = 1.0;

   if (err2 < err1) err1 = err2;

   x = 0.5 * sum2;

   if (x<=y) err2 = 1.0;
   else err2 = y/x;

   if (err1<err2) return err1;
   else return err2;

} /* truncation() */
//Same as trunction() but assumes dofs are 2.0 and n.c. values are 0.0
REAL8 truncation_twospect(qfvars *vars, REAL8 u, REAL8 tausq)
{
   REAL8 sum2, prod1, prod2, prod3, x, y, err1, err2;
   INT4 ii, s;

   (vars->count)++;           //Increase counter

   prod2 = 0.0;
   prod3 = 0.0;
   s = 0;

   sum2 = (vars->sigsq + tausq) * u*u;
   prod1 = 2.0 * sum2;
   u *= 2.0;

   for (ii=0; ii<(INT4)vars->weights->length; ii++) {
      x = (u * vars->weights->data[ii])*(u * vars->weights->data[ii]);
      if (x > 1.0) {
         prod2 += 2 * log(x);
         prod3 += 2 * log1p(x);
         s += 2;
      }
      else prod1 += 2 * log1p(x);
   } /* for ii < vars->weights->length */

   prod2 += prod1;
   prod3 += prod1;
   x = exp1(-0.25*prod2)*LAL_1_PI;
   y = exp1(-0.25*prod3)*LAL_1_PI;

   err1 = 2.0*x/s;

   if (prod3>1.0) err2 = 2.5*y;
   else err2 = 1.0;

   if (err2 < err1) err1 = err2;

   x = 0.5 * sum2;

   if (x<=y) err2 = 1.0;
   else err2 = y/x;

   if (err1<err2) return err1;
   else return err2;

} /* truncation_twospect() */

//find u such that truncation(u) < accx and truncation(u / 1.2) > accx
void findu(qfvars *vars, REAL8* utx, REAL8 accx)
{

   REAL8 u, ut;
   INT4 ii;
   REAL8 divis[] = {2.0, 1.4, 1.2, 1.1};

   ut = *utx;
   u = 0.25*ut;
   if ( truncation(vars, u, 0.0)>accx ) {
      //for ( u=ut; truncation(vars, u, 0.0)>accx; u=ut) ut *= 4.0;
      u = ut;
      while (truncation(vars, u, 0.0)>accx) {
         ut *= 4.0;
         u = ut;
      }
   } else {
      ut = u;
      //for ( u=0.25*u; truncation(vars, u, 0.0) <=  accx; u=0.25*u ) ut = u;
      u *= 0.25;
      while (truncation(vars, u, 0.0) <=  accx) {
         ut = u;
         u *= 0.25;
      }
   }

   for (ii=0; ii<4; ii++) {
      u = ut/divis[ii];
      if ( truncation(vars, u, 0.0)<=accx )  ut = u;
   }

   *utx = ut;

} /* findu() */
//Same as findu() but uses *_twospect() functions
void findu_twospect(qfvars *vars, REAL8* utx, REAL8 accx)
{

   REAL8 u, ut;
   INT4 ii;
   REAL8 divis[] = {2.0, 1.4, 1.2, 1.1};

   ut = *utx;
   u = 0.25*ut;
   if ( truncation_twospect(vars, u, 0.0)>accx ) {
      u = ut;
      while (truncation_twospect(vars, u, 0.0)>accx) {
         ut *= 4.0;
         u = ut;
      }
   } else {
      ut = u;
      u *= 0.25;
      while (truncation_twospect(vars, u, 0.0) <=  accx) {
         ut = u;
         u *= 0.25;
      }
   }

   for (ii=0; ii<4; ii++) {
      u = ut/divis[ii];
      if ( truncation_twospect(vars, u, 0.0)<=accx )  ut = u;
   }

   *utx = ut;

} /* findu_twospect() */


//carry out integration with nterm terms, at stepsize interv.  if (! mainx) multiply integrand by 1.0-exp(-0.5*tausq*u^2)
void integrate(qfvars *vars, INT4 nterm, REAL8 interv, REAL8 tausq, INT4 mainx)
{

   REAL8 inpi, u, sum1, sum2, sum3, x, y, z;
   INT4 ii, jj;

   inpi = interv*LAL_1_PI;    //inpi = pi*(k + 1/2)

   for (ii=nterm; ii>=0; ii--) {
      u = (ii + 0.5)*interv;     //First part of eq 3 in Davies 1980, eq 9 in Davies 1973
      sum1 = -2.0*u*vars->c;     //Third sum, eq 13 of Davies 1980, the u*c term, will divide by 2 at the end
      sum2 = fabs(sum1);         //Davies 1980 says that the sine term can be replaced by the sum of abs vals of the arguement
      sum3 = -0.5*vars->sigsq * u*u;   //First part of eq 13 Davies 1980 in the exponential

      for (jj=(INT4)vars->weights->length-1; jj>=0; jj--) {
         x = 2.0 * vars->weights->data[jj] * u;    //2 * lambda_j * u
         y = x*x;    //4 * lambda_j**2 * u**2
         sum3 -= 0.25 * vars->dofs->data[jj] * gsl_sf_log_1plusx(y);    //product in eq 13 of Davies 1980
         y = vars->noncentrality->data[jj] * x / (1.0 + y);    //First sum argument in eq 13 of Davies 1980
         z = vars->dofs->data[jj] * atan(x) + y;      //Third sum argument in eq 13 of Davies 1980
         sum1 += z;        //Third sum in eq 13
         sum2 += fabs(z);
         sum3 -= 0.5 * x * y;    //Product
      } /* for jj=vars->weights->length-1 --> 0 */

      x = inpi * exp1(sum3) / u;
      if ( !mainx ) x *= (1.0 - exp1(-0.5 * tausq * u*u));  //For auxillary integration, we multiply by this factor)
      sum1 = sin(0.5 * sum1) * x;   //Now compute the sine
      sum2 *= 0.5*x;
      vars->integrationValue += sum1;     //integration value
      vars->integrationError += sum2;     //error on integration
   } /* for ii=nterm --> 0 */

} /* integrate() */
//Same as integrate() but assumes dofs are 2.0 and n.c. values are 0.0
void integrate_twospect(qfvars *vars, INT4 nterm, REAL8 interv, REAL8 tausq, INT4 mainx)
{

   REAL8 inpi, u, sum1, sum2, sum3, x, z;
   INT4 ii, jj;

   inpi = interv*LAL_1_PI;
   REAL8 neg2timesc = -2.0*vars->c, neghalftimessigsq = -0.5*vars->sigsq, neghalftimestausq = -0.5*tausq;

   for (ii=nterm; ii>=0; ii--) {
      u = (ii + 0.5)*interv;
      sum1 = neg2timesc*u;
      sum2 = fabs(sum1);
      sum3 = neghalftimessigsq * u*u;

      REAL8 twotimesu = 2.0*u;
      for (jj=(INT4)vars->weights->length-1; jj>=0; jj--) {
         x = twotimesu * vars->weights->data[jj];
         sum3 -= 0.5 * log1p((x*x));
         z = 2.0 * atan(x);
         sum1 += z;
         sum2 += fabs(z);
      } /* for jj=vars->weights->length-1 --> 0 */

      x = inpi * exp1(sum3) / u;
      if ( !mainx ) x *= (1.0 - exp1(neghalftimestausq * u*u));
      sum1 = sin(0.5*sum1) * x;
      sum2 *= 0.5*x;
      vars->integrationValue += sum1;
      vars->integrationError += sum2;
   } /* for ii=nterm --> 0 */

} /* integrate_twospect() */
//This is from eq 13 of Davies 1980 and makes more sense than the integrate() function while giving nearly identical results (last digits of double precision are only slightly different)
void integrate_eg(qfvars *vars, INT4 nterm, REAL8 interv, REAL8 tausq, INT4 mainx)
{

   INT4 ii, jj;

   for (ii=nterm; ii>=0; ii--) {
      REAL8 u = (ii + 0.5)*interv;

      REAL8 exptermarguementsum = 0.0, logofproductterm = 0.0, sinetermargumentsum = 0.0, sumofabssinesumargs = 0.0;
      for (jj=(INT4)vars->weights->length-1; jj>=0; jj--) {
         exptermarguementsum += (vars->weights->data[jj]*vars->weights->data[jj])*(vars->noncentrality->data[jj]*vars->noncentrality->data[jj])/(1.0 + 4.0*(u*u)*(vars->weights->data[jj]*vars->weights->data[jj]));

         logofproductterm += -0.25*vars->dofs->data[jj]*log1p(4.0*(u*u)*(vars->weights->data[jj]*vars->weights->data[jj]));

         sinetermargumentsum += 0.5*vars->dofs->data[jj]*atan(2.0*u*vars->weights->data[jj]) + (vars->noncentrality->data[jj]*vars->noncentrality->data[jj])*u*vars->weights->data[jj]/(1.0 + 4.0*(u*u)*(vars->weights->data[jj]*vars->weights->data[jj]));

         sumofabssinesumargs += fabs(0.5*(2.0)*atan(2.0*u*vars->weights->data[jj]) + (vars->noncentrality->data[jj]*vars->noncentrality->data[jj])*u*vars->weights->data[jj]/(1.0 + 4.0*(u*u)*(vars->weights->data[jj]*vars->weights->data[jj])));
      }
      REAL8 firstterm = exp1(-2.0*(u*u)*exptermarguementsum - 0.5*(u*u)*vars->sigsq);
      REAL8 secondterm = exp1(logofproductterm);
      REAL8 thirdterm = sin(sinetermargumentsum - u*vars->c);
      REAL8 together = firstterm * secondterm * thirdterm/(LAL_PI*(ii+0.5));
      REAL8 together2 = firstterm * secondterm * (sumofabssinesumargs + fabs(u*vars->c)) / (LAL_PI*(ii+0.5));
      if ( !mainx ) {
         together *= (1.0 - exp1(-0.5 * tausq * u*u));
         together2 *= (1.0 - exp1(-0.5 * tausq * u*u));
      }

      vars->integrationValue += together;
      vars->integrationError += together2;
   }

}
//Rewrite of integrate_eg() to make it fast
void integrate_twospect2(qfvars *vars, INT4 nterm, REAL8 interv, REAL8 tausq, INT4 mainx)
{

   INT4 ii, jj;

   for (ii=nterm; ii>=0; ii--) {
      REAL8 u = (ii + 0.5)*interv;
      REAL8 oneoverPiTimesiiPlusHalf = 1.0/(LAL_PI*(ii+0.5));

      REAL8 exptermarguementsum = 0.0, logofproductterm = 0.0, sinetermargumentsum = 0.0, sumofabssinesumargs = 0.0;

      for (jj=(INT4)vars->weights->length-1; jj>=0; jj--) {
         REAL8 twoUtimesWeight = 2.0*u*vars->weights->data[jj];
         REAL8 atanTwoUtimesWeight = atan(twoUtimesWeight);

         logofproductterm += -0.5*log1p(twoUtimesWeight*twoUtimesWeight);
         sinetermargumentsum += atanTwoUtimesWeight;
         sumofabssinesumargs += fabs(atanTwoUtimesWeight);
      }
      REAL8 firstterm = exp1(-2.0*(u*u)*exptermarguementsum);
      REAL8 secondterm = exp1(logofproductterm);
      REAL8 thirdterm = sin(sinetermargumentsum - u*vars->c);
      REAL8 together = firstterm * secondterm * thirdterm * oneoverPiTimesiiPlusHalf;
      REAL8 together2 = firstterm * secondterm * (sumofabssinesumargs + fabs(u*vars->c)) * oneoverPiTimesiiPlusHalf;
      if ( !mainx ) {
         REAL8 scalingfactor = (1.0 - exp1(-0.5 * tausq * u*u));
         together *= scalingfactor;
         together2 *= scalingfactor;
      }

      vars->integrationValue += together;
      vars->integrationError += together2;
   }

}
//Use fast/SSE/AVX to make the integration even faster
INT4 fast_integrate_twospect2(qfvars *vars, INT4 nterm, REAL8 interv, REAL8 tausq, INT4 mainx)
{

   INT4 ii, jj;

   REAL8Vector *scaledweightvector = NULL;
   XLAL_CHECK( (scaledweightvector = XLALCreateREAL8Vector(vars->weights->length)) != NULL, XLAL_EFUNC );

   for (ii=nterm; ii>=0; ii--) {
      REAL8 u = (ii + 0.5)*interv;
      REAL8 oneoverPiTimesiiPlusHalf = 1.0/(LAL_PI*(ii+0.5));
      REAL8 twoU = 2.0*u;

      REAL8 exptermarguementsum = 0.0, logofproductterm = 0.0, sinetermargumentsum = 0.0, sumofabssinesumargs = 0.0;
      if (vars->useSSE) XLAL_CHECK( sseScaleREAL8Vector(scaledweightvector, vars->weights, twoU) == XLAL_SUCCESS, XLAL_EFUNC );
      else if (vars->useAVX) XLAL_CHECK( avxScaleREAL8Vector(scaledweightvector, vars->weights, twoU) == XLAL_SUCCESS, XLAL_EFUNC );
      else for (jj=0; jj<(INT4)scaledweightvector->length; jj++) scaledweightvector->data[jj] = vars->weights->data[jj]*twoU;

      for (jj=(INT4)vars->weights->length-1; jj>=0; jj--) {
         REAL8 atanTwoUtimesWeight = atan(scaledweightvector->data[jj]);

         logofproductterm += -0.5*log1p(scaledweightvector->data[jj]*scaledweightvector->data[jj]);
         sinetermargumentsum += atanTwoUtimesWeight;
         sumofabssinesumargs += fabs(atanTwoUtimesWeight);
      }
      REAL8 firstterm = exp1(-twoU*u*exptermarguementsum);
      REAL8 secondterm = exp1(logofproductterm);
      REAL8 thirdterm = sin(sinetermargumentsum - u*vars->c);
      REAL8 together = firstterm * secondterm * thirdterm * oneoverPiTimesiiPlusHalf;
      REAL8 together2 = firstterm * secondterm * (sumofabssinesumargs + fabs(u*vars->c)) * oneoverPiTimesiiPlusHalf;
      if ( !mainx ) {
         REAL8 scalingfactor = (1.0 - exp1(-0.5 * tausq * u*u));
         together *= scalingfactor;
         together2 *= scalingfactor;
      }

      vars->integrationValue += together;
      vars->integrationError += together2;
   }

   XLALDestroyREAL8Vector(scaledweightvector);

   return XLAL_SUCCESS;

}


//Coefficient of tausq in error when convergence factor of exp1(-0.5*tausq*u^2) is used when df is evaluated at x
//Eq. 10 of Davies 1980
REAL8 coeff(qfvars *vars, REAL8 x)
{

   REAL8 axl, axl1, axl2, sxl, sum1, lj;
   INT4 ii, jj, t;

   (vars->count)++;

   //If the sort hasn't been done, then do it now!
   if (vars->arrayNotSorted) {
      order(vars);
      if (vars->arrayNotSorted) {
         fprintf(stderr,"%s: order() failed\n.", __func__);
         vars->fail = 1;
         return 1.0;
      }
   }
   axl = fabs(x);    //absolute value of the value of c

   if (x>0.0) sxl = 1.0;
   else sxl = -1.0;

   sum1 = 0.0;
   for ( ii=vars->weights->length-1; ii>=0; ii-- ) {
      t = vars->sorting->data[ii];
      if ( vars->weights->data[t] * sxl > 0.0 ) {
         lj = fabs(vars->weights->data[t]);
         axl1 = axl - lj*(vars->dofs->data[t] + vars->noncentrality->data[t]);
         axl2 = 8.0*lj/LAL_LN2;
         if ( axl1 > axl2 )  axl = axl1;
         else {
            if ( axl > axl2 ) axl = axl2;
            sum1 = (axl - axl1) / lj;
            for ( jj = ii-1; jj>=0; jj--) sum1 += (vars->dofs->data[vars->sorting->data[jj]] + vars->noncentrality->data[vars->sorting->data[jj]]);

            if (sum1 > 100.0) {
               vars->fail = 1;
               return 1.0;
            } else {
               return exp2(0.25*sum1)*LAL_1_PI/(axl*axl);
            }
         }
      }
   } /* for ii=vars->weights->length-1 --> 0 */

   if (sum1 > 100.0) {
      vars->fail = 1;
      return 1.0;
   } else {
      return exp2(0.25*sum1)*LAL_1_PI/(axl*axl);
   }

} /* coeff() */
//Same as coeff() but assuming dofs are 2.0 and n.c. values are 0.0
REAL8 coeff_twospect(qfvars *vars, REAL8 x)
{

   REAL8 axl, axl1, axl2, sxl, sum1, lj;
   INT4 ii, jj, t;

   (vars->count)++;

   axl = fabs(x);

   if (x>0.0) sxl = 1.0;
   else sxl = -1.0;

   sum1 = 0.0;
   for ( ii=vars->weights->length-1; ii>=0; ii-- ) {
      t = vars->sorting->data[ii];
      if ( vars->weights->data[t] * sxl > 0.0 ) {
         lj = fabs(vars->weights->data[t]);
         axl1 = axl - 2*lj;
         axl2 = 8.0*lj/LAL_LN2;
         if ( axl1 > axl2 )  axl = axl1;
         else {
            if ( axl > axl2 ) axl = axl2;
            sum1 = (axl - axl1) / lj;
            for ( jj = ii-1; jj>=0; jj--) sum1 += 2;

            if (sum1 > 100.0) {
               vars->fail = 1;
               return 1.0;
            } else {
               return exp2(0.25*sum1)*LAL_1_PI/(axl*axl);
            }
         }
      }
   } /* for ii=vars->weights->length-1 --> 0 */

   if (sum1 > 100.0) {
      vars->fail = 1;
      return 1.0;
   } else {
      return exp2(0.25*sum1)*LAL_1_PI/(axl*axl);
   }

} /* coeff_twospect() */


/*  distribution function of a linear combination of non-central
   chi-squared random variables :

input:
   vars             cdfwchisq parameters structure
   sigma            coefficient of standard normal variable
   acc              maximum error

output:
   ifault = 1       required accuracy NOT achieved
            2       round-off error possibly significant
            3       invalid parameters
            4       unable to locate integration parameters
            5       out of memory     */

REAL8 cdfwchisq(qfvars *vars, REAL8 sigma, REAL8 acc, INT4 *ifault)
{

   INT4 ii, nt, ntm;
   REAL8 acc1, almx, xlim, xnt, xntm;
   REAL8 utx, tausq, wnstd, intv, intv1, x, up, un, d1, d2;
   REAL8 qfval;
   INT4 rats[] = {1, 2, 4, 8};

   *ifault = 0;
   vars->count = 0;
   vars->integrationValue = 0.0;
   vars->integrationError = 0.0;
   qfval = -1.0;
   acc1 = acc;
   vars->arrayNotSorted = 1;
   vars->fail = 0;
   xlim = (REAL8)vars->lim;

   /* find wnmean, wnstd, wnmax and wnmin of weights, check that parameter values are valid */
   vars->sigsq = sigma*sigma;    //Sigma squared
   wnstd = vars->sigsq;          //weights*noise standard deviation initial value
   vars->wnmax = 0.0;            //Initial value for weights*noise maximum
   vars->wnmin = 0.0;            //Initial value for weights*noise minimum
   vars->wnmean = 0.0;           //Initial value for weights*noise 'mean'
   for (ii=0; ii<(INT4)vars->weights->length; ii++ ) {
      if ( vars->dofs->data[ii] < 0  ||  vars->noncentrality->data[ii] < 0.0 ) {
         *ifault = 3;
         return qfval;
      } /* return error if any degrees of freedom is less than 0 or noncentrality parameter is less than 0.0 */

      wnstd += vars->weights->data[ii]*vars->weights->data[ii] * (2 * vars->dofs->data[ii] + 4.0 * vars->noncentrality->data[ii]);
      vars->wnmean += vars->weights->data[ii] * (vars->dofs->data[ii] + vars->noncentrality->data[ii]);

      //Find maximum and minimum values
      if (vars->wnmax < vars->weights->data[ii]) {
         vars->wnmax = vars->weights->data[ii];
      } else if (vars->wnmin > vars->weights->data[ii]) {
         vars->wnmin = vars->weights->data[ii];
      }
   }

   if ( wnstd == 0.0  ) {
      if (vars->c>0.0) qfval = 1.0;
      else qfval = 0.0;
      return qfval;
   }

   if ( vars->wnmin == 0.0 && vars->wnmax == 0.0 && sigma == 0.0 ) {
      *ifault = 3;
      return qfval;
   }

   wnstd = sqrt(wnstd);

   //almx is absolute value maximum of weights
   if (vars->wnmax < -vars->wnmin) almx = -vars->wnmin;
   else almx = vars->wnmax;

   /* starting values for findu, cutoff */
   utx = 16.0/wnstd;
   up = 4.5/wnstd;
   un = -up;

   /* truncation point with no convergence factor */
   findu(vars, &utx, 0.5*acc1);

   /* does convergence factor help */
   if (vars->c!=0.0  && almx>0.07*wnstd) {
      tausq = 0.25*acc1/coeff(vars, vars->c);
      if (vars->fail) vars->fail = 0;
      else if (truncation(vars, utx, tausq) < 0.2*acc1) {
         vars->sigsq += tausq;
         findu(vars, &utx, 0.25*acc1);
      }
   }
   acc1 *= 0.5;

      /* find RANGE of distribution, quit if outside this */
   l1:
      d1 = cutoff(vars, acc1, &up) - vars->c;
      if (d1 < 0.0) {
         qfval = 1.0;
         return qfval;
      }
      d2 = vars->c - cutoff(vars, acc1, &un);
      if (d2 < 0.0) {
         qfval = 0.0;
         return qfval;
      }

      /* find integration interval */
      if (d1>d2) intv = LAL_TWOPI/d1;
      else intv = LAL_TWOPI/d2;

      /* calculate number of terms required for main and auxillary integrations */
      xnt = utx/intv;
      xntm = 3.0/sqrt(acc1);

      if (xnt>xntm*1.5) {
         //parameters for auxillary integration
         if (xntm>xlim) {
            *ifault = 1;
            return qfval;
         }
         ntm = (INT4)round(xntm);
         intv1 = utx/ntm;
         x = LAL_TWOPI/intv1;
         if (x<=fabs(vars->c)) goto l2;

         //calculate convergence factor
         REAL8 coeffvalplusx = coeff(vars, vars->c+x);
         REAL8 coeffvalminusx = coeff(vars, vars->c-x);
         tausq = (1.0/3.0)*acc1/(1.1*(coeffvalminusx + coeffvalplusx));
         if (vars->fail) goto l2;
         acc1 = (2.0/3.0)*acc1;

         //auxillary integration
         //fprintf(stderr,"Num terms in auxillary integration %d\n", ntm);
         integrate(vars, ntm, intv1, tausq, 0);
         xlim -= xntm;
         vars->sigsq += tausq;

         //find truncation point with new convergence factor
         findu(vars, &utx, 0.25*acc1);
         acc1 *= 0.75;
         goto l1;
      }

      /* main integration */
   l2:
      if (xnt > xlim) {
         *ifault = 1;
         return qfval;
      }
      nt = (INT4)round(xnt);
      //fprintf(stderr,"Num terms in main integration %d\n", nt);
      integrate(vars, nt, intv, 0.0, 1);
      qfval = 0.5 - vars->integrationValue;

      /* test whether round-off error could be significant allow for radix 8 or 16 machines */
      up = vars->integrationError;
      x = up + 0.1*acc;
      for (ii=0; ii<4; ii++) {
         if (rats[ii] * x == rats[ii] * up) *ifault = 2;
      }

      return qfval;
} /* cdfwchisq() */
//Re-write of the cdfwchisq() function, but in a better way, without go-to's and to be more LAL compliant (though not totally!)
REAL8 cdfwchisq_twospect(qfvars *vars, REAL8 sigma, REAL8 acc, INT4 *ifault)
{

   INT4 ii, nt, ntm;
   REAL8 acc1, almx, xlim, xnt, xntm;
   REAL8 utx, tausq, wnstd, intv, intv1, x, up, un, d1, d2;
   REAL8 qfval;
   INT4 rats[] = {1, 2, 4, 8};

   //Initialize values
   *ifault = 0;
   vars->count = 0;
   vars->integrationValue = 0.0;
   vars->integrationError = 0.0;
   qfval = -1.0;
   acc1 = acc;
   vars->fail = 0;
   xlim = (REAL8)vars->lim;

   /* find wnmean, wnstd, wnmax and wnmin of weights, check that parameter values are valid */
   vars->sigsq = sigma*sigma;    //Sigma squared
   wnstd = vars->sigsq;          //weights*noise standard deviation initial value
   vars->wnmax = 0.0;            //Initial value for weights*noise maximum
   vars->wnmin = 0.0;            //Initial value for weights*noise minimum
   vars->wnmean = 0.0;           //Initial value for weights*noise 'mean'
   for (ii=0; ii<(INT4)vars->weights->length; ii++ ) {

      wnstd += 4*vars->weights->data[ii]*vars->weights->data[ii];  //weight_i^2 * 2 * 2
      vars->wnmean += 2*vars->weights->data[ii];                   //2*weight_i

      //Find maximum and minimum values of the weights
      if (vars->wnmax < vars->weights->data[ii]) vars->wnmax = vars->weights->data[ii];
      else if (vars->wnmin > vars->weights->data[ii]) vars->wnmin = vars->weights->data[ii];

   } /* for ii < vars->weights->length */

   //If somehow the wnstd value was 0, then output either 1 or 0
   if ( wnstd == 0.0 ) {
      if (vars->c>0.0) qfval = 1.0;
      else qfval = 0.0;
      return qfval;
   } /* if wnstd==0 */

   //If the min and max weight are zero, then there needs to be an error
   if ( vars->wnmin == 0.0 && vars->wnmax == 0.0 ) {
      *ifault = 3;
      return qfval;
   }

   //Do the square root of wnstd
   wnstd = sqrt(wnstd);

   //almx is absolute value maximum of weights
   if (vars->wnmax < -vars->wnmin) almx = -vars->wnmin;
   else almx = vars->wnmax;

   /* starting values for findu, cutoff */
   utx = 16.0/wnstd;
   up = 4.5/wnstd;
   un = -up;

   /* truncation point with no convergence factor */
   findu_twospect(vars, &utx, 0.5*acc1);

   /* does convergence factor help? */
   if (vars->c!=0.0  && almx>0.07*wnstd) {
      tausq = 0.25*acc1/coeff_twospect(vars, vars->c);
      if (vars->fail) vars->fail = 0;
      else if (truncation_twospect(vars, utx, tausq) < 0.2*acc1) {
         vars->sigsq += tausq;
         findu_twospect(vars, &utx, 0.25*acc1);
      }
   }
   acc1 *= 0.5;

   BOOLEAN contin = 1;

   /* find RANGE of distribution, quit if outside this */
   while (contin) {
      d1 = cutoff_twospect(vars, acc1, &up) - vars->c;
      if (d1 < 0.0) {
         qfval = 1.0;
         return qfval;
      }
      d2 = vars->c - cutoff_twospect(vars, acc1, &un);
      if (d2 < 0.0) {
         qfval = 0.0;
         return qfval;
      }

      /* find integration interval */
      if (d1>d2) intv = LAL_TWOPI/d1;
      else intv = LAL_TWOPI/d2;

      /* calculate number of terms required for main and auxillary integrations */
      xnt = utx/intv;
      xntm = 3.0/sqrt(acc1);

      if (xnt>xntm*1.5) {
         //parameters for auxillary integration
         if (xntm>xlim) {
            *ifault = 1;
            return qfval;
         }
         ntm = (INT4)round(xntm);
         intv1 = utx/ntm;
         x = LAL_TWOPI/intv1;
         if (x<=fabs(vars->c)) {
            contin = 0;
         } else {
            //calculate convergence factor
            REAL8 coeffvalplusx = coeff_twospect(vars, vars->c+x);
            REAL8 coeffvalminusx = coeff_twospect(vars, vars->c-x);
            tausq = acc1/(1.1*(coeffvalminusx + coeffvalplusx)*3.0);
            if (vars->fail) {
               contin = 0;
            } else {
               acc1 = (2.0/3.0)*acc1;

               //auxillary integration
               //fprintf(stderr,"Num terms in auxillary integration %d\n", ntm);
               XLAL_CHECK_REAL8( fast_integrate_twospect2(vars, ntm, intv1, tausq, 0) == XLAL_SUCCESS, XLAL_EFUNC );
               xlim -= xntm;
               vars->sigsq += tausq;

               //find truncation point with new convergence factor
               findu_twospect(vars, &utx, 0.25*acc1);
               acc1 *= 0.75;
            }
         }
      } else {
         contin = 0;
      }
   }

   /* main integration */
   if (xnt > xlim) {
      *ifault = 1;
      return qfval;
   }
   nt = (INT4)round(xnt);  //number of terms in main integration
   XLAL_CHECK_REAL8( fast_integrate_twospect2(vars, nt, intv, 0.0, 1) == XLAL_SUCCESS, XLAL_EFUNC );
   qfval = 0.5 - vars->integrationValue;

   /* test whether round-off error could be significant allow for radix 8 or 16 machines */
   up = vars->integrationError;
   x = up + 0.1*acc;
   for (ii=0; ii<4; ii++) {
      if (rats[ii] * x == rats[ii] * up) *ifault = 2;
   }

   return qfval;
} /* cdfwchisq_twospect() */
