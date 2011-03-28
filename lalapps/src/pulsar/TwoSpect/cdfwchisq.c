/*
*  Copyright (C) 2010 Evan Goetz
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



//Exp function to avoid underflows
REAL8 exp1(REAL8 x)
{
   //if (x<-50.0) return 0.0;
   if (x<-700.0) return 0.0;
   else return exp(x);
} /* exp1() */

//Next special function routines based on the gsl functions
REAL8 twospect_log_1plusx(REAL8 x)
{
   
   if (fabs(x)<GSL_ROOT6_DBL_EPSILON) {
      return x*(1.0+x*(-.5+x*(1.0/3.0+x*(-.25+x*(.2+x*(-1.0/6.0+x*(1.0/7.0+x*(-.125+x*(1.0/9.0-.1*x)))))))));
   } else if (fabs(x)<0.5) {
      return twospect_log_1plusx_chebapprox(x);
   } else {
      return log(1.0 + x);
   }
   
}
REAL8 twospect_log_1plusx_mx(REAL8 x)
{
   
   if(fabs(x)<GSL_ROOT5_DBL_EPSILON) {
      return x*x*(-.5+x*(1.0/3.0+x*(-.25+x*(.2+x*(-1.0/6.0+x*(1.0/7.0+x*(-.125+x*(1.0/9.0-.1*x))))))));
   } else if (fabs(x)<0.5) {
      return twospect_log_1plusx_mx_chebapprox(x);
   } else {
      return log(1.0 + x) - x;
   }
   
}
REAL8 twospect_log_1plusx_chebapprox(REAL8 x)
{
   /* Chebyshev expansion for log(1 + x(t))/x(t) 
    x(t) = (4t-1)/(2(4-t))
    t(x) = (8x+1)/(2(x+2))
    -1/2 < x < 1/2
    -1 < t < 1
    */
   static REAL8 lopx_data[21] = {2.16647910664395270521272590407,
      -0.28565398551049742084877469679,
      0.01517767255690553732382488171,
      -0.00200215904941415466274422081,
      0.00019211375164056698287947962,
      -0.00002553258886105542567601400,
      2.9004512660400621301999384544e-06,
      -3.8873813517057343800270917900e-07,
      4.7743678729400456026672697926e-08,
      -6.4501969776090319441714445454e-09,
      8.2751976628812389601561347296e-10,
      -1.1260499376492049411710290413e-10,
      1.4844576692270934446023686322e-11,
      -2.0328515972462118942821556033e-12,
      2.7291231220549214896095654769e-13,
      -3.7581977830387938294437434651e-14,
      5.1107345870861673561462339876e-15,
      -7.0722150011433276578323272272e-16,
      9.7089758328248469219003866867e-17,
      -1.3492637457521938883731579510e-17,
      1.8657327910677296608121390705e-18};
   
   REAL8 t = 0.5*(8.0*x + 1.0)/(x+2.0);
   
   INT4 j;
   REAL8 d  = 0.0;
   REAL8 dd = 0.0;
   
   REAL8 y  = t;
   REAL8 y2 = 2.0 * y;
   
   REAL8 temp;
   for(j = 20; j>=1; j--) {
      temp = d;
      d = y2*d - dd + lopx_data[j];
      dd = temp;
   }
   
   d = x*(y*d - dd + 0.5 * lopx_data[0]);
   
   return d;
   
}
REAL8 twospect_log_1plusx_mx_chebapprox(REAL8 x)
{
   /* Chebyshev expansion for (log(1 + x(t)) - x(t))/x(t)^2
    * x(t) = (4t-1)/(2(4-t))
    * t(x) = (8x+1)/(2(x+2))
    * -1/2 < x < 1/2
    * -1 < t < 1
   */
   static REAL8 lopxmx_data[20] = {-1.12100231323744103373737274541,
      0.19553462773379386241549597019,
      -0.01467470453808083971825344956,
      0.00166678250474365477643629067,
      -0.00018543356147700369785746902,
      0.00002280154021771635036301071,
      -2.8031253116633521699214134172e-06,
      3.5936568872522162983669541401e-07,
      -4.6241857041062060284381167925e-08,
      6.0822637459403991012451054971e-09,
      -8.0339824424815790302621320732e-10,
      1.0751718277499375044851551587e-10,
      -1.4445310914224613448759230882e-11,
      1.9573912180610336168921438426e-12,
      -2.6614436796793061741564104510e-13,
      3.6402634315269586532158344584e-14,
      -4.9937495922755006545809120531e-15,
      6.8802890218846809524646902703e-16,
      -9.5034129794804273611403251480e-17,
      1.3170135013050997157326965813e-17};
   
   REAL8 t = 0.5*(8.0*x + 1.0)/(x+2.0);
   
   INT4 j;
   REAL8 d  = 0.0;
   REAL8 dd = 0.0;
   
   REAL8 y  = t;
   REAL8 y2 = 2.0 * y;
   
   REAL8 temp;
   for(j = 19; j>=1; j--) {
      temp = d;
      d = y2*d - dd + lopxmx_data[j];
      dd = temp;
   }
   
   d = x*x*(y*d - dd + 0.5 * lopxmx_data[0]);
   
   return d;
   
}


//count number of calls to errbound, truncation, coeff
void counter(qfvars *vars)
{
   
   vars->count++;
   
} /* counter() */


//find order of absolute values of weights
void order(qfvars *vars)
{
   
   const CHAR *fn = __func__;
   
   INT4 ascend = 1;     //To sort descending, set ascend to zero
   if ( XLALHeapIndex(vars->sorting->data, vars->weights->data, vars->weights->length, sizeof(REAL8), &ascend, compar) != 0) {
      fprintf(stderr,"%s: XLALHeapIndex() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   vars->ndtsrt = 0; //Signify that we have done the sorting
   
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
   
   counter(vars);
   
   xconst = u * vars->sigsq;
   sum1 = u * xconst;
   u *= 2.0;
   for (ii=vars->weights->length-1; ii>=0; ii--) {
      x = u * vars->weights->data[ii];
      y = 1.0 - x;
      xconst += vars->weights->data[ii] * (vars->noncentrality->data[ii] / y + vars->dofs->data[ii]) / y;
      sum1 += vars->noncentrality->data[ii] * (x*x/(y*y)) + vars->dofs->data[ii] * (x*x / y + gsl_sf_log_1plusx_mx(-x));
   }
   *cx = xconst;
   
   return exp1(-0.5 * sum1);
   
}/* errbound() */
REAL8 errbound_twospect(qfvars *vars, REAL8 u, REAL8* cx)
{
   
   REAL8 sum1, x, y, xconst;
   INT4 ii;
   
   counter(vars);
   
   xconst = u * vars->sigsq;
   sum1 = u * xconst;
   u *= 2.0;
   for (ii=vars->weights->length-1; ii>=0; ii--) {
      x = u * vars->weights->data[ii];
      y = 1.0 - x;
      xconst += 2.0*vars->weights->data[ii]/y;
      sum1 += 2 * (x*x / y + twospect_log_1plusx_mx(-x));
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
REAL8 truncation(qfvars *vars, REAL8 u, REAL8 tausq)
{
   REAL8 sum1, sum2, prod1, prod2, prod3, x, y, err1, err2;
   INT4 ii, s;

   counter(vars);
   
   sum1  = 0.0;
   prod2 = 0.0;
   prod3 = 0.0;
   s = 0;
   
   sum2 = (vars->sigsq + tausq) * u*u;
   prod1 = 2.0 * sum2;
   u *= 2.0;
   
   for (ii=0; ii<(INT4)vars->weights->length; ii++ ) {
      x = (u * vars->weights->data[ii])*(u * vars->weights->data[ii]);
      sum1 += vars->noncentrality->data[ii] * x / (1.0 + x);
      if (x > 1.0) {
         prod2 += vars->dofs->data[ii] * log(x);
         prod3 += vars->dofs->data[ii] * gsl_sf_log_1plusx(x);
         s += vars->dofs->data[ii];
      }
      else prod1 += vars->dofs->data[ii] * gsl_sf_log_1plusx(x);
   } /* for ii < vars->weights->length */
   
   sum1 *= 0.5;
   prod2 += prod1;
   prod3 += prod1;
   x = exp1(-sum1 - 0.25*prod2)*LAL_1_PI;
   y = exp1(-sum1 - 0.25*prod3)*LAL_1_PI;
   
   if (s==0) err1 = 1.0;
   else err1 = 2.0*x/s;
   
   if (prod3>1.0) err2 = 2.5*y;
   else err2 = 1.0;
   
   if (err2 < err1) err1 = err2;
   
   x = 0.5 * sum2;
   
   if (x<=y) err2 = 1.0;
   else err2 = y/x;
   
   if (err1<err2) return err1;
   else return err2;
   
} /* truncation() */
REAL8 truncation_twospect(qfvars *vars, REAL8 u, REAL8 tausq)
{
   REAL8 sum2, prod1, prod2, prod3, x, y, err1, err2;
   INT4 ii, s;
   
   counter(vars);
   
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
         prod3 += 2 * twospect_log_1plusx(x);
         s += 2;
      }
      else prod1 += 2 * twospect_log_1plusx(x);
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
   
   inpi = interv*LAL_1_PI;
   
   for (ii=nterm; ii>=0; ii--) {
      u = (ii + 0.5)*interv;
      sum1 = - 2.0*u*vars->c;
      sum2 = fabs(sum1);
      sum3 = -0.5*vars->sigsq * u*u;
      
      for (jj=(INT4)vars->weights->length-1; jj>=0; jj--) {
         x = 2.0 * vars->weights->data[jj] * u;
         y = x*x;
         sum3 -= 0.25 * vars->dofs->data[jj] * gsl_sf_log_1plusx(y);
         y = vars->noncentrality->data[jj] * x / (1.0 + y);
         z = vars->dofs->data[jj] * atan(x) + y;
         sum1 += z;
         sum2 += fabs(z);
         sum3 -= 0.5 * x * y;
      } /* for jj=vars->weights->length-1 --> 0 */
      
      x = inpi * exp1(sum3) / u;
      if ( !mainx ) x *= (1.0 - exp1(-0.5 * tausq * u*u));
      sum1 = sin(0.5 * sum1) * x;
      sum2 *= 0.5*x;
      vars->intl += sum1;
      vars->ersm += sum2;
   } /* for ii=nterm --> 0 */
   
} /* integrate() */
void integrate_twospect(qfvars *vars, INT4 nterm, REAL8 interv, REAL8 tausq, INT4 mainx)
{
   
   REAL8 inpi, u, sum1, sum2, sum3, x, z;
   INT4 ii, jj;
   
   inpi = interv*LAL_1_PI;
   
   for (ii=nterm; ii>=0; ii--) {
      u = (ii + 0.5)*interv;
      sum1 = - 2.0*u*vars->c;
      sum2 = fabs(sum1);
      sum3 = -0.5*vars->sigsq * u*u;
      
      for (jj=(INT4)vars->weights->length-1; jj>=0; jj--) {
         x = 2.0 * vars->weights->data[jj] * u;
         sum3 -= 0.5 * twospect_log_1plusx((x*x));
         z = 2.0 * atan(x);
         sum1 += z;
         sum2 += fabs(z);
      } /* for jj=vars->weights->length-1 --> 0 */
      
      x = inpi * exp1(sum3) / u;
      if ( !mainx ) x *= (1.0 - exp1(-0.5 * tausq * u*u));
      sum1 = sin(0.5 * sum1) * x;
      sum2 *= 0.5*x;
      vars->intl += sum1;
      vars->ersm += sum2;
   } /* for ii=nterm --> 0 */
   
} /* integrate_twospect() */

//Coefficient of tausq in error when convergence factor of exp1(-0.5*tausq*u^2) is used when df is evaluated at x
REAL8 coeff(qfvars *vars, REAL8 x)
{
   
   const CHAR *fn = __func__;
   
   REAL8 axl, axl1, axl2, sxl, sum1, lj;
   INT4 ii, jj, t;
   
   counter(vars);
   
   if (vars->ndtsrt) {
      order(vars);
      if (vars->ndtsrt) {
         fprintf(stderr,"%s: order() failed\n.", fn);
         vars->fail = 1;
         return 1.0;
      }
   }
   axl = fabs(x);
   
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
               return pow(2.0, 0.25*sum1)*LAL_1_PI/(axl*axl);
            }
         }
      }
   } /* for ii=vars->weights->length-1 --> 0 */
   
   if (sum1 > 100.0) { 
      vars->fail = 1; 
      return 1.0; 
   } else {
      return pow(2.0, 0.25*sum1)*LAL_1_PI/(axl*axl);
   }
   
} /* coeff() */
REAL8 coeff_twospect(qfvars *vars, REAL8 x)
{
   
   REAL8 axl, axl1, axl2, sxl, sum1, lj;
   INT4 ii, jj, t;
   
   counter(vars);
   
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
               return pow(2.0, 0.25*sum1)*LAL_1_PI/(axl*axl);
            }
         }
      }
   } /* for ii=vars->weights->length-1 --> 0 */
   
   if (sum1 > 100.0) { 
      vars->fail = 1; 
      return 1.0; 
   } else {
      return pow(2.0, 0.25*sum1)*LAL_1_PI/(axl*axl);
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
   vars->intl = 0.0;
   vars->ersm = 0.0;
   qfval = -1.0;
   acc1 = acc;
   vars->ndtsrt = 1;
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
   acc1 = 0.5*acc1;

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
      qfval = 0.5 - vars->intl;

      /* test whether round-off error could be significant allow for radix 8 or 16 machines */
      up = vars->ersm;
      x = up + 0.1*acc;
      for (ii=0; ii<4; ii++) {
         if (rats[ii] * x == rats[ii] * up) *ifault = 2;
      }
   
      return qfval;
} /* cdfwchisq() */
REAL8 cdfwchisq_twospect(qfvars *vars, REAL8 sigma, REAL8 acc, INT4 *ifault)
{
   
   INT4 ii, nt, ntm;
   REAL8 acc1, almx, xlim, xnt, xntm;
   REAL8 utx, tausq, wnstd, intv, intv1, x, up, un, d1, d2;
   REAL8 qfval;
   INT4 rats[] = {1, 2, 4, 8};
   
   *ifault = 0;
   vars->count = 0;
   vars->intl = 0.0;
   vars->ersm = 0.0;
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
      
      wnstd += 4*vars->weights->data[ii]*vars->weights->data[ii];
      vars->wnmean += 2*vars->weights->data[ii];
      
      //Find maximum and minimum values
      if (vars->wnmax < vars->weights->data[ii]) vars->wnmax = vars->weights->data[ii];
      else if (vars->wnmin > vars->weights->data[ii]) vars->wnmin = vars->weights->data[ii];
      
   } /* for ii < vars->weights->length */
   
   if ( wnstd == 0.0  ) {  
      if (vars->c>0.0) qfval = 1.0;
      else qfval = 0.0;
      return qfval;
   } /* if wnstd==0 */
   
   if ( vars->wnmin == 0.0 && vars->wnmax == 0.0 ) {
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
   findu_twospect(vars, &utx, 0.5*acc1);
   
   /* does convergence factor help */
   if (vars->c!=0.0  && almx>0.07*wnstd) {
      tausq = 0.25*acc1/coeff_twospect(vars, vars->c);
      if (vars->fail) vars->fail = 0;
      else if (truncation_twospect(vars, utx, tausq) < 0.2*acc1) {
         vars->sigsq += tausq;
         findu_twospect(vars, &utx, 0.25*acc1);
      }
   }
   acc1 = 0.5*acc1;
   
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
               integrate_twospect(vars, ntm, intv1, tausq, 0);
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
   integrate_twospect(vars, nt, intv, 0.0, 1);
   qfval = 0.5 - vars->intl;
   
   /* test whether round-off error could be significant allow for radix 8 or 16 machines */
   up = vars->ersm;
   x = up + 0.1*acc;
   for (ii=0; ii<4; ii++) {
      if (rats[ii] * x == rats[ii] * up) *ifault = 2;
   }
   
   return qfval;
} /* cdfwchisq_twospect() */

