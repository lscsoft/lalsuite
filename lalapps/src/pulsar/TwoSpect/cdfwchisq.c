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

#include "cdfwchisq.h"



//Exp function to avoid underflows
REAL8 exp1(REAL8 x)
{
   //if (x<-50.0) return 0.0;
   if (x<-700.0) return 0.0;
   else return exp(x);
}


//count number of calls to errbound, truncation, coeff
void counter(qfvars *vars)
{
   
   vars->count++;
   
}


//If first is 1, then log(1 + x) ; else  log(1 + x) - x
REAL8 log1(REAL8 x, INT4 first)
{
   
   /* if (fabs(x) > 0.1) {
      if (first) return log(1.0 + x);
      else return (log(1.0 + x) - x);
   } else {
      REAL8 s, s1, term, y, k;
      y = x / (2.0 + x);
      term = 2.0 * y*y*y;
      k = 3.0;
      if (first) s = 2.0*y;
      else s = -x*y;
      y = y*y;
      for (s1=s+term/k; s1!=s; s1=s+term/k) { 
         k += 2.0; 
         term *= y; 
         s = s1;
      }
      return s;
   } */
   
   if (first) return gsl_sf_log_1plusx(x);
   else return gsl_sf_log_1plusx_mx(x);
   
}

//find order of absolute values of weights
void order(qfvars *vars)
{
   
   const CHAR *fn = __func__;
   
   //INT4 ii, jj;
   
   //Determine which values are largest to and place element numbers in th[]
   /* for (ii=0; ii<(INT4)vars->weights->length; ii++) {
      INT4 insertionpoint = ii;
      if (ii==0) {
         vars->sorting->data[insertionpoint] = 0;
      } else {
         while (insertionpoint>0 && fabs(vars->weights->data[ii])>fabs(vars->weights->data[vars->sorting->data[insertionpoint-1]])) insertionpoint--;
         
         for (jj=ii; jj>insertionpoint; jj--) {
            vars->sorting->data[jj] = vars->sorting->data[jj-1];
         }
         vars->sorting->data[insertionpoint] = ii;
      }
   } */
   
   INT4 ascend = 0;     //To sort descending, set ascend to zero
   if ( XLALHeapIndex(vars->sorting->data, vars->weights->data, vars->weights->length, sizeof(REAL8), &ascend, compar) != 0) {
      fprintf(stderr,"%s: XLALHeapIndex() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   vars->ndtsrt = 0; //Signify that we have done the sorting
   
}

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
   
}

//find bound on tail probability using mgf, cutoff point returned to *cx
REAL8 errbound(qfvars *vars, REAL8 u, REAL8* cx)
{
   
   REAL8 sum1, x, y, xconst;
   INT4 ii;
   
   counter(vars);
   
   xconst = u * vars->sigsq;
   sum1 = u * xconst;
   u = 2.0 * u;
   for (ii=vars->weights->length-1; ii>=0; ii--) {
      x = u * vars->weights->data[ii];
      y = 1.0 - x;
      xconst += vars->weights->data[ii] * (vars->noncentrality->data[ii] / y + vars->dofs->data[ii]) / y;
      sum1 += vars->noncentrality->data[ii] * (x*x/(y*y)) + vars->dofs->data[ii] * (x*x / y + log1(-x, 0 ));
   }
   *cx = xconst;
   
   return exp1(-0.5 * sum1);
   
}

//find cutoff so that p(qf > cutoff) < accx  if (upn > 0, p(qf < cutoff) < accx otherwise
REAL8 cutoff(qfvars *vars, REAL8 accx, REAL8* upn)
{
   
   REAL8 u1, u2, u, rb, xconst, c1, c2;
   
   u2 = *upn;
   u1 = 0.0;
   c1 = vars->wnmean;
   if (u2>0.0) rb = 2.0*vars->wnmax;
   else rb = 2.0*vars->wnmin;
   
   for (u=u2/(1.0+u2*rb); errbound(vars, u, &c2)>accx; u=u2/(1.0+u2*rb)) {
      u1 = u2;
      c1 = c2;
      u2 *= 2.0;
   }
   
   for (u=(c1-vars->wnmean)/(c2-vars->wnmean); u<0.9; u=(c1-vars->wnmean)/(c2-vars->wnmean)) {
      u = 0.5*(u1 + u2);
      if (errbound(vars, u/(1.0+u*rb), &xconst)>accx) {
         u1 = u; 
         c1 = xconst;
      } else {
         u2 = u;
         c2 = xconst;
      }
   }
   *upn = u2;
   
   return c2;
   
}

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
   u = 2.0 * u;
   for (ii=0; ii<(INT4)vars->weights->length; ii++ ) {
      x = (u * vars->weights->data[ii])*(u * vars->weights->data[ii]);
      sum1 += vars->noncentrality->data[ii] * x / (1.0 + x);
      if (x > 1.0) {
         prod2 += + vars->dofs->data[ii] * log(x);
         prod3 += + vars->dofs->data[ii] * log1(x, 1 );
         s += vars->dofs->data[ii];
      }
      else  prod1 += vars->dofs->data[ii] * log1(x, 1 );
   }
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
   
}

//find u such that truncation(u) < accx and truncation(u / 1.2) > accx
void findu(qfvars *vars, REAL8* utx, REAL8 accx)
{
   
   REAL8 u, ut;
   INT4 ii;
   REAL8 divis[] = {2.0, 1.4, 1.2, 1.1};
   
   ut = *utx;
   u = 0.25*ut;
   if ( truncation(vars, u, 0.0)>accx ) {
      for ( u=ut; truncation(vars, u, 0.0)>accx; u=ut) ut *= 4.0;
   } else {
      ut = u;
      for ( u=u/4.0; truncation(vars, u, 0.0) <=  accx; u=u/4.0 ) ut = u;
   }
   for (ii=0; ii<4; ii++) {
      u = ut/divis[ii]; 
      if ( truncation(vars, u, 0.0)<=accx )  ut = u; 
   }
   *utx = ut;
   
}


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
         sum3 = sum3 - 0.25 * vars->dofs->data[jj] * log1(y, 1 );
         y = vars->noncentrality->data[jj] * x / (1.0 + y);
         z = vars->dofs->data[jj] * atan(x) + y;
         sum1 += z;
         sum2 += fabs(z);
         sum3 -= 0.5 * x * y;
      }
      x = inpi * exp1(sum3) / u;
      if ( !  mainx ) x *= (1.0 - exp1(-0.5 * tausq * u*u));
      sum1 = sin(0.5 * sum1) * x;
      sum2 *= 0.5*x;
      vars->intl += sum1;
      vars->ersm += sum2;
   }
   
}

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
   }
   
   if (sum1 > 100.0) { 
      vars->fail = 1; 
      return 1.0; 
   } else {
      return pow(2.0, 0.25*sum1)*LAL_1_PI/(axl*axl);
   }
   
}


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
   
   const CHAR *fn = __func__;
   
      INT4 ii, nt, ntm;
      REAL8 acc1, almx, xlim, xnt, xntm;
      REAL8 utx, tausq, wnstd, intv, intv1, x, up, un, d1, d2;
      REAL8 qfval;
      INT4 rats[] = {1, 2, 4, 8};
      
      //for ( ii = 0; ii<7; ii++ ) trace[ii] = 0.0;
      
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
      vars->sigsq = sigma*sigma;
      wnstd = vars->sigsq;
      vars->wnmax = 0.0;
      vars->wnmin = 0.0;
      vars->wnmean = 0.0;
      for (ii=0; ii<(INT4)vars->weights->length; ii++ ) {
         if ( vars->dofs->data[ii] < 0  ||  vars->noncentrality->data[ii] < 0.0 ) { 
            *ifault = 3;
            //trace[6] = (REAL8)vars->count;
            return qfval;
         }
         wnstd += vars->weights->data[ii]*vars->weights->data[ii] * (2 * vars->dofs->data[ii] + 4.0 * vars->noncentrality->data[ii]);
         vars->wnmean += vars->weights->data[ii] * (vars->dofs->data[ii] + vars->noncentrality->data[ii]);
         if (vars->wnmax < vars->weights->data[ii]) vars->wnmax = vars->weights->data[ii];
         else if (vars->wnmin > vars->weights->data[ii]) vars->wnmin = vars->weights->data[ii];
      }
      
      if ( wnstd == 0.0  ) {  
         if (vars->c>0.0) qfval = 1.0;
         else qfval = 0.0;
         //trace[6] = (REAL8)vars->count;
         return qfval;
      }
      
      if ( vars->wnmin == 0.0 && vars->wnmax == 0.0 && sigma == 0.0 ) {
         *ifault = 3;
         //trace[6] = (REAL8)vars->count;
         return qfval;
      }
      
      wnstd = sqrt(wnstd);
      if (vars->wnmax < - vars->wnmin) almx = -vars->wnmin;
      else almx = vars->wnmax;

      /* starting values for findu, cutoff */
      utx = 16.0/wnstd;
      up = 4.5/wnstd;
      un = -up;
      
      /* truncation point with no convergence factor */
      findu(vars, &utx, 0.5*acc1);
      
      /* does convergence factor help */
      if (vars->c!=0.0  && almx>0.07*wnstd) {
         REAL8 coeffval = coeff(vars, vars->c);
         if (coeffval == 1.0) {
            fprintf(stderr,"%s: coeff() failed.\n", fn);
            XLAL_ERROR_REAL8(fn, XLAL_REAL8_FAIL_NAN);
         }
         tausq = 0.25*acc1/coeffval;
         if (vars->fail) vars->fail = 0 ;
         else if (truncation(vars, utx, tausq) < 0.2*acc1) {
            vars->sigsq += tausq;
            findu(vars, &utx, 0.25*acc1);
            //trace[5] = sqrt(tausq);
         }
      }
      //trace[4] = utx;
      acc1 = 0.5*acc1;

      /* find RANGE of distribution, quit if outside this */
   l1:
      d1 = cutoff(vars, acc1, &up) - vars->c;
      if (d1 < 0.0) {
         qfval = 1.0;
         //trace[6] = (REAL8)vars->count;
         return qfval;
      }
      d2 = vars->c - cutoff(vars, acc1, &un);
      if (d2 < 0.0) {
         qfval = 0.0;
         //trace[6] = (REAL8)vars->count;
         return qfval;
      }
      
      /* find integration interval */
      if (d1>d2) intv = LAL_TWOPI/d1;
      else intv = LAL_TWOPI/d2;
      
      /* calculate number of terms required for main and auxillary integrations */
      xnt = utx/intv;
      xntm = 3.0/sqrt(acc1);
      
      /* if (xntm>xlim) {
         *ifault = 1;
         trace[6] = (REAL8)vars->count;
         return qfval;
      }
      ntm = (INT4)round(xntm);
      intv1 = utx/ntm;
      x = LAL_TWOPI/intv1;
      
      if (x>fabs(vars->c)) tausq = 0.33*acc1/(1.1*(coeff(vars, vars->c-x) + coeff(vars, vars->c+x)));
      
      while (xnt>1.5*xntm && x>fabs(vars->c) && vars->fail!=1) {
         acc1 = 0.67*acc1;
         integrate(vars, ntm, intv1, tausq, 0 );
         xlim -= xntm;
         vars->sigsq += tausq;
         trace[2] = trace[2] + 1; 
         trace[1] = trace[1] + ntm + 1;
         findu(vars, &utx, 0.25*acc1);
         acc1 = 0.75*acc1;
         
         d1 = cutoff(vars, acc1, &up) - vars->c;
         if (d1 < 0.0) {
            qfval = 1.0;
            trace[6] = (REAL8)vars->count;
            return qfval;
         }
         d2 = vars->c - cutoff(vars, acc1, &un);
         if (d2 < 0.0) {
            qfval = 0.0;
            trace[6] = (REAL8)vars->count;
            return qfval;
         }
         if (d1>d2) intv = LAL_TWOPI/d1;
         else intv = LAL_TWOPI/d2;
         xnt = utx/intv;
         xntm = 3.0/sqrt(acc1);
         if (xntm>xlim) {
            *ifault = 1;
            trace[6] = (REAL8)vars->count;
            return qfval;
         }
         ntm = (INT4)round(xntm);
         intv1 = utx/ntm;
         x = LAL_TWOPI/intv1;
         if (x>fabs(vars->c)) tausq = 0.33*acc1/(1.1*(coeff(vars, vars->c-x) + coeff(vars, vars->c+x)));
      } */
      
      if (xnt>xntm*1.5) {
         //parameters for auxillary integration
         if (xntm>xlim) {
            *ifault = 1;
            //trace[6] = (REAL8)vars->count;
            return qfval;
         }
         ntm = (INT4)round(xntm);
         intv1 = utx/ntm;
         x = LAL_TWOPI/intv1;
         if (x<=fabs(vars->c)) goto l2;
         
         //calculate convergence factor
         REAL8 coeffvalplusx = coeff(vars, vars->c+x);
         if (coeffvalplusx == 1.0) {
            fprintf(stderr,"%s: coeff() failed.\n", fn);
            XLAL_ERROR_REAL8(fn, XLAL_REAL8_FAIL_NAN);
         }
         REAL8 coeffvalminusx = coeff(vars, vars->c-x);
         if (coeffvalminusx == 1.0) {
            fprintf(stderr,"%s: coeff() failed.\n", fn);
            XLAL_ERROR_REAL8(fn, XLAL_REAL8_FAIL_NAN);
         }
         tausq = (1.0/3.0)*acc1/(1.1*(coeffvalminusx + coeffvalplusx));
         if (vars->fail) goto l2;
         acc1 = (2.0/3.0)*acc1;
         
         //auxillary integration
         integrate(vars, ntm, intv1, tausq, 0 );
         xlim -= xntm;
         vars->sigsq += tausq;
         //trace[2] = trace[2] + 1; 
         //trace[1] = trace[1] + ntm + 1;
         
         //find truncation point with new convergence factor
         findu(vars, &utx, 0.25*acc1);
         acc1 = 0.75*acc1;
         goto l1;
      }

      /* main integration */
   l2:
      //trace[3] = intv;
      if (xnt > xlim) {
         *ifault = 1;
         //trace[6] = (REAL8)vars->count;
         return qfval;
      }
      nt = (INT4)round(xnt);
      integrate(vars, nt, intv, 0.0, 1 );
      //trace[2] = trace[2] + 1;
      //trace[1] = trace[1] + nt + 1;
      qfval = 0.5 - vars->intl;
      //trace[0] = vars->ersm;

      /* test whether round-off error could be significant allow for radix 8 or 16 machines */
      up = vars->ersm;
      x = up + 0.1*acc;
      for (ii=0; ii<4; ii++) {
         if (rats[ii] * x == rats[ii] * up) *ifault = 2;
      }

      //trace[6] = (REAL8)vars->count;
      return qfval;
}

