/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Resample.h>

/* Macro for rounding a number.
REAL8 globalNum;
#define ROUND( num )                                                 \
  ( globalNum = (num),                                               \
    ( globalNum >= 0 ) ?                                             \
    (INT4)( globalNum + 0.5 ) :                                      \
    (INT4)( globalNum - 0.5 ) ) */

/** Simpler macro for rounding a number. */
#define ROUND( num )                                                 \
  (INT4)( num )

/* Test of another rounding macro.
#define ROUND( num )                                                 \
  (INT4)( num + 1 ) */

/* Prototypes for local functions. */
static INT2
ConstantRules( INT4 **tempRules, INT4 *nRules, REAL8 *startDiff,
	       REAL8 *stopDiff, INT4 d, REAL8 t, REAL8 tStop,
	       REAL8 dt, REAL4 *tBound, REAL4 *t0, REAL4 *poly,
	       INT4 nPoly, INT4 nSqrt );

static INT2
LinearRules( INT4 **tempRules, INT4 *nRules, REAL8 *startDiff,
	     REAL8 *stopDiff, INT4 d, REAL8 t, REAL8 tStop,
	     REAL8 dt, REAL4 *tBound, REAL4 *t0, REAL4 *poly,
	     INT4 nPoly, INT4 nSqrt );

static INT2
QuadraticRules( INT4 **tempRules, INT4 *nRules, REAL8 *startDiff,
		REAL8 *stopDiff, INT4 d, REAL8 t, REAL8 tStop,
		REAL8 dt, REAL4 *tBound, REAL4 *t0, REAL4 *poly,
		INT4 nPoly, INT4 nSqrt );

static void
FreeTempRules( INT4 **tempRules, INT4 nSqrt );


/**
 * \author Creighton, T. D.
 * \ingroup Resample_h
 * \brief Creates an object of type ResampleRules according to a
 * piecewise polynomial fit to the canonical time \f$\tau(t)\f$.
 *
 * This function creates an object <tt>**rules</tt> of type
 * \c ResampleRules, according to the polynomial fit to the canonical
 * time function found in <tt>*polyco</tt> and the sampling parameters
 * found in <tt>*params</tt>.  See the header Resample.h for a
 * description of these datatypes.  Initially the output handle must be a
 * valid handle (\c rules\f$\neq\f$\c NULL) but should not point to
 * an existing object (<tt>*rules</tt>=\c NULL).
 *
 * ### Algorithm ###
 *
 *
 * ### Formulae for computing \c ResampleRules: ###
 *
 * Since the resampling rules can be expected to be applied to datasets
 * with a huge number of sample times, the primary concern when designing
 * the algorithm was to maintain a low operation count per sample.  For
 * instance, a simple and easy approach would be to step through the time
 * domain at the rate of one resampled time interval \f$d\Delta t\f$ per
 * step, and record how many steps it takes for \f$\tau-t\f$ to change by one
 * (unresampled) interval \f$\Delta t\f$.  However, this approach requires us
 * to recalculate \f$\tau-t\f$ at \e each sample time, when it is only
 * likely to change significantly every thousand or so sample times.  A
 * much better approach is to invert the polynomial fit to give \f$t\f$ as a
 * function of \f$\tau-t\f$; that is, we change \f$\tau-t\f$ by a discrete amount
 * \f$\Delta t\f$, and find the corresponding change in \f$t\f$.  This function
 * need only be evaluated once for each distinct value of \f$\tau-t\f$,
 * rather than once for each distinct value of \f$t\f$.
 *
 * At present, LALCreateResampleRules() considers polynomial fits to
 * \f$\tau(t)\f$ only up to quadratic order.  This is almost certainly
 * sufficient for most conventional pulsar searches, where the dominant
 * high-order variation in \f$\tau\f$ over short periods (under one hour)
 * arises from the Earth's motion: \f$\tau-t\sim
 * (R_\mathrm{Earth}/c)\sin(2\pi t/P_\mathrm{Earth})\f$.  Since one
 * normally applies polynomial fits to \f$\tau(t)\f$ over hour-long time
 * stretches, \f$t\f$ is never more than half an hour
 * (\f$1/48P_\mathrm{Earth}\f$) from a fitting point, and the maximum error
 * in \f$\tau-t\f$ due to cubic and higher terms is \f$\sim
 * (R_\mathrm{Earth}/c)(2\pi/48)^3/3! \sim 8\mu\mathrm{s}\f$.  This is
 * insignificant at the 32~kHz sampling rate expected for LIGO.
 *
 * We do most of the computational work in dimensionless quantities
 * normalized by the resampled time interval \f$d\Delta t\f$ (where \f$\Delta
 * t\f$ is the unresampled time interval and \f$d\f$ is the decimation factor).
 * In this way, dimensionless times can be converted into numbers of
 * samples simply by rounding.  In any given fitting region
 * \f$[t_{\mathrm{bound}(i-1)},t_{\mathrm{bound}(i)})\f$ we define
 * dimensionless parameters:
 * \f{eqnarray}{
 * n  & = & \frac{\tau - t}{\Delta t}   \; , \nonumber\\
 * T  & = & \frac{t-t_{(i)}}{d\Delta t} \; , \nonumber\\
 * A_k & = & a_{k(i)}(d\Delta t)^{k-1}   \; . \nonumber
 * \tag{eq:dimensionless-params}
 * \f}
 * Eq.\eqref{eq_delta-tau} therefore
 * gives us, for a quadratic fit:
 * \f{equation}{
 * \tag{eq:dimensionless-tau}
 * \frac{n}{d} = A_0 + A_1 T + A_2 T^2 \; .
 * \f}
 * If one wants to know whether \f$n\f$ is increasing or decreasing with \f$T\f$,
 * that is given by the sign of \f$A_1+2A_2T\f$.
 *
 * To compute the resample rules, we want to find the (integral)
 * intervals in \f$T\f$ that correspond to (integral) shifts in \f$n\f$.
 * Inverting the quadratic formula, we have the obvious solutions:
 * \f{equation}{
 * \tag{eq:t-quadratic}
 * T = \left(-\frac{A_1}{2A_2}\right) +
 * s\times\sqrt{\left[\left(-\frac{A_1}{2A_2}\right)^2
 * -\frac{A_0}{A_2}\right] + \left(\frac{1}{dA_2}\right)\times n}
 * \; ,
 * \f}
 * where \f$s=\pm1\f$, depending on the sign of \f$A_2\f$ and whether \f$n\f$ is
 * increasing or decreasing with \f$T\f$ (\f$s\f$ has the same sign as \f$A_2\f$ if
 * \f$n\f$ is increasing, and the opposite if \f$n\f$ is decreasing).  Clearly
 * most of these coefficients need only be computed once for any given
 * polynomial fit, reducing the operation cost per evaluation of \f$T\f$.
 *
 * The procedure for computing a resample rule is to keep track of the
 * current values of \f$n\f$ and \f$T\f$, and whether \f$n\f$ is increasing or
 * decreasing with \f$T\f$.  One then increments/decrements \f$n\f$ by 1 (the
 * \e shift), computes a new value of \f$T\f$, and finds the difference
 * from the old value (the \e interval).  If the interval is less
 * than a whole number, then increase the size of the shift by 1 and try
 * again.  If the argument of the square root goes negative, then
 * \f$\tau-t\f$ has reached a turning point: reverse the direction of the
 * shift and the value of \f$s\f$.  If \f$T\f$ moves outside the current fitting
 * region, then get the new polynomial fit and recompute \f$A_k\f$ and \f$T\f$
 * (making sure to translate the old value of \f$T\f$, since the origin \f$T=0\f$
 * has moved), and determine anew whether \f$n\f$ is increasing or
 * decreasing.  Continue until \f$T\f$ moves out of the timespan for which we
 * want to compute resample rules.
 *
 * If the <tt>*polyco</tt> structure contains only two polynomial
 * coefficients per fitting interval, or if \f$A_2\f$ is dangerously small in
 * a given fitting interval, then the LALCreateResampleRules()
 * routine reverts to a linear fit:
 * \f{equation}{
 * \tag{eq:t-linear}
 * T = \left(-\frac{A_0}{A_1}\right) + \left(\frac{1}{dA_1}\right)
 * \times n \; .
 * \f}
 *
 * If the <tt>*polyco</tt> structure contains only \e one polynomial
 * coefficients per fitting interval, or if \f$A_1\f$ is also dangerously
 * small in a given fitting interval, then the
 * LALCreateResampleRules() routine reverts to a constant fit: a
 * shift is computed from the value of \f$n=dA_0\f$ at middle of the fitting
 * region (\f$T=0\f$).
 *
 * ### Computational details: ###
 *
 * The basic structure of the algorithm is an inner loop and an outer
 * loop.  The inner loop is iterated once each time that \f$n\f$ is
 * incremented or decremented and a new \f$T\f$ is computed; it terminates
 * when \f$T\f$ moves out of the current region of fit, or past the end of
 * the desired timespan.  The outer loop is iterated once for each set of
 * polynomial fitting parameters that cover the desired timespan, and
 * normally terminates only when \f$T\f$ leaves that timespan.  Loop
 * termination is done with \c break and \c return commands, so
 * as to avoid unnecessary repetition of tests.
 *
 * Once some initial setup is done by LALCreateResampleRules(), this
 * function calls distinct subroutines to perform quadratic, linear, or
 * constant fits to \f$\tau-t\f$, depending on the number of polynomial
 * coefficients per fitting region.  This avoids repeatedly querying the
 * number of coefficients.
 *
 * Since the length of the arrays in the \c ResampleRules object are
 * determined in the process of computing their contents, one must
 * allocate temporary storage with care.  On the one hand, one doesn't
 * want to allocate vastly more memory than will be required; on the
 * other hand, allocating nodes to a linked list every time one evaluates
 * \f$T\f$ would needlessly slow execution.  Unfortunately, the only hard
 * upper limit on the size of the \c ResampleRules arrays is the
 * total number \f$N_\mathrm{max}\f$ of resampled data in the timespan
 * covered by the rules, and this probably overestimates the actual
 * requirements by a factor of a thousand or so.
 *
 * The solution taken here is to allocate a two-dimensional array of
 * integers, with \f$2\sqrt{N_\mathrm{max}}\f$ rows of
 * \f$\sqrt{N_\mathrm{max}}\f$ elements.  Each pair of rows represents a
 * continued sequence of intervals and shifts, and a new pair is
 * allocated only when the previous one is full.  This requires the data
 * to be copied into the output array at the end of execution, and adds
 * an extra test to the inner loop to determine when new space is needed,
 * but the actual memory allocation is done infrequently and in
 * reasonably-sized blocks.  For instance, a week-long chunk of data
 * sampled at 2~kHz would have \f$N_\mathrm{max}=1.2\times10^9\f$, but might
 * actually require only \f$\sim2\f$~million integers.  These would end up
 * being stored in \f$\sim60\f$ arrays of length 34780.
 */
void
LALCreateResampleRules( LALStatus          *stat,
			ResampleRules      **rules,
			PolycoStruc        *polyco,
			ResampleParamStruc *params )
{
  INT4 nPoly;    /* Number of polynomial coefficients */
  REAL4 *tBound; /* Pointer to bounds of polynomial fitting regions */
  REAL4 *t0;     /* Pointer to polynomial fitting times */
  REAL4 *poly;   /* Pointer to polynomial coefficients */

  REAL8 tRuleStart; /* Start time for the resampling rules */
  REAL8 tRuleStop;  /* Stop time for the resampling rules */
  REAL8 tPolyStart; /* Start time for the polynomial fits */
  REAL8 startDiff;  /* Offset between tau and t at tRuleStart */
  REAL8 stopDiff=0; /* Offset between tau and t at tRuleStop */
  REAL8 dt;         /* Resampled time interval */

  INT4 n;           /* An index */
  INT4 nSqrt;       /* Square root of number of resampled data */
  INT4 nRules=0;    /* Number of resampling rules computed */
  INT4 **tempRules; /* Array of resampling rules */
  INT4 **row;       /* Pointer to a row in tempRules */
  INT4 *tempInterval; /* Pointer to interval elements in tempRules */
  INT4 *tempShift;    /* Pointer to shift elements in tempRules */
  INT4 *interval;     /* Pointer to interval elements in rules */
  INT2 *shift;        /* Pointer to shift elements in rules */

  INITSTATUS(stat);

  /* Check that the inputs all exist. */
  ASSERT(rules,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(polyco,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(polyco->t0,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(polyco->t0->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(polyco->polyco,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(polyco->polyco->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(params,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);

  /* Check that the output has not already been allocated. */
  ASSERT(!*rules,stat,RESAMPLEH_EOUT,RESAMPLEH_MSGEOUT);

  /* This should be obvious, but make sure that there are polynomial
     fitting times and coefficients defined for each polynomial
     fitting region. */
  ASSERT(polyco->tBound->length==polyco->t0->length,stat,
	 RESAMPLEH_ELENGTH,RESAMPLEH_MSGELENGTH);
  ASSERT(polyco->tBound->length==polyco->polyco->length,stat,
	 RESAMPLEH_ELENGTH,RESAMPLEH_MSGELENGTH);

  /* Check that the timespan being resampled is a subset of the
     timespan for which we have polynomial coefficients. */
  tRuleStart=params->start.gpsSeconds
    +(1.0e-9)*params->start.gpsNanoSeconds;
  tRuleStop=params->stop.gpsSeconds
    +(1.0e-9)*params->stop.gpsNanoSeconds;
  tPolyStart=polyco->start.gpsSeconds
    +(1.0e-9)*polyco->start.gpsNanoSeconds;
#ifndef LAL_NDEBUG
  REAL8 tPolyStop;  /* Stop time for the polynomial fits */
  tPolyStop=tPolyStart+polyco->tBound->data[polyco->tBound->length-1];
  ASSERT(tRuleStop<tPolyStop,stat,RESAMPLEH_ETIME,
	 RESAMPLEH_MSGETIME);
#endif
  ASSERT(tRuleStart>tPolyStart,stat,RESAMPLEH_ETIME,
	 RESAMPLEH_MSGETIME);

  /* Start assigning some local computational variables. */
  nPoly=polyco->polyco->vectorLength;
  tBound=polyco->tBound->data;
  t0=polyco->t0->data;
  poly=polyco->polyco->data;
  dt=params->decimate*params->deltaT;

  /* The number of resample rules can be no greater than the total
     number of resampled data in the timespan covered, and in general
     the number of rules will be much less.  To achieve computational
     efficiency without overly inflating our memory requirements, we
     store the data in a square array, allocating new rows only as
     necessary.  (By ``row'' I in fact mean two vectors of data, one
     storing intervals and the other shifts, so the array is actually
     not square but 2x1.)  Compute the length of each row or column in
     the square array. */
  nSqrt=(INT4)(sqrt((tRuleStop-tRuleStart)/dt));

  /* Now allocate the handle to the array, and the first row of
     intervals and shifts. */
  if(!(row=tempRules=(INT4 **)LALMalloc(2*nSqrt*sizeof(INT4 *)))){
    ABORT(stat,RESAMPLEH_EMEM,RESAMPLEH_MSGEMEM);
  }
  memset(tempRules,0,2*nSqrt*sizeof(INT4 *));
  if(!(*row=(INT4 *)LALMalloc(nSqrt*sizeof(INT4)))){
    FreeTempRules(tempRules,nSqrt);
    ABORT(stat,RESAMPLEH_EMEM,RESAMPLEH_MSGEMEM);
  }
  if(!(*(++row)=(INT4 *)LALMalloc(nSqrt*sizeof(INT4)))){
    FreeTempRules(tempRules,nSqrt);
    ABORT(stat,RESAMPLEH_EMEM,RESAMPLEH_MSGEMEM);
  }

  /* Determine the fitting region in which the resample rules start
     time lies. */
  tRuleStart-=tPolyStart;
  for(;tRuleStart>*(tBound++);poly+=nPoly,t0++)
    ;
  tBound--;

  XLALPrintError("tStop = %f\n",tRuleStop);

  /* Convert the start time into units of dt, relative to the initial
     fitting point. */
  tRuleStart=(tRuleStart-*t0)/dt;

  /* Convert the stop time into units of dt, relative to the start of
     polyco. */
  tRuleStop=(tRuleStop-tPolyStart)/dt;

  /* At present I know how to compute resample rules for constant,
     linear and quadratic polynomial fits.  I put these in two
     separate subroutines, and call one or the other based on the
     number of polynomial coefficients.  Thus I only have to test the
     number of coefficients once.  Note that the only error checking
     done within each subroutine is for memory allocation errors. */
  if(nPoly==1){
    if(ConstantRules(tempRules,&nRules,&startDiff,&stopDiff,
		     params->decimate,tRuleStart,tRuleStop,dt,
		     tBound,t0,poly,nPoly,nSqrt)){
      FreeTempRules(tempRules,nSqrt);
      ABORT(stat,RESAMPLEH_EMEM,RESAMPLEH_MSGEMEM);
    }
  }else if(nPoly==2){
    if(LinearRules(tempRules,&nRules,&startDiff,&stopDiff,
		   params->decimate,tRuleStart,tRuleStop,dt,
		   tBound,t0,poly,nPoly,nSqrt)){
      FreeTempRules(tempRules,nSqrt);
      ABORT(stat,RESAMPLEH_EMEM,RESAMPLEH_MSGEMEM);
    }
  }else{
    if(QuadraticRules(tempRules,&nRules,&startDiff,&stopDiff,
		      params->decimate,tRuleStart,tRuleStop,dt,
		      tBound,t0,poly,nPoly,nSqrt)){
      FreeTempRules(tempRules,nSqrt);
      ABORT(stat,RESAMPLEH_EMEM,RESAMPLEH_MSGEMEM);
    }
  }

  /* The resampling rules have been successfully computed!  Allocate
     the output array. */
  if(!(*rules=(ResampleRules *)LALMalloc(sizeof(ResampleRules)))){
    FreeTempRules(tempRules,nSqrt);
    ABORT(stat,RESAMPLEH_EMEM,RESAMPLEH_MSGEMEM);
  }
  memset(*rules,0,sizeof(ResampleRules));
  (*rules)->length=nRules;
  (*rules)->startDiff=startDiff*dt;
  (*rules)->stopDiff=stopDiff*dt;
  (*rules)->decimate=params->decimate;
  (*rules)->deltaT=params->deltaT;
  (*rules)->start.gpsSeconds=params->start.gpsSeconds;
  (*rules)->start.gpsNanoSeconds=params->start.gpsNanoSeconds;
  (*rules)->stop.gpsSeconds=params->stop.gpsSeconds;
  (*rules)->stop.gpsNanoSeconds=params->stop.gpsNanoSeconds;
  if(!((*rules)->interval=(INT4 *)LALMalloc(nRules*sizeof(INT4)))){
    FreeTempRules(tempRules,nSqrt);
    LALFree(*rules);
    *rules=NULL;
    ABORT(stat,RESAMPLEH_EMEM,RESAMPLEH_MSGEMEM);
  }
  if(!((*rules)->shift=(INT2 *)LALMalloc(nRules*sizeof(INT2)))){
    FreeTempRules(tempRules,nSqrt);
    LALFree((*rules)->interval);
    LALFree(*rules);
    *rules=NULL;
    ABORT(stat,RESAMPLEH_EMEM,RESAMPLEH_MSGEMEM);
  }

  /* Now move the resample rules into the output array. */
  row=tempRules;
  interval=(*rules)->interval;
  shift=(*rules)->shift;
  tempInterval=*(row++);
  tempShift=*(row++);
  n=nSqrt;
  while(nRules--){
    *(interval++)=*(tempInterval++);
    *(shift++)=*(tempShift++);
    if(!(--n)){
      tempInterval=*(row++);
      tempShift=*(row++);
      n=nSqrt;
    }
  }

  /* Free temporary memory and exit. */
  FreeTempRules(tempRules,nSqrt);
  RETURN(stat);
}


static INT2
ConstantRules( INT4 **tempRules, INT4 *nRules, REAL8 *startDiff,
	       REAL8 *stopDiff, INT4 d, REAL8 t, REAL8 tStop,
	       REAL8 dt, REAL4 *tBound, REAL4 *t0, REAL4 *poly,
	       INT4 nPoly, INT4 nSqrt )
{
  INT4 nRule=0; /* Running count of the number of shifts. */
  INT4 j=nSqrt; /* Countdown of space remaining on a tempRules row. */
  INT4 **row=tempRules;     /* Pointer to a tempRules row. */
  INT4 *intervals=*(row++); /* Pointer to a tempRules interval. */
  INT4 *shifts=*(row++);    /* Pointer to a tempRules shift. */
  REAL8 dInv=1.0/d;         /* Inverse of decimation factor. */

  INT4 interval; /* Distance to the next correction point. */
  INT4 shift=0;  /* Size of the correction. */
  INT4 n;        /* The current offset between tau and t. */

  REAL8 a0; /* The current dimensionless coefficient of (t-t0)^0. */
  REAL8 b1; /* Computational variable. */

  /* Get the initial dimensionless polynomial coefficients. */
  a0=poly[0]/dt;

  /* Compute the initial offset, and whether it is increasing or
     decreasing. */
  *startDiff=a0;
  n=ROUND(*startDiff*d);

  /* The outer loop steps through the time regions covered by each
     polynomial fit.  The calling routine should already have
     confirmed that there are enough such regions to cover the
     required timespan, so no testing need be done.  This loop does
     not terminate; rather, when the total resample rules timespan has
     been covered, the function will return.  */
  while(1){

    /* For a constant fit, just see whether a0 is significantly
       different from the current value of n.  If so, apply that shift
       at the fitting time. */
    shift=ROUND(a0)-n;
    if(shift){
      *(intervals++)=interval=ROUND(-t);
      *(shifts++)=shift;
      t+=interval+shift*dInv;
      n+=shift;
      shift=0;
      nRule++;

      /* See if we need to allocate more space for tempRules.  Return
	 an error if allocation fails. */
      if(!(--j)){
	XLALPrintError("t = %f\n",t*dt+*t0);
	if((intervals=*(row++)=(INT4 *)
	    LALMalloc(nSqrt*sizeof(INT4)))&&
	   (shifts=*(row++)=(INT4 *)
	    LALMalloc(nSqrt*sizeof(INT4))))
	  j=nSqrt;
	else
	  return 1;
      }
    }

    /* See whether the current region extends beyond the timespan for
       which we need to compute resample rules.  If so, then we're
       done; assign the remaining return values, and return to the
       main routine. */
    if(*tBound/dt>=tStop){
      /* Give a final interval indicating that the next shift is
         beyond the domain of the resampling rules.  The value of the
         shift is irrelevant; we set it to zero. */
      *intervals=(INT4)(tStop-*t0/dt-t)+1;
      *shifts=0;
      *nRules=nRule+1;
      *stopDiff=a0;
      return 0;
    }

    /* Otherwise, transform our times into the next fitting region,
       and get the new coefficients.  (This parallels the computations
       before the start of the loop.) */
    b1=*t0++;
    b1-=*t0;
    b1/=dt;
    t+=b1;
    tBound++;
    poly+=nPoly;
    a0=poly[0]/dt;
  }
}


static INT2
LinearRules( INT4 **tempRules, INT4 *nRules, REAL8 *startDiff,
	     REAL8 *stopDiff, INT4 d, REAL8 t, REAL8 tStop,
	     REAL8 dt, REAL4 *tBound, REAL4 *t0, REAL4 *poly,
	     INT4 nPoly, INT4 nSqrt )
{
  INT4 nRule=0; /* Running count of the number of shifts. */
  INT4 j=nSqrt; /* Countdown of space remaining on a tempRules row. */
  INT4 **row=tempRules;     /* Pointer to a tempRules row. */
  INT4 *intervals=*(row++); /* Pointer to a tempRules interval. */
  INT4 *shifts=*(row++);    /* Pointer to a tempRules shift. */
  REAL8 dInv=1.0/d;         /* Inverse of decimation factor. */

  INT2 sgn;      /* The sign of the resampling shift. */
  INT4 interval; /* Distance to the next correction point. */
  INT4 shift=0;  /* Size of the correction. */
  INT4 n;        /* The current offset between tau and t. */
  REAL8 tNext=t; /* The next correction point. */
  REAL8 tMax;    /* The start of the next fitting region. */

  REAL8 a0; /* The current dimensionless coefficient of (t-t0)^0. */
  REAL8 a1; /* The current dimensionless coefficient of (t-t0)^1. */

  REAL8 b1; /* Computational variable, normally -a0/a1. */
  REAL8 b2; /* Computational variable, normally 1/(d*a1). */

  /* Get the initial dimensionless polynomial coefficients. */
  a0=poly[0]/dt;
  a1=poly[1];

  /* Compute the initial offset, and whether it is increasing or
     decreasing. */
  *startDiff=a0 + a1*t;
  n=ROUND(*startDiff*d);
  sgn=(a1>=0.0) ? 1 : -1;

  /* The outer loop steps through the time regions covered by each
     polynomial fit.  The calling routine should already have
     confirmed that there are enough such regions to cover the
     required timespan, so no testing need be done.  This loop does
     not terminate; rather, when the total resample rules timespan has
     been covered, the function will return.  */
  while(1){

    /* Find the end of the current region ... */
    tMax=*tBound/dt;
    if(tStop<=tMax)
      tMax=tStop;
    /* ... relative to the current fitting point. */
    tMax-=*t0/dt;

    /* If a1 is effectively zero, revert to a constant fit. */
    if((fabs(a0)>=fabs(a1)*LAL_REAL8_MAX)||
       (tMax>=d*fabs(a1)*LAL_REAL8_MAX)){

      /* For a constant fit, just see whether a0 is significantly
	 different from the current value of n.  If so, apply that
	 shift at the fitting time. */
      shift=ROUND(a0)-n;
      if(shift){
	*(intervals++)=interval=ROUND(-t);
	*(shifts++)=shift;
	t+=interval+shift*dInv;
	tNext=tMax;
	n+=shift;
	shift=0;
	nRule++;

	/* See if we need to allocate more space for tempRules.
	   Return an error if allocation fails. */
	if(!(--j)){
	  XLALPrintError("t = %f\n",t*dt+*t0);
	  if((intervals=*(row++)=(INT4 *)
	      LALMalloc(nSqrt*sizeof(INT4)))&&
	     (shifts=*(row++)=(INT4 *)
	      LALMalloc(nSqrt*sizeof(INT4))))
	    j=nSqrt;
	  else
	    return 1;
	}
      }
    }else{

      /* This is the case of a linear fit over the current region.
	 Compute the linear fit coefficients. */
      b1 = -a0/a1;
      b2 = 1.0/(d*a1);

      /* Compute the shifts until tNext steps outside the current
	 fitting region.  We have to test this anyway within the loop,
	 so I'll use a break to exit the loop (so as to avoid an
	 unnecessary floating-point test). */
      while(1){
	n+=sgn;
	shift+=sgn;

	/* Compute the next shift point, and, if it's still in the
	   fitting region, compute the appropriate resample rule. */
	if((tNext=b1+b2*n)<tMax){
	  interval=ROUND(tNext-t);

	  /* However, if the interval is not at least one sample, do
	     nothing this iteration; keep incrementing n until the
	     interval is at least one datum. */
	  if(interval){
	    *(intervals++)=interval;
	    *(shifts++)=shift;
	    t+=interval+shift*dInv;
	    shift=0;
	    nRule++;

	    /* See if we need to allocate more space for tempRules.
	       Return an error if allocation fails. */
	    if(!(--j)){
	      XLALPrintError("t = %f\n",t*dt+*t0);
	      if((intervals=*(row++)=(INT4 *)
		  LALMalloc(nSqrt*sizeof(INT4)))&&
		 (shifts=*(row++)=(INT4 *)
		  LALMalloc(nSqrt*sizeof(INT4))))
		j=nSqrt;
	      else
		return 1;
	    }
	  }
	}else{

	  /* tNext has stepped into the next fitting region, so break
	     out of this loop.  Since we're not using this value of
	     tNext, we should decrement n and shift; they'll be
	     re-incremented before we try to compute tNext again. */
	  n-=sgn;
	  shift-=sgn;
	  break;
	}
      }
    }

    /* Now tNext is greater than tMax.  Either it has stepped into the
       next fitting region, or it has stepped past the timespan for
       which we need to compute resample rules.  If the latter, then
       we're done; assign the remaining return values, and return to
       the main routine. */
    if(tNext>=tStop-*t0/dt){
      /* Give a final interval indicating that the next shift is
         beyond the domain of the resampling rules.  The value of the
         shift is irrelevant, but for consistency we make it what it
         would have been from extrapolating the current fit. */
      tNext=tStop-*t0/dt;
      *intervals=(INT4)(tNext-t)+1;
      *shifts=sgn;
      *nRules=nRule+1;
      *stopDiff=a0 + a1*tNext;
      return 0;
    }

    /* Otherwise, transform our times into the next fitting region,
       and get the new coefficients.  (This parallels the computations
       before the start of the loop.) */
    b1=*t0++;
    b1-=*t0;
    b1/=dt;
    t+=b1;
    tNext+=b1;
    tBound++;
    poly+=nPoly;

    a0=poly[0]/dt;
    a1=poly[1];
    sgn=(a1>=0.0) ? 1 : -1;
  }
}


static INT2
QuadraticRules( INT4 **tempRules, INT4 *nRules, REAL8 *startDiff,
		REAL8 *stopDiff, INT4 d, REAL8 t, REAL8 tStop,
		REAL8 dt, REAL4 *tBound, REAL4 *t0, REAL4 *poly,
		INT4 nPoly, INT4 nSqrt )
{
  INT4 nRule=0; /* Running count of the number of shifts. */
  INT4 j=nSqrt; /* Countdown of space remaining on a tempRules row. */
  INT4 **row=tempRules;     /* Pointer to a tempRules row. */
  INT4 *intervals=*(row++); /* Pointer to a tempRules interval. */
  INT4 *shifts=*(row++);    /* Pointer to a tempRules shift. */
  REAL8 dInv=1.0/d;         /* Inverse of decimation factor. */

  INT2 sgn;      /* The sign of the resampling shift. */
  INT2 sgnRoot;  /* The sign of the root in the quadratic formula. */
  INT4 interval; /* Distance to the next correction point. */
  INT4 shift=0;  /* Size of the correction. */
  INT4 n;        /* The current offset between tau and t. */
  REAL8 tNext=t; /* The next correction point. */
  REAL8 tMax;    /* The start of the next fitting region. */

  REAL8 a0; /* The current dimensionless coefficient of (t-t0)^0. */
  REAL8 a1; /* The current dimensionless coefficient of (t-t0)^1. */
  REAL8 a2; /* The current dimensionless coefficient of (t-t0)^2. */
  REAL8 b1; /* Computational variable, normally -a1/(2*a2). */
  REAL8 b2; /* Computational variable, normally b1*b1-a0/a2. */
  REAL8 b3; /* Computational variable, normally 1/(d*a2). */
  REAL8 disc; /* Discriminant of the quadratic formula. */

  REAL8 sqrtMax=sqrt(LAL_REAL8_MAX);

  /* Get the initial dimensionless polynomial coefficients. */
  a0=poly[0]/dt;
  a1=poly[1];
  a2=poly[2]*dt;

  /* Compute the initial offset, whether it is increasing or
     decreasing, and which sign of the root to use. */
  b1=a1 + a2*t;
  *startDiff=a0 + b1*t;
  n=ROUND(*startDiff*d);
  b1+=a2*t;
  sgn=(b1>=0.0) ? 1 : -1;
  sgnRoot=(a2>=0.0) ? sgn : -sgn;

  /* The outer loop steps through the time regions covered by each
     polynomial fit.  The calling routine should already have
     confirmed that there are enough such regions to cover the
     required timespan, so no testing need be done.  This loop does
     not terminate; rather, when the total resample rules timespan has
     been covered, the function will return.  */
  while(1){

    /* Find the end of the current region ... */
    tMax=*tBound/dt;
    if(tStop<=tMax)
      tMax=tStop;
    /* ... relative to the current fitting point. */
    tMax-=*t0/dt;

    /* If a2 is effectively zero, revert to a linear fit. */
    if((fabs(a1)>=fabs(a2)*sqrtMax)||
       (fabs(a0)>=fabs(a2)*LAL_REAL8_MAX)||
       (tMax>=d*fabs(a2)*LAL_REAL8_MAX)){

      /* If a1 is effectively zero, revert to a constant fit. */
      if((fabs(a0)>=fabs(a1)*LAL_REAL8_MAX)||
	 (tMax>=d*fabs(a1)*LAL_REAL8_MAX)){

	/* For a constant fit, just see whether a0 is significantly
           different from the current value of n.  If so, apply that
           shift at the fitting time. */
	shift=ROUND(a0)-n;
	if(shift){
	  *(intervals++)=interval=ROUND(-t);
	  *(shifts++)=shift;
	  t+=interval+shift*dInv;
	  tNext=tMax;
	  n+=shift;
	  shift=0;
	  nRule++;

	  /* See if we need to allocate more space for tempRules.
	     Return an error if allocation fails. */
	  if(!(--j)){
	    XLALPrintError("t = %f\n",t*dt+*t0);
	    if((intervals=*(row++)=(INT4 *)
		LALMalloc(nSqrt*sizeof(INT4)))&&
	       (shifts=*(row++)=(INT4 *)
		LALMalloc(nSqrt*sizeof(INT4))))
	      j=nSqrt;
	    else
	      return 1;
	  }
	}
      }else{

	/* This is the case of a linear fit over the current region.
           Compute the linear fit coefficients. */
	b1 = -a0/a1;
	b2 = 1.0/(d*a1);

	/* Compute the shifts until tNext steps outside the current
	   fitting region.  We have to test this anyway within the
	   loop, so I'll use a break to exit the loop (so as to avoid
	   an unnecessary floating-point test). */
	while(1){
	  n+=sgn;
	  shift+=sgn;

	  /* Compute the next shift point, and, if it's still in the
	     fitting region, compute the appropriate resample rule. */
	  if((tNext=b1+b2*n)<tMax){
	    interval=ROUND(tNext-t);

	    /* However, if the interval is not at least one sample, do
	       nothing this iteration; keep incrementing n until the
	       interval is at least one datum. */
	    if(interval){
	      *(intervals++)=interval;
	      *(shifts++)=shift;
	      t+=interval+shift*dInv;
	      shift=0;
	      nRule++;

	      /* See if we need to allocate more space for tempRules.
		 Return an error if allocation fails. */
	      if(!(--j)){
		XLALPrintError("t = %f\n",t*dt+*t0);
		if((intervals=*(row++)=(INT4 *)
		    LALMalloc(nSqrt*sizeof(INT4)))&&
		   (shifts=*(row++)=(INT4 *)
		    LALMalloc(nSqrt*sizeof(INT4))))
		  j=nSqrt;
		else
		  return 1;
	      }
	    }
	  }else{

	    /* tNext has stepped into the next fitting region, so
	       break out of this loop.  Since we're not using this
	       value of tNext, we should decrement n and shift;
	       they'll be re-incremented before we try to compute
	       tNext again. */
	    n-=sgn;
	    shift-=sgn;
	    break;
	  }
	}
      }
    }else{

      /* Compute quadratic fit coefficients. */
      b1 = -0.5*a1/a2;
      b2 = b1*b1 - a0/a2;
      b3 = 1.0/(d*a2);

      /* Compute the shifts until tNext steps outside the current
         fitting region.  We have to test this anyway within the loop,
         so I'll use a break to exit the loop (so as to avoid an
         unnecessary floating-point test). */
      while(1){
	n+=sgn;
	shift+=sgn;

	/* If the discriminant goes negative, start incrementing n in
           the opposite direction. */
	if((disc=b2+b3*n)<0){
	  sgn*=-1;
	  sgnRoot=1;
	  n+=2*sgn;
	  shift+=2*sgn;
	  disc=b2+b3*n;
	}

	/* Compute the next shift point, and, if it's still in the
           fitting region, compute the appropriate resample rule. */
	if((tNext=b1+sgnRoot*sqrt(disc))<tMax){
	  interval=ROUND(tNext-t);

	  /* However, only do something if the interval is at least
             one sample; otherwise, keep incrementing n until the
             interval is at least one datum. */
	  if(interval){
	    *(intervals++)=interval;
	    *(shifts++)=shift;
	    t+=interval+shift*dInv;
	    shift=0;
	    nRule++;

	    /* See if we need to allocate more space for tempRules.
               Return an error if allocation fails. */
	    if(!(--j)){
	      XLALPrintError("t = %f\n",t*dt+*t0);
	      if((intervals=*(row++)=(INT4 *)
		  LALMalloc(nSqrt*sizeof(INT4)))&&
		 (shifts=*(row++)=(INT4 *)
		  LALMalloc(nSqrt*sizeof(INT4))))
		j=nSqrt;
	      else
		return 1;
	    }
	  }
	}else{
	  /* tNext has stepped into the next fitting region, so break
             out of this loop.  Since we're not using this value of
             tNext, we should decrement n and shift; they'll be
             re-incremented before we try to compute tNext again. */
	  n-=sgn;
	  shift-=sgn;
	  break;
	}
      }
    }

    /* Now tNext is greater than tMax.  Either it has stepped into the
       next fitting region, or it has stepped past the timespan for
       which we need to compute resample rules.  If the latter, then
       we're done; assign the remaining return values, and return to
       the main routine. */
    if(tNext>=tStop-*t0/dt){
      /* Give a final interval indicating that the next shift is
         beyond the domain of the resampling rules.  The value of the
         shift is irrelevant, but for consistency we make it what it
         would have been from extrapolating the current fit. */
      tNext=tStop-*t0/dt;
      *intervals=(INT4)(tNext-t)+1;
      *shifts=sgn;
      *nRules=nRule+1;
      b1=a1 + a2*tNext;
      *stopDiff=a0 + b1*tNext;
      return 0;
    }

    /* Otherwise, transform our times into the next fitting region,
       and get the new coefficients.  (This parallels the computations
       before the start of the loop.) */
    b1=*t0++;
    b1-=*t0;
    b1/=dt;
    t+=b1;
    tNext+=b1;
    tBound++;
    poly+=nPoly;

    a0=poly[0]/dt;
    a1=poly[1];
    a2=poly[2]*dt;

    b1=a1 + 2.0*a2*t;
    sgn=(b1>=0.0) ? 1 : -1;
    sgnRoot=(a2>=0.0) ? sgn : -sgn;
  }
}


static void
FreeTempRules(INT4 **tempRules, INT4 nSqrt)
     /* This routine frees the memory currently allocated to the array
	tempRules, which has up to 2*nSqrt rows of nSqrt columns.  The
	pointer tempRules is passed by value, so on return it will
	still point to the deallocated memory address. */
{
  INT4 **rows=tempRules;
  nSqrt*=2;
  while((nSqrt--)&&(*rows))
    LALFree(*(rows++));
  LALFree(tempRules);
  return;
}
