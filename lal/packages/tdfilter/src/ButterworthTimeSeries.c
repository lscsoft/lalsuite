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

#include <complex.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <math.h>
#include <lal/IIRFilter.h>
#include <lal/BandPassTimeSeries.h>

/**
 * \defgroup ButterworthTimeSeries_c Module ButterworthTimeSeries.c
 * \ingroup BandPassTimeSeries_h
 *
 * \author Creighton, T. D.
 *
 * \brief Applies a low- or high-pass Butterworth filter to a time series.
 *
 * ### Description ###
 *
 * These routines perform an in-place time-domain band-pass filtering of
 * a data sequence <tt>*series</tt>, using a Butterworth filter generated
 * from parameters <tt>*params</tt>.  The routines construct a filter with
 * the square root of the desired amplitude response, which it then
 * applied to the data once forward and once in reverse.  This gives the
 * full amplitude response with little or no frequency-dependent phase
 * shift.
 *
 * The routine <tt>LALDButterworthREAL4TimeSeries()</tt> applies a
 * double-precision filter to single-precision data, using
 * <tt>LALDIIRFilterREAL4Vector()</tt> and
 * <tt>LALDIIRFilterREAL4VectorR()</tt>.
 *
 * ### Algorithm ###
 *
 * The frequency response of a Butterworth low-pass filter is easiest to
 * express in terms of the transformed frequency variable \f$w=\tan(\pi
 * f\Delta t)\f$, where \f$\Delta t\f$ is the sampling interval (i.e.
 * <tt>series-\>deltaT</tt>).  In this parameter, then, the \e power
 * response (attenuation) of the filter is:
 * \f[
 * |R|^2 = \sqrt{a} = \frac{1}{1+(w/w_c)^{2n}} \; ,
 * \f]
 * where \f$n\f$ is the filter order and \f$w_c\f$ is the characteristic
 * frequency.  We have written the attenuation as \f$\sqrt{a}\f$ to emphasize
 * that the full attenuation \f$a\f$ is achieved only after filtering twice
 * (once forward, once in reverse).  Similarly, a Butterworth high-pass
 * filter is given by
 * \f[
 * |R|^2 = \sqrt{a} = \frac{1}{1+(w_c/w)^{2n}} \; .
 * \f]
 * If one is given a filter order \f$n\f$, then the characteristic frequency
 * can be determined from the attenuation at some any given frequency.
 * Alternatively, \f$n\f$ and \f$w_c\f$ can both be computed given attenuations
 * at two different frequencies.
 *
 * Frequencies in <tt>*params</tt> are assumed to be real frequencies \f$f\f$
 * given in the inverse of the units used for the sampling interval
 * <tt>series-\>deltaT</tt>.  In order to be used, the pass band parameters
 * must lie in the ranges given below; if a parameter lies outside of its
 * range, then it is ignored and the filter is calculated from the
 * remaining parameters.  If too many parameters are missing, the routine
 * will fail.  The acceptable parameter ranges are:
 *
 * <dl>
 * <dt>params-\>nMax</dt><dd>   = 1, 2, \f$\ldots\f$</dd>
 * <dt>params-\>f1, f2</dt><dd> \f$\in(0,\{2\times<tt>series-\>deltaT</tt>\}^{-1}) \f$</dd>
 * <dt>params-\>a1, a2</dt><dd> \f$\in (0,1) \f$</dd>
 * </dl>
 *
 * If both pairs of frequencies and amplitudes are given, then \c a1,
 * \c a2 specify the minimal requirements on the attenuation of the
 * filter at frequencies \c f1, \c f2.  Whether the filter is a
 * low- or high-pass filter is determined from the relative sizes of
 * these parameters.  In this case the \c nMax parameter is optional;
 * if given, it specifies an upper limit on the filter order.  If the
 * desired attenuations would require a higher order, then the routine
 * will sacrifice performance in the stop band in order to remain within
 * the specified \c nMax.
 *
 * If one of the frequency/attenuation pairs is missing, then the filter
 * is computed using the remaining pair and \c nMax (which must be
 * given).  The filter is taken to be a low-pass filter if \c f1,
 * \c a1 are given, and high-pass if \c f2, \c a2 are given.
 * If only one frequency and no corresponding attenuation is specified,
 * then it is taken to be the characteristic frequency (i.e. the
 * corresponding attenuation is assumed to be \f$\sqrt{a}=1/2\f$).  If none
 * of these conditions are met, the routine will return an error.
 *
 * The following table summarizes the decision algorithm.  A \f$\bullet\f$
 * symbol indicates that the parameter is specified in the range given
 * above.  A \f$\circ\f$ symbol indicates that the parameter is <i>not given</i>, i.e. not specified in the valid range.
 *
 * <table>
 * <tr><th>nMax</th><th>f1</th><th>a1</th><th>f2</th><th>a2</th><th>Procedure</th></tr>
 * <tr><td>\f$\circ\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>
 * Type of filter (low- or high-pass), \f$w_c\f$, and \f$n\f$ are
 * computed from all four transition-band parameters.</td></tr>
 * <tr><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>
 * Ditto, but if the resulting \f$n>\f$ \c nMax, \f$w_c\f$ is
 * computed from \c nMax and the (\c f,\c a)
 * pair with the \e larger \c a.</td></tr>
 * <tr><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\circ\f$</td><td>\f$\circ\f$</td><td>
 * Low-pass filter; \f$w_c\f$ is computed from \c nMax,
 * \c f1, and \c a1.</td></tr>
 * <tr><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\circ\f$</td><td>\f$\bullet\f$</td><td>
 * Ditto; \c a2 is ignored.</td></tr>
 * <tr><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\circ\f$</td><td>
 * Ditto; \c f2 is ignored.</td></tr>
 * <tr><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\circ\f$</td><td>\f$\circ\f$</td><td>\f$\circ\f$</td><td>
 * Low-pass filter; \f$w_c\f$ is computed as above with \c a1
 * treated as 1/4.</td></tr>
 * <tr><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\circ\f$</td><td>\f$\circ\f$</td><td>\f$\bullet\f$</td><td>
 * Ditto; \c a2 is ignored.</td></tr>
 * <tr><td>\f$\bullet\f$</td><td>\f$\circ\f$</td><td>\f$\circ\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>
 * High-pass filter; \f$w_c\f$ is computed from \c nMax,
 * \c f2, and \c a2.</td></tr>
 * <tr><td>\f$\bullet\f$</td><td>\f$\circ\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>
 * Ditto; \c a1 is ignored.</td></tr>
 * <tr><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\circ\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>
 * Ditto; \c f1 is ignored.</td></tr>
 * <tr><td>\f$\bullet\f$</td><td>\f$\circ\f$</td><td>\f$\circ\f$</td><td>\f$\bullet\f$</td><td>\f$\circ\f$</td><td>
 * High-pass filter; \f$w_c\f$ is computed as above with \c a2
 * treated as 1/4.</td></tr>
 * <tr><td>\f$\bullet\f$</td><td>\f$\circ\f$</td><td>\f$\bullet\f$</td><td>\f$\bullet\f$</td><td>\f$\circ\f$</td><td>
 * Ditto; \c a1 is ignored.</td></tr>
 * <tr><td colspan=5>Other</td><td>Subroutine returns an error.</td></tr>
 * </table>
 *
 * Once an order \f$n\f$ and characteristic frequency \f$w_c\f$ are known, the
 * zeros and poles of a ZPG filter are readily determined.  A stable,
 * physically realizable Butterworth filter will have \f$n\f$ poles evenly
 * spaced on the upper half of a circle of radius \f$w_c\f$; that is,
 * \f[
 * R = \frac{(-iw_c)^n}{\prod_{k=0}^{n-1}(w - w_c e^{2\pi i(k+1/2)/n})}
 * \f]
 * for a low-pass filter, and
 * \f[
 * R = \frac{w^n}{\prod_{k=0}^{n-1}(w - w_c e^{2\pi i(k+1/2)/n})}
 * \f]
 * for a high-pass filter.  By choosing only poles on the upper-half
 * plane, one ensures that after transforming to \f$z\f$ the poles will have
 * \f$|z|<1\f$.  Furthermore, the phase factor \f$(-i)^n\f$ in the numerator of
 * the low-pass filter is chosen so that the DC response is purely real;
 * this ensures that the response function in the \f$z\f$-plane will have a
 * real gain factor, and the resulting IIR filter will be physically
 * realizable.  The high-pass filter has a purely real response at
 * Nyquist (\f$w\rightarrow\infty\f$), which similarly gives a physical IIR
 * filter.
 *
 * Although higher orders \f$n\f$ would appear to produce better (i.e.
 * sharper) filter responses, one rapidly runs into numerical errors, as
 * one ends up adding and subtracting \f$n\f$ large numbers to obtain small
 * filter responses.  One way around this is to break the filter up into
 * several lower-order filters.  The routines in this module do just
 * that.  Poles are paired up across the imaginary axis, (and combined
 * with pairs of zeros at \f$w=0\f$ for high-pass filters,) to form \f$[n/2]\f$
 * second-order filters.  If \f$n\f$ is odd, there will be an additional
 * first-order filter, with one pole at \f$w=iw_c\f$ (and one zero at \f$w=0\f$
 * for a high-pass filter).
 *
 * Each ZPG filter in the \f$w\f$-plane is first transformed to the \f$z\f$-plane
 * by a bilinear transformation, and is then used to construct a
 * time-domain IIR filter.  Each filter is then applied to the time
 * series.  As mentioned in the description above, the filters are
 * designed to give an overall amplitude response that is the square root
 * of the desired attenuation; however, each time-domain filter is
 * applied to the data stream twice: once in the normal sense, and once
 * in the time-reversed sense.  This gives the full attenuation with very
 * little frequency-dependent phase shift.
 *
 */
/*@{*/

/* Prototype for a local input parsing routine. */
static int
XLALParsePassBandParamStruc( PassBandParamStruc *params,
			 INT4               *n,
			 REAL8              *wc,
			 REAL8              deltaT );


/** \see See \ref ButterworthTimeSeries_c for documentation */
int XLALButterworthREAL4TimeSeries( REAL4TimeSeries *series, PassBandParamStruc *params )
{
  INT4 n;    /* The filter order. */
  INT4 type; /* The pass-band type: high, low, or undeterminable. */
  INT4 i;    /* An index. */
  INT4 j;    /* Another index. */
  REAL8 wc;  /* The filter's transformed frequency. */

  /* Make sure the input pointers are non-null. */
  if ( ! params || ! series || ! series->data || ! series->data->data )
    XLAL_ERROR( XLAL_EFAULT );

  /* Parse the pass-band parameter structure.  I separate this into a
     local static subroutine because it's an icky mess of conditionals
     that would clutter the logic of the main routine. */
  type=XLALParsePassBandParamStruc(params,&n,&wc,series->deltaT);
  if(type<0)
    XLAL_ERROR( XLAL_EINVAL );

  /* An order n Butterworth filter has n poles spaced evenly along a
     semicircle in the upper complex w-plane.  By pairing up poles
     symmetric across the imaginary axis, the filter gan be decomposed
     into [n/2] filters of order 2, plus perhaps an additional order 1
     filter.  The following loop pairs up poles and applies the
     filters with order 2. */
  for(i=0,j=n-1;i<j;i++,j--){
    REAL8 theta=LAL_PI*(i+0.5)/n;
    REAL8 ar=wc*cos(theta);
    REAL8 ai=wc*sin(theta);
    REAL8IIRFilter *iirFilter=NULL;
    COMPLEX16ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      zpgFilter = XLALCreateCOMPLEX16ZPGFilter(2,2);
      if ( ! zpgFilter )
        XLAL_ERROR( XLAL_EFUNC );
      zpgFilter->zeros->data[0]=0.0;
      zpgFilter->zeros->data[1]=0.0;
      zpgFilter->gain=1.0;
    }else{
      zpgFilter = XLALCreateCOMPLEX16ZPGFilter(0,2);
      if ( ! zpgFilter )
        XLAL_ERROR( XLAL_EFUNC );
      zpgFilter->gain=-wc*wc;
    }
    zpgFilter->poles->data[0]=ar;
    zpgFilter->poles->data[0]+=ai*I;
    zpgFilter->poles->data[1]=-ar;
    zpgFilter->poles->data[1]+=ai*I;

    /* Transform to the z-plane and create the IIR filter. */
    if (XLALWToZCOMPLEX16ZPGFilter(zpgFilter)<0)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLAL_ERROR( XLAL_EFUNC );
    }
    iirFilter = XLALCreateREAL8IIRFilter(zpgFilter);
    if (!iirFilter)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLAL_ERROR( XLAL_EFUNC );
    }

    /* Filter the data, once each way. */
    if (XLALIIRFilterREAL4Vector(series->data,iirFilter)<0
        || XLALIIRFilterReverseREAL4Vector(series->data,iirFilter)<0)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLALDestroyREAL8IIRFilter(iirFilter);
      XLAL_ERROR( XLAL_EFUNC );
    }

    /* Free the filters. */
    XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
    XLALDestroyREAL8IIRFilter(iirFilter);
  }

  /* Next, this conditional applies the possible order 1 filter
     corresponding to an unpaired pole on the imaginary w axis. */
  if(i==j){
    REAL8IIRFilter *iirFilter=NULL;
    COMPLEX16ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      zpgFilter=XLALCreateCOMPLEX16ZPGFilter(1,1);
      if(!zpgFilter)
        XLAL_ERROR(XLAL_EFUNC);
      *zpgFilter->zeros->data=0.0;
      zpgFilter->gain=1.0;
    }else{
      zpgFilter=XLALCreateCOMPLEX16ZPGFilter(0,1);
      if(!zpgFilter)
        XLAL_ERROR(XLAL_EFUNC);
      zpgFilter->gain=-wc*I;
    }
    *zpgFilter->poles->data=wc*I;

    /* Transform to the z-plane and create the IIR filter. */
    if (XLALWToZCOMPLEX16ZPGFilter(zpgFilter)<0)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLAL_ERROR(XLAL_EFUNC);
    }
    iirFilter=XLALCreateREAL8IIRFilter(zpgFilter);
    if (!iirFilter)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLAL_ERROR(XLAL_EFUNC);
    }

    /* Filter the data, once each way. */
    if (XLALIIRFilterREAL4Vector(series->data,iirFilter)<0
        || XLALIIRFilterReverseREAL4Vector(series->data,iirFilter)<0)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLALDestroyREAL8IIRFilter(iirFilter);
      XLAL_ERROR( XLAL_EFUNC );
    }

    /* Free the filters. */
    XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
    XLALDestroyREAL8IIRFilter(iirFilter);
  }

  return 0;
}

/** \see See \ref ButterworthTimeSeries_c for documentation */
int XLALButterworthREAL8TimeSeries( REAL8TimeSeries *series, PassBandParamStruc *params )
{
  INT4 n;    /* The filter order. */
  INT4 type; /* The pass-band type: high, low, or undeterminable. */
  INT4 i;    /* An index. */
  INT4 j;    /* Another index. */
  REAL8 wc;  /* The filter's transformed frequency. */

  /* Make sure the input pointers are non-null. */
  if ( ! params || ! series || ! series->data || ! series->data->data )
    XLAL_ERROR( XLAL_EFAULT );

  /* Parse the pass-band parameter structure.  I separate this into a
     local static subroutine because it's an icky mess of conditionals
     that would clutter the logic of the main routine. */
  type=XLALParsePassBandParamStruc(params,&n,&wc,series->deltaT);
  if(type<0)
    XLAL_ERROR( XLAL_EINVAL );

  /* An order n Butterworth filter has n poles spaced evenly along a
     semicircle in the upper complex w-plane.  By pairing up poles
     symmetric across the imaginary axis, the filter gan be decomposed
     into [n/2] filters of order 2, plus perhaps an additional order 1
     filter.  The following loop pairs up poles and applies the
     filters with order 2. */
  for(i=0,j=n-1;i<j;i++,j--){
    REAL8 theta=LAL_PI*(i+0.5)/n;
    REAL8 ar=wc*cos(theta);
    REAL8 ai=wc*sin(theta);
    REAL8IIRFilter *iirFilter=NULL;
    COMPLEX16ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      zpgFilter = XLALCreateCOMPLEX16ZPGFilter(2,2);
      if ( ! zpgFilter )
        XLAL_ERROR( XLAL_EFUNC );
      zpgFilter->zeros->data[0]=0.0;
      zpgFilter->zeros->data[1]=0.0;
      zpgFilter->gain=1.0;
    }else{
      zpgFilter = XLALCreateCOMPLEX16ZPGFilter(0,2);
      if ( ! zpgFilter )
        XLAL_ERROR( XLAL_EFUNC );
      zpgFilter->gain=-wc*wc;
    }
    zpgFilter->poles->data[0]=ar;
    zpgFilter->poles->data[0]+=ai*I;
    zpgFilter->poles->data[1]=-ar;
    zpgFilter->poles->data[1]+=ai*I;

    /* Transform to the z-plane and create the IIR filter. */
    if (XLALWToZCOMPLEX16ZPGFilter(zpgFilter)<0)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLAL_ERROR( XLAL_EFUNC );
    }
    iirFilter = XLALCreateREAL8IIRFilter(zpgFilter);
    if (!iirFilter)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLAL_ERROR( XLAL_EFUNC );
    }

    /* Filter the data, once each way. */
    if (XLALIIRFilterREAL8Vector(series->data,iirFilter)<0
        || XLALIIRFilterReverseREAL8Vector(series->data,iirFilter)<0)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLALDestroyREAL8IIRFilter(iirFilter);
      XLAL_ERROR( XLAL_EFUNC );
    }

    /* Free the filters. */
    XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
    XLALDestroyREAL8IIRFilter(iirFilter);
  }

  /* Next, this conditional applies the possible order 1 filter
     corresponding to an unpaired pole on the imaginary w axis. */
  if(i==j){
    REAL8IIRFilter *iirFilter=NULL;
    COMPLEX16ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      zpgFilter=XLALCreateCOMPLEX16ZPGFilter(1,1);
      if(!zpgFilter)
        XLAL_ERROR(XLAL_EFUNC);
      *zpgFilter->zeros->data=0.0;
      zpgFilter->gain=1.0;
    }else{
      zpgFilter=XLALCreateCOMPLEX16ZPGFilter(0,1);
      if(!zpgFilter)
        XLAL_ERROR(XLAL_EFUNC);
      zpgFilter->gain=-wc*I;
    }
    *zpgFilter->poles->data=wc*I;

    /* Transform to the z-plane and create the IIR filter. */
    if (XLALWToZCOMPLEX16ZPGFilter(zpgFilter)<0)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLAL_ERROR(XLAL_EFUNC);
    }
    iirFilter=XLALCreateREAL8IIRFilter(zpgFilter);
    if (!iirFilter)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLAL_ERROR(XLAL_EFUNC);
    }

    /* Filter the data, once each way. */
    if (XLALIIRFilterREAL8Vector(series->data,iirFilter)<0
        || XLALIIRFilterReverseREAL8Vector(series->data,iirFilter)<0)
    {
      XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
      XLALDestroyREAL8IIRFilter(iirFilter);
      XLAL_ERROR( XLAL_EFUNC );
    }

    /* Free the filters. */
    XLALDestroyCOMPLEX16ZPGFilter(zpgFilter);
    XLALDestroyREAL8IIRFilter(iirFilter);
  }

  return 0;
}


/** \see See \ref ButterworthTimeSeries_c for documentation */
int XLALLowPassREAL4TimeSeries( REAL4TimeSeries *series,
    REAL8 frequency, REAL8 amplitude, INT4 filtorder )
{
  PassBandParamStruc params;
  params.nMax = filtorder;
  params.f1   = frequency;
  params.a1   = amplitude;
  params.f2   = -1;
  params.a2   = -1;
  if ( XLALButterworthREAL4TimeSeries( series, &params ) < 0 )
    XLAL_ERROR( XLAL_EFUNC );
  return 0;
}

/** \see See \ref ButterworthTimeSeries_c for documentation */
int XLALLowPassREAL8TimeSeries( REAL8TimeSeries *series,
    REAL8 frequency, REAL8 amplitude, INT4 filtorder )
{
  PassBandParamStruc params;
  params.nMax = filtorder;
  params.f1   = frequency;
  params.a1   = amplitude;
  params.f2   = -1;
  params.a2   = -1;
  if ( XLALButterworthREAL8TimeSeries( series, &params ) < 0 )
    XLAL_ERROR( XLAL_EFUNC );
  return 0;
}

/** \see See \ref ButterworthTimeSeries_c for documentation */
int XLALHighPassREAL4TimeSeries( REAL4TimeSeries *series,
    REAL8 frequency, REAL8 amplitude, INT4 filtorder )
{
  PassBandParamStruc params;
  params.nMax = filtorder;
  params.f2   = frequency;
  params.a2   = amplitude;
  params.f1   = -1;
  params.a1   = -1;
  if ( XLALButterworthREAL4TimeSeries( series, &params ) < 0 )
    XLAL_ERROR( XLAL_EFUNC );
  return 0;
}

/** \see See \ref ButterworthTimeSeries_c for documentation */
int XLALHighPassREAL8TimeSeries( REAL8TimeSeries *series,
    REAL8 frequency, REAL8 amplitude, INT4 filtorder )
{
  PassBandParamStruc params;
  params.nMax = filtorder;
  params.f2   = frequency;
  params.a2   = amplitude;
  params.f1   = -1;
  params.a1   = -1;
  if ( XLALButterworthREAL8TimeSeries( series, &params ) < 0 )
    XLAL_ERROR( XLAL_EFUNC );
  return 0;
}

/**
 * Deprecated.
 * \deprecated Use XLALButterworthREAL4TimeSeries() instead.
 */
void
LALButterworthREAL4TimeSeries( LALStatus          *stat,
			       REAL4TimeSeries    *series,
			       PassBandParamStruc *params )
{
  INT4 n;    /* The filter order. */
  INT4 type; /* The pass-band type: high, low, or undeterminable. */
  INT4 i;    /* An index. */
  INT4 j;    /* Another index. */
  REAL8 wc;  /* The filter's transformed frequency. */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR(stat);

  /* Make sure the input pointers are non-null. */
  ASSERT(params,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);
  ASSERT(series,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);
  ASSERT(series->data,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);
  ASSERT(series->data->data,stat,BANDPASSTIMESERIESH_ENUL,
	 BANDPASSTIMESERIESH_MSGENUL);

  /* Parse the pass-band parameter structure.  I separate this into a
     local static subroutine because it's an icky mess of conditionals
     that would clutter the logic of the main routine. */
  type=XLALParsePassBandParamStruc(params,&n,&wc,series->deltaT);
  if(type<0) {
    XLALClearErrno();
    ABORT(stat,BANDPASSTIMESERIESH_EBAD,BANDPASSTIMESERIESH_MSGEBAD);
  }

  /* An order n Butterworth filter has n poles spaced evenly along a
     semicircle in the upper complex w-plane.  By pairing up poles
     symmetric across the imaginary axis, the filter gan be decomposed
     into [n/2] filters of order 2, plus perhaps an additional order 1
     filter.  The following loop pairs up poles and applies the
     filters with order 2. */
  for(i=0,j=n-1;i<j;i++,j--){
    REAL4 theta=LAL_PI*(i+0.5)/n;
    REAL4 ar=wc*cos(theta);
    REAL4 ai=wc*sin(theta);
    REAL4IIRFilter *iirFilter=NULL;
    COMPLEX8ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      TRY(LALCreateCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter,2,2),
	  stat);
      zpgFilter->zeros->data[0]=0.0;
      zpgFilter->zeros->data[1]=0.0;
      zpgFilter->gain=1.0;
    }else{
      TRY(LALCreateCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter,0,2),
	  stat);
      zpgFilter->gain=-wc*wc;
    }
    zpgFilter->poles->data[0]=ar;
    zpgFilter->poles->data[0]+=ai*I;
    zpgFilter->poles->data[1]=-ar;
    zpgFilter->poles->data[1]+=ai*I;

    /* Transform to the z-plane and create the IIR filter. */
    LALWToZCOMPLEX8ZPGFilter(stat->statusPtr,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);
    LALCreateREAL4IIRFilter(stat->statusPtr,&iirFilter,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);

    /* Filter the data, once each way. */
    LALIIRFilterREAL4Vector(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL4IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);
    LALIIRFilterREAL4VectorR(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL4IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);

    /* Free the filters. */
    TRY(LALDestroyREAL4IIRFilter(stat->statusPtr,&iirFilter),stat);
    TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),stat);
  }

  /* Next, this conditional applies the possible order 1 filter
     corresponding to an unpaired pole on the imaginary w axis. */
  if(i==j){
    REAL4IIRFilter *iirFilter=NULL;
    COMPLEX8ZPGFilter *zpgFilter=NULL;

    /* Generate the filter in the w-plane. */
    if(type==2){
      TRY(LALCreateCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter,1,1),
	  stat);
      *zpgFilter->zeros->data=0.0;
      zpgFilter->gain=1.0;
    }else{
      TRY(LALCreateCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter,0,1),
	  stat);
      zpgFilter->gain=-wc*I;
    }
    *zpgFilter->poles->data=wc*I;

    /* Transform to the z-plane and create the IIR filter. */
    LALWToZCOMPLEX8ZPGFilter(stat->statusPtr,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);
    LALCreateREAL4IIRFilter(stat->statusPtr,&iirFilter,zpgFilter);
    BEGINFAIL(stat)
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
    ENDFAIL(stat);

    /* Filter the data, once each way. */
    LALIIRFilterREAL4Vector(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL4IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);
    LALIIRFilterREAL4VectorR(stat->statusPtr,series->data,iirFilter);
    BEGINFAIL(stat) {
      TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),
	  stat);
      TRY(LALDestroyREAL4IIRFilter(stat->statusPtr,&iirFilter),stat);
    } ENDFAIL(stat);

    /* Free the filters. */
    TRY(LALDestroyREAL4IIRFilter(stat->statusPtr,&iirFilter),stat);
    TRY(LALDestroyCOMPLEX8ZPGFilter(stat->statusPtr,&zpgFilter),stat);
  }

  /* Normal exit. */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


/**
 * Deprecated.
 * \deprecated Use XLALButterworthREAL8TimeSeries() instead.
 */
void
LALButterworthREAL8TimeSeries( LALStatus          *stat,
			       REAL8TimeSeries    *series,
			       PassBandParamStruc *params )
{
  INITSTATUS(stat);

  if (XLALButterworthREAL8TimeSeries(series,params)<0)
  {
    int code=xlalErrno;
    XLALClearErrno();
    switch(code){
      case XLAL_EFAULT:
        ABORT(stat,BANDPASSTIMESERIESH_ENUL,BANDPASSTIMESERIESH_MSGENUL);
      case XLAL_EINVAL:
        ABORT(stat,BANDPASSTIMESERIESH_EBAD,BANDPASSTIMESERIESH_MSGEBAD);
      default:
        ABORTXLAL(stat);
    }
  }

  RETURN(stat);
}

/**
 * Deprecated.
 * \deprecated
 */
void
LALDButterworthREAL4TimeSeries( LALStatus          *stat,
				REAL4TimeSeries    *series,
				PassBandParamStruc *params )
{
  INITSTATUS(stat);

  if (XLALButterworthREAL4TimeSeries(series,params)<0)
  {
    int code=xlalErrno;
    XLALClearErrno();
    switch(code){
      case XLAL_EFAULT:
        ABORT(stat,BANDPASSTIMESERIESH_ENUL,BANDPASSTIMESERIESH_MSGENUL);
      case XLAL_EINVAL:
        ABORT(stat,BANDPASSTIMESERIESH_EBAD,BANDPASSTIMESERIESH_MSGEBAD);
      default:
        ABORTXLAL(stat);
    }
  }

  RETURN(stat);
}


static int
XLALParsePassBandParamStruc( PassBandParamStruc *params,
			 INT4               *n,
			 REAL8              *wc,
			 REAL8              deltaT )
     /* This local function parses the pass band parameters according
        to the rules given in the documentation to this module,
        computing the order and characteristic frequency (in the
        w-variable) of the desired filter.  The function returns 2 for
        a high-pass filter, 1 for a low-pass filter, and -1 if the
        parameters were poorly specified. */
{
  /* Okay, first define all the temporary variables, and figure out
     just which parameter combinations have been given. */
  REAL8 w1=deltaT*params->f1;    /* Dimensionless frequencies; */
  REAL8 w2=deltaT*params->f2;    /* later transformed to w-plane. */
  REAL8 a1=params->a1;           /* Attenuation factors; later a */
  REAL8 a2=params->a2;           /* square root will be taken. */
  BOOLEAN w1Given = (w1>0.0)&&(w1<0.5);
  BOOLEAN w2Given = (w2>0.0)&&(w2<0.5);
  BOOLEAN both1Given = w1Given && (a1>0.0)&&(a1<1.0);
  BOOLEAN both2Given = w2Given && (a2>0.0)&&(a2<1.0);

  /* If both frequencies and both attenuations have been given,
     compute the required filter order and characteristic frequency
     from them. */
  if(both1Given && both2Given){
    REAL8 aLow  = w1<w2 ? sqrt(a1) : sqrt(a2);
    REAL8 aHigh = w1>w2 ? sqrt(a1) : sqrt(a2);
    REAL8 wLow  = w1<w2 ? tan(LAL_PI*w1) : tan(LAL_PI*w2);
    REAL8 wHigh = w1>w2 ? tan(LAL_PI*w1) : tan(LAL_PI*w2);

    /* First make sure that two different frequencies and attenuations
       have been specified. */
    if(w1==w2){
      XLALPrintInfo("XLAL Info - %s: ", __func__ );
      XLALPrintInfo("The two frequencies defining the transition band"
	      " are the same.\n");
      XLAL_ERROR(XLAL_EINVAL);
    }
    if(a1==a2){
      XLALPrintInfo("XLAL Info - %s: ", __func__ );
      XLALPrintInfo("The two attenuations across the transition band"
	      " are the same.\n");
      XLAL_ERROR(XLAL_EINVAL);
    }

    /* Now compute the filter order. */
    if(aHigh>aLow) /* High-pass. */
      *n=(INT4)(0.5*log((1.0/aHigh-1.0)/(1.0/aLow-1.0))
		/log(wLow/wHigh))+1;
    else /* Low-pass. */
      *n=(INT4)(0.5*log((1.0/aHigh-1.0)/(1.0/aLow-1.0))
		/log(wHigh/wLow))+1;

    /* If a positive params->nMax less than *n has been specified,
       reduce to that order, with appropriate warnings. */
    if((params->nMax>0)&&(params->nMax<*n)){
      XLALPrintWarning("XLAL Warning - %s: ", __func__ );
      XLALPrintWarning("Filter order required to achieve requested"
		 " performance exceeds\n"
		 "\tspecified limit.\n");
      XLALPrintWarning("\tRequired: %i  Limit: %i\n",*n,params->nMax);
      *n=params->nMax;
    }

    /* Compute the characteristic frequency of the filter using the
       filter order and the minimal attenuation in the pass-band. */
    if(aHigh>aLow){ /* High-pass */
      *wc=wHigh*pow(1.0/aHigh - 1.0,-0.5/((REAL8)(*n)));
      return 2;
    }else{ /* Low-pass */
      *wc=wLow*pow(1.0/aLow - 1.0,0.5/((REAL8)(*n)));
      return 1;
    }
  }

  /* If one end of the transition band has both frequency and
     attenuation specified, and the other does not, use only the end
     that does.  The filter order must also be specified for this and
     all future cases. */
  else{
    if(params->nMax<=0){
      XLALPrintInfo("XLAL Info - %s: ", __func__ );
      XLALPrintInfo("The filter order must be given if one or more"
	      " frequencies or\n"
	      "\tattenuations of the transition band have not been"
	      " (validly) specified.");
      XLAL_ERROR(XLAL_EINVAL);
    }
    if(both1Given){
      REAL8 wLow=tan(LAL_PI*w1);
      REAL8 aLow=sqrt(a1);
      *n=params->nMax;
      *wc=wLow*pow(1.0/aLow - 1.0, -0.5/((REAL8)(*n)));
      return 1;
    }else if(both2Given){
      REAL8 wHigh=tan(LAL_PI*w2);
      REAL8 aHigh=sqrt(a2);
      *n=params->nMax;
      *wc=wHigh*pow(1.0/aHigh - 1.0, 0.5/((REAL8)(*n)));
      return 2;
    }

    /* If no attenuations are given, then there had better be only one
       frequency given, otherwise we don't know whether to make a low-
       or a high-pass filter. */
    else if(w1Given && w2Given){
      XLALPrintInfo("XLAL Info - %s: ", __func__ );
      XLALPrintInfo("Neither attenuation has been specified, so only"
	      " one frequency should be.");
      XLAL_ERROR(XLAL_EINVAL);
    }

    /* Treat the given frequency as the characteristic frequency and
       the specified nMax as the filter order. */
    else{
      *n=params->nMax;
      if(w2Given){
	*wc=tan(LAL_PI*w2);
	return 2;
      }else if(w1Given){
	*wc=tan(LAL_PI*w1);
	return 1;
      }else{
        XLALPrintInfo("XLAL Info - %s: ", __func__ );
	XLALPrintInfo("No frequencies within the Nyquist band have"
		" been specified!");
	XLAL_ERROR(XLAL_EINVAL);
      }
    }
  }
}
/*@}*/
