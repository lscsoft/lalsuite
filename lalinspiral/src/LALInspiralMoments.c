/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, Benjamin Owen, B.S. Sathyaprakash, Thomas Cokelaer
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

/** \defgroup LALInspiralMoments_c Module LALInspiralMoments.c
 * \ingroup LALInspiralBank_h
 * \brief Functions to calculate the moment of the noise power spectral density.
 * \author Brown, D. A., Cokelaer, T. and  Sathyaprakash, B. S.

The moments of the noise curve are defined as
\f{equation}{
I(q)  \equiv S_{h}(f_{0}) \int^{f_{c}/f_{0}}_{f_{s}/f_{0}}
\frac{x^{-q}}{S_{h}(x)} \, dx \,.
\f}
Because in practice we will always divide one of these moments by another, we
do not need to include the \f$S_{h}(f_{0})\f$ term, which always cancels.
This function calculates the integral
\f{equation}{
I = \int^{f_{c}/f_{0}}_{f_{s}/f_{0}} \frac{x^{-q}}{S_{h}(x)} \, dx \,.
\f}
It then divides this quantity by a normalisation constant which has been
passed to the function. In the case of calculating the components of the
metric for the signal manifold for the purpose of generating a template bank,
this constant is given by \f$I(7)\f$, because of the definition of the quantity
\f{equation}{
J(q) \equiv \frac{I(q)}{I(7/3)} \,.
\f}

\heading{Algorithm}
Given the exponent <tt>pars.ndx</tt> and limits of integration
<tt>pars.xmin</tt> and <tt>pars.xmax</tt> this function returns the moment of
the power spectral density specified by the frequency series
<tt>pars.shf</tt> according to
\f{equation}{
\mathtt{moment} = \int_{\mathtt{xmin}}^{\mathtt{xmax}}
\frac{x^{-\mathtt{ndx}}}{S_h(x)}\, dx \, .
\f}
*/

/*@{*/

#include <lal/LALInspiralBank.h>
#include <lal/Integrate.h>

/* Deprecation Warning */

/** \see See \ref LALInspiralMoments_c for documentation */
void
LALGetInspiralMoments (
    LALStatus            *status,
    InspiralMomentsEtc   *moments,
    REAL8FrequencySeries *psd,
    InspiralTemplate     *params
    )

{
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );
  XLALPrintDeprecationWarning("LALGetInspiralMoments", "XLALGetInspiralMoments");

  if (XLALGetInspiralMoments(moments, params->fLower, params->fCutoff, psd)!= XLAL_SUCCESS){
     ABORTXLAL( status );
  }
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

int
XLALGetInspiralMoments (
    InspiralMomentsEtc   *moments,
    REAL8 fLower,
    REAL8 fCutoff,
    REAL8FrequencySeries *psd
    )
{
  size_t k;
  REAL8 xmin;
  REAL8 xmax;
  REAL8 ndx;
  REAL8 norm;

  /* Check inputs */
  if (!moments){
    XLALPrintError("Moments is NULL\n");
    XLAL_ERROR(XLAL_EFAULT);
  }
  if (!psd){
    XLALPrintError("PSD is NULL\n");
    XLAL_ERROR(XLAL_EFAULT);
  }

  if (fLower <= 0 || fCutoff <= fLower){
    XLALPrintError("fLower must be between 0 and fCutoff\n");
    XLAL_ERROR(XLAL_EDOM);
  };

  /* Constants needed in computing the moments */
  moments->a01 = 3.L/5.L;
  moments->a21 = 11.L * LAL_PI/12.L;
  moments->a22 = 743.L/2016.L * cbrt(25.L/(2.L*LAL_PI*LAL_PI));
  moments->a31 = -3.L/2.L;
  moments->a41 = 617.L * LAL_PI * LAL_PI / 384.L;
  moments->a42 = 5429.L/5376.L * cbrt(25.L*LAL_PI/2.L);
  moments->a43 = 1.5293365L/1.0838016L * cbrt(5.L/(4.L*LAL_PI*LAL_PI*LAL_PI*LAL_PI));

  /* Divide all frequencies by fLower, a scaling that is used in solving */
  /* the moments integral                                                */
  psd->f0 /= fLower;
  psd->deltaF /= fLower;
  xmin = fLower / fLower;
  xmax = fCutoff / fLower;

  /* First compute the norm and print if requested */
  norm = 1.L;
  ndx = 7.L/3.L;
  moments->j[7]=XLALInspiralMoments(xmin, xmax, ndx, norm, psd);
  if (XLAL_IS_REAL8_FAIL_NAN(moments->j[7])){
    XLAL_ERROR(XLAL_EFUNC);
  }
  norm = moments->j[7];

  /* Then compute the normalised moments of the noise PSD from 1/3 to 17/3. */
  for ( k = 1; k <= 17; ++k )
  {
    ndx = (REAL8) k / 3.L;
    moments->j[k]=XLALInspiralMoments(xmin, xmax, ndx, norm, psd);
  }

  /* Moments are done: Rescale deltaF and f0 back to their original values */
  psd->deltaF *= fLower;
  psd->f0 *= fLower;

  return XLAL_SUCCESS;
}

/** \see See \ref LALInspiralMoments_c for documentation */
void
LALGetInspiralMomentsBCV (
    LALStatus               *status,
    InspiralMomentsEtcBCV   *moments,
    REAL8FrequencySeries    *psd,
    InspiralTemplate        *params
    )
{
  UINT4 k;
  InspiralMomentsIn in;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* doesn't seem to be needed. thomas, janvier 2004. I prefer to remove it for the moment.
   *  The factor is not important in the case of SPA approximation but is important in BCV
   *  case. Indeed on one hand we use quantity which are a ratio between two moments and
   *  consequently a factor 1 or 2 is not important. Howver in the case of BCV, we might
   *  use a moment alone. Thus a factor in the computation has an effect. */

  /*  for (i=0; i< psd->data->length ; i++)
  {
    psd->data->data[i] = psd->data->data[i] * 1e45;
  }
   */
  in.shf = psd;
  in.xmin = params->fLower;
  in.xmax = params->fCutoff;

  /* First compute the norm */
  in.norm = 1.L;
  for ( k = 0; k <= 22; ++k )
  {
    if (k <= 17)
    {
      /* positive value*/
      in.ndx = (REAL8)k / 3.L;
    }
    else
    {
      /* negative -1,-2 ...-6 */
      in.ndx = (17.- (REAL8)k) /3.L;
    }

    LALInspiralMoments( status->statusPtr, &moments->i[k], in );
    CHECKSTATUSPTR(status);
  }

  in.norm = moments->i[7] -2.*moments->alpha * moments->i[5] +
    moments->alpha * moments->alpha*moments->i[3];


  /* 17 */
  moments->M1[0][0] = (moments->i[17] -2.*moments->alpha * moments->i[15] +
      moments->alpha * moments->alpha*moments->i[13]) / in.norm;
  /* 14 */
  moments->M1[0][1] = (moments->i[14] -2.*moments->alpha * moments->i[12] +
      moments->alpha * moments->alpha*moments->i[10]) / in.norm;
  /* 11 */
  moments->M1[1][1] = (moments->i[11] -2.*moments->alpha * moments->i[9] +
      moments->alpha * moments->alpha*moments->i[7]) / in.norm;

  moments->M1[1][0]=moments->M1[0][1] ;

  /*  12 */
  moments->M2[0][0] = (moments->i[12] -2.*moments->alpha * moments->i[10] +
      moments->alpha * moments->alpha*moments->i[8]) / in.norm;
  /* 9 */

  moments->M2[0][1] = (moments->i[9] -2.*moments->alpha * moments->i[7] +
      moments->alpha * moments->alpha*moments->i[5]) / in.norm;
  /*  9 */

  moments->M2[1][0] = (moments->i[9] -2.*moments->alpha * moments->i[7] +
      moments->alpha * moments->alpha*moments->i[5]) / in.norm;
  /*  6 */
  moments->M2[1][1] = (moments->i[6] -2.*moments->alpha * moments->i[4] +
      moments->alpha * moments->alpha*moments->i[2]) / in.norm;

  /* 7 */
  moments->M3[0][0] = (moments->i[7] -2.*moments->alpha * moments->i[5] +
      moments->alpha * moments->alpha*moments->i[3]) / in.norm;
  /* 4 */
  moments->M3[0][1] = (moments->i[4] -2.*moments->alpha * moments->i[2] +
      moments->alpha * moments->alpha*moments->i[0]) / in.norm;
  /* 1 */
  moments->M3[1][1] = (moments->i[1] -2.*moments->alpha * moments->i[18] +
      moments->alpha * moments->alpha * moments->i[20]) / in.norm;

  moments->M3[1][0]=moments->M3[0][1] ;

  if ( lalDebugLevel & LALINFO )
  {
    LALPrintError( "#M1=\n");
    LALPrintError( "#%15.12lf %15.12lf \n# %15.12lf %15.12lf\n",
        moments->M1[0][0],
        moments->M1[0][1],
        moments->M1[1][0],
        moments->M1[1][1] );

    LALPrintError( "#M2=\n" );
    LALPrintError( "#%15.12lf %15.12lf \n# %15.12lf %15.12lf\n",
        moments->M2[0][0],
        moments->M2[0][1],

        moments->M2[1][0],
        moments->M2[1][1] );

    LALPrintError( "#M3=\n" );
    LALPrintError( "#%15.12lf %15.12lf \n# %15.12lf %15.12lf\n",
        moments->M3[0][0],
        moments->M3[0][1],
        moments->M3[1][0],
        moments->M3[1][1] );
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/** \see See \ref LALInspiralMoments_c for documentation */
void
LALInspiralMoments(
    LALStatus         *status,	/**< LAL status pointer */
    REAL8             *moment,	/**< [out] the value of the moment */
    InspiralMomentsIn  pars	/**< [in] input parameters */
    )

{
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALInspiralMoments", "XLALInspiralMoments");

  *moment = XLALInspiralMoments(pars.xmin, pars.xmax, pars.ndx, pars.norm, pars.shf);
  if (XLAL_IS_REAL8_FAIL_NAN(*moment)){
    ABORTXLAL( status );
  };
  RETURN (status);
}

REAL8
XLALInspiralMoments(
    REAL8 xmin,
    REAL8 xmax,
    REAL8 ndx,
    REAL8 norm,
    REAL8FrequencySeries *shf
    )

{
  REAL8 moment = 0;
  REAL8 f0, deltaF;
  size_t k, kMin, kMax;

  /* Check inputs */
  if (!shf || !(shf->data) || !(shf->data->data)) {
    XLALPrintError("PSD or its data are NULL\n");
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  }

  if (xmin <= 0 || xmax <= 0 || xmax <= xmin || norm <= 0) {
    XLALPrintError("xmin, xmax, and norm must be positive and xmax must be greater than xmin\n");
    XLAL_ERROR_REAL8(XLAL_EDOM);
  }

  /* set up and check domain of integration */
  /* NB: Although these are called f0 and deltaF, they are really supposed to
         be x0 and deltaX (x = f / f0). That is, you either need to have hacked
         the PSD's f0 and deltaF values before calling this function or be
         prepared to rescale the outputs. */
  f0 = shf->f0;
  deltaF = shf->deltaF;
  kMax = floor((xmax - f0) / deltaF);
  if ( (xmin < f0) || (kMax > shf->data->length) ) {
    XLALPrintError("PSD does not cover domain of integration\n");
    XLAL_ERROR_REAL8(XLAL_EDOM);
  }
  kMin = floor((xmin - f0) / deltaF);

  /* do the first point of the integral */
  if( shf->data->data[kMin] ) {
    const REAL8 f = f0 + kMin * deltaF;
    moment += pow( f, -(ndx) ) / ( 2.0 * shf->data->data[kMin] );
  }
  /* do the bulk of the integral */
  for ( k = kMin + 1; k < kMax; ++k ) {
    const REAL8 psd_val = shf->data->data[k];
    if ( psd_val ) {
      const REAL8 f = f0 + k * deltaF;
      moment += pow( f, -(ndx) ) / psd_val;
    }
  }
  /* Do the last point of the integral, but allow the integration domain
     to be open on the right if necessary. */
  if ( kMax < shf->data->length && shf->data->data[kMax] ) {
    const REAL8 f = f0 + kMax * deltaF;
    moment += pow( f, -(ndx) ) / ( 2.0 * shf->data->data[kMax] );
  }
  moment *= deltaF;

  /* now divide the moment by the user-specified norm */
  moment /= norm;

  return moment;
}
/*@}*/
