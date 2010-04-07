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

/*  <lalVerbatim file="LALInspiralMomentsCV">
Authors: Brown, D. A., Cokelaer, T. and  Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

#if 0
<lalLaTeX>

\subsection{Module \texttt{LALInspiralMoments.c}}

Module to calculate the moment of the noise power spectral density.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralMomentsCP}
\idx{LALInspiralMoments()}
\begin{itemize}
   \item \texttt{moment,} Output, the value of the moment
   \item \texttt{pars,} Input
\end{itemize}

\subsubsection*{Description}

The moments of the noise curve are defined as
\begin{equation}
I(q)  \equiv S_{h}(f_{0}) \int^{f_{c}/f_{0}}_{f_{s}/f_{0}}
\frac{x^{-q}}{S_{h}(x)} \, dx \,.
\end{equation}
Because in practice we will always divide one of these moments by another, we
do not need to include the $S_{h}(f_{0})$ term, which always cancels.
This function calculates the integral
\begin{equation}
I = \int^{f_{c}/f_{0}}_{f_{s}/f_{0}} \frac{x^{-q}}{S_{h}(x)} \, dx \,.
\end{equation}
It then divides this quantity by a normalisation constant which has been
passed to the function. In the case of calculating the components of the
metric for the signal manifold for the purpose of generating a template bank,
this constant is given by $I(7)$, because of the definition of the quantity
\begin{equation}
J(q) \equiv \frac{I(q)}{I(7/3)} \,.
\end{equation}

\subsubsection*{Algorithm}
Given the exponent \texttt{pars.ndx} and limits of integration
\texttt{pars.xmin} and \texttt{pars.xmax} this function returns the moment of
the power spectral density specified by the frequency series
\texttt{pars.shf} according to
\begin{equation}
\mathtt{moment} = \int_{\mathtt{xmin}}^{\mathtt{xmax}}
\frac{x^{-\mathtt{ndx}}}{S_h(x)}\, dx \, .
\end{equation}

\subsubsection*{Uses}
\begin{verbatim}
LALDRombergIntegrate
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralMomentsCV}}

</lalLaTeX>
#endif

#include <lal/LALInspiralBank.h>
#include <lal/Integrate.h>

NRCSID( LALINSPIRALMOMENTSC, "$Id$" );

/* <lalVerbatim file="LALInspiralMomentsCP"> */
void
LALGetInspiralMoments (
    LALStatus            *status,
    InspiralMomentsEtc   *moments,
    REAL8FrequencySeries *psd,
    InspiralTemplate     *params
    )
/* </lalVerbatim> */
{
  UINT4 k;
  InspiralMomentsIn in;

  INITSTATUS( status, "LALGetInspiralMoments", LALINSPIRALMOMENTSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( params, status,
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( params->fLower>0, status,
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( moments, status,
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( psd, status,
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );

  /* Constants needed in computing the moments */
  moments->a01 = 3.L/5.L;
  moments->a21 = 11.L * LAL_PI/12.L;
  moments->a22 = 743.L/2016.L * pow(25.L/(2.L*LAL_PI*LAL_PI),1.L/3.L);
  moments->a31 = -3.L/2.L;
  moments->a41 = 617.L * LAL_PI * LAL_PI / 384.L;
  moments->a42 = 5429.L/5376.L * pow(25.L*LAL_PI/2.L,1.L/3.L);
  moments->a43 = 1.5293365L/1.0838016L * pow(5.L/(4.L*pow(LAL_PI,4.L)),1.L/3.L);

  /* setup the input structure needed in the computation of the moments */
  in.shf = psd;

  /* Divide all frequencies by fLower, a scaling that is used in solving */
  /* the moments integral                                                */
  in.shf->f0 /= params->fLower;
  in.shf->deltaF /= params->fLower;
  in.xmin = params->fLower / params->fLower;
  in.xmax = params->fCutoff / params->fLower;

  /* First compute the norm and print if requested */
  in.norm = 1.L;
  in.ndx = 7.L/3.L;
  LALInspiralMoments( status->statusPtr, &moments->j[7], in );
  CHECKSTATUSPTR( status );
  in.norm = moments->j[7];

  if ( lalDebugLevel & LALINFO )
  {
    LALPrintError(
        "a01=%e a21=%e a22=%e a31=%e a41=%e a42=%e a43=%e j7=%e\n",
        moments->a01, moments->a21, moments->a22, moments->a31,
        moments->a41, moments->a42, moments->a43, moments->j[7] );
  }

  /* Then compute the normalised moments of the noise PSD from 1/3 to 17/3. */
  for ( k = 1; k <= 17; ++k )
  {
    in.ndx = (REAL8) k / 3.L;
    LALInspiralMoments( status->statusPtr, &moments->j[k], in );
    CHECKSTATUSPTR( status );

    if ( lalDebugLevel & LALINFO )
    {
      LALPrintError( "j%1i=%e\n", k, moments->j[k] );
    }
  }

  /* Moments are done: Rescale deltaF and f0 back to their original values */
  in.shf->deltaF *= params->fLower;
  in.shf->f0 *= params->fLower;

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


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
  double q;

  INITSTATUS( status, "LALGetInspiralMomentsBCV", LALINSPIRALMOMENTSC );
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
  q = 2* moments->n0; /*=10/3 */
  moments->M1[0][0] = (moments->i[17] -2.*moments->alpha * moments->i[15] +
      moments->alpha * moments->alpha*moments->i[13]) / in.norm;
  /* 14 */
  q = moments->n0 +moments->n15;
  moments->M1[0][1] = (moments->i[14] -2.*moments->alpha * moments->i[12] +
      moments->alpha * moments->alpha*moments->i[10]) / in.norm;
  /* 11 */
  q = 2 * moments->n15;
  moments->M1[1][1] = (moments->i[11] -2.*moments->alpha * moments->i[9] +
      moments->alpha * moments->alpha*moments->i[7]) / in.norm;

  moments->M1[1][0]=moments->M1[0][1] ;

  /*  12 */
  q = moments->n0;
  moments->M2[0][0] = (moments->i[12] -2.*moments->alpha * moments->i[10] +
      moments->alpha * moments->alpha*moments->i[8]) / in.norm;
  /* 9 */
  q = moments->n15;

  moments->M2[0][1] = (moments->i[9] -2.*moments->alpha * moments->i[7] +
      moments->alpha * moments->alpha*moments->i[5]) / in.norm;
  /*  9 */
  q = moments->n0-1;

  moments->M2[1][0] = (moments->i[9] -2.*moments->alpha * moments->i[7] +
      moments->alpha * moments->alpha*moments->i[5]) / in.norm;
  /*  6 */
  q = moments->n15-1;
  moments->M2[1][1] = (moments->i[6] -2.*moments->alpha * moments->i[4] +
      moments->alpha * moments->alpha*moments->i[2]) / in.norm;

  /* 7 */
  q = 0;
  moments->M3[0][0] = (moments->i[7] -2.*moments->alpha * moments->i[5] +
      moments->alpha * moments->alpha*moments->i[3]) / in.norm;
  /* 4 */
  q = -1;
  moments->M3[0][1] = (moments->i[4] -2.*moments->alpha * moments->i[2] +
      moments->alpha * moments->alpha*moments->i[0]) / in.norm;
  /* 1 */
  q = -2;
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


/* <lalVerbatim file="LALInspiralMomentsCP"> */
void
LALInspiralMoments(
    LALStatus         *status,
    REAL8             *moment,
    InspiralMomentsIn  pars
    )
/* </lalVerbatim> */
{
  REAL8 f;
  REAL8 momentTmp;
  REAL8 fMin;
  REAL8 fMax;
  REAL8 deltaF;
  UINT4 k;
  UINT4 kMin;
  UINT4 kMax;

  INITSTATUS( status, "LALInspiralMoments", LALINSPIRALMOMENTSC );

  ASSERT( pars.shf, status,
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( pars.shf->data, status,
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( pars.shf->data->data, status,
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );

  /* make sure that the minimum and maximum of the integral are within */
  /* the frequency series                                              */
  fMax = pars.shf->f0 + (REAL8) pars.shf->data->length * pars.shf->deltaF;
  if ( pars.xmin < pars.shf->f0 || pars.xmax > fMax+LAL_REAL4_EPS )
  {
    ABORT( status, LALINSPIRALBANKH_EFRANGE, LALINSPIRALBANKH_MSGEFRANGE );
  }

  /* the minimum and maximum frequency where we have four points */
  deltaF = pars.shf->deltaF;
  fMin = pars.shf->f0 + deltaF;
  fMax = pars.shf->f0 + ((REAL8) pars.shf->data->length - 2 ) * deltaF;

  if ( pars.xmin <= fMin )
  {
    kMin = 1;
  }
  else
  {
    kMin = (UINT8) floor( (pars.xmin - pars.shf->f0) / deltaF );
  }

  if ( pars.xmax >= fMax )
  {
    kMax = pars.shf->data->length - 1;
  }
  else
  {
    kMax = (UINT8) floor( (pars.xmax - pars.shf->f0) / deltaF );
  }

  /* the first and last points of the integral */
  momentTmp = 0.;
  f = pars.shf->f0 + (REAL8) kMin * deltaF;
  if( pars.shf->data->data[kMin] )
  {
    momentTmp = pow( f, -(pars.ndx) ) / ( 2.0 * pars.shf->data->data[kMin] );
  }
  *moment = momentTmp;

  momentTmp = 0.;
  f = pars.shf->f0 + (REAL8) kMax * deltaF;
  if( pars.shf->data->data[kMin] )
  {
    momentTmp = pow( f, -(pars.ndx) ) / ( 2.0 * pars.shf->data->data[kMin] );
  }
  *moment += momentTmp;
#if 0
  In the following line we should have kMax
  Changed by Sathya on June 30, 2002
  *moment += pow( f, -(pars.ndx) ) / ( 2.0 * pars.shf->data->data[kMin] );
#endif
  momentTmp = 0.;
  if ( pars.shf->data->data[kMax] )
  {
    momentTmp = pow( f, -(pars.ndx) ) / ( 2.0 * pars.shf->data->data[kMax] );
  }
  *moment += momentTmp;
  kMin++;
  kMax--;

  for ( k = kMin; k < kMax; ++k )
  {
    momentTmp = 0.;
    f = pars.shf->f0 + (REAL8) k * deltaF;
    if ( pars.shf->data->data[k] )
    {
      momentTmp = pow( f, -(pars.ndx) ) / pars.shf->data->data[k];
    }
    *moment += momentTmp;
  }

  *moment *= deltaF;

  /* now divide the moment by the specified norm */
  *moment /= pars.norm;

  RETURN (status);
}
