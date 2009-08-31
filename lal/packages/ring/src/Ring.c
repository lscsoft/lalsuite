/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton
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

/**** <lalVerbatim file="RingCV">
 * Author: Jolien Creighton
 * $Id$
 **** </lalVerbatim> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/Ring.h>

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{Ring.c}}
 *
 * Routines to generate ringdown waveforms and to make a ringdown template bank.
 *
 * \subsubsection*{Prototypes}
 * \input{RingCP}
 * \idx{XLALComputeRingTemplate()}
 * \idx{XLALComputeBlackHoleRing()}
 * \idx{XLALCreateRingTemplateBank()}
 * \idx{XLALDestroyRingTemplateBank()}
 *
 * \subsubsection*{Description}
 *
 * The routine \verb+LALComputeRingTemplate()+ computes the ringdown waveform
 * \begin{equation}
 *   r(t) = \left\{
 *   \begin{array}{ll}
 *     e^{-\pi ft/Q}\cos(2\pi ft) & \mbox{for $t\ge0$} \\
 *     0 & \mbox{for $t<0$}
 *   \end{array}
 *   \right.
 * \end{equation}
 * where the parameters $f$ and $Q$ are specified in the input structure.
 * The output must have an appropriate amount of memory allocated, and must
 * have the desired temporal spacing set.  Note: Ref.~\cite{JDECreighton}
 * used a different convention for the ringdown normlization: there the
 * ringdown waveform was taken to be $q(t)=(2\pi)^{1/2}r(t)$.
 *
 * The routine \verb+LALComputeBlackHoleRing()+ computes a waveform for a
 * black hole with the specified physical parameters (in the input structure).
 * The parameters are the black hole mass $M$ (in solar masses $M_\odot$), the
 * spin $S={\hat{a}}GM^2/c$ expressed in terms of the dimensionless mass
 * parameter ${\hat{a}}$, the fractional mass lost $\epsilon$ in ringdown
 * radiation expressed as a percent, and the distance to the typical source
 * (angle-averaged waveform) $r$ given in megaparsecs (Mpc).  The central
 * frequency and quality of the ringdown are approximated
 * as~\cite{EWLeaver,FEcheverria}:
 * \begin{equation}
 *   f \simeq 32\,\textrm{kHz}\times[1-0.63(1-{\hat{a}})^{3/10}](M_\odot/M)
 * \end{equation}
 * and
 * \begin{equation}
 *   Q \simeq 2(1-{\hat{a}})^{-9/20}.
 * \end{equation}
 * The strain waveform produced is $h(t)=A_q q(t)$ where the amplitude factor
 * is~\cite{JDECreighton}
 * \begin{equation}
 *   A_q = 2.415\times10^{-21}Q^{-1/2}[1-0.63(1-{\hat{a}})^{3/10}]^{-1/2}
 *   \left(\frac{\textrm{Mpc}}{r}\right)
 *   \left(\frac{M}{M_\odot}\right)
 *   \left(\frac{\epsilon}{0.01}\right)^{1/2}.
 * \end{equation}
 * Note that this is written $A_q$ to emphasize that it is the amplitude
 * factor for $q(t)$ rather than $r(t)$.
 *
 * The routine \verb+LALCreateRingTemplateBank()+ creates a bank of ringdown
 * templates that cover a set range in the parameters $f$ and $Q$.  The bank
 * is destroyed with \verb+LALDestroyRingTemplateBank()+.
 *
 * \subsubsection*{Algorithm}
 *
 * The waveform generation routines use recurrance relations for both the
 * exponentially-decaying envelope and for the co-sinusoid.
 *
 * The template placement algorithm is described above.
 *
 * \subsubsection*{Uses}
 *
 * %% List of any external functions called by this function.
 *
 * \subsubsection*{Notes}
 *
 * %% Any relevant notes.
 *
 * \vfill{\footnotesize\input{RingCV}}
 *
 **** </lalLaTeX> */


NRCSID( RINGC, "$Id$" );

static REAL4 ring_spin_factor( REAL4 a )
{
  return 1.0 - 0.63 * pow( 1.0 - a, 3.0/10.0 );
}

static REAL4 ring_quality_fn( REAL4 Q )
{
    return  1.0 + 7.0/(24.0*Q*Q);
}

/*
 *
 * XLAL Routines.
 *
 */


/* <lalVerbatim file="RingCP"> */
REAL4 XLALBlackHoleRingSpin( REAL4 Q )
/* </lalVerbatim> */
{
  return 1.0 - pow( 2.0/Q, 20.0/9.0 );
}

/* <lalVerbatim file="RingCP"> */
REAL4 XLALBlackHoleRingMass( REAL4 f, REAL4 Q )
/* </lalVerbatim> */
{
  const REAL4 c = LAL_C_SI;
  const REAL4 a = XLALBlackHoleRingSpin( Q );
  const REAL4 g = ring_spin_factor( a );
  return (c * c * c * g) / ( LAL_TWOPI * LAL_G_SI * LAL_MSUN_SI * f );
}

/* <lalVerbatim file="RingCP"> */
REAL4 XLALBlackHoleRingQuality( REAL4 a )
/* </lalVerbatim> */
{
  return 2.0 * pow( ( 1.0 - a ), -0.45 );
}

/* <lalVerbatim file="RingCP"> */
REAL4 XLALBlackHoleRingFrequency( REAL4 M, REAL4 a )
/* </lalVerbatim> */
{
  const REAL4 c = LAL_C_SI;
  const REAL4 g = ring_spin_factor( a );
  return (c * c * c * g) / ( LAL_TWOPI * LAL_G_SI * LAL_MSUN_SI * M );
}

/* Formulas for final mass and spin of a non-spinning binary */
/* Buonanno et al arxiv:0706.3732v3 */
/* <lalVerbatim file="RingCP"> */
REAL4 XLALNonSpinBinaryFinalBHSpin( REAL4 eta )
/* </lalVerbatim> */
{
  return sqrt(12.0) * eta - 2.9 * eta *eta;
}

/* <lalVerbatim file="RingCP"> */
REAL4 XLALNonSpinBinaryFinalBHMass( REAL4 eta, REAL4 mass1, REAL4 mass2 )
/* </lalVerbatim> */
{
  return ( 1 + ( sqrt(8.0/9.0) - 1) * eta - 0.498 * eta * eta) * (mass1 + mass2);
}


/* <lalVerbatim file="RingCP"> */
REAL4 XLALBlackHoleRingAmplitude( REAL4 f, REAL4 Q, REAL4 r, REAL4 epsilon )
/* </lalVerbatim> */
{
  const REAL4 c = LAL_C_SI;
  const REAL4 M = XLALBlackHoleRingMass( f, Q );
  const REAL4 a = XLALBlackHoleRingSpin( Q );
  const REAL4 g = ring_spin_factor( a );
  const REAL4 F = ring_quality_fn( Q );

  return sqrt(5.0/2.0 * epsilon) *
    ( (LAL_G_SI * M * LAL_MSUN_SI) / ( c * c * r * 1.0e6 * LAL_PC_SI) ) *
    (1.0 / sqrt( Q * F * g) );
}

/* <lalVerbatim file="RingCP"> */
REAL4 XLALBlackHoleRingEpsilon( REAL4 f, REAL4 Q, REAL4 r, REAL4 amplitude )
  /* </lalVerbatim> */
{
  const REAL4 c = LAL_C_SI;
  const REAL4 M = XLALBlackHoleRingMass( f, Q );
  const REAL4 a = XLALBlackHoleRingSpin( Q );
  const REAL4 g = ring_spin_factor( a );
  const REAL4 F = ring_quality_fn( Q );

  return (2.0/5.0 * amplitude * amplitude) *
    pow( ( c * c * r * 1.0e6 * LAL_PC_SI) / (LAL_G_SI * M * LAL_MSUN_SI), 2 ) *
    ( Q * F * g );
}

/* <lalVerbatim file="RingCP"> */
REAL8 XLAL2DRingMetricDistance( REAL8 fa, REAL8 fb, REAL8 Qa, REAL8 Qb )
/* </lalVerbatim> */
{
  REAL8 Q2 = Qa*Qa;
  REAL8 gQQ;
  REAL8 gff;
  REAL8 gQf;

  gQQ = ( 3.0 + 16.0 * Q2 * Q2) / ( Q2 * ( 1.0 + 4.0 * Q2 ) * ( 1.0 + 4.0 * Q2 ) );
  gff = ( 3.0 + 8.0 * Q2) / ( fa * fa);
  gQf = - 2.0 * ( 3.0 + 4.0 * Q2 ) / ( Qa * fa * ( 1.0 + 4.0 * Q2 ));

  return ( 1.0/8.0 * ( gQQ * pow(Qb-Qa,2) + gQf * (Qb-Qa) * (fb-fa) + gff * pow(fb-fa,2) ) );
}

/* <lalVerbatim file="RingCP"> */
REAL8 XLAL3DRingMetricDistance( REAL8 fa, REAL8 fb, REAL8 Qa, REAL8 Qb, REAL8 dt )
/* </lalVerbatim> */
{
  REAL8 gQQ, gff, gtt;
  REAL8 gQf, gtf, gtQ;
  REAL8 df, dQ, ds2;
  REAL8 f = (fa+fb)/2.;
  REAL8 Q = (Qa+Qb)/2.;
  REAL8 Q2 = Q*Q;

  gQQ = ( 1. + 28.*Q2*Q2 + 128.*Q2*Q2*Q2 + 64.*Q2*Q2*Q2*Q2) / ( 4. * Q2 * ( 1. + 6.*Q2 + 8.*Q2*Q2 )*( 1. + 6.*Q2 + 8.*Q2*Q2 ) );
  gff = ( 1. + 6.*Q2 + 16.*Q2*Q2) / ( 4. * f*f * ( 1. + 2.*Q2 ) );
  gtt = ( LAL_PI*LAL_PI * f*f ) * ( 1. + 4.*Q2 ) / ( Q2 );
  gQf = - ( 1. + 2.*Q2 + 8.*Q2*Q2 ) / ( 4.*Q*f * ( 1. + 6.*Q2 + 8.*Q2*Q2 ) );
  gtf = - ( LAL_PI * Q ) * ( 1. + 4.*Q2) / ( 1. + 2.*Q2 );
  gtQ = ( LAL_PI * f ) * ( 1. - 2.*Q2 ) / ( ( 1. + 2.*Q2 )*( 1. + 2.*Q2 ) );

  df = fb - fa;
  dQ = Qb - Qa;

  ds2 = ( gQQ * dQ*dQ + gff * df*df + gtt * dt*dt + gQf * 2.*dQ*df + gtf * 2.*dt*df + gtQ * 2.*dt*dQ );

  return ( ds2 );
}

REAL8 XLALRingdownTimeError( const SnglRingdownTable *table,  REAL8 lal_ring_ds_sq )
{
  REAL8 gtt;
  REAL8 f = table->frequency;
  REAL8 Q = table->quality;
  REAL8 Q2 = Q*Q;

  gtt = ( LAL_PI*LAL_PI * f*f ) * ( 1. + 4.*Q2 ) / ( Q2 );

  return ( sqrt( lal_ring_ds_sq / gtt ) );
}

/* <lalVerbatim file="RingCP"> */
int XLALComputeRingTemplate( REAL4TimeSeries *output, SnglRingdownTable *input )
/* </lalVerbatim> */
{
  static const char *func = "XLALComputeRingTemplate";
  const REAL8 efolds = 10;
  REAL8 amp = 1.0;
  REAL8 fac;
  REAL8 a;
  REAL8 y;
  REAL8 yy;
  UINT4 i;
  UINT4 n;

  if ( ! output || ! output->data || ! input )
    XLAL_ERROR( func, XLAL_EFAULT );

  if ( ! output->data->length )
    XLAL_ERROR( func, XLAL_EBADLEN );

  if ( output->deltaT <= 0 || input->quality <= 0 || input->frequency <= 0 )
    XLAL_ERROR( func, XLAL_EINVAL );

  /* exponential decay variables */
  fac = exp( - LAL_PI * input->frequency * output->deltaT / input->quality );
  n = ceil( - efolds / log( fac ) );

  /* oscillator variables */
  a = 2 * cos( 2 * LAL_PI * input->frequency * output->deltaT );
  y = sin( -2 * LAL_PI * input->frequency * output->deltaT +
      0.5 * LAL_PI + input->phase );
  yy = sin( -4 * LAL_PI * input->frequency * output->deltaT +
      0.5 * LAL_PI + input->phase );

  if ( n < output->data->length )
    memset( output->data->data + n, 0,
        ( output->data->length - n ) * sizeof( *output->data->data ) );
  else
    n = output->data->length;

  for ( i = 0; i < n; ++i )
  {
    REAL4 tmp = a * y - yy;
    yy = y;
    output->data->data[i] = amp * ( y = tmp );
    amp *= fac;
  }

  return 0;
}


/* <lalVerbatim file="RingCP"> */
int XLALComputeBlackHoleRing(
    REAL4TimeSeries     *output,
    SnglRingdownTable   *input,
    REAL4                dynRange
    )
/* </lalVerbatim> */
{
  static const char *func = "XLALComputeBlackHoleRing";
  const REAL4 amp = dynRange *
    XLALBlackHoleRingAmplitude(
        input->frequency, input->quality, input->eff_dist, input->epsilon );
  UINT4 i;

  if ( XLALComputeRingTemplate( output, input ) < 0 )
    XLAL_ERROR( func, XLAL_EFUNC );

  for ( i = 0; i < output->data->length; ++i )
    output->data->data[i] *= amp;

  return 0;
}


static int MakeBank( SnglRingdownTable *tmplt, RingTemplateBankInput *input )
{
  UINT4 count = 0;
  REAL4 dseff = 4 * sqrt( input->maxMismatch );
  REAL4 minlogf = log( input->minFrequency );
  REAL4 maxlogf = log( input->maxFrequency );
  REAL4 q = input->minQuality;

  while ( q < input->maxQuality )
  {
    REAL4 q2 = q * q;
    REAL4 logfreq = minlogf;

    while ( logfreq < maxlogf )
    {
      if ( tmplt )
      {
        tmplt[count].quality   = q;
        tmplt[count].frequency = exp( logfreq );
        tmplt[count].phase     = input->templatePhase;
        tmplt[count].epsilon   = input->templateEpsilon;
        tmplt[count].eff_dist  = input->templateDistance;
      }
      ++count;
      logfreq += dseff / sqrt( 3 + 8 * q2 );
    }

    q += dseff * q * ( 1 + 4 * q2 ) / sqrt( 3 + 16 * q2 * q2 );
  }

  return count;
}


/* <lalVerbatim file="RingCP"> */
RingTemplateBank *XLALCreateRingTemplateBank( RingTemplateBankInput *input )
/* </lalVerbatim> */
{
  static const char *func = "XLALCreateRingTemplateBank";
  UINT4 i;
  RingTemplateBank *bank;

  if ( ! input )
    XLAL_ERROR_NULL( func, XLAL_EFAULT );

  bank = LALCalloc( 1, sizeof( *bank ) );
  if ( ! bank )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );

  bank->numTmplt = MakeBank( NULL, input );
  bank->tmplt = LALCalloc( bank->numTmplt, sizeof( *bank->tmplt ) );
  if ( ! bank->tmplt )
  {
    LALFree( bank );
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  }
  for ( i = 0; i < bank->numTmplt - 1; ++i )
    bank->tmplt[i].next = bank->tmplt + i + 1;

  MakeBank( bank->tmplt, input );
  return bank;
}


/* <lalVerbatim file="RingCP"> */
void XLALDestroyRingTemplateBank( RingTemplateBank *bank )
/* </lalVerbatim> */
{
  if ( bank )
  {
    if ( bank->tmplt )
      LALFree( bank->tmplt );
    LALFree( bank );
  }
  return;
}
