/**** <lalVerbatim file="RingCV">
 * Author: Jolien Creighton
 * $Id$
 **** </lalVerbatim> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Ring.h>

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{Ring.c}}
 *
 * Routines to generate ringdown waveforms and to make a ringdown template bank.
 *
 * \subsubsection*{Prototypes}
 * \input{RingCP}
 * \idx{LALComputeRingTemplate()}
 * \idx{LALComputeBlackHoleRing()}
 * \idx{LALCreateRingTemplateBank()}
 * \idx{LALDestroyRingTemplateBank()}
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


/*
 *
 * XLAL Routines.
 *
 */


int XLALComputeRingTemplate( REAL4TimeSeries *output, RingTemplateInput *input )
{
  static const char *func = "XLALComputeRingTemplate";
  const REAL8 efolds = 10;
  REAL8 amp;
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
  /* amp = sqrt( 2 * LAL_PI ); */ /* OLD conventions of PRD 022001 (1999) */
  amp = 1; /* NEW conventions */
  fac = exp( - LAL_PI * input->frequency * output->deltaT / input->quality );
  n = ceil( - efolds / log( fac ) );

  /* oscillator variables */
  a = 2 * cos( 2 * LAL_PI * input->frequency * output->deltaT );
  y = sin( -2 * LAL_PI * input->frequency * output->deltaT + 0.5 * LAL_PI + input->phase );
  yy = sin( -4 * LAL_PI * input->frequency * output->deltaT + 0.5 * LAL_PI + input->phase );

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


int XLALComputeBlackHoleRing( REAL4TimeSeries *output, BlackHoleRingInput *input )
{
  static const char *func = "XLALComputeBlackHoleRing";
  RingTemplateInput tmplt;
  REAL4 ffac;
  REAL4 amp;
  UINT4 i;

  if ( ! input )
    XLAL_ERROR( func, XLAL_EFAULT );

  ffac = 1 - 0.63 * pow( 1 - input->dimensionlessSpin, 0.3 );
  tmplt.frequency = 32000 * ffac / input->solarMasses;
  tmplt.quality = 2 * pow( 1 - input->dimensionlessSpin, -0.45 );
  tmplt.phase = input->initialPhase;

  amp  = 2.415e-21; /* factor given in PRD 022001 (1999) */
  amp *= sqrt( 2 * LAL_PI ); /* convert NEW conventions to OLD conventions */
  amp *= sqrt( input->percentMassLoss / ( tmplt.quality * ffac ) );
  amp *= input->solarMasses / input->distanceMpc;

  if ( XLALComputeRingTemplate( output, &tmplt ) < 0 )
    XLAL_ERROR( func, XLAL_EFUNC );

  for ( i = 0; i < output->data->length; ++i )
    output->data->data[i] *= amp;
  
  return 0;
}


static int MakeBank( RingTemplateInput *tmplt, RingTemplateBankInput *input )
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
      }
      ++count;
      logfreq += dseff / sqrt( 3 + 8 * q2 );
    }

    q += dseff * q * ( 1 + 4 * q2 ) / sqrt( 3 + 16 * q2 * q2 );
  }

  return count;
}


RingTemplateBank *XLALCreateRingTemplateBank( RingTemplateBankInput *input )
{
  static const char *func = "XLALCreateRingTemplateBank";
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

  MakeBank( bank->tmplt, input );
  return bank;
}


void XLALDestroyRingTemplateBank( RingTemplateBank *bank )
{
  if ( bank )
  {
    if ( bank->tmplt )
      LALFree( bank->tmplt );
    LALFree( bank );
  }
  return;
}


/*
 *
 * LAL Routines.
 *
 */


/* <lalVerbatim file="RingCP"> */
void
LALComputeRingTemplate(
    LALStatus         *status,
    REAL4TimeSeries   *output,
    RingTemplateInput *input
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALComputeRingTemplate", RINGC );

  ASSERT( input, status, RINGH_ENULL, RINGH_MSGENULL );
  ASSERT( output, status, RINGH_ENULL, RINGH_MSGENULL );
  ASSERT( output->data, status, RINGH_ENULL, RINGH_MSGENULL );

  if ( XLALComputeRingTemplate( output, input ) < 0 )
  {
    int errnum = xlalErrno;
    XLALClearErrno();
    switch ( errnum )
    {
      case XLAL_EFAULT:
        ABORT( status, RINGH_ENULL, RINGH_MSGENULL );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="RingCP"> */
void
LALComputeBlackHoleRing(
    LALStatus          *status,
    REAL4TimeSeries    *output,
    BlackHoleRingInput *input
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALComputeBlackHoleRing", RINGC );

  ASSERT( input, status, RINGH_ENULL, RINGH_MSGENULL );
  ASSERT( output, status, RINGH_ENULL, RINGH_MSGENULL );
  ASSERT( output->data, status, RINGH_ENULL, RINGH_MSGENULL );

  if ( XLALComputeBlackHoleRing( output, input ) < 0 )
  {
    int errnum = XLALGetBaseErrno();
    XLALClearErrno();
    switch ( errnum )
    {
      case XLAL_EFAULT:
        ABORT( status, RINGH_ENULL, RINGH_MSGENULL );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="RingCP"> */
void
LALCreateRingTemplateBank(
    LALStatus              *status,
    RingTemplateBank      **output,
    RingTemplateBankInput  *input
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALCreateRingTemplateBank", RINGC );

  ASSERT( input, status, RINGH_ENULL, RINGH_MSGENULL );
  ASSERT( output, status, RINGH_ENULL, RINGH_MSGENULL );
  ASSERT( ! *output, status, RINGH_ENNUL, RINGH_MSGENNUL );

  *output = XLALCreateRingTemplateBank( input );
  if ( ! *output )
  {
    int errnum = xlalErrno;
    XLALClearErrno();
    switch ( errnum )
    {
      case XLAL_EFAULT:
        ABORT( status, RINGH_ENNUL, RINGH_MSGENNUL );
      case XLAL_ENOMEM:
        ABORT( status, RINGH_EALOC, RINGH_MSGEALOC );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="RingCP"> */
void
LALDestroyRingTemplateBank(
    LALStatus         *status,
    RingTemplateBank **bank
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALDestroyRingTemplateBank", RINGC );
  ASSERT( bank, status, RINGH_ENULL, RINGH_MSGENULL );
  ASSERT( *bank, status, RINGH_ENULL, RINGH_MSGENULL );
  XLALDestroyRingTemplateBank( *bank );
  *bank = NULL;
  RETURN( status );
}
