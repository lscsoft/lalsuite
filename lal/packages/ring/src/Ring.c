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
 *   q(t) = \left\{
 *   \begin{array}{ll}
 *     (2\pi)^{1/2}e^{-\pi ft/Q}\cos(2\pi ft) & \mbox{for $t\ge0$} \\
 *     0 & \mbox{for $t<0$}
 *   \end{array}
 *   \right.
 * \end{equation}
 * where the parameters $f$ and $Q$ are specified in the input structure.
 * The output must have an appropriate amount of memory allocated, and must
 * have the desired temporal spacing set.
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
 * The strain waveform produced is $h(t)=Aq(t)$ where the amplitude factor
 * is~\cite{JDECreighton}
 * \begin{equation}
 *   A = 2.415\times10^{-21}Q^{-1/2}[1-0.63(1-{\hat{a}})^{3/10}]^{-1/2}
 *   \left(\frac{\textrm{Mpc}}{r}\right)
 *   \left(\frac{M}{M_\odot}\right)
 *   \left(\frac{\epsilon}{0.01}\right)^{1/2}.
 * \end{equation}
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

/* <lalVerbatim file="RingCP"> */
void
LALComputeRingTemplate(
    LALStatus         *status,
    REAL4TimeSeries   *output,
    RingTemplateInput *input
    )
{ /* </lalVerbatim> */
  const REAL4 efolds = 10;
  REAL4 amp;
  REAL4 fac;
  REAL4 a;
  REAL4 y;
  REAL4 yy;
  UINT4 i;
  UINT4 n;

  INITSTATUS( status, "LALComputeRingTemplate", RINGC );

  ASSERT( input, status, RINGH_ENULL, RINGH_MSGENULL );
  ASSERT( output, status, RINGH_ENULL, RINGH_MSGENULL );
  ASSERT( output->data, status, RINGH_ENULL, RINGH_MSGENULL );

  /* exponential decay variables */
  amp = sqrt( 2 * LAL_PI );
  fac = exp( - LAL_PI * input->frequency * output->deltaT / input->quality );
  n = ceil( - efolds / log( fac ) );

  /* oscillator variables */
  a = 2 * cos( 2 * LAL_PI * input->frequency * output->deltaT );
  y = sin( -2 * LAL_PI * input->frequency * output->deltaT + 0.5 * LAL_PI );
  yy = sin( -4 * LAL_PI * input->frequency * output->deltaT + 0.5 * LAL_PI );

  if ( n < output->data->length )
  {
    memset( output->data->data + n, 0,
        ( output->data->length - n ) * sizeof( *output->data->data ) );
  }
  else
  {
    n = output->data->length;
  }

  for ( i = 0; i < n; ++i )
  {
    REAL4 tmp = a * y - yy;
    yy = y;
    output->data->data[i] = amp * ( y = tmp );
    amp *= fac;
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
  RingTemplateInput tmplt;
  REAL4 ffac;
  REAL4 amp;
  UINT4 i;

  INITSTATUS( status, "LALComputeBlackHoleRing", RINGC );
  ATTATCHSTATUSPTR( status );

  ASSERT( input, status, RINGH_ENULL, RINGH_MSGENULL );
  ASSERT( output, status, RINGH_ENULL, RINGH_MSGENULL );
  ASSERT( output->data, status, RINGH_ENULL, RINGH_MSGENULL );

  ffac = 1 - 0.63 * pow( 1 - input->dimensionlessSpin, 0.3 );
  tmplt.frequency = 32000 * ffac / input->solarMasses;
  tmplt.quality = 2 * pow( 1 - input->dimensionlessSpin, -0.45 );

  amp  = 2.415e-21;
  amp *= sqrt( input->percentMassLoss / ( tmplt.quality * ffac ) );
  amp *= input->solarMasses / input->distanceMpc;

  TRY( LALComputeRingTemplate( status->statusPtr, output, &tmplt ), status );
  for ( i = 0; i < output->data->length; ++i )
  {
    output->data->data[i] *= amp;
  }
  
  DETATCHSTATUSPTR( status );
  RETURN( status );
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
    REAL4 q2   = q * q;
    REAL4 logf = minlogf;

    while ( logf < maxlogf )
    {
      if ( tmplt )
      {
        tmplt[count].quality   = q;
        tmplt[count].frequency = exp( logf );
      }
      ++count;
      logf += dseff / sqrt( 3 + 8 * q2 );
    }

    q += dseff * q * ( 1 + 4 * q2 ) / sqrt( 3 + 16 * q2 * q2 );
  }

  return count;
}


/* <lalVerbatim file="RingCP"> */
void
LALCreateRingTemplateBank(
    LALStatus              *status,
    RingTemplateBank      **output,
    RingTemplateBankInput  *input
    )
{ /* </lalVerbatim> */
  RingTemplateBank *bank;

  INITSTATUS( status, "LALCreateRingTemplateBank", RINGC );

  ASSERT( input, status, RINGH_ENULL, RINGH_MSGENULL );
  ASSERT( output, status, RINGH_ENULL, RINGH_MSGENULL );
  ASSERT( ! *output, status, RINGH_ENNUL, RINGH_MSGENNUL );

  bank = *output = LALCalloc( 1, sizeof( **output ) );
  if ( ! bank )
  {
    ABORT( status, RINGH_EALOC, RINGH_MSGEALOC );
  }
  bank->numTmplt = MakeBank( NULL, input );
  bank->tmplt = LALMalloc( bank->numTmplt * sizeof( *bank->tmplt ) );
  if ( ! bank->tmplt )
  {
    LALFree( *output );
    *output = NULL;
    ABORT( status, RINGH_EALOC, RINGH_MSGEALOC );
  }

  MakeBank( bank->tmplt, input );
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
  LALFree( (*bank)->tmplt );
  LALFree( *bank );
  *bank = NULL;
  RETURN( status );
}
