/************************** <lalVerbatim file="GenerateBurstCV">
Author: Brady, P. B.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\subsection{Module \texttt{GenerateBurst.c}}
\label{ss:GenerateBurst.c}

Computes one of the standard burst waveforms with specified $h_{rss}$.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{GenerateBurstCP}
\idx{LALGenerateBurst()}

\subsubsection*{Description}

This function computes one of the following burst waveforms:
\begin{description}
\item[Sine-Gaussian]:  a linearly polarized sine-Gaussian with the specified
frequency and decay constant.
\end{description}

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()                   LALFree()
LALSCreateVectorSequence()    LALSDestroyVectorSequence()
LALSCreateVector()            LALSDestroyVector()
LALDCreateVector()            LALDDestroyVector()
LALSnprintf()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GenerateBurstCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GenerateBurst.h>

NRCSID( GENERATEBURSTC, "$Id$" );

/* <lalVerbatim file="GenerateBurstCP"> */
void
LALGenerateBurst( 
    LALStatus          *stat, 
    CoherentGW         *output, 
    BurstParamStruc    *params 
    )
{ /* </lalVerbatim> */
  UINT4 n, i;          /* number of and index over samples */
  REAL8 t, dt;         /* time, interval */
  REAL8 t0, tau, gtime;  /* central time, decay time, gaussian time */
  REAL8 f0, phi0;      /* initial phase and frequency */
  REAL8 twopif0;       /* 2*pi*f0 */
  REAL8 f;             /* current value of frequency */
  REAL4 hrss;          /* root sum square strain for burst */
  REAL4 df = 0.0;      /* maximum difference between f */
  REAL8 phi;           /* current value of phase */
  REAL4 *fData;        /* pointer to frequency data */
  REAL8 *phiData;      /* pointer to phase data */
  REAL4 *aData;        /* pointer to frequency data */

  INITSTATUS( stat, "LALGenerateBurst", GENERATEBURSTC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter and output structures exist. */
  ASSERT( params, stat, GENERATEBURSTH_ENUL,
	  GENERATEBURSTH_MSGENUL );
  ASSERT( output, stat, GENERATEBURSTH_ENUL,
	  GENERATEBURSTH_MSGENUL );

  /* Make sure output fields don't exist. */
  ASSERT( !( output->a ), stat, GENERATEBURSTH_EOUT,
	  GENERATEBURSTH_MSGEOUT );
  ASSERT( !( output->f ), stat, GENERATEBURSTH_EOUT,
	  GENERATEBURSTH_MSGEOUT );
  ASSERT( !( output->phi ), stat, GENERATEBURSTH_EOUT,
	  GENERATEBURSTH_MSGEOUT );
  ASSERT( !( output->shift ), stat, GENERATEBURSTH_EOUT,
	  GENERATEBURSTH_MSGEOUT );

  /* Set up some other constants, to avoid repeated dereferencing. */
  n = params->length;
  dt = params->deltaT;
  f0 = params->f0;
  twopif0 = f0*LAL_TWOPI;

  /* Allocate output structures. */
  if ( ( output->a = (REAL4TimeVectorSeries *)
	 LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
    ABORT( stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM );
  }
  memset( output->a, 0, sizeof(REAL4TimeVectorSeries) );
  if ( ( output->f = (REAL4TimeSeries *)
	 LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
    LALFree( output->a ); output->a = NULL;
    ABORT( stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM );
  }
  memset( output->f, 0, sizeof(REAL4TimeSeries) );
  if ( ( output->phi = (REAL8TimeSeries *)
	 LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL ) {
    LALFree( output->a ); output->a = NULL;
    LALFree( output->f ); output->f = NULL;
    ABORT( stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM );
  }
  memset( output->phi, 0, sizeof(REAL8TimeSeries) );

  /* Set output structure metadata fields. */
  output->position = params->position;
  output->psi = params->psi;
  output->a->epoch = output->f->epoch = output->phi->epoch
    = params->epoch;
  output->a->deltaT = params->deltaT;
  output->f->deltaT = output->phi->deltaT = params->deltaT;
  output->a->sampleUnits = lalStrainUnit;
  output->f->sampleUnits = lalHertzUnit;
  output->phi->sampleUnits = lalDimensionlessUnit;
  LALSnprintf( output->a->name, LALNameLength, "Burst amplitudes" );
  LALSnprintf( output->a->name, LALNameLength, "Burst frequency" );
  LALSnprintf( output->a->name, LALNameLength, "Burst phase" );

  /* Allocate phase and frequency arrays. */
  LALSCreateVector( stat->statusPtr, &( output->f->data ), n );
  BEGINFAIL( stat ) {
    LALFree( output->a );   output->a = NULL;
    LALFree( output->f );   output->f = NULL;
    LALFree( output->phi ); output->phi = NULL;
  } ENDFAIL( stat );
  LALDCreateVector( stat->statusPtr, &( output->phi->data ), n );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr, &( output->f->data ) ),
	 stat );
    LALFree( output->a );   output->a = NULL;
    LALFree( output->f );   output->f = NULL;
    LALFree( output->phi ); output->phi = NULL;
  } ENDFAIL( stat );

  /* Allocate amplitude array. */
  {
    CreateVectorSequenceIn in; /* input to create output->a */
    in.length = 2;
    in.vectorLength = n;
    LALSCreateVectorSequence( stat->statusPtr, &(output->a->data), &in );
    BEGINFAIL( stat ) {
      TRY( LALSDestroyVector( stat->statusPtr, &( output->f->data ) ),
	   stat );
      TRY( LALDDestroyVector( stat->statusPtr, &( output->phi->data ) ),
	   stat );
      LALFree( output->a );   output->a = NULL;
      LALFree( output->f );   output->f = NULL;
      LALFree( output->phi ); output->phi = NULL;
    } ENDFAIL( stat );
  }

  /* Fill frequency and phase arrays. */
  fData = output->f->data->data;
  phiData = output->phi->data->data;
  aData = output->a->data->data; 

  /* this depends on the waveform type */
  switch ( params->burstType )
  {
    /* sine-Gaussian burst */
    case sineGaussian:
      for ( i = 0; i < n; i++ ) {
        t = i*dt;
        gtime = (t-t0)/tau;
        *(fData++) = f0;
        *(phiData++) = twopif0 * t;
        *(aData++) = hrss * exp( - gtime * gtime );
        *(aData++) = 0.0;
      }
    /* Gaussian burst */
    case Gaussian:
      for ( i = 0; i < n; i++ ) {
        t = i*dt;
        gtime = (t-t0)/tau;
        *(fData++) = 0.0;
        *(phiData++) = 0.0;
        *(aData++) = hrss * exp( - gtime * gtime );
        *(aData++) = 0.0;
      }
  }

  /* Set output field and return. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
