/*************************** <lalVerbatim file="SimulateCoherentGWCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\subsection{Module \texttt{SimulateCoherentGW.c}}
\label{ss:SimulateCoherentGW.c}

Computes the response of a detector to a coherent gravitational wave.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SimulateCoherentGWCP}
\index{\texttt{LALSimulateCoherentGW()}}

\subsubsection*{Description}

This function takes an adiabatically-varying gravitational waveform
given in \verb@*signal@, and estimates the corresponding response of
the detector whose polarization response and transfer function are
specified in \verb@*detector@.  The result is stored in
\verb@*output@.

The fields \verb@output->epoch@, \verb@output->deltaT@, and
\verb@output->data@ must already be set, in order to specify the time
period and sampling rate for which the response is required.  For the
input signal, \verb@signal->h@ is ignored, and the signal is treated
as zero at any time for which either \verb@signal->a@ or
\verb@signal->phi@ is not defined.

\subsubsection*{Algorithm}

The routine first accounts for the time delay between the detector and
the centre of the Earth, based on the detector position information
stored in \verb@*detector@ and the propagation direction specified in
\verb@*signal@.  The propagation delay from the detector to the
Earth's centre is computed for every minute of the desired output, and
stored in an array; subsequently, the value of the delay for each
output sample is determined by interpolating the two nearest points
from the table.  This guarantees accurate timing information to within
$$
\Delta t \lessim \frac{\pi^2 R_\mathrm{Earth}}{2c}
	\left(\frac{\mathrm{1~minute}}{\mathrm{1~day}}\right)^2
	\sim 50\mathrm{ns} \; ,
$$
which is much less than the sampling interval of LIGO.

Next, the polarization response functions of the detector
$F_{+,\times}(\alpha,\delta)$ are computed for every minute of the
signal's duration, using the position of the source in \verb@*signal@,
the detector information in \verb@*detector@, and the response
functions in \verb@???@.  Subsequently, the polarization functions are
estimated for each output sample by interpolating these precomputed
values.  Again, this guarantees a fractional accuracy of
$\lessim(\pi^2/2)(1\mathrm{min}/1\mathrm{day})^2\sim2\times10^{-6}$.

Next, the frequency response of the detector is estimated in the
adiabatic limit as follows:
\begin{itemize}
\item At each sample point in \verb@*output@, the propagation delay is
computed and added to the sample time, and the instantaneous amplitude
$A_+$ and phase $\phi_+$ are found by interpolating the nearest values
in \verb@signal->a@ and \verb@signal->phi@, respectively.  The
instantaneous frequency is estimated by differencing the two nearest
values of $\phi_+$, as $f\approx\Delta\phi_+/2\pi\Delta t$.
\item The complex transfer function of the detector at this frequency
is found by interpolating \verb@detector->transfer@.  The amplitude of
the transfer function is multiplied with $A_+$, the phase of the
transfer function is added to $\phi_+$, and the plus-contribution
$o_+$ to the detector output is computed as $o_+=A_+\cos\phi_+$.
\item The response to the cross-mode polarization is computed
similarly, with the exception that if \verb@signal->a@ contains only
one polarization amplitude then it is assumed that $A_\times=A_+$, and
if \verb@signal->phi@ contains only one polarization phase then it is
assumed that $\phi_\times=\phi_+-\pi/2$.
\item The final detector response $o$ is computed as
$o=(o_+F_+)+(o_\times F_\times)$.
\end{itemize}

\paragraph{A note on interpolation:}
Much of the computational work in this routine involves interpolating
various time series to find their values at specific output times.
The algorithm is summarized below.

Let $A_j = A( t_A + j\Delta t_A )$ be a sampled time series, which we
want to resample at new (output) time intervals $t_k = t_0 + k\Delta
t$.  We first precompute the following quantities:
\begin{eqnarray}
t_\mathrm{off} & = & \frac{t_0-t_A}{\Delta t_A}  \; , \nonumber \\
            dt & = & \frac{\Delta t}{\Delta t_A} \; . \nonumber
\end{eqnarray}
Then, for each output sample time $t_k$, we compute:
\begin{eqnarray}
t & = & t_\mathrm{off} + k \times dt \; , \nonumber \\
j & = & t \bmod1                     \; , \nonumber \\
f & = & t - j                        \; . \nonumber
\end{eqnarray}
The time series sampled at the new time is then:
$$
A(t_k) = f \times A_{j+1} + (1-f) \times A_j \; .
$$

\subsubsection*{Uses}
\begin{verbatim}
LALSCreateVector()      LALSCreateVectorSequence()
LALSDestroyVector()     LALSDestroyVectorSequence()
LALGPStoU()
LALUTime()
LALLMST()
LALCVectorAbs()
LALCVectorAngle()
LALUnwrapREAL4Angle()
LALWarning()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{SimulateCoherentGWCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/VectorOps.h>
#include <lal/SimulateCoherentGW.h>

#define POLDT   60 /* (integral) number of seconds between table
		      entries for polarization response */
#define DELAYDT 60 /* (integral) number of seconds between table
		      entries for propagation delay */

/* A macro that takes a detector time (in units of output->deltaT from
   output->epoch) and adds the propagation time interpolated from the
   delays tabulated in delayData.  The following parameters are
   required to be defined outside the macro, and are set by it:

   REAL4 realIndex;  the interpolation point in delayData
   INT4 intIndex;    the index immediately preceding realIndex
   REAL4 indexFrac;  the value realIndex - intIndex */
#define TCENTRE( time )                                              \
(                                                                    \
 realIndex = delayOff + (time)*delayDt,                              \
 intIndex = (INT4)floor( realIndex ),                                \
 indexFrac = realIndex - intIndex,                                   \
 time + indexFrac*delayData[intIndex+1]                              \
 + (1.0-indexFrac)*delayData[intIndex]                               \
)


NRCSID( SIMULATECOHERENTGWC, "$Id$" );

/* <lalVerbatim file="SimulateCoherentGWCP"> */
void
LALSimulateCoherentGW( LALStatus        *stat,
		       REAL4TimeSeries  *output,
		       CoherentGW       *signal,
		       DetectorResponse *detector )
{ /* </lalVerbatim> */
  INT4 aDim, phiDim; /* dimension of amplitude and phase vectors */
  INT4 i, n;         /* index over output->data, and its final value */
  INT4 nMax;         /* used to store limits on index ranges */
  REAL4 *outData;    /* pointer to output data */
  REAL4 radius;      /* radius of the Earth in lightseconds */
  CreateVectorSequenceIn in; /* structure for creating vectors and
                                vector sequences */
  BOOLEAN fFlag = 0; /* 1 if frequency left detector->transfer range */
  LIGOTimeGPS gpsTime; /* current detector time */

  /* The amplitude, phase, polarization response, and propagation
     delay are stored in arrays that must be interpolated.  For a
     quantity x, we define a pointer xData to the data array.  At some
     time t measured in units of output->deltaT, the interpolation
     point in xData is given by ( xOff + t*xDt ), where xOff is an
     offset and xDt is a relative sampling rate. */
  REAL4VectorSequence *polResponse = NULL;
  REAL4Vector *delay = NULL;
  REAL4 *aData, *phiData, *polData, *delayData;
  REAL4 aOff, phiOff, polOff, delayOff;
  REAL4 aDt, phiDt, polDt, delayDt;

  /* Frequencies in the detector transfer function are interpolated
     similarly, except everything is normalized with respect to
     detector->transfer->deltaF. */
  REAL4Vector *aTransfer = NULL;
  REAL4Vector *phiTransfer = NULL;
  REAL4 *aTransData, *phiTransData;
  REAL4 fOff;
  REAL4 df;

  /* Variables required by the TCENTRE() macro, above. */
  REAL4 realIndex;
  INT4 intIndex;
  REAL4 indexFrac;

  INITSTATUS( stat, "LALSimulateCoherentGW", SIMULATECOHERENTGWC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structures and their fields exist. */
  ASSERT( signal, stat, SIMULATECOHERENTGWH_ENUL,
	  SIMULATECOHERENTGWH_MSGENUL );
  if ( !( signal->a ) ) {
    ABORT( stat, SIMULATECOHERENTGWH_ESIG,
	   SIMULATECOHERENTGWH_MSGESIG );
  }
  ASSERT( signal->a->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  ASSERT( signal->a->data->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  if ( !( signal->phi ) ) {
    ABORT( stat, SIMULATECOHERENTGWH_ESIG,
	   SIMULATECOHERENTGWH_MSGESIG );
  }
  ASSERT( signal->phi->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  ASSERT( signal->phi->data->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  ASSERT( detector, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  ASSERT( detector->transfer, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  ASSERT( detector->transfer->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  ASSERT( detector->transfer->data->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  ASSERT( output, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  ASSERT( output->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  ASSERT( output->data->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );

  /* Define temporary variables to access the dimensions and data of
     signal->a and signal->phi, as well as detector->transfer. */
  aDim = signal->a->data->vectorLength;
  aData = signal->a->data->data;
  phiDim = signal->phi->data->vectorLength;
  phiData = signal->phi->data->data;
  outData = output->data->data;
  ASSERT( ( aDim == 1 ) || ( aDim == 2 ), stat,
	  SIMULATECOHERENTGWH_EDIM, SIMULATECOHERENTGWH_MSGEDIM );
  ASSERT( ( phiDim == 1 ) || ( phiDim == 2 ), stat,
	  SIMULATECOHERENTGWH_EDIM, SIMULATECOHERENTGWH_MSGEDIM );

  /* Make sure we never divide by zero. */
  ASSERT( signal->a->deltaT != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  ASSERT( signal->phi->deltaT != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  ASSERT( detector->transfer->deltaF != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  aDt = output->deltaT / signal->a->deltaT;
  phiDt = output->deltaT / signal->a->deltaT;
  df = 1.0 /
    ( LAL_TWOPI*signal->phi->deltaT*detector->transfer->deltaF );
  fOff = detector->transfer->f0/detector->transfer->deltaF;
  ASSERT( aDt != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  ASSERT( phiDt != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );

  /* Generate the table of propagation delays. */
  delayDt = output->deltaT/DELAYDT;
  in.length = (INT4)( output->data->length*output->deltaT/DELAYDT ) + 1;
  TRY( LALSCreateVector( stat->statusPtr, &delay, in.length ), stat );
  delayData = delay->data;
  radius = ( LAL_REARTH_SI/LAL_C_SI )*cos( detector->latitude ) /
    output->deltaT;
  gpsTime.gpsSeconds = output->epoch.gpsSeconds;
  gpsTime.gpsNanoSeconds = output->epoch.gpsNanoSeconds;
  for ( i = 0; i < (INT4)( in.length ); i++ ) {
    LIGOTimeUnix unixTime; /* GPS time at current point */
    LALDate date;          /* LALDate time at current point */
    REAL8 lmst;            /* siderial time at current point (radians) */

    LALGPStoU( stat->statusPtr, &unixTime, &gpsTime );
    BEGINFAIL( stat )
      TRY( LALSDestroyVector( stat->statusPtr, &delay ), stat );
    ENDFAIL( stat );
    LALUtime( stat->statusPtr, &date, &unixTime );
    BEGINFAIL( stat )
      TRY( LALSDestroyVector( stat->statusPtr, &delay ), stat );
    ENDFAIL( stat );
    LALLMST1( stat->statusPtr, &lmst, &date, detector->longitude,
	      MST_RAD );
    BEGINFAIL( stat )
      TRY( LALSDestroyVector( stat->statusPtr, &delay ), stat );
    ENDFAIL( stat );
    delayData[i] = radius*cos( lmst - signal->ra );
    gpsTime.gpsSeconds += DELAYDT;
  }

  /* Generate the table of polarization response functions. */
  polDt = output->deltaT/POLDT;
  in.length = (INT4)( output->data->length*output->deltaT/POLDT ) + 1;
  in.vectorLength = 2;
  LALSCreateVectorSequence( stat->statusPtr, &polResponse, &in );
  BEGINFAIL( stat )
    TRY( LALSDestroyVector( stat->statusPtr, &delay ), stat );
  ENDFAIL( stat );
  polData = polResponse->data;
  for ( i = 0; i < 2*(INT4)( in.length ); i+=2 ) {
    /* This is a temporary kludge... */
    polData[i] = 1.0;
    polData[i+1] = 0.0;
  }

  /* Decompose the transfer function into an amplitude and phase
     response. */
  in.length = detector->transfer->data->length;
  LALSCreateVector( stat->statusPtr, &phiTransfer, in.length );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr, &delay ), stat );
    TRY( LALSDestroyVectorSequence( stat->statusPtr, &polResponse ),
	 stat );
  } ENDFAIL( stat );
  LALCVectorAngle( stat->statusPtr, phiTransfer,
		   detector->transfer->data );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr, &delay ), stat );
    TRY( LALSDestroyVector( stat->statusPtr, &phiTransfer ), stat );
    TRY( LALSDestroyVectorSequence( stat->statusPtr, &polResponse ),
	 stat );
  } ENDFAIL( stat );
  LALUnwrapREAL4Angle( stat->statusPtr, phiTransfer, phiTransfer );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr, &delay ), stat );
    TRY( LALSDestroyVector( stat->statusPtr, &phiTransfer ), stat );
    TRY( LALSDestroyVectorSequence( stat->statusPtr, &polResponse ),
	 stat );
  } ENDFAIL( stat );
  LALSCreateVector( stat->statusPtr, &aTransfer, in.length );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr, &delay ), stat );
    TRY( LALSDestroyVector( stat->statusPtr, &phiTransfer ), stat );
    TRY( LALSDestroyVectorSequence( stat->statusPtr, &polResponse ),
	 stat );
  } ENDFAIL( stat );
  LALCVectorAbs( stat->statusPtr, aTransfer,
		 detector->transfer->data );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr, &delay ), stat );
    TRY( LALSDestroyVector( stat->statusPtr, &phiTransfer ), stat );
    TRY( LALSDestroyVector( stat->statusPtr, &aTransfer ), stat );
    TRY( LALSDestroyVectorSequence( stat->statusPtr, &polResponse ),
	 stat );
  } ENDFAIL( stat );
  phiTransData = phiTransfer->data;
  aTransData = aTransfer->data;

  /* Compute offsets for interpolating the signal, delay, and
     response functions. */
  aOff = ( output->epoch.gpsSeconds -
	   signal->a->epoch.gpsSeconds ) /
    signal->a->deltaT;
  aOff += ( output->epoch.gpsNanoSeconds -
	    signal->a->epoch.gpsNanoSeconds ) * 1.0e-9 /
    signal->a->deltaT;
  phiOff = ( output->epoch.gpsSeconds -
	   signal->phi->epoch.gpsSeconds ) /
    signal->phi->deltaT;
  phiOff += ( output->epoch.gpsNanoSeconds -
	      signal->phi->epoch.gpsNanoSeconds ) * 1.0e-9 /
    signal->phi->deltaT;
  polOff = 0.0;
  delayOff = 0.0;

  /* Compute initial value of i, ensuring that we will never index
     signal->a or signal->phi below their range. */
  i = 0;
  if ( aOff + ( i - radius )*aDt < 0.0 ) {
    i = (INT4)( -aOff/aDt - radius );
    while ( ( i < (INT4)( output->data->length ) ) &&
	    ( aOff + TCENTRE( i )*aDt < 0.0 ) )
      i++;
  }
  if ( phiOff + ( i - radius )*phiDt < 0.0 ) {
    i = (INT4)( -phiOff/phiDt - radius );
    while ( ( i < (INT4)( output->data->length ) ) && 
	    ( phiOff + TCENTRE( i )*phiDt < 0.0 ) )
      i++;
  }
  if ( i >= (INT4)( output->data->length ) )
    LALWarning( stat, "Signal starts after the end of the output"
		" time series." );
  if ( i < 0 )
    i = 0;

  /* Compute final value of i, ensuring that we will never index
     signal->a or signal->phi above their range. */
  n = output->data->length - 1;
  nMax = signal->a->data->length - 1;
  if ( aOff + ( n + radius )*aDt > nMax ) {
    n = (INT4)( ( nMax - aOff )/aDt + radius + 1.0 );
    while ( ( n >= 0 ) &&
	    ( aOff + TCENTRE( n )*aDt >= nMax ) )
      n--;
  }
  nMax = signal->phi->data->length - 1;
  if ( phiOff + ( n + radius )*phiDt > nMax ) {
    n = (INT4)( ( nMax - phiOff )/phiDt + radius + 1.0 );
    while ( ( n >= 0 ) &&
	    ( phiOff + TCENTRE( n )*phiDt >= nMax ) )
      n--;
  }
  if ( n < 0 )
    LALWarning( stat, "Signal ends before the start of the output"
		" time series." );
  if ( n >= (INT4)( output->data->length ) )
    n = output->data->length - 1;

  /* Set output to zero where the signal is not defined. */
  if ( i > 0 )
    memset( output->data->data, 0, i*sizeof(REAL4) );
  if ( ( nMax = output->data->length - n ) > 0 )
    memset( output->data->data + n + 1, 0, nMax*sizeof(REAL4) );

  /* Keep track of the frequency range of the transfer function, so
     that we don't try to interpolate it out of its range. */
  nMax = detector->transfer->data->length - 1;

  /* Start computing responses. */
  for ( ; i <= n; i++ ) {
    REAL8 iCentre = TCENTRE( i );  /* value of i + propagation delays */
    REAL8 x;                       /* interpolation point in arrays */
    INT4 j;                        /* array index preceding x */
    REAL4 frac;                    /* value of x - j */
    REAL4 aPlus, aCross;           /* current signal amplitudes */
    REAL4 phiPlus, phiCross = 0.0; /* current signal phases */
    REAL4 fPlus, fCross = 0.0;     /* current signal frequencies */
    REAL4 aTrans, phiTrans;        /* values of the transfer function */

    /* Interpolate the signal amplitude. */
    x = aOff + iCentre*aDt;
    j = (INT4)floor( x );
    frac = (REAL4)( x - j );
    if ( aDim == 2 ) {
      j *= 2;
      aPlus = frac*aData[j+2] + ( 1.0 - frac )*aData[j];
      aCross = frac*aData[j+3] + ( 1.0 - frac )*aData[j+1];
    } else {
      aPlus = frac*aData[j+1] + ( 1.0 - frac )*aData[j];
      aCross = aPlus;
    }

    /* Interpolate the signal phase, and find the frequency. */
    x = phiOff + iCentre*phiDt;
    j = (INT4)floor( x );
    frac = (REAL4)( x - j );
    if ( phiDim == 2 ) {
      j *= 2;
      phiPlus = frac*phiData[j+2] + ( 1.0 - frac )*phiData[j];
      fPlus = ( phiData[j+2] - phiData[j] )*df;
      phiCross = frac*phiData[j+3] + ( 1.0 - frac )*phiData[j+1];
      fCross = ( phiData[j+3] - phiData[j+1] )*df;
    } else {
      phiPlus = frac*phiData[j+1] + ( 1.0 - frac )*phiData[j];
      fPlus = ( phiData[j+1] - phiData[j] )*df;
    }

    /* Interpolate the polarization response. */
    x = polOff + iCentre*polDt;
    j = (INT4)floor( x );
    frac = (REAL4)( x - j );
    j *= 2;
    aPlus *= frac*polData[j+2] + ( 1.0 - frac )*polData[j];
    aCross *= frac*polData[j+3] + ( 1.0 - frac )*polData[j+1];

    /* Interpolate the transfer function and compute the output. */
    x = fPlus - fOff;
    if ( ( x < 0.0 ) || ( x >= nMax ) ) {
      aTrans = 0.0;
      phiTrans = 0.0;
      fFlag = 1;
    } else {
      j = (INT4)floor( x );
      frac = (REAL4)( x - j );
      aTrans = frac*aTransData[j+1] + ( 1.0 - frac )*aTransData[j];
      phiTrans = frac*phiTransData[j+1] + ( 1.0 - frac )*phiTransData[j];
    }
    aPlus *= aTrans;
    phiPlus += phiTrans;
    if ( phiDim == 2 ) {
      x = fCross - fOff;
      if ( ( x < 0.0 ) || ( x >= nMax ) ) {
	aTrans = 0.0;
	phiTrans = 0.0;
	fFlag = 1;
      } else {
	j = (INT4)floor( x );
	frac = (REAL4)( x - j );
	aTrans = frac*aTransData[j+1] + ( 1.0 - frac )*aTransData[j];
	phiTrans = frac*phiTransData[j+1] + ( 1.0 - frac )*phiTransData[j];
      }
      aCross *= aTrans;
      phiCross += phiTrans;
      outData[i] = aPlus*cos( phiPlus ) + aCross*cos( phiCross );
    } else {
      aCross *= aTrans;
      outData[i] = aPlus*cos( phiPlus ) + aCross*sin( phiPlus );
    }
  }

  /* Warn if we ever stepped outside of the frequency domain of the
     transfer function. */
  if ( fFlag )
    LALWarning( stat, "Signal passed outside of the frequency domain"
		" of the transfer function (transfer function is"
		" treated as zero outside its specified domain)" );

  /* Cleanup and exit. */
  TRY( LALSDestroyVector( stat->statusPtr, &delay ), stat );
  TRY( LALSDestroyVector( stat->statusPtr, &phiTransfer ), stat );
  TRY( LALSDestroyVector( stat->statusPtr, &aTransfer ), stat );
  TRY( LALSDestroyVectorSequence( stat->statusPtr, &polResponse ),
       stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
