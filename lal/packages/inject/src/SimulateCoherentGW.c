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

This function takes a quasiperiodic gravitational waveform given in
\verb@*signal@, and estimates the corresponding response of the
detector whose polarization response and transfer function are
specified in \verb@*detector@.  The result is stored in
\verb@*output@.

The fields \verb@output->epoch@, \verb@output->deltaT@, and
\verb@output->data@ must already be set, in order to specify the time
period and sampling rate for which the response is required.  For the
input signal, \verb@signal->h@ is ignored, and the signal is treated
as zero at any time for which either \verb@signal->a@ or
\verb@signal->phi@ is not defined.

This routine will convert \verb@signal->position@ to equatorial
coordinates, if necessary.

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
quasiperiodic limit as follows:
\begin{itemize}
\item At each sample point in \verb@*output@, the propagation delay is
computed and added to the sample time, and the instantaneous
amplitudes $A_1$, $A_2$, frequency $f$, phase $\phi$, and polarization
shift $\Phi$ are found by interpolating the nearest values in
\verb@signal->a@, \verb@signal->f@, \verb@signal->phi@, and
\verb@signal->shift@, respectively.  If \verb@signal->f@ is not
defined at that point in time, then $f$ is estimated by differencing
the two nearest values of $\phi$, as $f\approx\Delta\phi/2\pi\Delta
t$.  If \verb@signal->shift@ is not defined, then $\Phi$ is treated as
zero.
\item The complex transfer function of the detector at this frequency
is found by interpolating \verb@detector->transfer@.  The amplitude of
the transfer function is multiplied with $A_1$ and $A_2$, and the phase of the
transfer function is added to $\phi$,
\item The plus and cross contributions $o_+$, $o_\times$ to the
detector output are computed as in Eqs.~\ref{eq:quasiperiodic-hplus}
and~\ref{eq:quasiperiodic-hcross} of \verb@SimulateCoherentGW.h@, but
using the response-adjusted amplitudes and phase.
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
LALSCreateVector()              LALSCreateVectorSequence()
LALSDestroyVector()             LALSDestroyVectorSequence()
LALConvertSkyCoordinates()      LALGeodeticToGeocentric()
LALCVectorAbs()                 LALCVectorAngle()
LALUnwrapREAL4Angle()
LALWarning()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{SimulateCoherentGWCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/DetectorSite.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>
#include <lal/VectorOps.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/SkyCoordinates.h>

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
  INT4 i, n;          /* index over output->data, and its final value */
  INT4 nMax;          /* used to store limits on index ranges */
  INT4 fInit, fFinal; /* index range for which signal->f is defined */
  INT4 shiftInit, shiftFinal; /* ditto for signal->shift */
  REAL4 *outData;             /* pointer to output data */
  REAL8 delayMin, delayMax;   /* min and max values of time delay */
  CreateVectorSequenceIn in;  /* structure for creating vectors and
				 vector sequences */
  BOOLEAN fFlag = 0; /* 1 if frequency left detector->transfer range */
  BOOLEAN pFlag = 0; /* 1 if frequency was estimated from phase */

  /* The amplitude, frequency, phase, polarization shift, polarization
     response, and propagation delay are stored in arrays that must be
     interpolated.  For a quantity x, we define a pointer xData to the
     data array.  At some time t measured in units of output->deltaT,
     the interpolation point in xData is given by ( xOff + t*xDt ),
     where xOff is an offset and xDt is a relative sampling rate. */
  REAL4VectorSequence *polResponse = NULL;
  REAL8Vector *delay = NULL;
  REAL4 *aData, *fData, *phiData, *shiftData, *polData;
  REAL8 *delayData;
  REAL8 aOff, fOff, phiOff, shiftOff, polOff, delayOff;
  REAL8 aDt, fDt, phiDt, shiftDt, polDt, delayDt;

  /* Frequencies in the detector transfer function are interpolated
     similarly, except everything is normalized with respect to
     detector->transfer->deltaF. */
  REAL4Vector *aTransfer = NULL;
  REAL4Vector *phiTransfer = NULL;
  REAL4Vector *phiTemp = NULL;
  REAL4 *aTransData, *phiTransData;
  REAL8 f0;
  REAL8 phiFac, fFac;

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
  if ( signal->f ) {
    ASSERT( signal->f->data, stat,
	    SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
    ASSERT( signal->f->data->data, stat,
	    SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  }
  if ( signal->shift ) {
    ASSERT( signal->shift->data, stat,
	    SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
    ASSERT( signal->shift->data->data, stat,
	    SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  }
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

  /* Check dimensions of amplitude array. */
  ASSERT( signal->a->data->vectorLength == 2, stat,
	  SIMULATECOHERENTGWH_EDIM, SIMULATECOHERENTGWH_MSGEDIM );

  /* Make sure we never divide by zero. */
  ASSERT( signal->a->deltaT != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  ASSERT( signal->phi->deltaT != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  aDt = output->deltaT / signal->a->deltaT;
  phiDt = output->deltaT / signal->a->deltaT;
  ASSERT( aDt != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  ASSERT( phiDt != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  if ( signal->f ) {
    ASSERT( signal->f->deltaT != 0.0, stat,
	    SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
    fDt = output->deltaT / signal->f->deltaT;
    ASSERT( fDt != 0.0, stat,
	    SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  } else
    fDt = 0.0;
  if ( signal->shift ) {
    ASSERT( signal->shift->deltaT != 0.0, stat,
	    SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
    shiftDt = output->deltaT / signal->shift->deltaT;
    ASSERT( shiftDt != 0.0, stat,
	    SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  } else
    shiftDt = 0.0;
  ASSERT( detector->transfer->deltaF != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  fFac = 1.0 / detector->transfer->deltaF;
  phiFac = fFac / ( LAL_TWOPI*signal->phi->deltaT );
  f0 = detector->transfer->f0/detector->transfer->deltaF;

  /* Define temporary variables to access the data of signal->a,
     signal->f, and signal->phi, as well as detector->transfer. */
  aData = signal->a->data->data;
  phiData = signal->phi->data->data;
  outData = output->data->data;
  if ( signal->f )
    fData = signal->f->data->data;
  else
    fData = NULL;
  if ( signal->shift )
    shiftData = signal->shift->data->data;
  else
    shiftData = NULL;

  /* Generate the table of propagation delays. */
  delayDt = output->deltaT/DELAYDT;
  in.length = (UINT4)( output->data->length*output->deltaT/DELAYDT )
    + 3;
  TRY( LALDCreateVector( stat->statusPtr, &delay, in.length ), stat );
  delayData = delay->data;
  if ( detector->site ) {
    SkyPosition source;      /* position of source on sky */
    LIGOTimeGPS gpsTime;     /* detector time when we compute delay */
    LALPlaceAndGPS event;    /* spacetime point where we compute delay */
    DetTimeAndASource input; /* input to time delay function */

    /* Arrange nested pointers, and set initial values. */
    event.p_detector = detector->site;
    event.p_gps = &gpsTime;
    input.p_det_and_time = &event;
    input.p_source = &source;
    gpsTime = output->epoch;
    gpsTime.gpsSeconds -= DELAYDT/2;
    delayMin = delayMax = LAL_REARTH_SI / ( LAL_C_SI*output->deltaT );
    delayMin *= -1;

    /* Convert source position to equatorial coordinates. */
    source = signal->position;
    if ( source.system != COORDINATESYSTEM_EQUATORIAL ) {
      ConvertSkyParams params; /* parameters for conversion */
      EarthPosition location;  /* location of detector */
      params.gpsTime = &gpsTime;
      params.system = COORDINATESYSTEM_EQUATORIAL;
      if ( source.system == COORDINATESYSTEM_HORIZON ) {
	params.zenith = &( location.geodetic );
	location.x = detector->site->location[0];
	location.y = detector->site->location[1];
	location.z = detector->site->location[2];
	LALGeodeticToGeocentric( stat->statusPtr, &location );
	BEGINFAIL( stat )
	  TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
	ENDFAIL( stat );
      }
      LALConvertSkyCoordinates( stat->statusPtr, &source, &source,
				&params );
      BEGINFAIL( stat )
	TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
      ENDFAIL( stat );
    }

    /* Compute table. */
    for ( i = 0; (UINT4)( i ) < in.length; i++ ) {
      REAL8 tDelay; /* propagation time */
      LALTimeDelayFromEarthCenter( stat->statusPtr, &tDelay, &input );
      BEGINFAIL( stat )
	TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
      ENDFAIL( stat );
      /* TimeDelayFromEarthCenter() measures propagation delay from
         geocentre to detector, which is the opposite of what we want.
         We also want it normalized. */
      tDelay /= -output->deltaT;
      delayData[i] = tDelay;
      if ( tDelay < delayMin )
	delayMin = tDelay;
      if ( tDelay > delayMax )
	delayMax = tDelay;
      gpsTime.gpsSeconds += DELAYDT;
    }
  } else {
    LALInfo( stat, "Detector site absent; simulating h_plus with no"
	     " propagation delays." );
    memset( delayData, 0, in.length*sizeof(REAL8) );
    delayMin = delayMax = 0.0;
  }

  /* Generate the table of polarization response functions. */
  polDt = output->deltaT/POLDT;
  in.length = (UINT4)( output->data->length*output->deltaT/POLDT )
    + 3;
  in.vectorLength = 2;
  LALSCreateVectorSequence( stat->statusPtr, &polResponse, &in );
  BEGINFAIL( stat )
    TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
  ENDFAIL( stat );
  polData = polResponse->data;
  for ( i = 0; (UINT4)( i ) < 2*in.length; i+=2 ) {
    /* This is a temporary kludge... */
    polData[i] = 1.0;
    polData[i+1] = 0.0;
  }

  /* Decompose the transfer function into an amplitude and phase
     response. */
  in.length = detector->transfer->data->length;
  LALSCreateVector( stat->statusPtr, &phiTemp, in.length );
  BEGINFAIL( stat ) {
    TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
    TRY( LALSDestroyVectorSequence( stat->statusPtr, &polResponse ),
	 stat );
  } ENDFAIL( stat );
  LALCVectorAngle( stat->statusPtr, phiTemp,
		   detector->transfer->data );
  BEGINFAIL( stat ) {
    TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
    TRY( LALSDestroyVector( stat->statusPtr, &phiTemp ), stat );
    TRY( LALSDestroyVectorSequence( stat->statusPtr, &polResponse ),
	 stat );
  } ENDFAIL( stat );
  LALSCreateVector( stat->statusPtr, &phiTransfer, in.length );
  BEGINFAIL( stat ) {
    TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
    TRY( LALSDestroyVector( stat->statusPtr, &phiTemp ), stat );
    TRY( LALSDestroyVectorSequence( stat->statusPtr, &polResponse ),
	 stat );
  } ENDFAIL( stat );
  LALUnwrapREAL4Angle( stat->statusPtr, phiTransfer, phiTemp );
  BEGINFAIL( stat ) {
    TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
    TRY( LALSDestroyVector( stat->statusPtr, &phiTemp ), stat );
    TRY( LALSDestroyVector( stat->statusPtr, &phiTransfer ), stat );
    TRY( LALSDestroyVectorSequence( stat->statusPtr, &polResponse ),
	 stat );
  } ENDFAIL( stat );
  TRY( LALSDestroyVector( stat->statusPtr, &phiTemp ), stat );
  LALSCreateVector( stat->statusPtr, &aTransfer, in.length );
  BEGINFAIL( stat ) {
    TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
    TRY( LALSDestroyVector( stat->statusPtr, &phiTransfer ), stat );
    TRY( LALSDestroyVectorSequence( stat->statusPtr, &polResponse ),
	 stat );
  } ENDFAIL( stat );
  LALCVectorAbs( stat->statusPtr, aTransfer,
		 detector->transfer->data );
  BEGINFAIL( stat ) {
    TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
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
  if ( signal->f ) {
    fOff = ( output->epoch.gpsSeconds -
	     signal->f->epoch.gpsSeconds ) /
      signal->f->deltaT;
    fOff += ( output->epoch.gpsNanoSeconds -
	      signal->f->epoch.gpsNanoSeconds ) * 1.0e-9 /
      signal->f->deltaT;
  } else
    fOff = 0.0;
  if ( signal->shift ) {
    shiftOff = ( output->epoch.gpsSeconds -
		 signal->shift->epoch.gpsSeconds ) /
      signal->shift->deltaT;
    shiftOff += ( output->epoch.gpsNanoSeconds -
		  signal->shift->epoch.gpsNanoSeconds ) * 1.0e-9 /
      signal->shift->deltaT;
  } else
    shiftOff = 0.0;
  polOff = 0.5;
  delayOff = 0.5;

  /* Compute initial value of i, ensuring that we will never index
     signal->a or signal->phi below their range. */
  i = 0;
  if ( aOff + ( i + delayMin )*aDt < 0.0 ) {
    i = (INT4)( -aOff/aDt + delayMin );
    while ( ( i < (INT4)( output->data->length ) ) &&
	    ( aOff + TCENTRE( i )*aDt < 0.0 ) )
      i++;
  }
  if ( phiOff + ( i + delayMin )*phiDt < 0.0 ) {
    i = (INT4)( -phiOff/phiDt + delayMin );
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
  if ( aOff + ( n + delayMax )*aDt > nMax ) {
    n = (INT4)( ( nMax - aOff )/aDt + delayMax + 1.0 );
    while ( ( n >= 0 ) &&
	    ( aOff + TCENTRE( n )*aDt >= nMax ) )
      n--;
  }
  nMax = signal->phi->data->length - 1;
  if ( phiOff + ( n + delayMax )*phiDt > nMax ) {
    n = (INT4)( ( nMax - phiOff )/phiDt + delayMax + 1.0 );
    while ( ( n >= 0 ) &&
	    ( phiOff + TCENTRE( n )*phiDt >= nMax ) )
      n--;
  }
  if ( n < 0 )
    LALWarning( stat, "Signal ends before the start of the output"
		" time series." );
  if ( n >= (INT4)( output->data->length ) )
    n = output->data->length - 1;

  /* Compute the values of i for which signal->f is given. */
  if ( signal->f ) {
    fInit = (INT4)( -fOff/fDt + delayMin );
    while ( ( fInit < (INT4)( output->data->length ) ) &&
	    ( fOff + TCENTRE( fInit )*fDt < 0.0 ) )
      fInit++;
    fFinal = (INT4)( ( nMax - fOff )/fDt + delayMax + 1.0 );
    while ( ( fFinal >= 0 ) &&
	    ( fOff + TCENTRE( fFinal )*fDt >= nMax ) )
      fFinal--;
  } else {
    fInit = n + 1;
    fFinal = i - 1;
  }

  /* Compute the values of i for which signal->shift is given. */
  if ( signal->shift ) {
    shiftInit = (INT4)( -shiftOff/shiftDt + delayMin );
    while ( ( shiftInit < (INT4)( output->data->length ) ) &&
	    ( shiftOff + TCENTRE( shiftInit )*shiftDt < 0.0 ) )
      shiftInit++;
    shiftFinal = (INT4)( ( nMax - shiftOff )/shiftDt + delayMax + 1.0 );
    while ( ( shiftFinal >= 0 ) &&
	    ( shiftOff + TCENTRE( shiftFinal )*shiftDt >= nMax ) )
      shiftFinal--;
  } else {
    shiftInit = n + 1;
    shiftFinal = i - 1;
  }

  /* Set output to zero where the signal is not defined. */
  if ( i > 0 )
    memset( output->data->data, 0, i*sizeof(REAL4) );
  if ( ( nMax = output->data->length - n - 1 ) > 0 )
    memset( output->data->data + n + 1, 0, nMax*sizeof(REAL4) );

  /* Keep track of the frequency range of the transfer function, so
     that we don't try to interpolate it out of its range. */
  nMax = detector->transfer->data->length - 1;

  /* Start computing responses. */
  for ( ; i <= n; i++ ) {
    REAL8 iCentre = TCENTRE( i );  /* value of i + propagation delays */
    REAL8 x;                /* interpolation point in arrays */
    INT4 j;                 /* array index preceding x */
    REAL4 frac;             /* value of x - j */
    REAL4 a1, a2;           /* current signal amplitudes */
    REAL4 phi = 0.0;        /* current signal phase */
    REAL4 f = 0.0;          /* current signal frequency */
    REAL4 shift = 0.0;      /* current signal polarization shift */
    REAL4 aTrans, phiTrans; /* current values of the transfer function */
    REAL4 oPlus, oCross;    /* current output amplitudes */

    /* Interpolate the signal amplitude. */
    x = aOff + iCentre*aDt;
    j = (INT4)floor( x );
    frac = (REAL4)( x - j );
    j *= 2;
    a1 = frac*aData[j+2] + ( 1.0 - frac )*aData[j];
    a2 = frac*aData[j+3] + ( 1.0 - frac )*aData[j+1];

    /* Interpolate the polarization shift. */
    if ( ( i < shiftInit ) || ( i > shiftFinal ) )
      shift = 0.0;
    else {
      x = shiftOff + iCentre*shiftDt;
      j = (INT4)floor( x );
      frac = (REAL4)( x - j );
      shift = frac*shiftData[j+1] + ( 1.0 - frac )*shiftData[j];
    }

    /* Interpolate the signal phase, and find the frequency. */
    x = phiOff + iCentre*phiDt;
    j = (INT4)floor( x );
    frac = (REAL4)( x - j );
    phi = frac*phiData[j+1] + ( 1.0 - frac )*phiData[j];
    if ( ( i < fInit ) || ( i > fFinal ) ) {
      f = ( phiData[j+1] - phiData[j] )*phiFac;
      pFlag = 1;
    } else {
      x = fOff + iCentre*fDt;
      j = (INT4)floor( x );
      frac = (REAL4)( x - j );
      f = frac*fData[j+1] + ( 1.0 - frac )*fData[j];
      f *= fFac;
    }

    /* Interpolate the transfer function. */
    x = f - f0;
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
    a1 *= aTrans;
    a2 *= aTrans;
    phi += phiTrans;

    /* Compute components of output. */
    oPlus = a1*cos( shift )*cos( phi ) - a2*sin( shift )*sin( phi );
    oCross = a1*sin( shift )*cos( phi ) + a2*cos( shift )*sin( phi );

    /* Interpolate the polarization response, and compute output. */
    x = polOff + iCentre*polDt;
    j = (INT4)floor( x );
    frac = (REAL4)( x - j );
    j *= 2;
    oPlus *= frac*polData[j+2] + ( 1.0 - frac )*polData[j];
    oCross *= frac*polData[j+3] + ( 1.0 - frac )*polData[j+1];
    outData[i] = oPlus + oCross;
  }

  /* Warn if we ever stepped outside of the frequency domain of the
     transfer function, or if we had to estimate f from phi. */
  if ( fFlag )
    LALWarning( stat, "Signal passed outside of the frequency domain"
		" of the transfer function (transfer function is"
		" treated as zero outside its specified domain)" );
  if ( pFlag )
    LALInfo( stat, "Signal frequency was estimated by differencing"
	     " the signal phase" );

  /* Cleanup and exit. */
  TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
  TRY( LALSDestroyVector( stat->statusPtr, &phiTransfer ), stat );
  TRY( LALSDestroyVector( stat->statusPtr, &aTransfer ), stat );
  TRY( LALSDestroyVectorSequence( stat->statusPtr, &polResponse ),
       stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
