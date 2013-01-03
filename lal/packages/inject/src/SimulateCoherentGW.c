/*
*  Copyright (C) 2007 Jolien Creighton, Reinhard Prix, Stephen Fairhurst, Teviet Creighton
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

/* Use more efficient trig routines for solaris, if available and
   requested. */
#include <config.h>
#ifdef HAVE_SUNMATH_H
#include <sunmath.h>
#if defined HAVE_LIBSUNMATH && defined ONLINE
#define USE_SINCOSP 1
#endif
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/TimeDelay.h>
#include <lal/LALBarycenter.h>
#include <lal/VectorOps.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/SkyCoordinates.h>

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

/**
\author Creighton, T. D.

\brief Computes the response of a detector to a coherent gravitational wave.

This function takes a quasiperiodic gravitational waveform given in
<tt>*signal</tt>, and estimates the corresponding response of the
detector whose position, orientation, and transfer function are
specified in <tt>*detector</tt>.  The result is stored in
<tt>*output</tt>.

The fields <tt>output-\>epoch</tt>, <tt>output->deltaT</tt>, and
<tt>output-\>data</tt> must already be set, in order to specify the time
period and sampling rate for which the response is required.  If
<tt>output-\>f0</tt> is nonzero, idealized heterodyning is performed (an
amount \f$2\pi f_0(t-t_0)\f$ is subtracted from the phase before computing
the sinusoid, where \f$t_0\f$ is the heterodyning epoch defined in
\c detector).  For the input signal, <tt>signal-\>h</tt> is ignored,
and the signal is treated as zero at any time for which either
<tt>signal-\>a</tt> or <tt>signal-\>phi</tt> is not defined.

This routine will convert <tt>signal-\>position</tt> to equatorial
coordinates, if necessary.

\heading{Algorithm}

The routine first accounts for the time delay between the detector and
the solar system barycentre, based on the detector position
information stored in <tt>*detector</tt> and the propagation direction
specified in <tt>*signal</tt>.  Values of the propagation delay are
precomuted at fixed intervals and stored in a table, with the
intervals \f$\Delta T_\mathrm{delay}\f$ chosen such that the value
interpolated from adjacent table entries will never differ from the
true value by more than some timing error \f$\sigma_T\f$.  This implies
that:
\f[
\Delta T_\mathrm{delay} \leq \sqrt{
	\frac{8\sigma_T}{\max\{a/c\}} } \; ,
\f]
where \f$\max\{a/c\}=1.32\times10^{-10}\mathrm{s}^{-1}\f$ is the maximum
acceleration of an Earth-based detector in the barycentric frame.  The
total propagation delay also includes Einstein and Shapiro delay, but
these are more slowly varying and thus do not constrain the table
spacing.  At present, a 400s table spacing is hardwired into the code,
implying \f$\sigma_T\approx3\mu\f$s, comparable to the stated accuracy of
<tt>LALBarycenter()</tt>.

Next, the polarization response functions of the detector
\f$F_{+,\times}(\alpha,\delta)\f$ are computed for every 10 minutes of the
signal's duration, using the position of the source in <tt>*signal</tt>,
the detector information in <tt>*detector</tt>, and the function
<tt>LALComputeDetAMResponseSeries()</tt>.  Subsequently, the
polarization functions are estimated for each output sample by
interpolating these precomputed values.  This guarantees that the
interpolated value is accurate to \f$\sim0.1\%\f$.

Next, the frequency response of the detector is estimated in the
quasiperiodic limit as follows:
<ul>
<li> At each sample point in <tt>*output</tt>, the propagation delay is
computed and added to the sample time, and the instantaneous
amplitudes \f$A_1\f$, \f$A_2\f$, frequency \f$f\f$, phase \f$\phi\f$, and polarization
shift \f$\Phi\f$ are found by interpolating the nearest values in
<tt>signal-\>a</tt>, <tt>signal-\>f</tt>, <tt>signal-\>phi</tt>, and
<tt>signal-\>shift</tt>, respectively.  If <tt>signal-\>f</tt> is not
defined at that point in time, then \f$f\f$ is estimated by differencing
the two nearest values of \f$\phi\f$, as \f$f\approx\Delta\phi/2\pi\Delta
t\f$.  If <tt>signal-\>shift</tt> is not defined, then \f$\Phi\f$ is treated as
zero.</li>
<li> The complex transfer function of the detector the frequency \f$f\f$
is found by interpolating <tt>detector-\>transfer</tt>.  The amplitude of
the transfer function is multiplied with \f$A_1\f$ and \f$A_2\f$, and the
phase of the transfer function is added to \f$\phi\f$,</li>
<li> The plus and cross contributions \f$o_+\f$, \f$o_\times\f$ to the
detector output are computed as in Eqs.\eqref{eq_quasiperiodic_hpluscross}
of \ref SimulateCoherentGW_h, but
using the response-adjusted amplitudes and phase.</li>
<li> The final detector response \f$o\f$ is computed as
\f$o=(o_+F_+)+(o_\times F_\times)\f$.</li>
</ul>

\heading{A note on interpolation:}
Much of the computational work in this routine involves interpolating
various time series to find their values at specific output times.
The algorithm is summarized below.

Let \f$A_j = A( t_A + j\Delta t_A )\f$ be a sampled time series, which we
want to resample at new (output) time intervals \f$t_k = t_0 + k\Delta
t\f$.  We first precompute the following quantities:
\f{eqnarray}{
t_\mathrm{off} & = & \frac{t_0-t_A}{\Delta t_A}  \; , \nonumber \\
            dt & = & \frac{\Delta t}{\Delta t_A} \; . \nonumber
\f}
Then, for each output sample time \f$t_k\f$, we compute:
\f{eqnarray}{
t & = & t_\mathrm{off} + k \times dt \; , \nonumber \\
j & = & \lfloor t \rfloor            \; , \nonumber \\
f & = & t - j                        \; , \nonumber
\f}
where \f$\lfloor x\rfloor\f$ is the "floor" function; i.e.\ the largest
integer \f$\leq x\f$.  The time series sampled at the new time is then:
\f[
A(t_k) = f \times A_{j+1} + (1-f) \times A_j \; .
\f]

\heading{Notes}

The major computational hit in this routine comes from computing the
sine and cosine of the phase angle in
Eqs.\eqref{eq_quasiperiodic_hpluscross} of
\ref SimulateCoherentGW_h.  For better online performance, these can
be replaced by other (approximate) trig functions.  Presently the code
uses the native \c libm functions by default, or the function
<tt>sincosp()</tt> in \c libsunmath \e if this function is
available \e and the constant \c ONLINE is defined.
Differences at the level of 0.01 begin to appear only for phase
arguments greater than \f$10^{14}\f$ or so (corresponding to over 500
years between phase epoch and observation time for frequencies of
around 1kHz).

To activate this feature, be sure that \ref sunmath.h and
\c libsunmath are on your system, and add <tt>-DONLINE</tt> to the
<tt>--with-extra-cppflags</tt> configuration argument.  In future this
flag may be used to turn on other efficient trig algorithms on other
(non-Solaris) platforms.

*/
void
LALSimulateCoherentGW( LALStatus        *stat,
		       REAL4TimeSeries  *output,
		       CoherentGW       *CWsignal,
		       DetectorResponse *detector )
{
  INT4 i, n;          /* index over output->data, and its final value */
  INT4 nMax;          /* used to store limits on index ranges */
  INT4 fInit, fFinal; /* index range for which CWsignal->f is defined */
  INT4 shiftInit, shiftFinal; /* ditto for CWsignal->shift */
  UINT4 dtDelayBy2 = 400;     /* delay table half-interval (s) */
  UINT4 dtPolBy2 = 300;       /* polarization table half-interval (s) */
  REAL4 *outData;             /* pointer to output data */
  REAL8 delayMin, delayMax;   /* min and max values of time delay */
  SkyPosition source;         /* source sky position */
  BOOLEAN transfer;  /* 1 if transfer function is specified */
  BOOLEAN fFlag = 0; /* 1 if frequency left detector->transfer range */
  BOOLEAN pFlag = 0; /* 1 if frequency was estimated from phase */

  /* The amplitude, frequency, phase, polarization shift, polarization
     response, and propagation delay are stored in arrays that must be
     interpolated.  For a quantity x, we define a pointer xData to the
     data array.  At some time t measured in units of output->deltaT,
     the interpolation point in xData is given by ( xOff + t*xDt ),
     where xOff is an offset and xDt is a relative sampling rate. */
  LALDetAMResponseSeries polResponse;
  REAL8Vector *delay = NULL;
  REAL4 *aData, *fData, *shiftData, *plusData, *crossData;
  REAL8 *phiData, *delayData;
  REAL8 aOff, fOff, phiOff, shiftOff, polOff, delayOff;
  REAL8 aDt, fDt, phiDt, shiftDt, polDt, delayDt;

  /* Frequencies in the detector transfer function are interpolated
     similarly, except everything is normalized with respect to
     detector->transfer->deltaF. */
  REAL4Vector *aTransfer = NULL;
  REAL4Vector *phiTransfer = NULL;
  REAL4Vector *phiTemp = NULL;
  REAL4 *aTransData = NULL, *phiTransData = NULL;
  REAL8 f0 = 1.0;
  REAL8 phiFac = 1.0, fFac = 1.0;

  /* Heterodyning phase factor LAL_TWOPI*output->f0*output->deltaT,
     and phase offset at the start of the series
     LAL_TWOPI*output->f0*(time offset). */
  REAL8 heteroFac, phi0;

  /* Variables required by the TCENTRE() macro, above. */
  REAL8 realIndex;
  INT4 intIndex;
  REAL8 indexFrac;

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structures and their fields exist. */
  ASSERT( CWsignal, stat, SIMULATECOHERENTGWH_ENUL,
	  SIMULATECOHERENTGWH_MSGENUL );
  if ( !( CWsignal->a ) ) {
    ABORT( stat, SIMULATECOHERENTGWH_ESIG,
	   SIMULATECOHERENTGWH_MSGESIG );
  }
  ASSERT( CWsignal->a->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  ASSERT( CWsignal->a->data->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  if ( !( CWsignal->phi ) ) {
    ABORT( stat, SIMULATECOHERENTGWH_ESIG,
	   SIMULATECOHERENTGWH_MSGESIG );
  }
  ASSERT( CWsignal->phi->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  ASSERT( CWsignal->phi->data->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  if ( CWsignal->f ) {
    ASSERT( CWsignal->f->data, stat,
	    SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
    ASSERT( CWsignal->f->data->data, stat,
	    SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  }
  if ( CWsignal->shift ) {
    ASSERT( CWsignal->shift->data, stat,
	    SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
    ASSERT( CWsignal->shift->data->data, stat,
	    SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  }
  ASSERT( detector, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  if ( ( transfer = ( detector->transfer != NULL ) ) ) {
    ASSERT( detector->transfer->data, stat,
	    SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
    ASSERT( detector->transfer->data->data, stat,
	    SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  }
  ASSERT( output, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  ASSERT( output->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );
  ASSERT( output->data->data, stat,
	  SIMULATECOHERENTGWH_ENUL, SIMULATECOHERENTGWH_MSGENUL );

  /* Check dimensions of amplitude array. */
  ASSERT( CWsignal->a->data->vectorLength == 2, stat,
	  SIMULATECOHERENTGWH_EDIM, SIMULATECOHERENTGWH_MSGEDIM );

  /* Make sure we never divide by zero. */
  ASSERT( CWsignal->a->deltaT != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  ASSERT( CWsignal->phi->deltaT != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  aDt = output->deltaT / CWsignal->a->deltaT;
  phiDt = output->deltaT / CWsignal->phi->deltaT;
  ASSERT( aDt != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  ASSERT( phiDt != 0.0, stat,
	  SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  if ( CWsignal->f ) {
    ASSERT( CWsignal->f->deltaT != 0.0, stat,
	    SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
    fDt = output->deltaT / CWsignal->f->deltaT;
    ASSERT( fDt != 0.0, stat,
	    SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  } else
    fDt = 0.0;
  if ( CWsignal->shift ) {
    ASSERT( CWsignal->shift->deltaT != 0.0, stat,
	    SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
    shiftDt = output->deltaT / CWsignal->shift->deltaT;
    ASSERT( shiftDt != 0.0, stat,
	    SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
  } else
    shiftDt = 0.0;
  if ( transfer ) {
    ASSERT( detector->transfer->deltaF != 0.0, stat,
	    SIMULATECOHERENTGWH_EBAD, SIMULATECOHERENTGWH_MSGEBAD );
    fFac = 1.0 / detector->transfer->deltaF;
    phiFac = fFac / ( LAL_TWOPI*CWsignal->phi->deltaT );
    f0 = detector->transfer->f0/detector->transfer->deltaF;
  }
  heteroFac = LAL_TWOPI*output->f0*output->deltaT;
  phi0 = (REAL8)( output->epoch.gpsSeconds -
		  detector->heterodyneEpoch.gpsSeconds );
  phi0 += 0.000000001*(REAL8)( output->epoch.gpsNanoSeconds -
			       detector->heterodyneEpoch.gpsNanoSeconds );
  phi0 *= LAL_TWOPI*output->f0;
  if ( phi0 > 1.0/LAL_REAL8_EPS ) {
    LALWarning( stat, "REAL8 arithmetic is not sufficient to maintain"
		" heterodyne phase to within a radian." );
  }

  /* Check units on input, and set units on output. */
  {
    BOOLEAN unitsOK;
    LALUnitPair pair;

    pair.unitOne = &(CWsignal->f->sampleUnits);
    pair.unitTwo = &lalHertzUnit;
    TRY( LALUnitCompare( stat->statusPtr, &unitsOK, &pair ), stat );
    ASSERT( unitsOK, stat, SIMULATECOHERENTGWH_EUNIT,
	    SIMULATECOHERENTGWH_MSGEUNIT );
    pair.unitOne = &(CWsignal->phi->sampleUnits);
    pair.unitTwo = &lalDimensionlessUnit;
    TRY( LALUnitCompare( stat->statusPtr, &unitsOK, &pair ), stat );
    ASSERT( unitsOK, stat, SIMULATECOHERENTGWH_EUNIT,
	    SIMULATECOHERENTGWH_MSGEUNIT );
    if( CWsignal->shift ) {
      pair.unitOne = &(CWsignal->shift->sampleUnits);
      TRY( LALUnitCompare( stat->statusPtr, &unitsOK, &pair ), stat );
      ASSERT( unitsOK, stat, SIMULATECOHERENTGWH_EUNIT,
	      SIMULATECOHERENTGWH_MSGEUNIT );
    }
    if ( transfer ) {
      pair.unitOne = &(CWsignal->a->sampleUnits);
      pair.unitTwo = &(detector->transfer->sampleUnits);
      TRY( LALUnitMultiply( stat->statusPtr, &(output->sampleUnits),
			    &pair ), stat );
    } else
      output->sampleUnits = CWsignal->a->sampleUnits;
    snprintf( output->name, LALNameLength, "response to %s",
		 CWsignal->a->name );
  }

  /* Define temporary variables to access the data of CWsignal->a,
     CWsignal->f, and CWsignal->phi. */
  aData = CWsignal->a->data->data;
  phiData = CWsignal->phi->data->data;
  outData = output->data->data;
  if ( CWsignal->f )
    fData = CWsignal->f->data->data;
  else
    fData = NULL;
  if ( CWsignal->shift )
    shiftData = CWsignal->shift->data->data;
  else
    shiftData = NULL;

  /* Convert source position to equatorial coordinates, if
     required. */
  if ( detector->site ) {
    source = CWsignal->position;
    if ( source.system != COORDINATESYSTEM_EQUATORIAL ) {
      ConvertSkyParams params; /* parameters for conversion */
      EarthPosition location;  /* location of detector */
      params.gpsTime = &( output->epoch );
      params.system = COORDINATESYSTEM_EQUATORIAL;
      if ( source.system == COORDINATESYSTEM_HORIZON ) {
	params.zenith = &( location.geodetic );
	location.x = detector->site->location[0];
	location.y = detector->site->location[1];
	location.z = detector->site->location[2];
	TRY( LALGeocentricToGeodetic( stat->statusPtr, &location ),
	     stat );
      }
      TRY( LALConvertSkyCoordinates( stat->statusPtr, &source,
				     &source, &params ), stat );
    }
  }

  /* Generate the table of propagation delays.
  dtDelayBy2 = (UINT4)( 38924.9/sqrt( output->f0 +
				      1.0/output->deltaT ) ); */
  delayDt = output->deltaT/( 2.0*dtDelayBy2 );
  nMax = (UINT4)( output->data->length*delayDt ) + 3;
  TRY( LALDCreateVector( stat->statusPtr, &delay, nMax ), stat );
  delayData = delay->data;

  /* Compute delay from solar system barycentre. */
  if ( detector->site && detector->ephemerides ) {
    LIGOTimeGPS gpsTime;   /* detector time when we compute delay */
    EarthState state;      /* Earth position info at that time */
    BarycenterInput input; /* input structure to LALBarycenter() */
    EmissionTime emit;     /* output structure from LALBarycenter() */

    /* Arrange nested pointers, and set initial values. */
    gpsTime = input.tgps = output->epoch;
    gpsTime.gpsSeconds -= dtDelayBy2;
    input.tgps.gpsSeconds -= dtDelayBy2;
    input.site = *(detector->site);
    for ( i = 0; i < 3; i++ )
      input.site.location[i] /= LAL_C_SI;
    input.alpha = source.longitude;
    input.delta = source.latitude;
    input.dInv = 0.0;
    delayMin = delayMax = 1.1*LAL_AU_SI/( LAL_C_SI*output->deltaT );
    delayMax *= -1;

    /* Compute table. */
    for ( i = 0; i < nMax; i++ ) {
      REAL8 tDelay; /* propagation time */
      LALBarycenterEarth( stat->statusPtr, &state, &gpsTime,
			  detector->ephemerides );
      BEGINFAIL( stat )
	TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
      ENDFAIL( stat );
      LALBarycenter( stat->statusPtr, &emit, &input, &state );
      BEGINFAIL( stat )
	TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
      ENDFAIL( stat );
      delayData[i] = tDelay = emit.deltaT/output->deltaT;
      if ( tDelay < delayMin )
	delayMin = tDelay;
      if ( tDelay > delayMax )
	delayMax = tDelay;
      gpsTime.gpsSeconds += 2*dtDelayBy2;
      input.tgps.gpsSeconds += 2*dtDelayBy2;
    }
  }

  /* Compute delay from Earth centre. */
  else if ( detector->site ) {
    LIGOTimeGPS gpsTime;     /* detector time when we compute delay */

    LALInfo( stat, "Ephemeris field absent; computing propagation"
	     " delays from Earth centre" );

    /* Arrange nested pointers, and set initial values. */
    gpsTime = output->epoch;
    gpsTime.gpsSeconds -= dtDelayBy2;
    delayMin = delayMax = LAL_REARTH_SI / ( LAL_C_SI*output->deltaT );
    delayMin *= -1;

    /* Compute table. */
    for ( i = 0; i < nMax; i++ ) {
      REAL8 tDelay; /* propagation time */
      tDelay = XLALTimeDelayFromEarthCenter( detector->site->location, source.longitude, source.latitude, &gpsTime );
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
      gpsTime.gpsSeconds += 2*dtDelayBy2;
    }
  }

  /* No information from which to compute delays. */
  else {
    LALInfo( stat, "Detector site absent; simulating hplus with no"
	     " propagation delays" );
    memset( delayData, 0, nMax*sizeof(REAL8) );
    delayMin = delayMax = 0.0;
  }

  /* Generate the table of polarization response functions. */
  polDt = output->deltaT/( 2.0*dtPolBy2 );
  nMax = (UINT4)( output->data->length*polDt ) + 3;
  memset( &polResponse, 0, sizeof( LALDetAMResponseSeries ) );
  polResponse.pPlus = (REAL4TimeSeries *)
    LALMalloc( sizeof(REAL4TimeSeries) );
  polResponse.pCross = (REAL4TimeSeries *)
    LALMalloc( sizeof(REAL4TimeSeries) );
  polResponse.pScalar = (REAL4TimeSeries *)
    LALMalloc( sizeof(REAL4TimeSeries) );
  if ( !polResponse.pPlus || !polResponse.pCross ||
       !polResponse.pScalar ) {
    if ( polResponse.pPlus )
      LALFree( polResponse.pPlus );
    if ( polResponse.pCross )
      LALFree( polResponse.pCross );
    if ( polResponse.pScalar )
      LALFree( polResponse.pScalar );
    TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
    ABORT( stat, SIMULATECOHERENTGWH_EMEM,
	   SIMULATECOHERENTGWH_MSGEMEM );
  }
  memset( polResponse.pPlus, 0, sizeof(REAL4TimeSeries) );
  memset( polResponse.pCross, 0, sizeof(REAL4TimeSeries) );
  memset( polResponse.pScalar, 0, sizeof(REAL4TimeSeries) );
  LALSCreateVector( stat->statusPtr, &( polResponse.pPlus->data ),
		    nMax );
  BEGINFAIL( stat ) {
    LALFree( polResponse.pPlus );
    LALFree( polResponse.pCross );
    LALFree( polResponse.pScalar );
    TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
  } ENDFAIL( stat );
  LALSCreateVector( stat->statusPtr, &( polResponse.pCross->data ),
		    nMax );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr,
			    &( polResponse.pPlus->data ) ), stat );
    LALFree( polResponse.pPlus );
    LALFree( polResponse.pCross );
    LALFree( polResponse.pScalar );
    TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
  } ENDFAIL( stat );
  LALSCreateVector( stat->statusPtr, &( polResponse.pScalar->data ),
		    nMax );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr,
			    &( polResponse.pPlus->data ) ), stat );
    TRY( LALSDestroyVector( stat->statusPtr,
			    &( polResponse.pCross->data ) ), stat );
    LALFree( polResponse.pPlus );
    LALFree( polResponse.pCross );
    LALFree( polResponse.pScalar );
    TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
  } ENDFAIL( stat );
  plusData = polResponse.pPlus->data->data;
  crossData = polResponse.pCross->data->data;
  if ( detector->site ) {
    LALSource polSource;     /* position and polarization angle */
    LALDetAndSource input;            /* response input structure */
    LALTimeIntervalAndNSample params; /* response parameter structure */

    /* Arrange nested pointers, and set initial values. */
    polSource.equatorialCoords = source;
    polSource.orientation = (REAL8)( CWsignal->psi );
    input.pSource = &polSource;
    input.pDetector = detector->site;
    params.epoch = output->epoch;
    params.epoch.gpsSeconds -= dtPolBy2;
    params.deltaT = 2.0*dtPolBy2;
    params.nSample = nMax;

    /* Compute table of responses. */
    LALComputeDetAMResponseSeries( stat->statusPtr, &polResponse,
				   &input, &params );
    BEGINFAIL( stat ) {
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pPlus->data ) ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pCross->data ) ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pScalar->data ) ), stat );
      LALFree( polResponse.pPlus );
      LALFree( polResponse.pCross );
      LALFree( polResponse.pScalar );
      TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
    } ENDFAIL( stat );
  } else {
    /* No detector site, so just simulate response to hplus. */
    for ( i = 0; i < nMax; i++ ) {
      plusData[i] = 1.0;
      crossData[i] = 0.0;
    }
  }
  /* Free memory for the unused scalar mode. */
  TRY( LALSDestroyVector( stat->statusPtr,
			  &( polResponse.pScalar->data ) ), stat );
  LALFree( polResponse.pScalar );

  /* Decompose the transfer function into an amplitude and phase
     response. */
  if ( transfer ) {
    nMax = detector->transfer->data->length;
    LALSCreateVector( stat->statusPtr, &phiTemp, nMax );
    BEGINFAIL( stat ) {
      TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pPlus->data ) ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pCross->data ) ), stat );
      LALFree( polResponse.pPlus );
      LALFree( polResponse.pCross );
    } ENDFAIL( stat );
    LALCVectorAngle( stat->statusPtr, phiTemp,
		     detector->transfer->data );
    BEGINFAIL( stat ) {
      TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
      TRY( LALSDestroyVector( stat->statusPtr, &phiTemp ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pPlus->data ) ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pCross->data ) ), stat );
      LALFree( polResponse.pPlus );
      LALFree( polResponse.pCross );
    } ENDFAIL( stat );
    LALSCreateVector( stat->statusPtr, &phiTransfer, nMax );
    BEGINFAIL( stat ) {
      TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
      TRY( LALSDestroyVector( stat->statusPtr, &phiTemp ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pPlus->data ) ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pCross->data ) ), stat );
      LALFree( polResponse.pPlus );
      LALFree( polResponse.pCross );
    } ENDFAIL( stat );
    LALUnwrapREAL4Angle( stat->statusPtr, phiTransfer, phiTemp );
    BEGINFAIL( stat ) {
      TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
      TRY( LALSDestroyVector( stat->statusPtr, &phiTemp ), stat );
      TRY( LALSDestroyVector( stat->statusPtr, &phiTransfer ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pPlus->data ) ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pCross->data ) ), stat );
      LALFree( polResponse.pPlus );
      LALFree( polResponse.pCross );
    } ENDFAIL( stat );
    TRY( LALSDestroyVector( stat->statusPtr, &phiTemp ), stat );
    LALSCreateVector( stat->statusPtr, &aTransfer, nMax );
    BEGINFAIL( stat ) {
      TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
      TRY( LALSDestroyVector( stat->statusPtr, &phiTransfer ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pPlus->data ) ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pCross->data ) ), stat );
      LALFree( polResponse.pPlus );
      LALFree( polResponse.pCross );
    } ENDFAIL( stat );
    LALCVectorAbs( stat->statusPtr, aTransfer,
		   detector->transfer->data );
    BEGINFAIL( stat ) {
      TRY( LALDDestroyVector( stat->statusPtr, &delay ), stat );
      TRY( LALSDestroyVector( stat->statusPtr, &phiTransfer ), stat );
      TRY( LALSDestroyVector( stat->statusPtr, &aTransfer ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pPlus->data ) ), stat );
      TRY( LALSDestroyVector( stat->statusPtr,
			      &( polResponse.pCross->data ) ), stat );
      LALFree( polResponse.pPlus );
      LALFree( polResponse.pCross );
    } ENDFAIL( stat );
    phiTransData = phiTransfer->data;
    aTransData = aTransfer->data;
  }

  /* Compute offsets for interpolating the signal, delay, and
     response functions. */
  aOff = ( output->epoch.gpsSeconds -
	   CWsignal->a->epoch.gpsSeconds ) /
    CWsignal->a->deltaT;
  aOff += ( output->epoch.gpsNanoSeconds -
	    CWsignal->a->epoch.gpsNanoSeconds ) * 1.0e-9 /
    CWsignal->a->deltaT;
  phiOff = ( output->epoch.gpsSeconds -
	   CWsignal->phi->epoch.gpsSeconds ) /
    CWsignal->phi->deltaT;
  phiOff += ( output->epoch.gpsNanoSeconds -
	      CWsignal->phi->epoch.gpsNanoSeconds ) * 1.0e-9 /
    CWsignal->phi->deltaT;
  if ( CWsignal->f ) {
    fOff = ( output->epoch.gpsSeconds -
	     CWsignal->f->epoch.gpsSeconds ) /
      CWsignal->f->deltaT;
    fOff += ( output->epoch.gpsNanoSeconds -
	      CWsignal->f->epoch.gpsNanoSeconds ) * 1.0e-9 /
      CWsignal->f->deltaT;
  } else
    fOff = 0.0;
  if ( CWsignal->shift ) {
    shiftOff = ( output->epoch.gpsSeconds -
		 CWsignal->shift->epoch.gpsSeconds ) /
      CWsignal->shift->deltaT;
    shiftOff += ( output->epoch.gpsNanoSeconds -
		  CWsignal->shift->epoch.gpsNanoSeconds ) * 1.0e-9 /
      CWsignal->shift->deltaT;
  } else
    shiftOff = 0.0;
  polOff = 0.5;
  delayOff = 0.5;

  /* Compute initial value of i, ensuring that we will never index
     CWsignal->a or CWsignal->phi below their range. */
  i = 0;
  if ( aOff + ( i + delayMin )*aDt < 0.0 ) {
    INT4 j = (INT4)floor( -aOff/aDt - delayMax );
    if ( i < j )
      i = j;
    while ( ( i < (INT4)( output->data->length ) ) &&
	    ( aOff + TCENTRE( i )*aDt < 0.0 ) )
      i++;
  }
  if ( phiOff + ( i + delayMin )*phiDt < 0.0 ) {
    INT4 j = (INT4)( -phiOff/phiDt - delayMax );
    if ( i < j )
      i = j;
    while ( ( i < (INT4)( output->data->length ) ) &&
	    ( phiOff + TCENTRE( i )*phiDt < 0.0 ) )
      i++;
  }
  if ( i >= (INT4)( output->data->length ) ) {
    LALWarning( stat, "Signal starts after the end of the output"
		" time series." );
    i = (INT4)( output->data->length );
  }

  /* Compute final value of i, ensuring that we will never index
     CWsignal->a or CWsignal->phi above their range. */
  n = output->data->length - 1;
  nMax = CWsignal->a->data->length - 1;
  if ( aOff + ( n + delayMax )*aDt > nMax ) {
    INT4 j = (INT4)( ( nMax - aOff )/aDt - delayMin + 1.0 );
    if ( n > j )
      n = j;
    while ( ( n >= 0 ) &&
	    ( (INT4)floor(aOff + TCENTRE( n )*aDt) > nMax ) )
      n--;
  }
  nMax = CWsignal->phi->data->length - 1;
  if ( phiOff + ( n + delayMax )*phiDt > nMax ) {
    INT4 j = (INT4)( ( nMax - phiOff )/phiDt - delayMin + 1.0 );
    if ( n > j )
      n = j;
    while ( ( n >= 0 ) &&
	    ( (INT4)floor(phiOff + TCENTRE( n )*phiDt) > nMax ) )
      n--;
  }
  if ( n < 0 ) {
    LALWarning( stat, "CWsignal ends before the start of the output"
		" time series." );
    n = -1;
  }
  
  /* Compute the values of i for which CWsignal->f is given. */
  if ( CWsignal->f ) {
    fInit = i;
    if ( fOff + ( fInit + delayMin )*fDt < 0.0 ) {
      INT4 j = (INT4)floor( -fOff/fDt - delayMax );
      if ( fInit < j )
	fInit = j;
      while ( ( fInit <= n ) &&
	      ( fOff + TCENTRE( fInit )*fDt < 0.0 ) )
	fInit++;
    }
    fFinal = n;
    nMax = CWsignal->f->data->length - 1;
    if ( fOff + ( fFinal + delayMax )*fDt > nMax ) {
      INT4 j = (INT4)( ( nMax - fOff )/fDt - delayMin + 1.0 );
      if ( fFinal > j )
	fFinal = j;
      while ( ( fFinal >= i ) &&
	      ( (INT4)floor(fOff + TCENTRE( fFinal )*fDt) > nMax ) )
	fFinal--;
    }
  } else {
    fInit = n + 1;
    fFinal = i - 1;
  }

  /* Compute the values of i for which CWsignal->shift is given. */
  if ( CWsignal->shift ) {
    shiftInit = i;
    if ( shiftOff + ( shiftInit + delayMin )*shiftDt < 0.0 ) {
      INT4 j = (INT4)floor( -shiftOff/shiftDt - delayMax );
      if ( shiftInit < j )
	shiftInit = j;
      while ( ( shiftInit <= n ) &&
	      ( shiftOff + TCENTRE( shiftInit )*shiftDt < 0.0 ) )
	shiftInit++;
    }
    shiftFinal = n;
    nMax = CWsignal->shift->data->length - 1;
    if ( shiftOff + ( shiftFinal + delayMax )*shiftDt > nMax ) {
      INT4 j = (INT4)( ( nMax - shiftOff )/shiftDt - delayMin + 1.0 );
      if ( shiftFinal > j )
	shiftFinal = j;
      while ( ( shiftFinal >= i ) &&
	      ( (INT4)floor(shiftOff + TCENTRE( shiftFinal )*shiftDt) > nMax ) )
	shiftFinal--;
    }
  } else {
    shiftInit = n + 1;
    shiftFinal = i - 1;
  }

  /* Set output to zero where the CWsignal is not defined. */
  if ( i > 0 )
    memset( output->data->data, 0, i*sizeof(REAL4) );
  if ( ( nMax = output->data->length - n - 1 ) > 0 )
    memset( output->data->data + n + 1, 0, nMax*sizeof(REAL4) );

  /* Keep track of the frequency range of the transfer function, so
     that we don't try to interpolate it out of its range. */
  if ( transfer )
    nMax = detector->transfer->data->length - 1;

  /* Start computing responses. */
  for ( ; i <= n; i++ ) {
    REAL8 iCentre = TCENTRE( i );  /* value of i + propagation delays */
    REAL8 x;                /* interpolation point in arrays */
    INT4 j;                 /* array index preceding x */
    REAL8 frac;             /* value of x - j */
    REAL4 a1, a2;           /* current signal amplitudes */
    REAL8 phi = 0.0;        /* current signal phase */
    REAL4 f = 0.0;          /* current signal frequency */
    REAL4 shift = 0.0;      /* current signal polarization shift */
    REAL4 aTrans, phiTrans; /* current values of the transfer function */
    REAL4 oPlus, oCross;    /* current output amplitudes */
#if USE_SINCOSP
    REAL8 sp, cp, ss, cs;   /* sine and cosine of shift and phase */
#endif

    /* Interpolate the signal amplitude. */
    x = aOff + iCentre*aDt;
    j = (INT4)floor( x );
    frac = (REAL8)( x - j );
    j *= 2;
    
    /* Handle special case where output lands on final sample - no interpolation */
    if(i==n){
      a1=aData[j];
      a2=aData[j+1];
    }
    else
    {
      a1 = frac*aData[j+2] + ( 1.0 - frac )*aData[j];
      a2 = frac*aData[j+3] + ( 1.0 - frac )*aData[j+1];
    }
    
    /* Interpolate the polarization shift. */
    if ( ( i < shiftInit ) || ( i > shiftFinal ) )
      shift = 0.0;
    else {
      x = shiftOff + iCentre*shiftDt;
      j = (INT4)floor( x );
      frac = (REAL8)( x - j );
      if(i==n) shift=shiftData[j];
      else     shift = frac*shiftData[j+1] + ( 1.0 - frac )*shiftData[j];
    }

    /* Interpolate the signal phase, and apply any heterodyning. */
    x = phiOff + iCentre*phiDt;
    j = (INT4)floor( x );
    frac = (REAL8)( x - j );
    phi = frac*phiData[j+1] + ( 1.0 - frac )*phiData[j];
    phi -= heteroFac*i + phi0;

    /* Compute the frequency and apply the transfer function. */
    if ( transfer ) {
      if ( ( i < fInit ) || ( i > fFinal ) ) {
	f = ( phiData[j+1] - phiData[j] )*phiFac;
	pFlag = 1;
      } else {
	x = fOff + iCentre*fDt;
	j = (INT4)floor( x );
	frac = (REAL8)( x - j );
	if(i==n) f=fData[j];
	else     f = frac*fData[j+1] + ( 1.0 - frac )*fData[j];
	f *= fFac;
      }
      x = f - f0;
      if ( ( x < 0.0 ) || ( x > nMax ) ) {
	aTrans = 0.0;
	phiTrans = 0.0;
	fFlag = 1;
      } else  {
	j = (INT4)floor( x );
	frac = (REAL8)( x - j );
	if(i==n)
	{
	  aTrans=aTransData[j];
	  phiTrans=phiTransData[j];
	} else
	{
	  aTrans = frac*aTransData[j+1] + ( 1.0 - frac )*aTransData[j];
	  phiTrans = frac*phiTransData[j+1] + ( 1.0 - frac )*phiTransData[j];
	}
      }
      a1 *= aTrans;
      a2 *= aTrans;
      phi += phiTrans;
    }

    /* Compute components of output. */
#if USE_SINCOSP
    sincosp( shift, &ss, &cs );
    sincosp( phi, &sp, &cp );
    oPlus  = a1*cs*cp - a2*ss*sp;
    oCross = a1*ss*cp + a2*cs*sp;
#else
    oPlus = a1*cos( shift )*cos( phi ) - a2*sin( shift )*sin( phi );
    oCross = a1*sin( shift )*cos( phi ) + a2*cos( shift )*sin( phi );
#endif

    /* Interpolate the polarization response, and compute output. */
    x = polOff + i*polDt;
    j = (INT4)floor( x );
    frac = (REAL8)( x - j );
    if(i==n)
    {
      oPlus*=plusData[j];
      oCross*=crossData[j];      
    } else {
      oPlus *= frac*plusData[j+1] + ( 1.0 - frac )*plusData[j];
      oCross *= frac*crossData[j+1] + ( 1.0 - frac )*crossData[j];
    }
 
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
  if ( transfer ) {
    TRY( LALSDestroyVector( stat->statusPtr, &phiTransfer ), stat );
    TRY( LALSDestroyVector( stat->statusPtr, &aTransfer ), stat );
  }
  TRY( LALSDestroyVector( stat->statusPtr,
			  &( polResponse.pPlus->data ) ), stat );
  TRY( LALSDestroyVector( stat->statusPtr,
			  &( polResponse.pCross->data ) ), stat );
  LALFree( polResponse.pPlus );
  LALFree( polResponse.pCross );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );

} /* LALSimulateCoherentGW() */
