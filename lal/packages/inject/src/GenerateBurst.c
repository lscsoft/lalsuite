/*
 * Copyright (C) 2007 Jolien Creighton, Patrick Brady, Saikat Ray-Majumder,
 * Xavier Siemens, Teviet Creighton, Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>
#include <gsl/gsl_randist.h>


/* FIXME:  which of these are still needed? */
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/VectorOps.h>
#include <lal/Inject.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/TimeFreqFFT.h>
#include <lal/GenerateBurst.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/RealFFT.h>


NRCSID(GENERATEBURSTC, "$Id$");


/*
 * ============================================================================
 *
 *                    Fill a time series with white noise
 *
 * ============================================================================
 */


static void gaussian_noise(REAL8TimeSeries *series, REAL8 rms, gsl_rng *rng)
{
	unsigned i;

	for(i = 0; i < series->data->length; i++)
		series->data->data[i] = gsl_ran_gaussian(rng, rms);
}


/*
 * ============================================================================
 *
 *            Construct a Band- and Time-Limited White Noise Burst
 *
 * ============================================================================
 */


/*
 * Given the Fourier transform of a real-valued function h(t), compute and
 * return the integral of the square if its derivative:
 *
 * \int \dot{h}^{2} \diff t.
 *
 * Uses Kahan's compensated summation algorithm, and does summation from
 * lowest to highest frequency assuming that high frequency components tend
 * to add more to the magnitude of the derivative.
 *
 * The normalization factors in this function assume that
 * XLALREAL8FreqTimeFFT() will be used to convert the frequency series to
 * the time domain.
 */


static double integrate_hdot_squared_dt(COMPLEX16FrequencySeries *fseries)
{
	unsigned i;
	double e = 0.0;
	double sum = 0.0;

	for(i = 0; i < fseries->data->length; i++) {
		double tmp = sum;
		/* what we want to add = f^{2} |\tilde{s}(f)|^{2} + "error
		 * from last iteration" */
		double x = (fseries->f0 + i * fseries->deltaF) * (fseries->f0 + i * fseries->deltaF) * (fseries->data->data[i].re * fseries->data->data[i].re + fseries->data->data[i].im * fseries->data->data[i].im) + e;
		/* add */
		sum += x;
		/* negative of what was actually added */
		e = tmp - sum;
		/* what didn't get added, add next time */
		e += x;
	}

	/* because we've only summed the positive frequency components */

	sum *= 2;

	/* 4 \pi^{2} \delta f */

	sum *= LAL_TWOPI * LAL_TWOPI * fseries->deltaF;

	return sum;
}


/*
 * Parameters:
 *
 * duration
 * 	duration of waveform in seconds
 * frequency
 * 	centre frequency of waveform in Hertz
 * bandwidth
 * 	bandwidth of waveform in Hertz
 * int_hdot_squared
 * 	waveform is normalized so that \int \dot{h}^{2} \diff t equals this
 * delta_t
 * 	the sample rate of the time series to construct
 * rng
 * 	a GSL random number generator to be used to produce Gaussian random
 * 	variables
 */


REAL8TimeSeries *XLALBandAndTimeLimitedWhiteNoiseBurst(REAL8 duration, REAL8 frequency, REAL8 bandwidth, REAL8 int_hdot_squared, REAL8 delta_t, gsl_rng *rng)
{
	static const char func[] = "XLALBandAndTimeLimitedWhiteNoiseBurst";
	int length;
	LIGOTimeGPS epoch;
	REAL8TimeSeries *series;
	COMPLEX16FrequencySeries *fseries;
	REAL8Window *window;
	REAL8FFTPlan *plan;
	REAL8 norm_factor;
	unsigned i;

	/* check input */

	if(duration < 0 || bandwidth < 0 || duration * bandwidth < LAL_2_PI || frequency < bandwidth / 2 || int_hdot_squared < 0 || delta_t <= 0)
		XLAL_ERROR_NULL(func, XLAL_EINVAL);

	/* length of the injection time series is 10 * duration, rounded to
	 * the nearest odd integer */

	length = (int) (10.0 * duration / delta_t / 2.0);
	length = 2 * length + 1;

	/* the middle sample is t = 0 */

	XLALGPSSetREAL8(&epoch, -(length - 1) / 2 * delta_t);

	/* allocate the time series */

	series = XLALCreateREAL8TimeSeries("BTLWNB", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	if(!series)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	/* fill with independent zero-mean unit variance Gaussian random
	 * numbers */

	gaussian_noise(series, 1, rng);

	/* apply the time-domain Gaussian window.  the window function's
	 * shape parameter is ((length - 1) / 2) * delta_t / \sigma_{t} where
	 *
	 * \sigma_{t} = \sqrt{duration^{2} / 4 - 1 / (\pi^{2} bandwidth^{2})}
	 *
	 * is the compensated time-domain window duration */

	window = XLALCreateGaussREAL8Window(series->data->length, ((series->data->length - 1) / 2) * delta_t / sqrt(duration * duration / 4.0 - 1.0 / (LAL_PI * LAL_PI * bandwidth * bandwidth)));
	if(!window) {
		XLALDestroyREAL8TimeSeries(series);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	for(i = 0; i < window->data->length; i++)
		series->data->data[i] *= window->data->data[i];
	XLALDestroyREAL8Window(window);

	/* apply an additional Hann window to taper the time series to 0 */

	window = XLALCreateHannREAL8Window(series->data->length);
	if(!window) {
		XLALDestroyREAL8TimeSeries(series);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	for(i = 0; i < window->data->length; i++)
		series->data->data[i] *= window->data->data[i];
	XLALDestroyREAL8Window(window);

	/* transform to the frequency domain */

	plan = XLALCreateForwardREAL8FFTPlan(series->data->length, 0);
	fseries = XLALCreateCOMPLEX16FrequencySeries(NULL, &epoch, 0.0, 0.0, &lalDimensionlessUnit, series->data->length / 2 + 1);
	if(!plan || !fseries) {
		XLALDestroyCOMPLEX16FrequencySeries(fseries);
		XLALDestroyREAL8FFTPlan(plan);
		XLALDestroyREAL8TimeSeries(series);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	i = XLALREAL8TimeFreqFFT(fseries, series, plan);
	XLALDestroyREAL8FFTPlan(plan);
	if(i) {
		XLALDestroyCOMPLEX16FrequencySeries(fseries);
		XLALDestroyREAL8TimeSeries(series);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/* apply the frequency-domain Gaussian window.  the window
	 * function's shape parameter is computed similarly to that of the
	 * time-domain window, with \sigma_{f} = \Delta f / 2.  the window
	 * is created with its peak on the middle sample, which we need to
	 * shift to the sample corresponding to the injection's centre
	 * frequency. */

	window = XLALCreateGaussREAL8Window(2 * fseries->data->length, fseries->data->length * fseries->deltaF / (bandwidth / 2.0));
	if(!window) {
		XLALDestroyCOMPLEX16FrequencySeries(fseries);
		XLALDestroyREAL8TimeSeries(series);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	/* FIXME:  it's possible the start index should have 0.5 added or
	 * subtracted because the peak of the Gaussian might be 1/2 bin
	 * away from where this expression considers it to be */
	XLALResizeREAL8Sequence(window->data, fseries->data->length - frequency / fseries->deltaF, fseries->data->length);
	for(i = 0; i < window->data->length; i++) {
		fseries->data->data[i].re *= window->data->data[i];
		fseries->data->data[i].im *= window->data->data[i];
	}
	XLALDestroyREAL8Window(window);

	/* normalize the waveform to achieve the desired \int \dot{h}^{2}
	 * dt */

	norm_factor = sqrt(int_hdot_squared / integrate_hdot_squared_dt(fseries));
	for(i = 0; i < fseries->data->length; i++) {
		fseries->data->data[i].re *= norm_factor;
		fseries->data->data[i].im *= norm_factor;
	}

	/* transform to the time domain */

	plan = XLALCreateReverseREAL8FFTPlan(length, 0);
	if(!plan) {
		XLALDestroyCOMPLEX16FrequencySeries(fseries);
		XLALDestroyREAL8TimeSeries(series);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	i = XLALREAL8FreqTimeFFT(series, fseries, plan);
	XLALDestroyREAL8FFTPlan(plan);
	XLALDestroyCOMPLEX16FrequencySeries(fseries);
	if(i) {
		XLALDestroyREAL8TimeSeries(series);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/* apply a Tukey window to ensure continuity at the start and end
	 * of the injection.  the window's shape parameter sets what
	 * fraction of the window is used by the tapers */

	window = XLALCreateTukeyREAL8Window(series->data->length, 0.5);
	if(!window) {
		XLALDestroyREAL8TimeSeries(series);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	for(i = 0; i < window->data->length; i++)
		series->data->data[i] *= window->data->data[i];
	XLALDestroyREAL8Window(window);

	/* done */

	return series;
}


/*
 * ============================================================================
 *
 *                   Legacy Code --- Please Update to XLAL!
 *
 * ============================================================================
 */


static void
LALGenerateBurst( 
    LALStatus          *stat, 
    CoherentGW         *output,
    SimBurstTable      *simBurst,
    BurstParamStruc    *params 
    )
{
  UINT4 n, i;          /* number of and index over samples */
  REAL8 t, dt, duration;         /* time, interval */
  REAL8 t0, tau, gtime;  /* central time, decay time, gaussian time */
  REAL8 f0/*, phi0*/;      /* initial phase and frequency */
  REAL8 twopif0;       /* 2*pi*f0 */
  /* REAL8 f; */            /* current value of frequency */
  REAL4 hpeak;         /* peak strain for burst */
  /* REAL4 df = 0.0;*/      /* maximum difference between f */
  /* REAL8 phi; */          /* current value of phase */
  REAL4 *fData;        /* pointer to frequency data */
  REAL8 *phiData;      /* pointer to phase data */
  REAL4 *aData;        /* pointer to frequency data */
  LIGOTimeGPS startTime;  /* start time of injection */
  LALTimeInterval dummyInterval;

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
  duration = (REAL8)(simBurst->dtplus + simBurst->dtminus);
  dt = params->deltaT;
  if ( ( n = (INT4) (2.0 * duration / dt) ) == 0 )
  {
    ABORT(stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM );
  }
  /* notice the factor of 2 in the definition of n confusingly makes injections 
     twice as long as the variable duration */

  /* start time of data is peak time duration */
  TRY( LALFloatToInterval( stat->statusPtr, &dummyInterval, &duration ), stat );
  TRY( LALDecrementGPS( stat->statusPtr, &startTime, 
        &(simBurst->geocent_peak_time), &dummyInterval), stat);

  /* Generic burst parameters */
  hpeak = simBurst->hpeak;
  tau = (REAL8)simBurst->tau;
  f0 = (REAL8)simBurst->freq;
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
  output->position.longitude = simBurst->longitude;
  output->position.latitude = simBurst->latitude;
  output->position.system = params->system;
  output->psi = simBurst->polarization;
  output->a->epoch = output->f->epoch = output->phi->epoch = startTime;
  output->a->deltaT = params->deltaT;
  output->f->deltaT = output->phi->deltaT = params->deltaT;
  output->a->sampleUnits = lalStrainUnit;
  output->f->sampleUnits = lalHertzUnit;
  output->phi->sampleUnits = lalDimensionlessUnit;
  LALSnprintf( output->a->name, LALNameLength, "Burst amplitudes" );
  LALSnprintf( output->f->name, LALNameLength, "Burst frequency" );
  LALSnprintf( output->phi->name, LALNameLength, "Burst phase" );

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
    in.length = n;
    in.vectorLength = 2;
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
  if( !( strcmp( simBurst->waveform, "SineGaussian" ) ) )
  {
    /* find the peak time as a REAL8 relative to start of segment */
    TRY( LALDeltaGPS( stat->statusPtr, &dummyInterval, 
          &(simBurst->geocent_peak_time), &startTime ), stat );
    TRY( LALIntervalToFloat( stat->statusPtr, &t0, 
          &dummyInterval ), stat );

    /* construct the signal */
    for ( i = 0; i < n; i++ ) {
      t = i*dt;
      gtime = (t-t0)/tau;
      *(fData++) = f0;
      *(phiData++) = twopif0 * (t-t0);
      *(aData++) = hpeak * exp( - gtime * gtime );
      *(aData++) = 0.0;
    }
  }
  else if ( !( strcmp( simBurst->waveform, "Gaussian" ) ) )
  {
    /* find the peak time as a REAL8 relative to start of segment */
    TRY( LALDeltaGPS( stat->statusPtr, &dummyInterval, 
          &(simBurst->geocent_peak_time), &startTime ), stat );
    TRY( LALIntervalToFloat( stat->statusPtr, &t0, 
          &dummyInterval ), stat );

    /* construct the signal */
    for ( i = 0; i < n; i++ ) {
      t = i*dt;
      gtime = (t-t0)/tau;
      *(fData++) = 0.0;
      *(phiData++) = 0.0;
      *(aData++) = hpeak * exp( - gtime * gtime );
      *(aData++) = 0.0;
    }
  }
  else if ( !( strcmp( simBurst->waveform, "Ringdown" ) ) )
  {
    TRY( LALDeltaGPS( stat->statusPtr, &dummyInterval, 
          &(simBurst->geocent_peak_time), &startTime ), stat );
    TRY( LALIntervalToFloat( stat->statusPtr, &t0, 
          &dummyInterval ), stat );
    for ( i = 0; i < n; i++ )
    {
      t = i * dt;
      gtime = ( t - t0 ) / tau;
      *fData++   = f0;
      *phiData++ = twopif0 * ( t - t0 );
      if ( gtime > 0 )
        *aData++ = hpeak * exp( - gtime );
      else
        *aData++ = 0;
      *aData++   = 0;
    }
  }
  else if ( !( strcmp( simBurst->waveform, "Ringup" ) ) )
  {
    TRY( LALDeltaGPS( stat->statusPtr, &dummyInterval, 
          &(simBurst->geocent_peak_time), &startTime ), stat );
    TRY( LALIntervalToFloat( stat->statusPtr, &t0, 
          &dummyInterval ), stat );
    for ( i = 0; i < n; i++ )
    {
      t = i * dt;
      gtime = ( t - t0 ) / tau;
      *fData++   = f0;
      *phiData++ = twopif0 * ( t - t0 );
      if ( gtime < 0 )
        *aData++ = hpeak * exp( gtime );
      else
        *aData++ = 0;
      *aData++   = 0;
    }
  }

  else if ( !( strcmp( simBurst->waveform, "StringCusp" ) ) )
  {
    TRY( LALDeltaGPS( stat->statusPtr, &dummyInterval, 
          &(simBurst->geocent_peak_time), &startTime ), stat );
    TRY( LALIntervalToFloat( stat->statusPtr, &t0, 
          &dummyInterval ), stat );

    /* I use the simburst table as follows: The duration is still
       dtplus+dtminus; hpeak is the amplitude of the cusp, not the
       value of the strain at the peak, t0 is the central time of the
       cusp, which introduces a phase into the waveform. The low
       frequency cutoff will be fixed at 1Hz there's nothing special
       about 1Hz except that it low compared to the ferquecny at which
       we should be high-passing the data; the high frequency cutoff
       is given by f0 */
    {
      REAL4Vector *vector = NULL;
      COMPLEX8Vector *vtilde = NULL;
      REAL4 dfreq=1/(2*duration); /* the factor of two here is becaus the length of injections 
				     is actually twice the value of the variable duration */
      RealFFTPlan *rplan=NULL;
      REAL4 flow=1;

      /* create vector that will hold frequency domain template */
      TRY( LALCCreateVector( stat->statusPtr, &vtilde, n / 2 + 1 ), stat );

      for (i=0; i < vtilde->length-1; i++)
	{
	  REAL4 freq=i*dfreq;
          /* Set the FD template */
	  vtilde->data[i].re = hpeak *  pow((sqrt(1+pow(flow,2)*pow(freq,-2))),-8) * pow(freq,-4.0/3.0); 

	  if(freq>=f0)
	    {
	      vtilde->data[i].re *= exp(1-freq/f0); 
	    }

	  vtilde->data[i].im = vtilde->data[i].re * sin(-LAL_TWOPI*freq*duration);
	  vtilde->data[i].re = vtilde->data[i].re * cos(-LAL_TWOPI*freq*duration);
	}

      /* set dc to zero */
      vtilde->data[0].re = 0;
      vtilde->data[0].im = 0;
      /* set nyquist to zero */
      vtilde->data[vtilde->length - 1].re = 0;
      vtilde->data[vtilde->length - 1].im = 0;

      /* Create vector to store h(t) */
      TRY( LALSCreateVector( stat->statusPtr, &vector, n ), stat );

      /* Create fft plan */
      TRY( LALCreateReverseRealFFTPlan( stat->statusPtr, &rplan, n, 0 ), stat );
      /* Reverse FFT */
      TRY( LALReverseRealFFT( stat->statusPtr, vector, vtilde,  rplan), stat );

      /* multiply times dfreq to make sure units are correct */
      for ( i = 0 ; i < vector->length; i++ )
	vector->data[i] *= dfreq;
 
      /* make sure injection starts precisely at 0 */
      for ( i = 0 ; i < vector->length; i++ )
	vector->data[i] -= vector->data[0]; 

      for ( i = 0; i < n; i++ )
	{
	  *fData++   = 0.0;
	  *phiData++ = 0.0;
	  *aData++ = vector->data[i];
	  *aData++ = 0;
	}

      /* free the data */
      TRY( LALSDestroyVector( stat->statusPtr, &vector ), stat );
      TRY( LALCDestroyVector( stat->statusPtr, &vtilde ), stat );
      TRY( LALDestroyRealFFTPlan( stat->statusPtr, &rplan ), stat );
    }
  }
  else if ( !( strcmp( simBurst->waveform, "warren" ) ) )
  {
    /* set everything to 0 */
    for ( i = 0; i < n; i++ ) {
      *(fData++) = 0.0;
      *(phiData++) = 0.0;
      *(aData++) = 0.0;
      *(aData++) = 0.0;
    }
  }
  else
  {
    ABORT( stat, GENERATEBURSTH_ETYP, GENERATEBURSTH_MSGETYP );
  }

  /* Set output field and return. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}



void
LALBurstInjectSignals( 
    LALStatus               *stat, 
    REAL4TimeSeries         *series, 
    SimBurstTable           *injections,
    COMPLEX8FrequencySeries *resp,
    INT4                     calType
    )
{
  UINT4              k;
  INT4               injStartTime;
  INT4               injStopTime;
  DetectorResponse   detector;
  COMPLEX8Vector    *unity = NULL;
  CoherentGW         waveform;
  BurstParamStruc    burstParam;
  REAL4TimeSeries    signal;
  SimBurstTable     *simBurst=NULL;
  LALDetector       *tmpDetector=NULL /*,*nullDetector=NULL*/;
  COMPLEX8FrequencySeries    *transfer = NULL;

  INITSTATUS( stat, "LALBurstInjectSignals", GENERATEBURSTC );
  ATTATCHSTATUSPTR( stat );

  /* set up start and end of injection zone TODO: fix this hardwired 10 */
  injStartTime = series->epoch.gpsSeconds - 10;
  injStopTime = series->epoch.gpsSeconds + 10 + (INT4)(series->data->length
      * series->deltaT);

  /* 
   *compute the transfer function 
   */

  /* allocate memory and copy the parameters describing the freq series */
  memset( &detector, 0, sizeof( DetectorResponse ) );
  transfer = (COMPLEX8FrequencySeries *)
    LALCalloc( 1, sizeof(COMPLEX8FrequencySeries) );
  if ( ! transfer ) 
  {
    ABORT( stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM );
  }
  memcpy( &(transfer->epoch), &(resp->epoch),
      sizeof(LIGOTimeGPS) );
  transfer->f0 = resp->f0;
  transfer->deltaF = resp->deltaF;

  tmpDetector = detector.site = (LALDetector *) LALMalloc( sizeof(LALDetector) );
  /* set the detector site */
  switch ( series->name[0] )
  {
    case 'H':
      *(detector.site) = lalCachedDetectors[LALDetectorIndexLHODIFF];
      LALWarning( stat, "computing waveform for Hanford." );
      break;
    case 'L':
      *(detector.site) = lalCachedDetectors[LALDetectorIndexLLODIFF];
      LALWarning( stat, "computing waveform for Livingston." );
      break;
    default:
      LALFree( detector.site );
      detector.site = NULL;
      tmpDetector = NULL;
      LALWarning( stat, "Unknown detector site, computing plus mode "
          "waveform with no time delay" );
      break;
  }

  /* set up units for the transfer function */
  {
    RAT4 negOne = { -1, 0 };
    LALUnit unit;
    LALUnitPair pair;
    pair.unitOne = &lalADCCountUnit;
    pair.unitTwo = &lalStrainUnit;
    LALUnitRaise( stat->statusPtr, &unit, pair.unitTwo, &negOne );
    CHECKSTATUSPTR( stat );
    pair.unitTwo = &unit;
    LALUnitMultiply( stat->statusPtr, &(transfer->sampleUnits),
        &pair );
    CHECKSTATUSPTR( stat );
  }

  /* invert the response function to get the transfer function */
  LALCCreateVector( stat->statusPtr, &( transfer->data ),
      resp->data->length );
  CHECKSTATUSPTR( stat );

  LALCCreateVector( stat->statusPtr, &unity, resp->data->length );
  CHECKSTATUSPTR( stat );
  for ( k = 0; k < resp->data->length; ++k ) 
  {
    unity->data[k].re = 1.0;
    unity->data[k].im = 0.0;
  }

  LALCCVectorDivide( stat->statusPtr, transfer->data, unity,
      resp->data );
  CHECKSTATUSPTR( stat );

  LALCDestroyVector( stat->statusPtr, &unity );
  CHECKSTATUSPTR( stat );

  /* Set up a time series to hold signal in ADC counts */
  signal.deltaT = series->deltaT;
  if ( ( signal.f0 = series->f0 ) != 0 )
  {
    ABORT( stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM );
  }
  signal.sampleUnits = lalADCCountUnit;

  signal.data=NULL;
  LALSCreateVector( stat->statusPtr, &(signal.data), 
      series->data->length );
  CHECKSTATUSPTR( stat );

  /* loop over list of waveforms and inject into data stream */
  for ( simBurst = injections; simBurst; simBurst = simBurst->next )
  {
    /* only do the work if the burst is in injection zone */
    if( (injStartTime - simBurst->geocent_peak_time.gpsSeconds) *
        (injStopTime - simBurst->geocent_peak_time.gpsSeconds) > 0 )
      continue;

    /* set the burt params */
    burstParam.deltaT = series->deltaT;
    if( !( strcmp( simBurst->coordinates, "HORIZON" ) ) )
    {
      burstParam.system = COORDINATESYSTEM_HORIZON;
    }
    else if ( !( strcmp( simBurst->coordinates, "ZENITH" ) ) )
    {
      /* set coordinate system for completeness */
      burstParam.system = COORDINATESYSTEM_EQUATORIAL;
      detector.site = NULL;
    }
    else if ( !( strcmp( simBurst->coordinates, "GEOGRAPHIC" ) ) )
    {
      burstParam.system = COORDINATESYSTEM_GEOGRAPHIC;
    }
    else if ( !( strcmp( simBurst->coordinates, "EQUATORIAL" ) ) )
    {
      burstParam.system = COORDINATESYSTEM_EQUATORIAL;
    }
    else if ( !( strcmp( simBurst->coordinates, "ECLIPTIC" ) ) )
    {
      burstParam.system = COORDINATESYSTEM_ECLIPTIC;
    }
    else if ( !( strcmp( simBurst->coordinates, "GALACTIC" ) ) )
    {
      burstParam.system = COORDINATESYSTEM_GALACTIC;
    }
    else
      burstParam.system = COORDINATESYSTEM_EQUATORIAL;

    /* generate the burst */
    memset( &waveform, 0, sizeof(CoherentGW) );
    LALGenerateBurst( stat->statusPtr, &waveform, simBurst, &burstParam );
    CHECKSTATUSPTR( stat );

    /* must set the epoch of signal since it's used by coherent GW */
    signal.epoch = waveform.a->epoch;
    memset( signal.data->data, 0, signal.data->length * sizeof(REAL4) );

    /* decide which way to calibrate the data; defaul to old way */
    if( calType )
      detector.transfer=NULL;
    else
      detector.transfer=transfer;
    
    /* convert this into an ADC signal */
    LALSimulateCoherentGW( stat->statusPtr, 
        &signal, &waveform, &detector );
    CHECKSTATUSPTR( stat );

    /* if calibration using RespFilt */
    if( calType == 1 )
      XLALRespFilt(&signal, transfer);

    /* inject the signal into the data channel */
    LALSSInjectTimeSeries( stat->statusPtr, series, &signal );
    CHECKSTATUSPTR( stat );

    /* free memory in coherent GW structure.  TODO:  fix this */
    LALSDestroyVectorSequence( stat->statusPtr, &( waveform.a->data ) );
    CHECKSTATUSPTR( stat );
    LALSDestroyVector( stat->statusPtr, &( waveform.f->data ) );
    CHECKSTATUSPTR( stat );
    LALDDestroyVector( stat->statusPtr, &( waveform.phi->data ) );
    CHECKSTATUSPTR( stat );
    LALFree( waveform.a );   waveform.a = NULL;
    LALFree( waveform.f );   waveform.f = NULL;
    LALFree( waveform.phi );  waveform.phi = NULL;

    /* reset the detector site information in case it changed */
    detector.site = tmpDetector;
  }

  /* destroy the signal */
  LALSDestroyVector( stat->statusPtr, &(signal.data) );
  CHECKSTATUSPTR( stat );

  LALCDestroyVector( stat->statusPtr, &( transfer->data ) );
  CHECKSTATUSPTR( stat );

  if ( detector.site ) LALFree( detector.site );
  LALFree( transfer );

  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
