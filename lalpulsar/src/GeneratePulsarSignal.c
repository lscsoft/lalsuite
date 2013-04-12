/*
 * Copyright (C) 2004, 2005 Reinhard Prix
 * Copyright (C) 2004, 2005 Greg Mendell
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

/* NOTES: */
/* 07/14/04 gam; add functions LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
/* 10/08/04 gam; fix indexing into trig lookup tables (LUTs) by having table go from -2*pi to 2*pi */
/* 10/12/04 gam; When computing fCross and fPlus need to use 2.0*psi. */
/* 09/07/05 gam; Add Dterms parameter to LALFastGeneratePulsarSFTs; use this to fill in SFT bins with fake data as per LALDemod else fill in bin with zero */

#include <math.h>
#include <gsl/gsl_math.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/AVFactories.h>
#include <lal/TimeSeries.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/SFTutils.h>
#include <lal/Window.h>

#include <lal/GeneratePulsarSignal.h>

/*----------------------------------------------------------------------*/
/* Internal helper functions */
static int XLALcheck_timestamp_bounds (const LIGOTimeGPSVector *timestamps, LIGOTimeGPS t0, LIGOTimeGPS t1);
static int XLALcheckNoiseSFTs ( const SFTVector *sfts, REAL8 f0, REAL8 f1, REAL8 deltaF );
static int XLALcorrect_phase ( SFTtype *sft, LIGOTimeGPS tHeterodyne );

/*----------------------------------------------------------------------*/

extern INT4 lalDebugLevel;

static REAL8 eps = 1.e-14;	/* maximal REAL8 roundoff-error (used for determining if some REAL8 frequency corresponds to an integer "bin-index" */

/* ----- DEFINES ----- */

/*---------- Global variables ----------*/
/* empty init-structs for the types defined in here */
static SpinOrbitCWParamStruc emptyCWParams;
static CoherentGW emptySignal;

const PulsarSignalParams empty_PulsarSignalParams;
const SFTParams empty_SFTParams;
const SFTandSignalParams empty_SFTandSignalParams;
static LALUnit emptyUnit;

/** Generate a time-series at the detector for a given pulsar.
 */
REAL4TimeSeries *
XLALGeneratePulsarSignal ( const PulsarSignalParams *params /**< input params */
                           )
{
  XLAL_CHECK_NULL ( params != NULL, XLAL_EINVAL, "Invalid NULL input 'params'\n" );

  int ret;
  /*----------------------------------------------------------------------
   *
   * First call GenerateSpinOrbitCW() to generate the source-signal
   *
   *----------------------------------------------------------------------*/
  SpinOrbitCWParamStruc sourceParams = emptyCWParams;
  sourceParams.psi = params->pulsar.psi;
  sourceParams.aPlus = params->pulsar.aPlus;
  sourceParams.aCross = params->pulsar.aCross;
  sourceParams.phi0 = params->pulsar.phi0;
  sourceParams.f0 = params->pulsar.f0;

  sourceParams.position = params->pulsar.position;
  /* set source position: make sure it's "normalized", i.e. [0<=alpha<2pi]x[-pi/2<=delta<=pi/2] */
  XLALNormalizeSkyPosition ( &(sourceParams.position.longitude), &(sourceParams.position.latitude) );

  /* if pulsar is in binary-orbit, set binary parameters */
  if (params->orbit)
    {
      /*------------------------------------------------------------ */
      /* temporary fix for comparison with Chris' code */
      /*
	TRY (LALConvertGPS2SSB (status->statusPtr, &tmpTime, params->orbit->orbitEpoch, params), status);
	sourceParams.orbitEpoch = tmpTime;
      */
      sourceParams.orbitEpoch =  params->orbit->tp;
      sourceParams.omega = params->orbit->argp;
      /* ------- here we do conversion to Teviets preferred variables -------*/
      sourceParams.rPeriNorm = params->orbit->asini*(1.0 - params->orbit->ecc);
      sourceParams.oneMinusEcc = 1.0 - params->orbit->ecc;
      sourceParams.angularSpeed = (LAL_TWOPI/params->orbit->period)*sqrt((1.0 +params->orbit->ecc)/pow((1.0 - params->orbit->ecc),3.0));
    }
  else
    {
      sourceParams.rPeriNorm = 0.0;		/* this defines an isolated pulsar */
    }

  if ( params->pulsar.refTime.gpsSeconds != 0)
    {
      sourceParams.spinEpoch = params->pulsar.refTime;   /* pulsar reference-time in SSB frame (TDB) */
    }
  else	/* if not given: use startTime converted to SSB as tRef ! */
    {
      LIGOTimeGPS tmpTime;
      ret = XLALConvertGPS2SSB ( &tmpTime, params->startTimeGPS, params );
      XLAL_CHECK_NULL ( ret == XLAL_SUCCESS, XLAL_EFUNC );
      sourceParams.spinEpoch = tmpTime;
    }

  /* sampling-timestep and length for source-parameters */
  /* in seconds; hardcoded; was 60s in makefakedata_v2,
   * but for fast binaries (e.g. SCO-X1) we need faster sampling
   * This does not seem to affect performance a lot (~4% in makefakedata),
   * but we'll nevertheless make this sampling faster for binaries and slower
   * for isolated pulsars */
  if (params->orbit)
    sourceParams.deltaT = 5;	/* for binaries */
  else
    sourceParams.deltaT = 60;	/* for isolated pulsars */

  /* start-time in SSB time */
  LIGOTimeGPS t0;
  ret = XLALConvertGPS2SSB ( &t0, params->startTimeGPS, params );
  XLAL_CHECK_NULL ( ret == XLAL_SUCCESS, XLAL_EFUNC );

  t0.gpsSeconds -= (UINT4)sourceParams.deltaT; /* start one time-step earlier to be safe */

  /* end time in SSB */
  LIGOTimeGPS t1 = params->startTimeGPS;
  XLALGPSAdd ( &t1, params->duration );
  ret = XLALConvertGPS2SSB ( &t1, t1, params );
  XLAL_CHECK_NULL ( ret == XLAL_SUCCESS, XLAL_EFUNC );

  /* get duration of source-signal */
  REAL8 SSBduration = XLALGPSDiff ( &t1, &t0 );
  SSBduration += 2.0 * sourceParams.deltaT; /* add two time-steps to be safe */

  sourceParams.epoch = t0;
  sourceParams.length = (UINT4) ceil( SSBduration / sourceParams.deltaT );

  /* we use frequency-spindowns, but GenerateSpinOrbitCW wants f_k = fkdot / (f0 * k!) */
  if ( params->pulsar.spindown )
    {
      UINT4 numSpindowns = params->pulsar.spindown->length;
      sourceParams.f = XLALCreateREAL8Vector ( numSpindowns );
      XLAL_CHECK_NULL ( sourceParams.f != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector(%d) failed.", numSpindowns );

      UINT4 kFact = 1;
      for ( UINT4 i = 0; i < numSpindowns; i++ )
	{
	  sourceParams.f->data[i] = params->pulsar.spindown->data[i] / (kFact * params->pulsar.f0);
	  kFact *= i + 2;
	} // i < numSpindowns
    } // if pulsar.spindown

  /* finally, call the function to generate the source waveform */
  CoherentGW sourceSignal = emptySignal;

  XLAL_CHECK_NULL ( XLALGenerateSpinOrbitCW ( &sourceSignal, &sourceParams ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* free spindown-vector right away, so we don't forget */
  if (sourceParams.f) {
    XLALDestroyREAL8Vector ( sourceParams.f );
  }

  /* check that sampling interval was short enough */
  XLAL_CHECK_NULL ( sourceParams.dfdt <= 2.0, XLAL_ETOL, "GenerateSpinOrbitCW() returned df*dt = %f > 2.0\n", sourceParams.dfdt );

  /*----------------------------------------------------------------------
   *
   * Now call the function to translate the source-signal into a (heterodyned)
   * signal at the detector
   *
   *----------------------------------------------------------------------*/
  /* first set up the detector-response */
  DetectorResponse detector;
  detector.transfer = params->transfer;
  detector.site = params->site;
  detector.ephemerides = params->ephemerides;

  /* *contrary* to makefakedata_v2, we use the GPS start-time of the timeseries as
   * the heterodyning epoch, in order to make sure that the global phase of the
   * SFTs it correct: this is necessary to allow correct parameter-estimation on the
   * pulsar signal-phase in heterodyned SFTs
   */
  detector.heterodyneEpoch = params->startTimeGPS;

  /* NOTE: a timeseries of length N*dT has no timestep at N*dT !! (convention) */
  UINT4 numSteps = (UINT4) ceil( params->samplingRate * params->duration );
  REAL8 dt = 1.0 / params->samplingRate;
  REAL8 fHet = params->fHeterodyne;

  /* ok, we  need to prepare the output time-series */
  REAL4TimeSeries *output = XLALCreateREAL4TimeSeries ( "", &(params->startTimeGPS), fHet, dt, &emptyUnit, numSteps );
  XLAL_CHECK_NULL ( output != NULL, XLAL_EFUNC, "XLALCreateREAL4TimeSeries() failed with xlalErrno = %d\n", xlalErrno );

  // internal interpolation parameters for LALSimulateCoherentGW()
  sourceSignal.dtDelayBy2 = params->dtDelayBy2;
  sourceSignal.dtPolBy2   = params->dtPolBy2;

  XLAL_CHECK_NULL ( XLALSimulateCoherentGW ( output, &sourceSignal, &detector ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* set 'name'-field of timeseries to contain the right "channel prefix" for the detector */
  CHAR *name = XLALGetChannelPrefix ( params->site->frDetector.name );
  XLAL_CHECK_NULL ( name != NULL, XLAL_EFUNC );
  strcpy ( output->name, name );
  XLALFree ( name );

  /*----------------------------------------------------------------------*/
  /* Free all allocated memory that is not returned */
  XLALDestroyREAL4VectorSequence ( sourceSignal.a->data );
  XLALFree ( sourceSignal.a );
  XLALDestroyREAL4TimeSeries ( sourceSignal.f );
  XLALDestroyREAL8TimeSeries ( sourceSignal.phi );

  return output;

} /* XLALGeneratePulsarSignal() */

/**
 * \deprecated Use XLALGeneratePulsarSignal() instead.
 */
void
LALGeneratePulsarSignal (LALStatus *status,		   /**< pointer to LALStatus structure */
			 REAL4TimeSeries **signalvec, 	   /**< output time-series */
			 const PulsarSignalParams *params) /**< input params */
{
  INITSTATUS(status);

  REAL4TimeSeries *output = XLALGeneratePulsarSignal ( params );
  if ( output == NULL ) {
    XLALPrintError ("XLALGeneratePulsarSignal() failed with xlalErrno = %d\n", xlalErrno );
    ABORTXLAL( status );
  }

  (*signalvec) = output;

  RETURN (status);

} /* LALGeneratePulsarSignal() */

/** Turn the given time-series into (v2-)SFTs and add noise if given.
 */
SFTVector *
XLALSignalToSFTs ( const REAL4TimeSeries *signalvec, 	/**< input time-series */
                   const SFTParams *params		/**< params for output-SFTs */
                  )
{
  XLAL_CHECK_NULL ( signalvec != NULL, XLAL_EINVAL, "Invalid NULL input 'signalvec'\n");
  XLAL_CHECK_NULL ( params != NULL, XLAL_EINVAL, "Invalid NULL input 'params'\n");

  REAL8 f0 = signalvec->f0;		/* lowest frequency */
  REAL8 dt = signalvec->deltaT;		/* timeseries timestep */
  REAL8 Band = 1.0 / (2.0 * dt);	/* NOTE: frequency-band is determined by sampling-rate! */
  REAL8 deltaF = 1.0 / params->Tsft;	/* frequency-resolution */

  int ret;
  /* if noiseSFTs are given: check they are consistent with signal! */
  if ( params->noiseSFTs )
    {
      ret = XLALcheckNoiseSFTs ( params->noiseSFTs, f0, f0 + Band, deltaF );
      XLAL_CHECK_NULL ( ret == XLAL_SUCCESS, XLAL_EFUNC );
    }

  /* make sure that number of timesamples/SFT is an integer (up to possible rounding errors) */
  REAL8 REALnumTimesteps = params->Tsft / dt;		/* this is a float!*/
  UINT4 numTimesteps = (UINT4) (REALnumTimesteps + 0.5);		/* number of time-samples in an Tsft, round to closest int */
  XLAL_CHECK_NULL ( fabs ( REALnumTimesteps - numTimesteps ) / REALnumTimesteps < eps, XLAL_ETOL,
                    "Inconsistent sampling-step (dt=%g) and Tsft=%g: must be integer multiple Tsft/dt = %g >= %g\n",
                    dt, params->Tsft, REALnumTimesteps, eps );

  /* Prepare FFT: compute plan for FFTW */
  RealFFTPlan *pfwd = XLALCreateForwardREAL4FFTPlan ( numTimesteps, 0 );
  XLAL_CHECK_NULL ( pfwd != NULL, XLAL_EFUNC, "XLALCreateForwardREAL4FFTPlan(%d,0) failed.\n", numTimesteps );

  /* get some info about time-series */
  LIGOTimeGPS tStart = signalvec->epoch;	/* start-time of time-series */

  /* get last possible start-time for an SFT */
  REAL8 duration =  (UINT4) (1.0* signalvec->data->length * dt + 0.5); /* total duration rounded to seconds */
  LIGOTimeGPS tLast = tStart;
  XLALGPSAdd( &tLast, duration - params->Tsft );
  XLAL_CHECK_NULL ( xlalErrno == XLAL_SUCCESS, XLAL_EFUNC );

  /* for simplicity we _always_ work with timestamps.
   * Therefore, we have to generate them now if none have been provided by the user. */
  LIGOTimeGPSVector *timestamps;
  if ( params->timestamps == NULL )
    {
      timestamps = XLALMakeTimestamps ( tStart, duration, params->Tsft );
      XLAL_CHECK_NULL ( timestamps != NULL, XLAL_EFUNC );
      /* see if the last timestamp is valid (can fit a full SFT in?), if not, drop it */
      LIGOTimeGPS lastTs = timestamps->data[timestamps->length-1];
      if ( XLALGPSDiff ( &lastTs, &tLast ) > 0 )	// if lastTs > tLast
	timestamps->length --;
    }
  else	/* if given, use those, and check they are valid */
    {
      timestamps = params->timestamps;
    }

  /* check that all timestamps lie within [tStart, tLast] */
  XLAL_CHECK_NULL ( XLALcheck_timestamp_bounds ( timestamps, tStart, tLast) == XLAL_SUCCESS, XLAL_EFUNC );

  UINT4 numSFTs = timestamps->length;			/* number of SFTs to produce */
  /* check that we have the right number of noise-SFTs */
  if ( params->noiseSFTs ) {
    XLAL_CHECK_NULL ( params->noiseSFTs->length == numSFTs, XLAL_EDOM, "Inconsistent number of SFTs in timestamps (%d) and noise-SFTs (%d)\n",
                      numSFTs, params->noiseSFTs->length );
  }

  /* check that if the user gave a window then the length should be correct */
  if ( params->window ) {
    XLAL_CHECK_NULL ( numTimesteps == params->window->data->length, XLAL_EDOM, "Inconsistent window-length =%d, differs from numTimesteps=%d\n",
                      params->window->data->length, numTimesteps );
  }

  /* prepare SFT-vector for return */
  UINT4 numBins = (UINT4)(numTimesteps/2) + 1;		/* number of frequency-bins per SFT */

  SFTVector *sftvect = XLALCreateSFTVector ( numSFTs, numBins );
  XLAL_CHECK_NULL ( sftvect != NULL, XLAL_EFUNC, "XLALCreateSFTVector(numSFTs=%d, numBins=%d) failed.\n", numSFTs, numBins );

  LIGOTimeGPS tPrev = tStart;	/* initialize */
  UINT4 totalIndex = 0;		/* timestep-index to start next FFT from */

  /* Assign memory to timeStretchCopy */
  REAL4Vector *timeStretchCopy = XLALCreateREAL4Vector ( numTimesteps );
  XLAL_CHECK_NULL ( timeStretchCopy != NULL, XLAL_EFUNC, "XLALCreateREAL4Vector(%d) failed.\n", numTimesteps );

  /* main loop: apply FFT the requested time-stretches */
  for (UINT4 iSFT = 0; iSFT < numSFTs; iSFT++ )
    {
      SFTtype *thisSFT = &(sftvect->data[iSFT]);	/* point to current SFT-slot */

      /* find the start-bin for this SFT in the time-series */
      REAL8 delay = XLALGPSDiff ( &(timestamps->data[iSFT]), &tPrev );

      /* round properly: picks *closest* timestep (==> "nudging") !!  */
      INT4 relIndexShift = (INT4) ( delay / signalvec->deltaT + 0.5 );
      totalIndex += relIndexShift;

      REAL4Vector timeStretch;
      timeStretch.length = numTimesteps;
      timeStretch.data = signalvec->data->data + totalIndex; /* point to the right sample-bin */
      memcpy ( timeStretchCopy->data, timeStretch.data, numTimesteps * sizeof(*timeStretch.data) );

      /* fill the header of the i'th output SFT */
      REAL8 realDelay = (REAL4)( relIndexShift * signalvec->deltaT );  /* cast to REAL4 to avoid rounding-errors*/
      LIGOTimeGPS tmpTime = tPrev;
      XLALGPSAdd ( &tmpTime, realDelay );

      strcpy ( thisSFT->name, signalvec->name );
      /* set the ACTUAL timestamp! (can be different from requested one ==> "nudging") */
      thisSFT->epoch = tmpTime;
      thisSFT->f0 = signalvec->f0;			/* minimum frequency */
      thisSFT->deltaF = 1.0 / params->Tsft;	/* frequency-spacing */

      tPrev = tmpTime;				/* prepare next loop */

      /* ok, issue at least a warning if we have "nudged" an SFT-timestamp */
      if ( lalDebugLevel > 0 )
	{
	  REAL8 diff = XLALGPSDiff ( &(timestamps->data[iSFT]), &tmpTime );
	  if (diff != 0)
	    {
	      XLALPrintError ("Warning: timestamp %d had to be 'nudged' by %e s to fit with time-series\n", iSFT, diff );
	      /* double check if magnitude of nudging seems reasonable .. */
	      XLAL_CHECK_NULL ( fabs(diff) < signalvec->deltaT, XLAL_ETOL, "Nudged by more (%g) than deltaT=%g ... this sounds wrong! (We better stop)\n",
                                fabs(diff), signalvec->deltaT );
            } // if nudging
        } /* if lalDebugLevel */

      /* Now window the current time series stretch, if necessary */
      if ( params->window )
        {
	  const float A = 1.0 / sqrt(params->window->sumofsquares / params->window->data->length);
	  for( UINT4 idatabin = 0; idatabin < timeStretchCopy->length; idatabin++ )
            {
              timeStretchCopy->data[idatabin] *= A * params->window->data->data[idatabin];
            }
        } // if window

      /* the central step: FFT the ith time-stretch into an SFT-slot */
      ret = XLALREAL4ForwardFFT ( thisSFT->data, timeStretchCopy, pfwd );
      XLAL_CHECK_NULL ( ret == XLAL_SUCCESS, XLAL_EFUNC, "XLALREAL4ForwardFFT() failed.\n");

      /* normalize DFT-data to conform to v2 ( ie. COMPLEX8FrequencySeries ) specification ==> multiply DFT by dt */
      COMPLEX8 *data = thisSFT->data->data;
      for ( UINT4 i = 0; i < numBins ; i ++ )
	{
	  data->realf_FIXME *= dt;
	  data->imagf_FIXME *= dt;
	  data ++;
	} /* for i < numBins */

      /* correct heterodyning-phase, IF NECESSARY */
      if ( ( (INT4)signalvec->f0 != signalvec->f0  ) || (signalvec->epoch.gpsNanoSeconds != 0) || (thisSFT->epoch.gpsNanoSeconds != 0) )
	{
	  /* theterodyne = signalvec->epoch!*/
	  ret = XLALcorrect_phase ( thisSFT, signalvec->epoch);
          XLAL_CHECK_NULL ( ret == XLAL_SUCCESS, XLAL_EFUNC, "XLALcorrect_phase() failed.\n");
	} /* if phase-correction necessary */

      /* Now add the noise-SFTs if given */
      if (params->noiseSFTs)
	{
	  SFTtype *thisNoiseSFT = &( params->noiseSFTs->data[iSFT] );
	  UINT4 index0n = round ( (thisSFT->f0 - thisNoiseSFT->f0) / thisSFT->deltaF );

	  data  = thisSFT->data->data;
	  COMPLEX8 *noise = &( thisNoiseSFT->data->data[index0n] );
	  for ( UINT4 j=0; j < numBins; j++ )
	    {
	      data->realf_FIXME += crealf(*noise);
	      data->imagf_FIXME += cimagf(*noise);
	      data++;
	      noise++;
	    } /* for j < numBins */

	} /* if noiseSFTs */

    } /* for iSFT < numSFTs */

  /* free stuff */
  XLALDestroyREAL4FFTPlan ( pfwd );
  XLALDestroyREAL4Vector ( timeStretchCopy );

  /* did we create timestamps ourselves? */
  if (params->timestamps == NULL) {
    XLALDestroyTimestampVector ( timestamps );	// if yes, free them
  }

  return sftvect;

} /* XLALSignalToSFTs() */

/**
 * \deprecated Use XLALSignalToSFTs() instead
 */
void
LALSignalToSFTs (LALStatus *status,		/**< pointer to LALStatus structure */
		 SFTVector **outputSFTs,	/**< [out] SFT-vector */
		 const REAL4TimeSeries *signalvec, /**< input time-series */
		 const SFTParams *params)	/**< params for output-SFTs */
{
  INITSTATUS(status);

  if ( outputSFTs == NULL ) {
    ABORT ( status, XLAL_EINVAL, "Invalid NULL input 'outputSFTs'");
  }
  if ( (*outputSFTs) != NULL ) {
    ABORT ( status, XLAL_EINVAL, "Input-pointer (*outputSFTs) must be NULL");
  }

  SFTVector *out = XLALSignalToSFTs ( signalvec, params );
  if ( out == NULL ) {
    XLALPrintError ("XLALSignalToSFTs() failed with xlalErrno = %d\n", xlalErrno );
    ABORTXLAL ( status );
  }

  (*outputSFTs) = out;

  RETURN (status);

} /* LALSignalToSFTs() */


/* 07/14/04 gam */
/**
 * Wrapper for LALComputeSky() and  LALComputeDetAMResponse() that finds the sky
 * constants and \f$F_+\f$ and \f$F_\times\f$ for use with LALFastGeneratePulsarSFTs().
 *  Uses the output of LALComputeSkyAndZeroPsiAMResponse() and the same inputs as
 * LALGeneratePulsarSignal() and LALSignalToSFTs().
 * This function used LALComputeSkyBinary() if params->pSigParams->orbit is not
 * NULL, else it uses LALComputeSky() to find the skyConsts.
 * NOTE THAT THIS FUNCTION COMPUTES \f$F_+\f$ and \f$F_x\f$ for ZERO Psi!!!
 * LALFastGeneratePulsarSFTs() used these to find \f$F_+\f$ and \f$F_x\f$ for NONZERO Psi.
 */
void
LALComputeSkyAndZeroPsiAMResponse (LALStatus *status,		/**< pointer to LALStatus structure */
                                   SkyConstAndZeroPsiAMResponse *output,	/**< output */
                                   const SFTandSignalParams *params)		/**< params */
{
  INT4 i;
  INT4 numSFTs;                      /* number of SFTs */
  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */
  CSParams *csParams   = NULL;       /* ComputeSky parameters */
  CSBParams *csbParams = NULL;       /* ComputeSkyBinary parameters */
  SkyPosition tmp;
  EarthState earth;
  EmissionTime emit;
  LALDetAMResponse response;  /* output of LALComputeDetAMResponse */
  LALDetAndSource      *das;  /* input for LALComputeDetAMResponse */
  REAL8 halfTsft;             /* half the time of one SFT */
  LIGOTimeGPS midTS;          /* midpoint time for an SFT */

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  numSFTs = params->pSFTParams->timestamps->length; /* number of SFTs */
  halfTsft = 0.5*params->pSFTParams->Tsft;          /* half the time of one SFT */

  /* setup baryinput for LALComputeSky */
  baryinput.site = *(params->pSigParams->site);
  /* account for a quirk in LALBarycenter(): -> see documentation of type BarycenterInput */
  baryinput.site.location[0] /= LAL_C_SI;
  baryinput.site.location[1] /= LAL_C_SI;
  baryinput.site.location[2] /= LAL_C_SI;
  if (params->pSigParams->pulsar.position.system != COORDINATESYSTEM_EQUATORIAL) {
      ABORT (status, GENERATEPULSARSIGNALH_EBADCOORDS, GENERATEPULSARSIGNALH_MSGEBADCOORDS);
  }
  TRY( LALNormalizeSkyPosition (status->statusPtr, &tmp, &(params->pSigParams->pulsar.position)), status);
  baryinput.alpha = tmp.longitude;
  baryinput.delta = tmp.latitude;
  baryinput.dInv = 0.e0;      /* following makefakedata_v2 */

  if (params->pSigParams->orbit) {
    /* LALComputeSkyBinary parameters */
    csbParams=(CSBParams *)LALMalloc(sizeof(CSBParams));
    csbParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
    if (params->pSigParams->pulsar.spindown) {
       csbParams->spinDwnOrder=params->pSigParams->pulsar.spindown->length;
    } else {
       csbParams->spinDwnOrder=0;
    }
    csbParams->mObsSFT=numSFTs;
    csbParams->tSFT=params->pSFTParams->Tsft;
    csbParams->tGPS=params->pSFTParams->timestamps->data;
    csbParams->skyPos[0]=params->pSigParams->pulsar.position.longitude;
    csbParams->skyPos[1]=params->pSigParams->pulsar.position.latitude;
    csbParams->OrbitalEccentricity = params->pSigParams->orbit->ecc; /* Orbital eccentricy */
    csbParams->ArgPeriapse = params->pSigParams->orbit->argp;       /* argument of periapsis (radians) */
    csbParams->TperiapseSSB = params->pSigParams->orbit->tp; /* time of periapsis passage (in SSB) */
    /* compute semi-major axis and orbital period */
    csbParams->SemiMajorAxis = params->pSigParams->orbit->asini;
    csbParams->OrbitalPeriod = params->pSigParams->orbit->period;
    csbParams->baryinput=&baryinput;
    csbParams->emit = &emit;
    csbParams->earth = &earth;
    csbParams->edat=params->pSigParams->ephemerides;

    /* Call LALComputeSkyBinary */
    TRY ( LALComputeSkyBinary (status->statusPtr, output->skyConst, 0, csbParams), status);
    LALFree(csbParams->skyPos);
    LALFree(csbParams);
  } else {
    /* LALComputeSky parameters */
    csParams=(CSParams *)LALMalloc(sizeof(CSParams));
    csParams->tGPS=params->pSFTParams->timestamps->data;
    csParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
    csParams->mObsSFT=numSFTs;
    csParams->tSFT=params->pSFTParams->Tsft;
    csParams->edat=params->pSigParams->ephemerides;
    csParams->baryinput=&baryinput;
    if (params->pSigParams->pulsar.spindown) {
       csParams->spinDwnOrder=params->pSigParams->pulsar.spindown->length;
    } else {
       csParams->spinDwnOrder=0;
    }
    csParams->skyPos[0]=params->pSigParams->pulsar.position.longitude;
    csParams->skyPos[1]=params->pSigParams->pulsar.position.latitude;
    csParams->earth = &earth;
    csParams->emit = &emit;

    /* Call LALComputeSky */
    TRY ( LALComputeSky (status->statusPtr, output->skyConst, 0, csParams), status);
    LALFree(csParams->skyPos);
    LALFree(csParams);
  }

  /* Set up das, the Detector and Source info */
  das = (LALDetAndSource *)LALMalloc(sizeof(LALDetAndSource));
  das->pSource = (LALSource *)LALMalloc(sizeof(LALSource));
  das->pDetector = params->pSigParams->site;
  das->pSource->equatorialCoords.latitude = params->pSigParams->pulsar.position.latitude;
  das->pSource->equatorialCoords.longitude = params->pSigParams->pulsar.position.longitude;
  das->pSource->orientation = 0.0;  /* NOTE THIS FUNCTION COMPUTE F_+ and F_x for ZERO Psi!!! */
  das->pSource->equatorialCoords.system = params->pSigParams->pulsar.position.system;

  /* loop that calls LALComputeDetAMResponse to find F_+ and F_x at the midpoint of each SFT for ZERO Psi */
  for(i=0; i<numSFTs; i++) {
      /* Find mid point from timestamp, half way through SFT. */
      midTS = params->pSFTParams->timestamps->data[i];
      XLALGPSAdd(&midTS, halfTsft);
      TRY ( LALComputeDetAMResponse(status->statusPtr, &response, das, &midTS), status);
      output->fPlusZeroPsi[i] = response.plus;
      output->fCrossZeroPsi[i] = response.cross;
  }
  LALFree(das->pSource);
  LALFree(das);

  DETATCHSTATUSPTR( status );
  RETURN (status);
} /* LALComputeSkyAndZeroPsiAMResponse */

/* 07/14/04 gam */
/**
 * Fast generation of Fake SFTs for a pure pulsar signal.
 * Uses the output of LALComputeSkyAndZeroPsiAMResponse and the same inputs
 * as LALGeneratePulsarSignal and LALSignalToSFTs.  The fake signal is
 * Taylor expanded to first order about the midpoint time of each SFT.
 * Analytic expressions are used to find each SFT
 */
void
LALFastGeneratePulsarSFTs (LALStatus *status,
                           SFTVector **outputSFTs,
                           const SkyConstAndZeroPsiAMResponse *input,
                           const SFTandSignalParams *params)
{
  INT4 numSFTs;                 /* number of SFTs */
  REAL4 N;                      /* N = number of time-samples that would have been used to generate SFTs directly */
  INT4 iSFT;                    /* index that gives which SFT in an SFTVector */
  INT4 SFTlen;                  /* number of frequency bins in an SFT */
  REAL8 tSFT, f0, band, f0Signal, deltaF;
  REAL4 fPlus, fCross, psi, phi0Signal;
  /* REAL4 halfAPlus, halfACross, cosPsi, sinPsi; */ /* 10/12/04 gam */
  REAL4 halfAPlus, halfACross, cos2Psi, sin2Psi;
  REAL8 realA, imagA, xSum, ySum, xTmp, yTmp; /* xSum, ySum and xTmp are the same as xSum, ySum, and x in LALDemod; yTmp is -y from LALDemod plus phi0Signal */
  REAL8 realQcc, imagQcc, realPcc, imagPcc, realTmp, imagTmp;  /* Pcc is the complex conjugate of P in LALDemod; Qcc is the complex conjugate of Q in LALDemod times exp(i*phi0Signal) */
  REAL8 kappa;   /* kappa = index of freq at midpoint of SFT which is usually not an integer */
  REAL8 real8TwoPi = (REAL8)LAL_TWOPI;
  REAL8 sin2PiKappa, oneMinusCos2PiKappa;
  SFTtype *thisSFT, *thisNoiseSFT;      /* SFT-pointers */
  SFTVector *sftvect = NULL;            /* return value. For better readability */
  INT4 j, k, k0, s, spOrder, tmpInt, index0n;
  INT4 jStart, jEnd, k1;
  BOOLEAN setToZero = 0;  /* 09/07/05 gam; flag that set whether to zero bins not within the Dterms loop */
  REAL8 smallX=0.000000001;
  /* Next are for LUT for trig calls */
  INT4 indexTrig;
  REAL8 halfResTrig = ((REAL8)params->resTrig)/2.0; /* 10/08/04 gam; fix indexing into trig lookup tables (LUTs) by having table go from -2*pi to 2*pi */
  REAL8 varTmp, dTmp, dTmp2, sinTmp, cosTmp;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* fprintf(stdout,"\n Hello from LALFastGeneratePulsarSFTs \n");
  fflush(stdout); */

  ASSERT (outputSFTs != NULL, status, GENERATEPULSARSIGNALH_ENULL, GENERATEPULSARSIGNALH_MSGENULL);
  /* ASSERT (*outputSFTs == NULL, status,  GENERATEPULSARSIGNALH_ENONULL,  GENERATEPULSARSIGNALH_MSGENONULL); */
  ASSERT (params != NULL, status, GENERATEPULSARSIGNALH_ENULL, GENERATEPULSARSIGNALH_MSGENULL);
  ASSERT (input != NULL, status, GENERATEPULSARSIGNALH_ENULL, GENERATEPULSARSIGNALH_MSGENULL);

  if ( params->pSFTParams->timestamps && params->pSFTParams->noiseSFTs) {
    ASSERT ( params->pSFTParams->timestamps->length == params->pSFTParams->noiseSFTs->length, status,
      GENERATEPULSARSIGNALH_ENUMSFTS,  GENERATEPULSARSIGNALH_MSGENUMSFTS);
  }

  /* 10/08/04 gam; fix indexing into trig lookup tables (LUTs) by having table go from -2*pi to 2*pi */
  if (params->resTrig > 0) {
     ASSERT ( fabs( ( params->trigArg[0] + ((REAL8)LAL_TWOPI) )/ ((REAL8)LAL_TWOPI) ) < ( 2.0e-6*((REAL8)LAL_TWOPI) / params->resTrig ), status,
           GENERATEPULSARSIGNALH_ELUTS,  GENERATEPULSARSIGNALH_MSGELUTS);
     ASSERT ( fabs( ( params->trigArg[params->resTrig] - ((REAL8)LAL_TWOPI) ) / ((REAL8)LAL_TWOPI)  ) < ( 2.0e-6*((REAL8)LAL_TWOPI) / params->resTrig ), status,
           GENERATEPULSARSIGNALH_ELUTS,  GENERATEPULSARSIGNALH_MSGELUTS);
  }

  /* SFT parameters */
  tSFT = params->pSFTParams->Tsft;                  /* SFT duration */
  deltaF = 1.0/tSFT;                                /* frequency resolution */
  f0 = params->pSigParams->fHeterodyne;             /* start frequency */
  k0 = (INT4)(f0*tSFT + 0.5);                       /* index of start frequency */
  band = 0.5*params->pSigParams->samplingRate;      /* frequency band */
  SFTlen = (INT4)(band*tSFT + 0.5);                 /* number of frequency-bins */
  numSFTs = params->pSFTParams->timestamps->length; /* number of SFTs */

  if ( (params->Dterms < 1) || (params->Dterms > SFTlen) ) {
     ABORT (status, GENERATEPULSARSIGNALH_EDTERMS, GENERATEPULSARSIGNALH_MSGEDTERMS);
  }

  /* prepare SFT-vector for return */
  if (*outputSFTs == NULL) {
    TRY (LALCreateSFTVector (status->statusPtr, &sftvect, numSFTs, SFTlen), status);
    setToZero = 1; /* 09/07/05 gam; allocated memory for the output SFTs, zero bins not within the Dterms loop */
  } else {
    sftvect = *outputSFTs;  /* Assume memory already allocated for SFTs */
    setToZero = 0; /* 09/07/05 gam; it's up to the user in this case to initialize bin to zero */
  }

  /* if noiseSFTs are given: check they are consistent with signal! */
  if (params->pSFTParams->noiseSFTs) {
    int ret = XLALcheckNoiseSFTs ( params->pSFTParams->noiseSFTs, f0, f0 + band, deltaF );
    if ( ret != XLAL_SUCCESS ) {
      ABORT ( status, XLAL_EFAILED, "XLALcheckNoiseSFTs() failed" );
    }
  }

  /* Signal parameters */
  N = (REAL4)(2*params->nSamples);   /* N = number of time-samples that would have been used to generate SFTs directly */
  halfAPlus = 0.5*N*params->pSigParams->pulsar.aPlus;
  halfACross = 0.5*N*params->pSigParams->pulsar.aCross;
  psi = params->pSigParams->pulsar.psi;
  /* cosPsi = (REAL4)cos(psi);
  sinPsi = (REAL4)sin(psi); */ /* 10/12/04 gam */
  cos2Psi = (REAL4)cos(2.0*psi);
  sin2Psi = (REAL4)sin(2.0*psi);
  f0Signal = params->pSigParams->pulsar.f0;
  phi0Signal = params->pSigParams->pulsar.phi0;
  if (params->pSigParams->pulsar.spindown) {
    spOrder = params->pSigParams->pulsar.spindown->length;
  } else {
    spOrder = 0;
  }

  /* loop that generates each SFT */
  for (iSFT = 0; iSFT < numSFTs; iSFT++)  {

      thisSFT = &(sftvect->data[iSFT]); /* select the SFT to work on */

      /* find fPlus, fCross, and the real and imaginary parts of the modulated amplitude, realA and imagA */
      /* fPlus = input->fPlusZeroPsi[iSFT]*cosPsi + input->fCrossZeroPsi[iSFT]*sinPsi;
      fCross = input->fCrossZeroPsi[iSFT]*cosPsi - input->fPlusZeroPsi[iSFT]*sinPsi; */ /* 10/12/04 gam */
      fPlus = input->fPlusZeroPsi[iSFT]*cos2Psi + input->fCrossZeroPsi[iSFT]*sin2Psi;
      fCross = input->fCrossZeroPsi[iSFT]*cos2Psi - input->fPlusZeroPsi[iSFT]*sin2Psi;
      realA = (REAL8)(halfAPlus*fPlus);
      imagA = (REAL8)(halfACross*fCross);

      /* Compute sums used to find the phase at the beginning of each SFT and kappa associated with fOneHalf*/
      /* xSum and ySum are the same as xSum and ySum in LALDemod */
      tmpInt = 2*iSFT*(spOrder+1)+1;
      xSum = 0.0;
      ySum = 0.0;
      for(s=0;s<spOrder;s++) {
        xSum += params->pSigParams->pulsar.spindown->data[s] * input->skyConst[tmpInt + 2 + 2*s];
        ySum += params->pSigParams->pulsar.spindown->data[s] * input->skyConst[tmpInt + 1 + 2*s];
      }

      /* find kappa associated with fOneHalf */
      /* fOneHalf = (f0Signal*input->skyConst[tmpInt] + xSum)/tSFT; */
      /* Do not need to actually compute this, just kappa */
      /* kappa = REAL8 index associated with fOneHalf; usually not an integer  */
      kappa = f0Signal*input->skyConst[tmpInt] + xSum;

      if (params->resTrig > 0) {
        /* if (params->resTrig > 0) use LUT for trig calls to find sin and cos, else will use standard sin and cos */

        /* Compute phase at the beginning of each SFT, called yTmp */
        /* yTmp is -y from LALDemod plus phi0Signal */
        /* Qcc is the complex conjugate of Q in LALDemod times exp(i*phi0Signal) */
        /* Using LUT to find cos(yTmp) and sin(yTmp) */
        yTmp = phi0Signal/real8TwoPi + f0Signal*input->skyConst[tmpInt-1] + ySum;
        varTmp = yTmp-(INT4)yTmp;
        /* indexTrig=(INT4)(varTmp*params->resTrig+0.5); */ /* 10/08/04 gam */
        indexTrig=(INT4)((varTmp + 1.0)*halfResTrig + 0.5);
        dTmp = real8TwoPi*varTmp - params->trigArg[indexTrig];
        dTmp2 = 0.5*dTmp*dTmp;
        sinTmp = params->sinVal[indexTrig];
        cosTmp = params->cosVal[indexTrig];
        imagQcc = sinTmp + dTmp*cosTmp - dTmp2*sinTmp;
        realQcc = cosTmp - dTmp*sinTmp - dTmp2*cosTmp;

        /* Find sin(2*pi*kappa) and 1 - cos(2*pi*kappa) */
        /* Using LUT to find sin(2*pi*kappa) and 1 - cos(2*pi*kappa) */
        varTmp = kappa-(INT4)kappa;
        /* indexTrig=(INT4)(varTmp*params->resTrig+0.5); */
        indexTrig=(INT4)((varTmp + 1.0)*halfResTrig + 0.5); /* 10/08/04 gam */
        dTmp = real8TwoPi*varTmp - params->trigArg[indexTrig];
        dTmp2 = 0.5*dTmp*dTmp;
        sinTmp = params->sinVal[indexTrig];
        cosTmp = params->cosVal[indexTrig];
        sin2PiKappa = sinTmp + dTmp*cosTmp - dTmp2*sinTmp;
        oneMinusCos2PiKappa = 1.0 - cosTmp + dTmp*sinTmp + dTmp2*cosTmp;

      } else {
        /* if (params->resTrig > 0) use LUT for trig calls to find sin and cos, else will use standard sin and cos */

        /* Compute phase at the beginning of each SFT, called yTmp */
        /* yTmp is -y from LALDemod plus phi0Signal */
        /* Qcc is the complex conjugate of Q in LALDemod times exp(i*phi0Signal) */
        yTmp = phi0Signal + real8TwoPi*(f0Signal*input->skyConst[tmpInt-1] + ySum);
        realQcc = cos(yTmp);
        imagQcc = sin(yTmp);

        /* Find sin(2*pi*kappa) and 1 - cos(2*pi*kappa) */
        /* use xTmp as temp storage for 2\pi\kappa; note xTmp = 2\pi(\kappa -k) is used in loop below */
        xTmp = real8TwoPi*kappa;
        sin2PiKappa = sin(xTmp);
        oneMinusCos2PiKappa = 1.0 - cos(xTmp);

      } /* END if (params->resTrig > 0) else ... */

      /* 09/07/05 gam; use Dterms to fill in SFT bins with fake data as per LALDemod else fill in bin with zero */
      k1=(INT4)kappa-params->Dterms+1; /* This is the same as k1 in LALDemod */
      jStart = k1 - k0;
      if (jStart < 0) jStart = 0;
      jEnd = k1 + 2*params->Dterms - k0;
      if (jEnd > SFTlen) jEnd = SFTlen;

      /* fill in the data */
      if (setToZero) {
        for (j=0; j<jStart; j++) {
          thisSFT->data->data[j].realf_FIXME = 0.0;
          thisSFT->data->data[j].imagf_FIXME = 0.0;
        }
      }
      /* This is the same as the inner most loop over k in LALDemod */
      for (j=jStart; j<jEnd; j++) {
          k = k0 + j;  /* k is the index of the frequency associated with index j */
          /* xTmp is the same as x in LALDemod */
          xTmp=real8TwoPi*(kappa - ((REAL8)k));
          /* Pcc is the complex conjugate of P in LALDemod */
          if (fabs(xTmp) < smallX) {
             /* If xTmp is small we need correct xTmp->0 limit */
             realPcc=1.0;
             imagPcc=0.0;
          } else {
             realPcc=sin2PiKappa/xTmp;
             imagPcc=oneMinusCos2PiKappa/xTmp;
          }
          realTmp = realQcc*realPcc - imagQcc*imagPcc;
          imagTmp = realQcc*imagPcc + imagQcc*realPcc;
          thisSFT->data->data[j].realf_FIXME = (REAL4)(realTmp*realA - imagTmp*imagA);
          thisSFT->data->data[j].imagf_FIXME = (REAL4)(realTmp*imagA + imagTmp*realA);
      } /* END for (j=jStart; j<jEnd; j++) */
      if (setToZero) {
        for (j=jEnd; j<SFTlen; j++) {
          thisSFT->data->data[j].realf_FIXME = 0.0;
          thisSFT->data->data[j].imagf_FIXME = 0.0;
        }
      }
      /* fill in SFT metadata */
      thisSFT->epoch = params->pSFTParams->timestamps->data[iSFT];
      thisSFT->f0 = f0;          /* start frequency */
      thisSFT->deltaF = deltaF;  /* frequency resolution */

      /* Now add the noise-SFTs if given */
      if (params->pSFTParams->noiseSFTs) {
        thisNoiseSFT = &(params->pSFTParams->noiseSFTs->data[iSFT]);
        index0n = (INT4)( (thisSFT->f0 - thisNoiseSFT->f0)*tSFT + 0.5 );
        for (j=0; j < SFTlen; j++)
        {
           thisSFT->data->data[j].realf_FIXME += crealf(thisNoiseSFT->data->data[index0n + j]);
           thisSFT->data->data[j].imagf_FIXME += cimagf(thisNoiseSFT->data->data[index0n + j]);
        } /* for j < SFTlen */
      }
  } /* for iSFT < numSFTs */

  /* prepare SFT-vector for return */
  if (*outputSFTs == NULL) {
    *outputSFTs = sftvect;
  } /* else sftvect already points to same memory as *outputSFTs */

  DETATCHSTATUSPTR( status );
  RETURN (status);
} /* LALFastGeneratePulsarSFTs () */



/*--------------- some useful helper-functions ---------------*/

/** Convert earth-frame GPS time into barycentric-frame SSB time for given source.
 * \note The only fields used in params are: \a site, \a pulsar.position
 * and \a ephemerides.
 */
int
XLALConvertGPS2SSB ( LIGOTimeGPS *SSBout, 		/**< [out] arrival-time in SSB */
                     LIGOTimeGPS GPSin, 		/**< [in]  GPS-arrival time at detector */
                     const PulsarSignalParams *params 	/**< define source-location and detector */
                     )
{
  XLAL_CHECK ( SSBout != NULL, XLAL_EINVAL, "Invalid NULL input 'SSBout'\n" );
  XLAL_CHECK ( params != NULL, XLAL_EINVAL, "Invalid NULL input 'params'\n" );

  BarycenterInput baryinput = empty_BarycenterInput;
  baryinput.site = *(params->site);
  /* account for a quirk in LALBarycenter(): -> see documentation of type BarycenterInput */
  baryinput.site.location[0] /= LAL_C_SI;
  baryinput.site.location[1] /= LAL_C_SI;
  baryinput.site.location[2] /= LAL_C_SI;
  XLAL_CHECK ( params->pulsar.position.system == COORDINATESYSTEM_EQUATORIAL, XLAL_EDOM, "Non-equatorial coords not implemented yet\n");

  SkyPosition tmp = params->pulsar.position;
  XLALNormalizeSkyPosition ( &(tmp.longitude), &(tmp.latitude) );
  baryinput.alpha = tmp.longitude;
  baryinput.delta = tmp.latitude;
  baryinput.dInv = 0.e0;	/* following makefakedata_v2 */

  baryinput.tgps = GPSin;


  EarthState earth;
  EmissionTime emit;

  int ret = XLALBarycenterEarth ( &earth, &GPSin, params->ephemerides );
  XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC );

  ret = XLALBarycenter ( &emit, &baryinput, &earth );
  XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC );

  (*SSBout) = emit.te;

  return XLAL_SUCCESS;

} /* XLALConvertGPS2SSB() */

/**
 * Convert  barycentric frame SSB time into earth-frame GPS time
 *
 * NOTE: this uses simply the inversion-routine used in the original
 *       makefakedata_v2
 */
int XLALConvertSSB2GPS ( LIGOTimeGPS *GPSout,			/**< [out] GPS-arrival-time at detector */
                         LIGOTimeGPS SSBin,			/**< [in] input: signal arrival time at SSB */
                         const PulsarSignalParams *params	/**< params defining source-location and detector */
                         )
{
  XLAL_CHECK ( GPSout != NULL, XLAL_EINVAL, "Invalid NULL output-pointer 'GPSout'\n");
  XLAL_CHECK ( params != NULL, XLAL_EINVAL, "Invalid NULL input 'params'\n");

  /*
   * To start root finding, use SSBpulsarparams as guess
   * (not off by more than 400 secs!
   */
  LIGOTimeGPS GPSguess = SSBin;
  UINT4 flip_flop_counter = 0;
  INT8 delta;

  /* now find GPS time corresponding to SSBin by iterations */
  UINT4 iterations;
  for ( iterations = 0; iterations < 100; iterations++ )
    {
      LIGOTimeGPS SSBofguess;

      /* find SSB time of guess */
      int ret = XLALConvertGPS2SSB ( &SSBofguess, GPSguess, params);
      XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC );

      /* compute difference between that and what we want */
      delta = XLALGPSToINT8NS( &SSBin ) - XLALGPSToINT8NS( &SSBofguess );

      /* if we are within 1ns of the result increment the flip-flop counter */
      if ( abs(delta) == 1) {
        flip_flop_counter ++;
      }

      /* break if we've converged: let's be strict to < 1 ns ! */
      /* also break if the flip-flop counter has reached 3 */
      if ( (delta == 0) || (flip_flop_counter >= 3) ) {
	break;
      }

      /* use delta to make next guess */
      INT8 guess = XLALGPSToINT8NS ( &GPSguess );
      guess += delta;

      XLALINT8NSToGPS( &GPSguess, guess );

    } /* for iterations < 100 */

  /* check for convergence of root finder */
  if ( iterations == 100 ) {
    XLAL_ERROR ( XLAL_EFAILED, "SSB->GPS iterative conversion failed to converge to <= 1ns within 100 iterations: delta = %d ns\n", delta );
  }

  /* if we exited because of flip-flop and final delta was +1 then round up to the higher value */
  /* otherwise we are already at the higher value and we do nothing */
  if ( (flip_flop_counter == 3) && (delta == +1) )
    {
      XLALINT8NSToGPS( &GPSguess, XLALGPSToINT8NS ( &GPSguess ) + 1 );
    }

  /* Now that we've found the GPS time that corresponds to the given SSB time */
  (*GPSout) = GPSguess;

  return XLAL_SUCCESS;

} /* XLALConvertSSB2GPS() */


/* ***********************************************************************
 * the following are INTERNAL FUNCTIONS not to be called outside of this
 * module
 ************************************************************************/

/** Check that all timestamps given lie within the range [t0, t1]
 *
 *  return: 0 if ok, ERROR if not
 */
int
XLALcheck_timestamp_bounds ( const LIGOTimeGPSVector *timestamps, LIGOTimeGPS t0, LIGOTimeGPS t1 )
{
  XLAL_CHECK ( timestamps != NULL, XLAL_EINVAL, "Invalid NULL input 'timestamps'\n");
  UINT4 numTimestamps = timestamps->length;
  XLAL_CHECK ( numTimestamps > 0, XLAL_EDOM, "Invalid zero-length vector 'timestamps'\n");
  REAL8 t0R = XLALGPSGetREAL8 ( &t0 );
  REAL8 t1R = XLALGPSGetREAL8 ( &t1 );
  XLAL_CHECK ( t1R >= t0R, XLAL_EDOM, "Invalid negative time range: t0=%f must be <= t1=%f]\n", t0R, t1R );

  for (UINT4 i = 0; i < numTimestamps; i ++)
    {
      LIGOTimeGPS *ti = &(timestamps->data[i]);

      REAL8 tiR = XLALGPSGetREAL8 ( ti );

      REAL8 diff0 = XLALGPSDiff ( ti, &t0 );
      XLAL_CHECK ( diff0 >= 0, XLAL_EDOM, "Timestamp i=%d  outside of bounds: t_i = %f < [%f,%f]\n", i, tiR, t0R, t1R );

      REAL8 diff1 = XLALGPSDiff ( &t1, ti );
      XLAL_CHECK ( diff1 >= 0, XLAL_EDOM, "Timestamp i=%d  outside of bounds: t_i = %f > [%f,%f]\n", i, tiR, t0R, t1R );

    } /* for i < numTimestamps */

  return XLAL_SUCCESS;

} /* XLALcheck_timestamp_bounds() */

/** Check if frequency-range and resolution of noiseSFTs is consistent with signal-band [f0, f1]
 * \note All frequencies f are required to correspond to integer *bins* f/dFreq, ABORT if not
 *
 * return XLAL_SUCCESS if everything fine, error-code otherwise
 */
static int
XLALcheckNoiseSFTs ( const SFTVector *sfts, REAL8 f0, REAL8 f1, REAL8 deltaF )
{
  XLAL_CHECK ( sfts != NULL, XLAL_EINVAL, "Invalid NULL input 'sfts'\n" );
  XLAL_CHECK ( (f0 >= 0) && (f1 > 0) && (deltaF > 0), XLAL_EDOM, "Invalid non-positive frequency input: f0 = %g, f1 = %g, deltaF = %g\n", f0, f1, deltaF );

  int ret;

  for ( UINT4 i = 0; i < sfts->length; i++ )
    {
      SFTtype *thisSFT = &(sfts->data[i]);
      REAL8 deltaFn    = thisSFT->deltaF;
      REAL8 fn0        = thisSFT->f0;
      REAL8 fn1        = f0 + thisSFT->data->length * deltaFn;

      ret = gsl_fcmp ( deltaFn, deltaF, eps );
      XLAL_CHECK ( ret == 0, XLAL_ETOL, "Time-base of noise-SFTs Tsft_n=%f differs from signal-SFTs Tsft=%f\n", 1.0/deltaFn, 1.0/deltaF );

      XLAL_CHECK ( (f0 >= fn0) && (f1 <= fn1), XLAL_EDOM, "Signal frequency-band [%f,%f] is not contained in noise SFTs [%f,%f]\n", f0, f1, fn0, fn1 );

      /* all frequencies here must correspond to exact integer frequency-bins (wrt dFreq = 1/TSFT) */
      REAL8 binReal    = f0 / deltaF;
      REAL8 binRounded = round ( binReal );
      ret = gsl_fcmp ( binReal, binRounded, eps );
      XLAL_CHECK ( ret == 0, XLAL_ETOL, "Signal-band frequency f0/deltaF = %.16g differs from integer bin by more than %g relative deviation.\n", binReal, eps );

      binReal = f1 / deltaF;
      binRounded = round ( binReal );
      ret = gsl_fcmp ( binReal, binRounded, eps );
      XLAL_CHECK ( ret == 0, XLAL_ETOL, "Signal-band frequency f1/deltaF = %.16g differs from integer bin by more than %g relative deviation.\n", binReal, eps );

      binReal = fn0 / deltaF;
      binRounded = round ( binReal );
      ret = gsl_fcmp ( binReal, binRounded, eps );
      XLAL_CHECK ( ret == 0, XLAL_ETOL, "Noise-SFT start frequency fn0/deltaF = %.16g differs from integer bin by more than %g relative deviation.\n", binReal, eps );

    } /* for i < numSFTs */

  return XLAL_SUCCESS;

} /* XLALcheckNoiseSFTs() */



/** Yousuke's phase-correction function, taken from makefakedata_v2
 */
int
XLALcorrect_phase ( SFTtype *sft, LIGOTimeGPS tHeterodyne )
{
  XLAL_CHECK ( sft != NULL, XLAL_EINVAL, "Invalid NULL input 'sft'\n");

  REAL8 deltaT = XLALGPSDiff ( &(sft->epoch), &tHeterodyne );
  REAL8 deltaFT = deltaT * sft->f0;	// freq * time

  /* check if we really need to do anything here? (i.e. is deltaT an integer?) */
  if ( fabs (deltaFT - (INT4) deltaFT ) > eps )
    {
      XLALPrintWarning ("XLALcorrect_phase(): we DO need to apply heterodyning phase-correction\n");

      REAL8 deltaPhase = deltaFT * LAL_TWOPI;	// 'phase' = freq * time * 2pi

      REAL8 cosx = cos ( deltaPhase );
      REAL8 sinx = sin ( deltaPhase );

      for (UINT4 i = 0; i < sft->data->length; i++ )
	{
	  COMPLEX8 fvec1 = sft->data->data[i];
	  sft->data->data[i].realf_FIXME = crealf(fvec1) * cosx - cimagf(fvec1) * sinx;
	  sft->data->data[i].imagf_FIXME = cimagf(fvec1) * cosx + crealf(fvec1) * sinx;
	} /* for i < length */

    } /* if deltaFT not integer */

  return XLAL_SUCCESS;

} /* XLALcorrect_phase() */
