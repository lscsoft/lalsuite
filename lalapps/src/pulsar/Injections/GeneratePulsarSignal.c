/************************************ <lalVerbatim file="GeneratePulsarSignalCV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{GeneratePulsarSignal}
\label{ss:GeneratePulsarSignal.c}

Module to generate subsequent sky-position.

\subsubsection*{Prototypes}
\input{GeneratePulsarSignalCP}
\idx{NextSkyPosition()}

\subsubsection*{Description}

This module generates a fake pulsar-signal, either for an isolated or a binary pulsar. 

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GeneratePulsarSignalCV}}

******************************************************* </lalLaTeX> */
#include <math.h>

#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>

#include "GeneratePulsarSignal.h"

/*----------------------------------------------------------------------*/
/* prototypes for internal functions */
void write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series);
static void write_timeSeriesR8 (FILE *fp, const REAL8TimeSeries *series);
static LIGOTimeGPSVector* make_timestamps (const LIGOTimeGPS *tStart, REAL8 duration, REAL8 Tsft);
static int check_timestamp_bounds (const LIGOTimeGPSVector *timestamps, LIGOTimeGPS t0, LIGOTimeGPS t1);
static void checkNoiseSFTs (LALStatus *stat, const SFTVector *sfts, REAL8 f0, REAL8 f1, REAL8 deltaF);
/*----------------------------------------------------------------------*/

NRCSID( PULSARSIGNALC, "$Id$");

extern INT4 lalDebugLevel;

static REAL8 eps = 1.e-14;	/* maximal REAL8 roundoff-error (used for determining if some number is an INT) */

#define LTT 1000
/* FIXME!!!!!: 
   -> we use this time-constant to make the signal start earlier and last longer than
   really required

   This is done according to makefakedata_v2, but it really points
   to a bug in SimulateCoherentGW(), which should be traced at some point!!
*/

/***********************************************************************
 * generate a time-series at the detector for a given pulsar
 ***********************************************************************/
void
LALGeneratePulsarSignal (LALStatus *stat, REAL4TimeSeries **signal, const PulsarSignalParams *params)
{
  static SpinOrbitCWParamStruc sourceParams;
  static CoherentGW sourceSignal;
  DetectorResponse detector;
  UINT4 SSBduration;
  LIGOTimeGPS time = {0,0};
  REAL4TimeSeries *output;
  UINT4 i;

  INITSTATUS( stat, "LALGeneratePulsarSignal", PULSARSIGNALC );
  ATTATCHSTATUSPTR (stat);

  ASSERT (signal != NULL, stat, PULSARSIGNALH_ENULL, PULSARSIGNALH_MSGENULL);
  ASSERT (*signal == NULL, stat,  PULSARSIGNALH_ENONULL,  PULSARSIGNALH_MSGENONULL);

  /*----------------------------------------------------------------------
   *
   * First call GenerateSpinOrbitCW() to generate the source-signal
   *
   *----------------------------------------------------------------------*/
  sourceParams.psi = params->pulsar.psi;
  sourceParams.aPlus = params->pulsar.aPlus;
  sourceParams.aCross = params->pulsar.aCross;
  sourceParams.phi0 = params->pulsar.phi0;
  sourceParams.f0 = params->pulsar.f0;
  /* set source position: make sure it's "normalized", i.e. [0<=alpha<2pi]x[-pi<=delta<=pi] */
  TRY( LALNormalizeSkyPosition (stat->statusPtr, &(sourceParams.position), &(params->pulsar.position)), stat);

  /* if pulsar is in binary-orbit, set binary parameters */
  if (params->orbit)
    {
      sourceParams.orbitEpoch = params->orbit->orbitEpoch;
      sourceParams.omega = params->orbit->omega;
      sourceParams.rPeriNorm = params->orbit->rPeriNorm;
      sourceParams.oneMinusEcc = params->orbit->oneMinusEcc;
      sourceParams.angularSpeed = params->orbit->angularSpeed;
    }
  else
    sourceParams.rPeriNorm = 0.0;		/* this defines an isolated pulsar */


  /* pulsar reference-time in SSB frame !*/
  if (params->pulsar.TRefSSB.gpsSeconds != 0)
    sourceParams.spinEpoch = params->pulsar.TRefSSB;
  else
    sourceParams.spinEpoch = time;	/* use start-time */

  /* sampling-timestep and length for source-parameters */
  sourceParams.deltaT = 60;	/* in seconds; hardcoded default from makefakedata_v2 */

  /* start-time in SSB time */
  TRY (ConvertGPS2SSB (stat->statusPtr, &time, params->startTimeGPS, params), stat);
  sourceParams.epoch = time;
  /* ----------
     FIXME: argh, circumvent some mysterious bug in SimulateCoherentGW()... */
  sourceParams.epoch.gpsSeconds -= 0.75*LTT;
  /*----------*/

  time = params->startTimeGPS;
  time.gpsSeconds += params->duration;
  TRY (ConvertGPS2SSB (stat->statusPtr, &time, time, params), stat);	 /* convert time to SSB */
  SSBduration = time.gpsSeconds - sourceParams.spinEpoch.gpsSeconds;
  sourceParams.length = (UINT4)( 1.0* SSBduration / sourceParams.deltaT );
  /* ----------
     FIXME: argh, circumvent some mysterious bug in SimulateCoherentGW()... */
  sourceParams.length   += (UINT4)(1.5*LTT/ sourceParams.deltaT);
  /*----------*/

  /* we use frequency-spindowns, but GenerateSpinOrbitCW wants it f0-normalized,
     so we have to do that here: */
  if (params->pulsar.spindown)
    {
      TRY ( LALDCreateVector (stat->statusPtr, &(sourceParams.f), params->pulsar.spindown->length), stat);
      for (i=0; i < sourceParams.f->length; i++)
	sourceParams.f->data[i] = params->pulsar.spindown->data[i] / params->pulsar.f0;
    }

  /* finally, call the function to generate the source waveform */
  TRY ( LALGenerateSpinOrbitCW (stat->statusPtr, &sourceSignal, &sourceParams), stat);

  /* free spindown-vector right away, so we don't forget */
  if (sourceParams.f) {
    TRY ( LALDDestroyVector (stat->statusPtr, &(sourceParams.f) ), stat);
  }

  /* check that sampling interval was short enough */
  if ( sourceParams.dfdt > 2.0 )  /* taken from makefakedata_v2 */
    {
      LALPrintError ("GenerateSpinOrbitCW() returned df*dt = %f > 2.0", sourceParams.dfdt);
      ABORT (stat, PULSARSIGNALH_ESAMPLING, PULSARSIGNALH_MSGESAMPLING);
    }

  if (lalDebugLevel >= 3) {
    TRY (PrintGWSignal (stat->statusPtr, &sourceSignal, "signal2.agr"), stat);
  }

  /*----------------------------------------------------------------------
   *
   * Now call the function to translate the source-signal into a (heterodyned)
   * signal at the detector 
   *
   *----------------------------------------------------------------------*/
  /* first set up the detector-response */
  detector.transfer = params->transferFunction;
  detector.site = params->site;
  detector.ephemerides = params->ephemerides;
  /* we need to set the heterodyne epoch in GPS time, but we want to use
   * the pulsar reference-time for that (which is in SSB), so we have to convert it first
   */
  TRY ( ConvertSSB2GPS (stat->statusPtr, &time, params->pulsar.TRefSSB, params), stat);
  detector.heterodyneEpoch = time;
  
  /* ok, we  need to prepare the output time-series */
  if ( (output = LALCalloc (1, sizeof (*output) )) == NULL) {
    ABORT (stat,  PULSARSIGNALH_EMEM,  PULSARSIGNALH_MSGEMEM);
  }
  LALCreateVector (stat->statusPtr, &(output->data), (UINT4)( params->samplingRate * params->duration) );
  BEGINFAIL(stat) {
    LALFree (output);
  } ENDFAIL(stat);

  output->deltaT = 1.0 / params->samplingRate;
  output->f0 = params->fHeterodyne;
  output->epoch = params->startTimeGPS;

  
  TRY ( LALSimulateCoherentGW (stat->statusPtr, output, &sourceSignal, &detector ), stat );

			       
  /*----------------------------------------------------------------------*/
  /* Free all allocated memory that is not returned */
  TRY (LALSDestroyVectorSequence ( stat->statusPtr, &(sourceSignal.a->data)), stat);
  LALFree ( sourceSignal.a );
  TRY (LALSDestroyVector (stat->statusPtr, &(sourceSignal.f->data ) ), stat);
  LALFree (sourceSignal.f);
  TRY (LALDDestroyVector (stat->statusPtr, &(sourceSignal.phi->data )), stat);
  LALFree (sourceSignal.phi);

  *signal = output;

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* LALGeneratePulsarSignal() */



/*----------------------------------------------------------------------
 * take time-series as input, convert it into SFTs and add noise if given
 *----------------------------------------------------------------------*/
void
LALSignalToSFTs (LALStatus *stat, SFTVector **outputSFTs, const REAL4TimeSeries *signal, const SFTParams *params)
{
  UINT4 numSFTs;			/* number of SFTs */
  UINT4 numSamples;			/* number of time-samples in each Tsft */
  UINT4 iSFT;
  REAL8 Band, samples, f0, deltaF;
  RealFFTPlan *pfwd = NULL;
  LIGOTimeGPSVector *timestamps = NULL;
  REAL4Vector timeStretch = {0,0};
  LIGOTimeGPS tStart;			/* start time of input time-series */
  LIGOTimeGPS tLast;			/* start-time of last _possible_ SFT */
  REAL8 duration, delay;
  UINT4 SFTlen;				/* number of samples in an SFT */
  UINT4 indexShift;
  UINT4 index0n;			/* first frequency-bin to use from noise-SFT */
  SFTtype *thisSFT, *thisNoiseSFT;	/* SFT-pointers */
  REAL4 renorm;				/* renormalization-factor of taking only part of an SFT */
  SFTVector *sftvect = NULL;		/* return value. For better readability */
  UINT4 j;


  INITSTATUS( stat, "AddSignalToSFTs", PULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );
  
  ASSERT (outputSFTs != NULL, stat, PULSARSIGNALH_ENULL, PULSARSIGNALH_MSGENULL);
  ASSERT (*outputSFTs == NULL, stat,  PULSARSIGNALH_ENONULL,  PULSARSIGNALH_MSGENONULL);
  ASSERT (signal != NULL, stat, PULSARSIGNALH_ENULL, PULSARSIGNALH_MSGENULL);
  ASSERT (params != NULL, stat, PULSARSIGNALH_ENULL, PULSARSIGNALH_MSGENULL);
  if ( params->timestamps && params->noiseSFTs) {
    ASSERT ( params->timestamps->length == params->noiseSFTs->length, stat,  
	     PULSARSIGNALH_ENUMSFTS,  PULSARSIGNALH_MSGENUMSFTS);
  }

  f0 = signal->f0;				/* lowest frequency */
  Band = 1.0 / (2.0 * signal->deltaT);		/* NOTE: frequency-band is determined by sampling-rate! */
  deltaF = 1.0 / params->Tsft;			/* frequency-resolution */

  /* if noiseSFTs are given: check they are consistent with signal! */
  if (params->noiseSFTs) {
    TRY (checkNoiseSFTs (stat->statusPtr, params->noiseSFTs, f0, f0 + Band, deltaF), stat);
  }
    
  /* make sure that number of samples is an integer (up to possible rounding errors */
  samples = 2.0 * Band * params->Tsft;		/* this is a float!*/
  numSamples = (UINT4) (samples + 0.5);		/* round to int */
  ASSERT ( fabs(samples - numSamples)/samples < eps, stat, 
	   PULSARSIGNALH_EINCONSBAND, PULSARSIGNALH_MSGEINCONSBAND);

  
  /* Prepare FFT: compute plan for FFTW */
  /* FIXME: put some buffering + measure best plan */
  TRY (LALCreateForwardRealFFTPlan(stat->statusPtr, &pfwd, numSamples, 0), stat); 	

  /* get some info about time-series */
  tStart = signal->epoch;					/* start-time of time-series */
  duration = 1.0* signal->data->length * signal->deltaT;	/* total duration in seconds */
  /* get last possible start-time for an SFT */
  TRY ( LALAddFloatToGPS (stat->statusPtr, &tLast, &tStart, duration - params->Tsft), stat);

  /* for simplicity we _always_ work with timestamps. 
   * Therefore, we have to generate them now if none have been provided by the user. */
  if (params->timestamps == NULL)
    {
      timestamps = make_timestamps (&tStart, duration, params->Tsft );
      if (timestamps == NULL) {
	ABORT (stat, PULSARSIGNALH_EMEM, PULSARSIGNALH_MSGEMEM);
      }
    }
  else	/* if given, nothing to do */
    timestamps = params->timestamps;

  /* check that all timestamps lie within [tStart, tLast] */
  if ( check_timestamp_bounds (timestamps, tStart, tLast) != 0) {
    ABORT (stat, PULSARSIGNALH_ETIMEBOUND, PULSARSIGNALH_MSGETIMEBOUND);
  }

  /* prepare SFT-vector for return */
  numSFTs = timestamps->length;			/* number of SFTs to produce */
  SFTlen = (UINT4)(numSamples/2) + 1;		/* number of frequency-bins per SFT */

  LALCreateSFTVector (stat->statusPtr, &sftvect, numSFTs, SFTlen);
  BEGINFAIL (stat) {
    if (params->timestamps == NULL)
      LALDestroyTimestampVector (stat->statusPtr, &timestamps);
  } ENDFAIL (stat);

  /* main loop: FFT the stuff */
  for (iSFT = 0; iSFT < numSFTs; iSFT++)
    {
      thisSFT = &(sftvect->data[iSFT]);	/* point to current SFT-slot */

      /* find the start-bin for this SFT in the time-series */
      TRY ( LALDeltaFloatGPS (stat->statusPtr, &delay, &(timestamps->data[iSFT]), &tStart), stat);
      indexShift = (UINT4) (delay / signal->deltaT + 0.5);	/* round properly */
      timeStretch.length = numSamples;    
      timeStretch.data = signal->data->data + indexShift;	/* point to the right sample-bin */
      /* FIXME: double-check that timestamps are consistent with time-series sampling! */

      /* the central step: FFT the ith time-stretch into an SFT-slot */
      LALForwardRealFFT (stat->statusPtr, thisSFT->data, &timeStretch, pfwd);
      BEGINFAIL(stat) {
	LALDestroySFTVector (stat->statusPtr, &sftvect);
      } ENDFAIL(stat);

      /* prepare the i'th output SFT with the result */
      thisSFT->epoch = timestamps->data[iSFT];	/* set the proper timestamp */
      thisSFT->f0 = signal->f0;			/* minimum frequency */
      thisSFT->deltaF = 1.0 / params->Tsft;	/* frequency-spacing */

      /* Now add the noise-SFTs if given */
      if (params->noiseSFTs)
	{
	  thisNoiseSFT = &(params->noiseSFTs->data[iSFT]);
	  index0n = (UINT4) ((thisSFT->f0 - thisNoiseSFT->f0) / thisSFT->deltaF);
	  /* The renormalization follows strictly makefakedata_v2, even if I admittedly
	     don't quite understand this... */
	  renorm = 1.0*SFTlen/(thisNoiseSFT->data->length);
	  for (j=0; j < SFTlen; j++)
	    {
	      thisSFT->data->data[j].re += renorm * thisNoiseSFT->data->data[index0n + j].re;
	      thisSFT->data->data[j].im += renorm * thisNoiseSFT->data->data[index0n + j].im;
	    } /* for j < SFTlen */
	}

    } /* for iSFT < numSFTs */ 



  /* clean up */ /* FIXME: buffering */
  LALDestroyRealFFTPlan(stat->statusPtr, &pfwd);

  /* did we get timestamps or did we make them? */
  if (params->timestamps == NULL)
    {
      LALFree (timestamps->data);
      LALFree (timestamps);
      timestamps = NULL;
    }

  *outputSFTs = sftvect;

  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* LALSignalToSFTs() */




/*----------------------------------------------------------------------
 * convert earth-frame GPS time into barycentric-frame SSB time for given source
 *----------------------------------------------------------------------*/
void
ConvertGPS2SSB (LALStatus* stat, LIGOTimeGPS *SSBout, LIGOTimeGPS GPSin, const PulsarSignalParams *params)
{
  EarthState earth;
  EmissionTime emit;
  BarycenterInput baryinput;
  SkyPosition tmp;

  INITSTATUS( stat, "ConvertGPS2SSB", PULSARSIGNALC );
  ATTATCHSTATUSPTR (stat);

  ASSERT (SSBout != NULL, stat,  PULSARSIGNALH_ENULL,  PULSARSIGNALH_MSGENULL);
  ASSERT (params != NULL, stat,  PULSARSIGNALH_ENULL,  PULSARSIGNALH_MSGENULL);

  baryinput.site = *(params->site);
  /* account for a quirk in LALBarycenter(): -> see documentation of type BarycenterInput */
  baryinput.site.location[0] /= LAL_C_SI;
  baryinput.site.location[1] /= LAL_C_SI;
  baryinput.site.location[2] /= LAL_C_SI;
  if (params->pulsar.position.system != COORDINATESYSTEM_EQUATORIAL)
    {
      /* FIXME: do proper conversion or error-reporting */
      printf ("\nARGH: non-equatorial coords not implemented here yet, sorry!\n");
      exit(-1);	/* suicide */
    }

  TRY( LALNormalizeSkyPosition (stat->statusPtr, &tmp, &(params->pulsar.position)), stat);
  baryinput.alpha = tmp.longitude;
  baryinput.delta = tmp.latitude;
  baryinput.dInv = 0.e0;	/* following makefakedata_v2 */

  baryinput.tgps = GPSin;

  TRY (LALBarycenterEarth(stat->statusPtr, &earth, &GPSin, params->ephemerides), stat);
  TRY (LALBarycenter(stat->statusPtr, &emit, &baryinput, &earth), stat);

  *SSBout = emit.te;

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* ConvertGPS2SSB() */


/*----------------------------------------------------------------------
 * convert  barycentric frame SSB time into earth-frame GPS time
 *
 * NOTE: this uses simply the inversion-routine used in the original
 *       makefakedata_v2
 *----------------------------------------------------------------------*/
void
ConvertSSB2GPS (LALStatus *stat, LIGOTimeGPS *GPSout, LIGOTimeGPS SSBin, const PulsarSignalParams *params)
{
  LIGOTimeGPS SSBofguess;
  LIGOTimeGPS GPSguess;
  INT4 iterations, E9=1000000000;
  INT8 delta, guess;

  INITSTATUS( stat, "ConvertSSB2GPS", PULSARSIGNALC );
  ATTATCHSTATUSPTR (stat);


  /* 
   * To start root finding, use SSBpulsarparams as guess 
   * (not off by more than 400 secs! 
   */
  GPSguess = SSBin;

  /* now find GPS time corresponding to SSBin by iterations */
  for (iterations = 0; iterations < 100; iterations++) 
    {
      /* find SSB time of guess */
      TRY ( ConvertGPS2SSB (stat->statusPtr, &SSBofguess, GPSguess, params), stat);

      /* compute difference between that and what we want */
      delta  = SSBin.gpsSeconds;
      delta -= SSBofguess.gpsSeconds;
      delta *= E9;
      delta += SSBin.gpsNanoSeconds;
      delta -= SSBofguess.gpsNanoSeconds;
      
      /* break if we've converged: let's be strict to < 1 ns ! */
      if (delta == 0)
	break;

      /* use delta to make next guess */
      guess  = GPSguess.gpsSeconds;
      guess *= E9;
      guess += GPSguess.gpsNanoSeconds;
      guess += delta;

      GPSguess.gpsSeconds = guess / E9;
      guess -= GPSguess.gpsSeconds * E9;	/* get ns remainder */
      GPSguess.gpsNanoSeconds = guess;

    } /* for iterations < 100 */

  /* check for convergence of root finder */
  if (iterations == 100) {
    ABORT ( stat, PULSARSIGNALH_ESSBCONVERT, PULSARSIGNALH_MSGESSBCONVERT);
  }

  /* Now that we've found the GPS time that corresponds to the given SSB time */
  *GPSout = GPSguess;
  
  
  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* ConvertSSB2GPS() */



/*----------------------------------------------------------------------
 * export a REAL4 time-series in an xmgrace graphics file 
 *----------------------------------------------------------------------*/
void
PrintR4TimeSeries (LALStatus *stat, const REAL4TimeSeries *series, const CHAR *fname)
{
  FILE *fp = NULL;
  UINT4 set;
  CHAR *xmgrHeader = 
    "@version 50103\n"
    "@xaxis label \"T (sec)\"\n"
    "@yaxis label \"A\"\n";

  /* Set up shop. */
  INITSTATUS( stat, "PrintR4TimeSeries", PULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );

  fp = fopen (fname, "w");

  if( !fp ) {
    ABORT (stat, PULSARSIGNALH_ESYS, PULSARSIGNALH_MSGESYS);
  }
  
  fprintf (fp, xmgrHeader);
  fprintf (fp, "@title \"REAL4 Time-Series: %s\"\n", series->name);
  fprintf (fp, "@subtitle \"GPS-start: %f, f0 = %e\"\n", 
	   1.0*series->epoch.gpsSeconds + series->epoch.gpsNanoSeconds / 1.0e9, series->f0);

  /* Print set header. */
  set = 0;
  fprintf( fp, "\n@target G0.S%d\n@type xy\n", set);
  fprintf( fp, "@s%d color (0,0,0)\n", set );

  write_timeSeriesR4 (fp, series);
      
  fclose(fp);

  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* PrintR4TimeSeries() */

/*----------------------------------------------------------------------
 * write a CoherentGW type signal into an xmgrace file
 *----------------------------------------------------------------------*/
void
PrintGWSignal (LALStatus *stat, const CoherentGW *signal, const CHAR *fname)
{
  FILE *fp;
  UINT4 set;
  CHAR *xmgrHeader = 
    "@version 50103\n"
    "@xaxis label \"T (sec)\"\n"
    "@yaxis label \"A\"\n";


  INITSTATUS( stat, "PrintGWSignal", PULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );

  fp = fopen (fname, "w");

  if( !fp ) {
    ABORT (stat, PULSARSIGNALH_ESYS, PULSARSIGNALH_MSGESYS);
  }

  fprintf (fp, xmgrHeader);
  fprintf (fp, "@title \"GW source signal\"\n");
  fprintf (fp, "@subtitle \"position: (alpha,delta) = (%f, %f), psi = %f\"\n", 
	   signal->position.longitude, signal->position.latitude, signal->psi);


  /* Print set header. */
  set = 0;
  fprintf( fp, "\n@target G0.S%d\n@type xy\n", set);
  fprintf( fp, "@s%d color (0,0,0)\n", set );
  write_timeSeriesR4 (fp, signal->f);

  /* Print set header. */
  set ++;
  fprintf( fp, "\n@target G0.S%d\n@type xy\n", set);
  fprintf( fp, "@s%d color (0,0,0)\n", set );
  write_timeSeriesR8 (fp, signal->phi);

  fclose (fp);

  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* PrintGWSignal() */


/* internal helper-function */
void
write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series)
{
  REAL8 timestamp; 
  UINT4 i;

  if (series == NULL)
    {
      printf ("\nempty input!\n");
      return; 
    }

  timestamp = 1.0*series->epoch.gpsSeconds + series->epoch.gpsNanoSeconds * 1.0e-9;

  for( i = 0; i < series->data->length; i++)
    {
      fprintf( fp, "%16.9f %e\n", timestamp, series->data->data[i] );
      timestamp += series->deltaT;
    }

  return;

} /* write_timeSeriesR4() */

void
write_timeSeriesR8 (FILE *fp, const REAL8TimeSeries *series)
{
  REAL8 timestamp; 
  UINT4 i;

  timestamp = 1.0*series->epoch.gpsSeconds + series->epoch.gpsNanoSeconds * 1.0e-9;

  for( i = 0; i < series->data->length; i++)
    {
      fprintf( fp, "%f %e\n", timestamp, series->data->data[i] );
      timestamp += series->deltaT;
    }

  return;

} /* write_timeSeriesR4() */




void 
write_SFT (LALStatus *stat, const SFTtype *sft, const CHAR *fname)
{
  FILE *fp;
  REAL4 rpw,ipw;
  INT4 i;
  REAL8 Tsft;
  

  struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
  } header;

  INITSTATUS( stat, "write_SFT", PULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );

  fp = fopen (fname, "w");

  if( !fp ) {
    ABORT (stat, PULSARSIGNALH_ESYS, PULSARSIGNALH_MSGESYS);
  }

  header.endian=1.0;
  header.gps_sec = sft->epoch.gpsSeconds;
  header.gps_nsec = sft->epoch.gpsNanoSeconds;

  Tsft = 1.0 / sft->deltaF;
  header.tbase = Tsft;
  header.firstfreqindex = (UINT4)( sft->f0 / sft->deltaF );
  header.nsamples=sft->data->length - 1;		/* following makefakedata_v2, but WHY -1 ?? */

  /* write header */
  if (fwrite( (void*)&header, sizeof(header), 1, fp) != 1) {
    ABORT (stat, PULSARSIGNALH_ESYS, PULSARSIGNALH_MSGESYS);
  }

  for (i=0; i < header.nsamples; i++)
    {
      rpw = sft->data->data[i].re;
      ipw = sft->data->data[i].im;

      if (fwrite((void*)&rpw, sizeof(REAL4), 1, fp) != 1) { 
	ABORT (stat, PULSARSIGNALH_ESYS, PULSARSIGNALH_MSGESYS); 
      }
      if (fwrite((void*)&ipw, sizeof(REAL4),1,fp) != 1) {  
	ABORT (stat, PULSARSIGNALH_ESYS, PULSARSIGNALH_MSGESYS); 
      }
        
    } /* for i < nsamples */
  
  fclose(fp);

  DETATCHSTATUSPTR( stat );
  RETURN (stat);
  
} /* write_SFT() */

REAL4 mymax (REAL4 x, REAL4 y)
{
  return (x > y ? x : y);
}
/* little debug-function: compare two sft's and print relative errors */
void
compare_SFTs (const SFTtype *sft1, const SFTtype *sft2)
{
  static LALStatus status;
  UINT4 i;
  REAL4 errpow= 0, errph = 0;
  REAL4 re1, re2, im1, im2, pow1, pow2, phase1, phase2;
  REAL8 Tdiff;

  if (sft1->data->length != sft2->data->length) 
    {
      printf ("\ncompare_SFTs(): lengths differ! %d != %d\n", sft1->data->length, sft2->data->length);
      return;
    }
  LALDeltaFloatGPS (&status, &Tdiff, &(sft1->epoch), &(sft2->epoch));
  if ( Tdiff != 0.0 ) 
    printf ("epochs differ: (%d s, %d ns)  vs (%d s, %d ns)\n", 
	    sft1->epoch.gpsSeconds, sft1->epoch.gpsNanoSeconds, sft2->epoch.gpsSeconds, sft2->epoch.gpsNanoSeconds);

  if ( sft1->f0 != sft2->f0)
    printf ("fmin differ: %fHz vs %fHz\n", sft1->f0, sft2->f0);

  if ( sft1->deltaF != sft2->deltaF )
    printf ("deltaF differs: %fHz vs %fHz\n", sft1->deltaF, sft2->deltaF);

  for (i=0; i < sft1->data->length; i++)
    {
      re1 = sft1->data->data[i].re;
      im1 = sft1->data->data[i].im;
      re2 = sft2->data->data[i].re;
      im2 = sft2->data->data[i].im;

      pow1 = sqrt(re1*re1 + im1*im1);
      pow2 = sqrt(re2*re2 + im2*im2);

      phase1 = atan2 (im1, re1);
      phase2 = atan2 (im2, re2);

      errpow = mymax (errpow, (pow1-pow2)/pow1 );
      errph = mymax (errph, phase1 - phase2 );

    } /* for i */


  printf ("maximal relative error in power: %e, max phaseerror: %e radians\n", errpow, errph);

  return;

} /* compare_SFTs() */



/* create a whole vector of <numSFT> SFTs with <SFTlen> frequency-bins */
void
LALCreateSFTVector (LALStatus *stat, SFTVector **output, UINT4 numSFTs, UINT4 SFTlen)
{
  UINT4 iSFT, j;
  SFTVector *vect;	/* vector to be returned */
  COMPLEX8Vector *data = NULL;

  INITSTATUS( stat, "CreateSFTVector", PULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );

  ASSERT (output != NULL, stat, PULSARSIGNALH_ENULL,  PULSARSIGNALH_MSGENULL);
  ASSERT (*output == NULL, stat, PULSARSIGNALH_ENONULL,  PULSARSIGNALH_MSGENONULL);

  vect = LALCalloc (1, sizeof(*vect) );
  if (vect == NULL) {
    ABORT (stat, PULSARSIGNALH_EMEM, PULSARSIGNALH_MSGEMEM);
  }

  vect->length = numSFTs;

  vect->data = LALCalloc (1, numSFTs * sizeof ( *vect->data ) );
  if (vect->data == NULL) {
    LALFree (vect);
    ABORT (stat, PULSARSIGNALH_EMEM, PULSARSIGNALH_MSGEMEM);
  }

  for (iSFT=0; iSFT < numSFTs; iSFT ++)
    {
      LALCCreateVector (stat->statusPtr, &data , SFTlen);
      BEGINFAIL (stat) { /* crap, we have to de-allocate as far as we got so far... */
	for (j=0; j<iSFT; j++)
	  LALCDestroyVector (stat->statusPtr, (COMPLEX8Vector**)&(vect->data[j].data) );
	LALFree (vect->data);
	LALFree (vect);
      } ENDFAIL (stat);

      vect->data[iSFT].data = data;
      data = NULL;

    } /* for iSFT < numSFTs */

  *output = vect;

  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* LALCreateSFTVector() */


void
LALDestroySFTVector (LALStatus *stat, SFTVector **vect)
{
  UINT4 i;

  INITSTATUS( stat, "DestroySFTVector", PULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );

  ASSERT (vect != NULL, stat, PULSARSIGNALH_ENULL,  PULSARSIGNALH_MSGENULL);
  ASSERT (*vect != NULL, stat, PULSARSIGNALH_ENULL,  PULSARSIGNALH_MSGENULL);

  
  for (i=0; i < (*vect)->length; i++)
    LALCDestroyVector (stat->statusPtr, &((*vect)->data[i].data) );

  LALFree ( (*vect)->data );
  LALFree ( *vect );

  *vect = NULL;

  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* LALDestroySFTVector() */


void
LALDestroyTimestampVector (LALStatus *stat, LIGOTimeGPSVector **vect)
{
  UINT4 i;

  INITSTATUS( stat, "DestroyTimestampVector", PULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );

  ASSERT (vect != NULL, stat, PULSARSIGNALH_ENULL,  PULSARSIGNALH_MSGENULL);
  ASSERT (*vect != NULL, stat, PULSARSIGNALH_ENULL,  PULSARSIGNALH_MSGENULL);

  for (i=0; i < (*vect)->length; i++)
    LALFree ( (*vect)->data);

  LALFree ( *vect );
  
  *vect = NULL;

  DETATCHSTATUSPTR( stat );
  RETURN (stat);
  
} /* LALDestroyTimestampVector() */




/***********************************************************************
 * given a start-time, duration and Tsft, returns a list of timestamps
 * covering this time-stretch.
 * returns NULL on out-of-memory
 ***********************************************************************/
LIGOTimeGPSVector* 
make_timestamps (const LIGOTimeGPS *tStart, REAL8 duration, REAL8 Tsft)
{
  LALStatus status;
  UINT4 i;
  UINT4 numSFTs;
  LIGOTimeGPS tt;
  LIGOTimeGPSVector *timestamps = NULL;

  numSFTs = (UINT4)( duration / Tsft );			/* floor */
  timestamps = LALCalloc (1, sizeof( *timestamps) );

  timestamps->length = numSFTs;
  if ( (timestamps->data = LALCalloc (1, numSFTs * sizeof (*timestamps->data) )) == NULL) 
    return (NULL);

  tt = *tStart;	/* initialize to start-time */
  for (i = 0; i < numSFTs; i++)
    {
      timestamps->data[i] = tt;
      /* get next time-stamp */
      /* NOTE: we add the interval Tsft successively instead of
       * via iSFT*Tsft, because this way we avoid possible ns-rounding problems
       * with a REAL8 interval, which becomes critial from about 100days on...
       */
      LALAddFloatToGPS( &status, &tt, &tt, Tsft);
      if (status.statusCode)	/* error */
	return (NULL);

    } /* for i < numSFTs */
  
  return (timestamps);
  
} /* make_timestamps() */

/* check that all timestamps given lie within the range [t0, t1] 
   return: 0 if ok, -1 if not
*/
int
check_timestamp_bounds (const LIGOTimeGPSVector *timestamps, LIGOTimeGPS t0, LIGOTimeGPS t1)
{
  static LALStatus status1, status2;
  REAL8 diff0, diff1;
  UINT4 i;

  for (i = 0; i < timestamps->length; i ++)
    {
      LALDeltaFloatGPS (&status1, &diff0, &(timestamps->data[i]), &t0);
      LALDeltaFloatGPS (&status2, &diff1, &t1, &(timestamps->data[i]));

      /* catch errors */
      if (status1.statusCode || status2.statusCode)
	return (-2);

      /* check timestamps-bounds */
      if ( (diff0 < 0) || (diff1 < 0) ) 
	return (-1);

    } /* for i < numSFTs */

  return (0);

} /* check_timestamp_bounds */

/* write an SFT in xmgrace-readable format */
void 
LALwriteSFTtoXMGR (LALStatus *stat, const SFTtype *sft, const CHAR *fname)
{
  FILE *fp;

  REAL4 val;
  UINT4 i;
  REAL8 Tsft, freqBand;
  REAL8 f0, df, ff;
  UINT4 set, nsamples;

  CHAR *xmgrHeader = 
    "@version 50103\n"
    "@xaxis label \"f (Hz)\"\n";

  INITSTATUS( stat, "LALwriteSFTtoXMGR", PULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );

  fp = fopen (fname, "w");

  if( !fp ) {
    ABORT (stat, PULSARSIGNALH_ESYS, PULSARSIGNALH_MSGESYS);
  }


  f0 = sft->f0;
  df = sft->deltaF;
  nsamples = sft->data->length;
  Tsft = 1.0 / sft->deltaF;
  freqBand = nsamples * df;

  fprintf (fp, xmgrHeader);
  fprintf (fp, "@subtitle \"epoch = (%d s, %d ns), Tsft = %f\"\n", 
	   sft->epoch.gpsSeconds, sft->epoch.gpsNanoSeconds, Tsft);

  set = 0;
  /* Print set header. */
  fprintf( fp, "\n@target G0.S%d\n@type xy\n", set);
  fprintf( fp, "@s%d color (0,0,1)\n", set );
  for (i=0; i < nsamples; i++)
    {
      ff = f0 + i*df;
      val = sft->data->data[i].re;

      fprintf(fp, "%f %e\n", ff, val);
        
    } /* for i < nsamples */

  set ++;
  /* Print set header. */
  fprintf( fp, "\n@target G0.S%d\n@type xy\n", set);
  fprintf( fp, "@s%d color (0,1,0)\n", set );
  for (i=0; i < nsamples; i++)
    {
      ff = f0 + i*df;
      val = sft->data->data[i].im;

      fprintf(fp, "%f %e\n", ff, val);
        
    } /* for i < nsamples */

  set ++;
  /* Print set header. */
  fprintf( fp, "\n@target G0.S%d\n@type xy\n", set);
  fprintf( fp, "@s%d color (1,0,0)\n", set );
  for (i=0; i < nsamples; i++)
    {
      ff = f0 + i*df;
      val = (sft->data->data[i].re * sft->data->data[i].re + sft->data->data[i].im*sft->data->data[i].im) / freqBand;
      
      fprintf(fp, "%f %e\n", ff, val);
        
    } /* for i < nsamples */

  
  fclose(fp);

  DETATCHSTATUSPTR( stat );
  RETURN (stat);
  
} /* write_SFT() */

/* if sky-position is not in "normal-range", normalize it correspondingly */
/* based on Alicia's function with some additional "unwinding" added  */
void
LALNormalizeSkyPosition (LALStatus *stat, SkyPosition *posOut, const SkyPosition *posIn)
{
  SkyPosition tmp;	/* allow posOut == posIn */

  INITSTATUS( stat, "NormalizeSkyPosition", PULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );
  
  ASSERT (posIn, stat, PULSARSIGNALH_ENULL ,  PULSARSIGNALH_MSGENULL );
  ASSERT (posOut, stat, PULSARSIGNALH_ENULL ,  PULSARSIGNALH_MSGENULL );

  tmp = *posIn;
  
  /* FIRST STEP: completely "unwind" positions, i.e. make sure that 
   * [0 <= alpha < 2pi] and [-pi < delta <= pi] */
  /* normalize longitude */
  while (tmp.longitude < 0)
    tmp.longitude += LAL_TWOPI;
  while (tmp.longitude >= LAL_TWOPI)
    tmp.longitude -= LAL_TWOPI;

  /* pre-normalize (unwind) latitude */
  while (tmp.latitude <= -LAL_PI)
    tmp.latitude += LAL_TWOPI;
  while (tmp.latitude > LAL_TWOPI)
    tmp.latitude -= LAL_TWOPI;

  /* SECOND STEP: get latitude into canonical interval [-pi/2 <= delta <= pi/2 ] */
  /* this requires also a change in longitude by adding/subtracting PI */
  if (tmp.latitude > LAL_PI_2)
    {
      tmp.latitude = LAL_PI - tmp.latitude;
      if (tmp.longitude < LAL_PI)
	{
	  tmp.longitude += LAL_PI;
	}
      else
	{
	  tmp.longitude -= LAL_PI;
	}
    }

  if (tmp.latitude < -LAL_PI_2)
    {
      tmp.latitude = -LAL_PI - tmp.latitude;
      if (tmp.longitude <= LAL_PI)
	{
	  tmp.longitude += LAL_PI;
	}
      else
	{
	  tmp.longitude -= LAL_PI;
	}
    }

  *posOut = tmp;

  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* LALNormalizeSkyPosition() */

/*----------------------------------------------------------------------
 * check if frequency-range and resolution of noiseSFTs is consistent with signal
 * ABORT if not
 *----------------------------------------------------------------------*/
void
checkNoiseSFTs (LALStatus *stat, const SFTVector *sfts, REAL8 f0, REAL8 f1, REAL8 deltaF)
{
  UINT4 i;
  SFTtype *thisSFT;
  REAL8 fn0, fn1, deltaFn, shift;
  UINT4 nshift;

  INITSTATUS( stat, "checkNoiseSFTs", PULSARSIGNALC);

  for (i=0; i < sfts->length; i++)
    {
      thisSFT = &(sfts->data[i]);
      deltaFn = thisSFT->deltaF;
      fn0 = thisSFT->f0;
      fn1 = f0 + thisSFT->data->length * deltaFn;
      
      if (deltaFn != deltaF) {
	ABORT (stat,  PULSARSIGNALH_ENOISEDELTAF,  PULSARSIGNALH_MSGENOISEDELTAF);
      }

      if ( (f0 < fn0) || (f1 > fn1) ) {
	ABORT (stat, PULSARSIGNALH_ENOISEBAND, PULSARSIGNALH_MSGENOISEBAND);
      }
      
      shift = (fn0 - f0) / deltaF;
      /* frequency bins have to coincide! ==> check that shift is integer!  */
      nshift = (UINT4)(shift+0.5);
      if ( fabs( nshift - shift) > eps ) {
	ABORT (stat, PULSARSIGNALH_ENOISEBINS, PULSARSIGNALH_MSGENOISEBINS);
      }

    } /* for i < numSFTs */

  RETURN (stat);

} /* checkNoiseSFTs() */

/* dump an SFT into a text-file */
void 
dump_SFT (LALStatus *stat, const SFTtype *sft, const CHAR *fname)
{
  FILE *fp;

  REAL4 valre, valim;
  UINT4 i;
  REAL8 Tsft, freqBand;
  REAL8 f0, df, ff;
  UINT4 nsamples;

  INITSTATUS( stat, "dump_SFT", PULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );

  fp = fopen (fname, "w");

  if( !fp ) {
    ABORT (stat, PULSARSIGNALH_ESYS, PULSARSIGNALH_MSGESYS);
  }


  f0 = sft->f0;
  df = sft->deltaF;
  nsamples = sft->data->length;
  Tsft = 1.0 / sft->deltaF;
  freqBand = nsamples * df;

  for (i=0; i < nsamples; i++)
    {
      ff = f0 + i*df;
      valre = sft->data->data[i].re;
      valim = sft->data->data[i].im;
      fprintf(fp, "%f %e %e\n", ff, valre, valim);
        
    } /* for i < nsamples */
  
  fclose(fp);

  DETATCHSTATUSPTR( stat );
  RETURN (stat);
  
} /* write_SFT() */
