/************************************ <lalVerbatim file="GeneratePulsarSignalCV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{GeneratePulsarSignal}
\label{ss:GeneratePulsarSignal.c}

Module to generate subsequent sky-positions.

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
static void write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series, UINT4 set);
static void write_timeSeriesR8 (FILE *fp, const REAL8TimeSeries *series, UINT4 set);
static LIGOTimeGPSVector* make_timestamps (const LIGOTimeGPS *tStart, REAL8 duration, REAL8 Tsft);
static int check_timestamp_bounds (const LIGOTimeGPSVector *timestamps, LIGOTimeGPS t0, LIGOTimeGPS t1);
/*----------------------------------------------------------------------*/

NRCSID( PULSARSIGNALC, "$Id$");

extern INT4 lalDebugLevel;

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
LALGeneratePulsarSignal (LALStatus *stat, REAL4TimeSeries *signal, const PulsarSignalParams *params)
{
  static SpinOrbitCWParamStruc sourceParams;
  static CoherentGW sourceSignal;
  DetectorResponse detector;
  UINT4 SSBduration;
  LIGOTimeGPS time = {0,0};

  INITSTATUS( stat, "LALGeneratePulsarSignal", PULSARSIGNALC );
  ATTATCHSTATUSPTR (stat);

  ASSERT (signal != NULL, stat, PULSARSIGNALH_ENULL, PULSARSIGNALH_MSGENULL);
  ASSERT (signal->data == NULL, stat,  PULSARSIGNALH_ENONULL,  PULSARSIGNALH_MSGENONULL);

  /*----------------------------------------------------------------------
   *
   * First call GenerateSpinOrbitCW() to generate the source-signal
   *
   *----------------------------------------------------------------------*/
  sourceParams.position.system = COORDINATESYSTEM_EQUATORIAL;
  sourceParams.position = params->pulsar.position;
  sourceParams.psi = params->pulsar.psi;
  sourceParams.aPlus = params->pulsar.aPlus;
  sourceParams.aCross = params->pulsar.aCross;
  sourceParams.phi0 = params->pulsar.phi0;
  sourceParams.f0 = params->pulsar.f0;

  /* we use frequency-spindowns, but GenerateSpinOrbitCW wants it f0-normalized,
     so we have to do that here: */
  if (params->pulsar.spindown)
    {
      UINT4 i;
      TRY ( LALDCreateVector (stat->statusPtr, &(sourceParams.f), params->pulsar.spindown->length), stat);
      for (i=0; i < sourceParams.f->length; i++)
	sourceParams.f->data[i] = params->pulsar.spindown->data[i] / params->pulsar.f0;
    }

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

  /*
   * finally, call the function to generate the source waveform 
   */
  TRY ( LALGenerateSpinOrbitCW (stat->statusPtr, &sourceSignal, &sourceParams), stat);
  /* check that sampling interval was short enough */
  if ( sourceParams.dfdt > 2.0 )  /* taken from makefakedata_v2 */
    {
      LALPrintError ("GenerateSpinOrbitCW() returned df*dt = %f > 2.0", sourceParams.dfdt);
      ABORT (stat, PULSARSIGNALH_ESAMPLING, PULSARSIGNALH_MSGESAMPLING);
    }

  TRY (PrintGWSignal (stat->statusPtr, &sourceSignal, "signal2.agr"), stat);

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
  
  /* ok, but we also need to prepare the output time-series */
  signal->data = LALMalloc (sizeof(*(signal->data)));
  signal->data->length = (UINT4)( params->samplingRate * params->duration );
  signal->data->data = LALMalloc ( signal->data->length * (sizeof(*(signal->data->data))) );
  signal->deltaT = 1.0 / params->samplingRate;
  signal->f0 = params->fHeterodyne;
  signal->epoch = params->startTimeGPS;


  printf ("\nLaenge time-series: %d samples\n", signal->data->length);
  
  TRY ( LALSimulateCoherentGW (stat->statusPtr, signal, &sourceSignal, &detector ), stat );

			       
  /*----------------------------------------------------------------------*/
  /* Free all allocated memory that is not returned */
  TRY (LALSDestroyVectorSequence ( stat->statusPtr, &(sourceSignal.a->data)), stat);
  LALFree ( sourceSignal.a );
  TRY (LALSDestroyVector (stat->statusPtr, &(sourceSignal.f->data ) ), stat);
  LALFree (sourceSignal.f);
  TRY (LALDDestroyVector (stat->statusPtr, &(sourceSignal.phi->data )), stat);
  LALFree (sourceSignal.phi);

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
  UINT4 numSamples;			/* number of samples in each Tsft */
  UINT4 iSFT;
  REAL8 Band, samples;
  RealFFTPlan *pfwd = NULL;
  LIGOTimeGPSVector *timestamps = NULL;
  REAL4Vector timeStretch = {0,0};
  LIGOTimeGPS tStart;			/* start time of input time-series */
  LIGOTimeGPS tLast;			/* start-time of last _possible_ SFT */
  REAL8 duration, delay;
  UINT4 SFTlen;				/* number of samples in an SFT */
  UINT4 indexShift;
  const REAL8 eps = 1.e-14;		
  SFTtype *thisSFT;
  SFTVector *sftvect = NULL;		/* return value. For better readability */

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

  Band = 1.0 / (2.0 * signal->deltaT);		/* NOTE: frequency-band is determined by sampling-rate! */
  samples = 2.0 * Band * params->Tsft;		/* this is a float!*/
  numSamples = (UINT4) (samples + 0.5);		/* round to int */
  /* now make sure that number of samples is an integer (up to possible rounding errors */
  ASSERT ( abs(samples - numSamples)/samples < eps, stat, 
	   PULSARSIGNALH_EINCONSBAND, PULSARSIGNALH_MSGEINCONSBAND);
  
  /* Prepare FFT: compute plan for FFTW */
  /* FIXME: put some buffering */
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

      /* now prepare the i'th output SFT with the result */
      thisSFT->epoch = timestamps->data[iSFT];	/* set the proper timestamp */
      thisSFT->f0 = signal->f0;			/* minimum frequency */
      thisSFT->deltaF = 1.0 / params->Tsft;	/* frequency-spacing */

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
      printf ("\nARGH: non-equatorial coords not implemented here yet!\n");
      exit(-1);	/* suicide */
    }
  baryinput.alpha = params->pulsar.position.longitude;
  baryinput.delta = params->pulsar.position.latitude;
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
LALPrintR4TimeSeries (LALStatus *stat, const REAL4TimeSeries *series, const CHAR *fname)
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

  set = 0;

  write_timeSeriesR4 (fp, series, set);
      
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

  set = 0;
  write_timeSeriesR4 (fp, signal->f, set++);
  write_timeSeriesR8 (fp, signal->phi, set++);

  fclose (fp);

  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* PrintGWSignal() */


/* internal helper-function */
void
write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series, UINT4 set)
{
  REAL8 timestamp; 
  UINT4 i;


  if (series == NULL)
    {
      printf ("\nset %d is empty\n", set);
      return; 
    }

  /* Print set header. */
  fprintf( fp, "\n@target G0.S%d\n@type xy\n", set);
  fprintf( fp, "@s%d color (0,0,0)\n", set );

  timestamp = 1.0*series->epoch.gpsSeconds + series->epoch.gpsNanoSeconds * 1.0e-9;

  for( i = 0; i < series->data->length; i++)
    {
      fprintf( fp, "%16.9f %e\n", timestamp, series->data->data[i] );
      timestamp += series->deltaT;
    }

  return;

} /* write_timeSeriesR4() */

void
write_timeSeriesR8 (FILE *fp, const REAL8TimeSeries *series, UINT4 set)
{
  REAL8 timestamp; 
  UINT4 i;
  /* Print set header. */
  fprintf( fp, "\n@target G0.S%d\n@type xy\n", set);
  fprintf( fp, "@s%d color (0,0,0)\n", set );

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
  REAL8 freqBand;

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

  freqBand = (sft->data->length -1 ) * sft->deltaF;
  header.tbase = (REAL4) ( (sft->data->length - 1.0) / freqBand);
  header.firstfreqindex = (INT4)(sft->f0 * header.tbase + 0.5);
  header.nsamples=sft->data->length - 1;

  printf ("\nwrite_SFT():\nfreqBand = %f\n", freqBand);
  printf ("f0 = %f\ndf = %f\nnsamples = %d\n", sft->f0, sft->deltaF, header.nsamples);
  printf ("Tsft = %f\n", header.tbase);  

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

/* little debug-function: compare two sft's and print relative errors */
void
compare_SFTs (const SFTtype *sft1, const SFTtype *sft2)
{
  UINT4 i;
  REAL4 maxdiff = 0;

  if (sft1->data->length != sft2->data->length) 
    {
      printf ("\ncompare_SFTs(): lengths differ! %d != %d\n", sft1->data->length, sft2->data->length);
      return;
    }
#define MAX(x,y) (x > y ? x : y)

  for (i=0; i < sft1->data->length; i++)
    {
      maxdiff = MAX ( maxdiff, abs ( 2.0*(sft1->data->data[i].re - sft2->data->data[i].re) 
				     / (sft1->data->data[i].re + sft2->data->data[i].re) ) );
      maxdiff = MAX ( maxdiff, abs ( 2.0*(sft1->data->data[i].im - sft2->data->data[i].im) 
				     / (sft1->data->data[i].im + sft2->data->data[i].im) ) );
    }

  printf ("\ncompare_SFTs(): maximal relative difference: %f\n", maxdiff);

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
  /*  fprintf (fp, "@title \"GW source signal\"\n"); */
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
