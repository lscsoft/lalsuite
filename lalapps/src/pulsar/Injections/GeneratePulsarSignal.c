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
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

#include "GeneratePulsarSignal.h"

/*----------------------------------------------------------------------*/
/* prototypes for internal functions */
static void write_timeSeriesR4 (FILE *fp, REAL4TimeSeries *series, UINT4 set);
static void write_timeSeriesR8 (FILE *fp, REAL8TimeSeries *series, UINT4 set);
/*----------------------------------------------------------------------*/

NRCSID( PULSARSIGNALC, "$Id$");

extern INT4 lalDebugLevel;

void
LALGeneratePulsarSignal (LALStatus *stat, REAL4TimeSeries *signal, PulsarSignalParams *params)
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
  sourceParams.position.longitude = params->pulsar.Alpha;
  sourceParams.position.latitude =  params->pulsar.Delta;
  sourceParams.psi = params->pulsar.psi;
  sourceParams.aPlus = params->pulsar.aPlus;
  sourceParams.aCross = params->pulsar.aCross;
  sourceParams.phi0 = params->pulsar.phi0;
  sourceParams.f0 = params->pulsar.f0;
  sourceParams.f = params->pulsar.f;

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

  /* start-time in SSB time */
  TRY (ConvertGPS2SSB (stat->statusPtr, &time, params->startTimeGPS, params), stat);
  sourceParams.epoch = time;

  /* pulsar reference-time in SSB frame !*/
  if (params->pulsar.TRefSSB.gpsSeconds != 0)
    sourceParams.spinEpoch = params->pulsar.TRefSSB;
  else
    sourceParams.spinEpoch = time;	/* use start-time */

  /* sampling-timestep and length for source-parameters */
  sourceParams.deltaT = 60;	/* in seconds; hardcoded default from makefakedata_v2 */
  time = params->startTimeGPS;
  time.gpsSeconds += params->duration;
  TRY (ConvertGPS2SSB (stat->statusPtr, &time, time, params), stat);	 /* convert time to SSB */
  SSBduration = time.gpsSeconds - sourceParams.spinEpoch.gpsSeconds +2 ;
  sourceParams.length = (UINT4)( 1.0* SSBduration / sourceParams.deltaT ) + 1;

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
 * convert earth-frame GPS time into barycentric-frame SSB time for given source
 *----------------------------------------------------------------------*/
void
ConvertGPS2SSB (LALStatus* stat, LIGOTimeGPS *SSBout, LIGOTimeGPS GPSin, PulsarSignalParams *params)
{

  LIGOTimeGPS ssb;
  EarthState earth;
  EmissionTime emit;
  BarycenterInput baryinput;
  /* just here for checking that we get the same as in makefakedata_v2 (remove eventually) */
  REAL8 doubleTime;
  REAL8 Ts=GPSin.gpsSeconds;
  REAL8 Tns=GPSin.gpsNanoSeconds;


  INITSTATUS( stat, "ConvertGPS2SSB", PULSARSIGNALC );
  ATTATCHSTATUSPTR (stat);

  ASSERT (SSBout != NULL, stat,  PULSARSIGNALH_ENULL,  PULSARSIGNALH_MSGENULL);
  ASSERT (params != NULL, stat,  PULSARSIGNALH_ENULL,  PULSARSIGNALH_MSGENULL);

  baryinput.tgps = GPSin;
  baryinput.site = *(params->site);
  /* account for a quirk in LALBarycenter(): -> see documentation of type BarycenterInput */
  baryinput.site.location[0] /= LAL_C_SI;
  baryinput.site.location[1] /= LAL_C_SI;
  baryinput.site.location[2] /= LAL_C_SI;
  baryinput.alpha = params->pulsar.Alpha;
  baryinput.delta = params->pulsar.Delta;
  baryinput.dInv = 0.e0;	/* following makefakedata_v2 */

  TRY (LALBarycenterEarth(stat->statusPtr, &earth, &GPSin, params->ephemerides), stat);
  TRY (LALBarycenter(stat->statusPtr, &emit, &baryinput, &earth), stat);

  /* RP: this was used in makefakedata_v2: check that we get the same! */
  doubleTime= emit.deltaT + Ts + Tns*1.E-9;
  TRY (LALFloatToGPS (stat->statusPtr, &ssb, &doubleTime), stat);

  printf ("\nDEBUG: checking  ConvertGPS2SSB(): difference = %e s\n", 
	  1.0*(ssb.gpsSeconds - emit.te.gpsSeconds)+ (ssb.gpsNanoSeconds - emit.te.gpsNanoSeconds)*1.0e-9);

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
ConvertSSB2GPS (LALStatus *stat, LIGOTimeGPS *GPSout, LIGOTimeGPS SSBin, PulsarSignalParams *params)
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
      
      /* break if we've converged: +/- 2 NanoSeconds seems good enough!  */
      if (delta > -2 && delta < 2)
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
LALPrintR4TimeSeries (LALStatus *stat, REAL4TimeSeries *series, const CHAR *fname)
{
  FILE *fp = NULL;
  UINT4 set;
  CHAR *xmgrHeader = 
    "@version 50103\n"
    "@xaxis label \"T (sec)\"\n"
    "@yaxis label \"A\"\n";

  /*
    "@world xmin -0.1\n"
    "@world xmax 6.4\n"
    "@world ymin -3.2\n"
    "@world ymax 3.2\n"
  */
  /*
  CHAR           name[LALNameLength];
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  LALUnit        sampleUnits;
  REAL4Sequence *data;
  */

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


void
PrintGWSignal (LALStatus *stat, CoherentGW *signal, const CHAR *fname)
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
write_timeSeriesR4 (FILE *fp, REAL4TimeSeries *series, UINT4 set)
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
write_timeSeriesR8 (FILE *fp, REAL8TimeSeries *series, UINT4 set)
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

/*----------------------------------------------------------------------
 * take time-series as input, convert it into SFTs and add noise if given
 *----------------------------------------------------------------------*/
void
AddSignalToSFTs (LALStatus *stat, COMPLEX8Vector **outputSFTs, REAL4TimeSeries *signal, SFTParams *params)
{
  UINT4 i, numSFTs, i0, iSFT;
  REAL8 duration;

  INITSTATUS( stat, "AddSignalToSFTs", PULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );
  
  ASSERT (outputSFTs != NULL, stat, PULSARSIGNALH_ENULL, PULSARSIGNALH_MSGENULL);
  ASSERT (signal != NULL, stat, PULSARSIGNALH_ENULL, PULSARSIGNALH_MSGENULL);
  ASSERT (params != NULL, stat, PULSARSIGNALH_ENULL, PULSARSIGNALH_MSGENULL);
  ASSERT (*outputSFTs == NULL, stat,  PULSARSIGNALH_ENONULL,  PULSARSIGNALH_MSGENONULL);

  /* determine number of SFT to be produced */
  duration = 1.0* signal->data->length * signal->deltaT;
  if (params->timestamps == NULL)
    numSFTs = (UINT4)( duration / params->Tsft );
  else
    {
      /* FIXME: stopped writing here!! */
      for (i=0; i<10; i++)
	;
    }


  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* AddSignalToSFTs() */
