/************************************ <lalVerbatim file="GeneratePulsarSignalCV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{GeneratePulsarSignal.c}}
\label{ss:GeneratePulsarSignal.c}

Module to generate simulated pulsar-signals.

\subsubsection*{Prototypes}
\idx{LALGeneratePulsarSignal()}
\idx{LALSignalToSFTs()}
\idx{LALConvertSSB2GPS()}
\idx{LALConvertGPS2SSB()}
\idx{LALCreateSFTtype()}
\idx{LALCreateSFTVector()}
\idx{LALDestroySFTtype()}
\idx{LALDestroySFTVector()}
\idx{LALDestroyTimestampVector()}

\input{GeneratePulsarSignalCP}

\newpage
\subsubsection*{Description}

\begin{itemize}
\item The main function \verb+LALGeneratePulsarSignal()+ generates a fake
pulsar-signal, either for an isolated or a binary pulsar. It returns a
time-series with the generated signal as received by the detector. 

\item The time-series can then be turned into a vector of short time-base FFTs
(so-called "SFTs") by using the function \verb+LALSignalToSFTs()+.
These SFTs are the data-format used by most frequency-domain pulsar codes,
therefore these functions can be directly used in a Monte-Carlo
injection driver. 
\end{itemize}


This module also contains a few more general-purpose helper-functions:

\begin{itemize}
\item Namely, \verb+LALConvertSSB2GPS()+ and \verb+LALConvertGPS2SSB()+
which convert arrival times for a given source (not necessarily a
pulsar!) the detector ("GPS") and the solar-system barycenter ("SSB"). 
NOTE: only the source-location (\verb+params->pulsar.position+), the
detector-site (\verb+params->site+) and the ephemeris-data
(\verb+params->ephemerides+)are used from the
\verb+PulsarSignalParams+-structure.  

\item The helper functions \verb+LALCreateSFTtype()+,
\verb+LALDestroySFTtype()+, \verb+LALCreateSFTVector()+
and\verb+LALDestroySFTVector()+  respectively allocate and free an
SFT-structs and SFT-vectors. 
Similarly, \verb+LALCreateTimestampVector()+ and
\verb+LALDestroyTimestampVector()+ allocate and free a bunch of  
GPS-timestamps respectively.

\subsubsection*{Algorithm}

\verb+LALGeneratePulsarSignal()+ is basically a wrapper for the two
LAL-functions \verb+GenerateSpinOrbitCW()+
(cf.~\ref{ss:GenerateSpinOrbitCW.c}) to produce the source-signal, and
\verb+LALSimulateCoherentGW()+ (cf.~\ref{ss:SimulateCoherentGW.c}) to
turn it into a time-series at the detector.


\verb+LALSignalToSFTs()+ uses \verb+LALForwardRealFFT()+
(cf.~\ref{ss:RealFFT.c}) appropriately on the input-timeseries to 
produce the required output-SFTs. 

\subsubsection*{Uses}
\begin{verbatim}
LALGenerateSpinOrbitCW()        LALSimulateCoherentGW()
LALBarycenterEarth()            LALBarycenter()
LALAddFloatToGPS()              LALDeltaFloatGPS()
LALCreateForwardRealFFTPlan()   LALDestroyRealFFTPlan()
LALForwardRealFFT()

LALCalloc()                     LALFree()
LALDCreateVector()              LALDDestroyVector()
LALCCreateVector()              LALCDestroyVector()
LALCreateVector()               LALSDestroyVector()
LALSDestroyVectorSequence()     LALPrintError()
\end{verbatim}

\subsubsection*{Notes}

\verb+LALSignalToSFTs()+ currently enforces the constraint of 
\verb+Tsft * Band+ being an integer, so that the number of
time-samples per SFT is even. This follows \verb+makefakedata_v2+.


Furthermore, if the timestamps for SFT-creation passed to
\verb+LALSignalToSFTs()+ do not fit exactly on a time-step of the
input time-series, it will be "nudged" to the closest one.
If \verb+lalDebugLevel>0+, a warning will be printed about this. 
The user could also see this effect in the actual timestamps of the
SFTs returned.


The FFTW-"plan" is currently created using the \verb+ESTIMATE+ flag,
which is fast but only yields an approximate plan. Better results
might be achieved by using \verb+MEASURE+ and an appropriate buffering
of the resulting plan (which doesnt change provided the SFT-length is
the same). Measuring the plan takes longer but might lead to
substantial speedups for many FFTs, which seems the most likely situation.

\vfill{\footnotesize\input{GeneratePulsarSignalCV}}

******************************************************* </lalLaTeX> */
#include <math.h>

#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>

#include <lal/GeneratePulsarSignal.h>

/*----------------------------------------------------------------------*/
static LIGOTimeGPSVector* make_timestamps (const LIGOTimeGPS *tStart, REAL8 duration, REAL8 Tsft);
static int check_timestamp_bounds (const LIGOTimeGPSVector *timestamps, LIGOTimeGPS t0, LIGOTimeGPS t1);
static void checkNoiseSFTs (LALStatus *stat, const SFTVector *sfts, REAL8 f0, REAL8 f1, REAL8 deltaF);
static void correct_phase (LALStatus* stat, SFTtype *sft, LIGOTimeGPS tHeterodyne);
/*----------------------------------------------------------------------*/

NRCSID( GENERATEPULSARSIGNALC, "$Id$");

extern INT4 lalDebugLevel;

static REAL8 eps = 1.e-14;	/* maximal REAL8 roundoff-error (used for determining if some number is an INT) */

/* some empty structs for initializing */
static LALStatus emptyStatus;	
static SpinOrbitCWParamStruc emptyCWParams;
static CoherentGW emptySignal;

/***********************************************************************
 * generate a time-series at the detector for a given pulsar
 ***********************************************************************/
/* <lalVerbatim file="GeneratePulsarSignalCP"> */
/* --------------- the central functions of this module --------------- */
void
LALGeneratePulsarSignal (LALStatus *stat, 
			 REAL4TimeSeries **signal, 	   /* output time-series */
			 const PulsarSignalParams *params) /* input params */
{ /* </lalVerbatim> */
  SpinOrbitCWParamStruc sourceParams = emptyCWParams; 
  CoherentGW sourceSignal = emptySignal;
  DetectorResponse detector;
  REAL8 SSBduration;
  LIGOTimeGPS t0, t1, tmpTime;
  REAL4TimeSeries *output;
  UINT4 i;

  INITSTATUS( stat, "LALGeneratePulsarSignal", GENERATEPULSARSIGNALC );
  ATTATCHSTATUSPTR (stat);

  ASSERT (signal != NULL, stat, GENERATEPULSARSIGNALH_ENULL, GENERATEPULSARSIGNALH_MSGENULL);
  ASSERT (*signal == NULL, stat,  GENERATEPULSARSIGNALH_ENONULL,  GENERATEPULSARSIGNALH_MSGENONULL);

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

  if ( params->pulsar.tRef.gpsSeconds != 0)
    sourceParams.spinEpoch = params->pulsar.tRef;   /* pulsar reference-time in SSB frame (TDB) */
  else	/* if not given: use startTime converted to SSB as tRef ! */
    {
      TRY (LALConvertGPS2SSB (stat->statusPtr, &tmpTime, params->startTimeGPS, params), stat);
      sourceParams.spinEpoch = tmpTime;
    }

  /* sampling-timestep and length for source-parameters */
  sourceParams.deltaT = 60;	/* in seconds; hardcoded default from makefakedata_v2 
				 * this should be largely enough, at it only concerns
				 * the amplitudes, and _phase_ of the signal, not h(t),
				 * which will then be used for interpolation */

  /* start-time in SSB time */
  TRY (LALConvertGPS2SSB (stat->statusPtr, &t0, params->startTimeGPS, params), stat);
  t0.gpsSeconds -= (UINT4)sourceParams.deltaT; /* start one time-step earlier to be safe */


  /* end time in SSB */
  t1 = params->startTimeGPS;
  TRY ( LALAddFloatToGPS (stat->statusPtr, &t1, &t1, params->duration), stat);
  TRY (LALConvertGPS2SSB (stat->statusPtr, &t1, t1, params), stat);	 /* convert time to SSB */

  /* get duration of source-signal */
  TRY (LALDeltaFloatGPS (stat->statusPtr, &SSBduration, &t1, &t0), stat);
  SSBduration += 2.0 * sourceParams.deltaT; /* add two time-steps to be safe */

  sourceParams.epoch = t0; 
  sourceParams.length = (UINT4) ceil( SSBduration / sourceParams.deltaT );

  /* we use frequency-spindowns, but GenerateSpinOrbitCW wants it f0-normalized,
     so we have to do that here: */
  sourceParams.f = NULL;
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
      ABORT (stat, GENERATEPULSARSIGNALH_ESAMPLING, GENERATEPULSARSIGNALH_MSGESAMPLING);
    }

  /*----------------------------------------------------------------------
   *
   * Now call the function to translate the source-signal into a (heterodyned)
   * signal at the detector 
   *
   *----------------------------------------------------------------------*/
  /* first set up the detector-response */
  detector.transfer = params->transfer;
  detector.site = params->site;
  detector.ephemerides = params->ephemerides;

  /* similar to the spirit of makefakedata_v2, we use the pulsar reference-time 
   * as t_h if given, so that several time-series generated for the same pulsar would be 
   * guaranteed to have consistent heterodyning. If tRef not given, use startTime.
   *
   * However, we round this to integer-seconds to allow simplifications concerning
   * the necessary phase-correction in the final SFTs. (see docu)
   */
  if ( params->pulsar.tRef.gpsSeconds != 0)
    detector.heterodyneEpoch = params->pulsar.tRef;
  else
    detector.heterodyneEpoch = params->startTimeGPS;

  detector.heterodyneEpoch.gpsNanoSeconds = 0;	/* make sure this only has seconds */
  
  /* ok, we  need to prepare the output time-series */
  if ( (output = LALCalloc (1, sizeof (*output) )) == NULL) {
    ABORT (stat,  GENERATEPULSARSIGNALH_EMEM,  GENERATEPULSARSIGNALH_MSGEMEM);
  }
  LALCreateVector (stat->statusPtr, &(output->data), (UINT4) ceil( params->samplingRate * params->duration) );
  BEGINFAIL(stat) {
    LALFree (output);
  } ENDFAIL(stat);

  output->deltaT = 1.0 / params->samplingRate;
  output->f0 = params->fHeterodyne;
  output->epoch = params->startTimeGPS;

  
  TRY ( LALSimulateCoherentGW_exp (stat->statusPtr, output, &sourceSignal, &detector ), stat );

			       
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
/* <lalVerbatim file="GeneratePulsarSignalCP"> */
void
LALSignalToSFTs (LALStatus *stat, 
		 SFTVector **outputSFTs,	/* output: SFT-vector */
		 const REAL4TimeSeries *signal, /* input: time-series */
		 const SFTParams *params)	/* params for output-SFTs */
{ /* </lalVerbatim> */
  UINT4 numSFTs;			/* number of SFTs */
  UINT4 numSamples;			/* number of time-samples in each Tsft */
  UINT4 iSFT;
  REAL8 Band, samplesb2, f0, deltaF;
  LIGOTimeGPSVector *timestamps = NULL;
  REAL4Vector timeStretch = {0,0};
  LIGOTimeGPS tStart;			/* start time of input time-series */
  LIGOTimeGPS tLast;			/* start-time of last _possible_ SFT */
  LIGOTimeGPS tmpTime;
  REAL8 duration, delay;
  UINT4 SFTlen;				/* number of samples in an SFT */
  UINT4 indexShift;
  UINT4 index0n;			/* first frequency-bin to use from noise-SFT */
  SFTtype *thisSFT, *thisNoiseSFT;	/* SFT-pointers */
  REAL4 renorm;				/* renormalization-factor of taking only part of an SFT */
  SFTVector *sftvect = NULL;		/* return value. For better readability */
  UINT4 j;
  RealFFTPlan *pfwd = NULL;		/* FFTW plan */

  INITSTATUS( stat, "LALSignalToSFTs", GENERATEPULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );
  
  ASSERT (outputSFTs != NULL, stat, GENERATEPULSARSIGNALH_ENULL, GENERATEPULSARSIGNALH_MSGENULL);
  ASSERT (*outputSFTs == NULL, stat,  GENERATEPULSARSIGNALH_ENONULL,  GENERATEPULSARSIGNALH_MSGENONULL);
  ASSERT (signal != NULL, stat, GENERATEPULSARSIGNALH_ENULL, GENERATEPULSARSIGNALH_MSGENULL);
  ASSERT (params != NULL, stat, GENERATEPULSARSIGNALH_ENULL, GENERATEPULSARSIGNALH_MSGENULL);
  if ( params->timestamps && params->noiseSFTs) {
    ASSERT ( params->timestamps->length == params->noiseSFTs->length, stat,  
	     GENERATEPULSARSIGNALH_ENUMSFTS,  GENERATEPULSARSIGNALH_MSGENUMSFTS);
  }

  f0 = signal->f0;				/* lowest frequency */
  Band = 1.0 / (2.0 * signal->deltaT);		/* NOTE: frequency-band is determined by sampling-rate! */
  deltaF = 1.0 / params->Tsft;			/* frequency-resolution */

  /* if noiseSFTs are given: check they are consistent with signal! */
  if (params->noiseSFTs) {
    TRY (checkNoiseSFTs (stat->statusPtr, params->noiseSFTs, f0, f0 + Band, deltaF), stat);
  }
    
  /* make sure that number of samples is an integer (up to possible rounding errors */
  samplesb2 = Band * params->Tsft;		/* this is a float!*/
  numSamples = 2 * (UINT4) (samplesb2 + 0.5);	/* round to int */
  ASSERT ( fabs(samplesb2 - numSamples/2)/samplesb2 < eps, stat, 
	   GENERATEPULSARSIGNALH_EINCONSBAND, GENERATEPULSARSIGNALH_MSGEINCONSBAND);
  
  /* Prepare FFT: compute plan for FFTW */
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
	ABORT (stat, GENERATEPULSARSIGNALH_EMEM, GENERATEPULSARSIGNALH_MSGEMEM);
      }
    }
  else	/* if given, nothing to do */
    timestamps = params->timestamps;

  /* check that all timestamps lie within [tStart, tLast] */
  if ( check_timestamp_bounds (timestamps, tStart, tLast) != 0) {
    ABORT (stat, GENERATEPULSARSIGNALH_ETIMEBOUND, GENERATEPULSARSIGNALH_MSGETIMEBOUND);
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


      /* fill the header of the i'th output SFT */
      TRY( LALAddFloatToGPS (stat->statusPtr, &tmpTime, &tStart, (indexShift * signal->deltaT)), stat);	
      thisSFT->epoch = tmpTime;			/* set the ACTUAL timestamp! */
      thisSFT->f0 = signal->f0;			/* minimum frequency */
      thisSFT->deltaF = 1.0 / params->Tsft;	/* frequency-spacing */

      /* ok, issue at least a warning if we have "nudged" an SFT-timestamp */
      if (lalDebugLevel >= 3)
	{
	  REAL8 diff;
	  TRY ( LALDeltaFloatGPS (stat->statusPtr, &diff, &(timestamps->data[iSFT]), &tmpTime), stat);
	  if (diff != 0) 
	    {
	      LALPrintError ("Warning: timestamp %d had to be 'nudged' by %e s to fit with time-series\n", iSFT, diff);
	      /* double check if magnitude of nudging seems reasonable .. */
	      if ( fabs(diff) >= signal->deltaT ) 
		{
		  LALPrintError ("WARNING: nudged by more than deltaT=%e... this sounds wrong! (We better stop)\n");
		  ABORT (stat, GENERATEPULSARSIGNALH_ENULL, GENERATEPULSARSIGNALH_MSGENULL );
		}
	    } /* if nudging */
	} /* if lalDebugLevel */

      /* the central step: FFT the ith time-stretch into an SFT-slot */
      LALForwardRealFFT (stat->statusPtr, thisSFT->data, &timeStretch, pfwd);
      BEGINFAIL(stat) {
	LALDestroySFTVector (stat->statusPtr, &sftvect);
      } ENDFAIL(stat);


      /* correct heterodyning-phase, IF NECESSARY */
      if ( ( (INT4)signal->f0 != signal->f0  )
	   || (signal->epoch.gpsNanoSeconds != 0)
	   || (thisSFT->epoch.gpsNanoSeconds != 0) )
	{
	  correct_phase (stat->statusPtr, thisSFT, signal->epoch);	/* theterodyne = signal->epoch!*/
	  BEGINFAIL (stat) {
	    LALDestroySFTVector (stat->statusPtr, &sftvect);
	  } ENDFAIL (stat);
	} /* if phase-correction necessary */

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

  /* free stuff */
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
/* <lalVerbatim file="GeneratePulsarSignalCP"> */
/*--------------- some useful helper-functions ---------------*/
void
LALConvertGPS2SSB (LALStatus* stat, 
		   LIGOTimeGPS *SSBout, 	/* output: arrival-time in SSB */
		   LIGOTimeGPS GPSin, 		/* input: GPS-arrival time at detector */
		   const PulsarSignalParams *params) /* define source-location and detector */
{ /* </lalVerbatim> */
  EarthState earth;
  EmissionTime emit;
  BarycenterInput baryinput;
  SkyPosition tmp;

  INITSTATUS( stat, "ConvertGPS2SSB", GENERATEPULSARSIGNALC );
  ATTATCHSTATUSPTR (stat);

  ASSERT (SSBout != NULL, stat,  GENERATEPULSARSIGNALH_ENULL,  GENERATEPULSARSIGNALH_MSGENULL);
  ASSERT (params != NULL, stat,  GENERATEPULSARSIGNALH_ENULL,  GENERATEPULSARSIGNALH_MSGENULL);

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

} /* LALConvertGPS2SSB() */


/*----------------------------------------------------------------------
 * convert  barycentric frame SSB time into earth-frame GPS time
 *
 * NOTE: this uses simply the inversion-routine used in the original
 *       makefakedata_v2
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="GeneratePulsarSignalCP"> */
void
LALConvertSSB2GPS (LALStatus *stat, 
		   LIGOTimeGPS *GPSout,		 /* output: GPS-arrival-time at detector */
		   LIGOTimeGPS SSBin, 		 /* input: signal arrival time at SSB */
		   const PulsarSignalParams *params) /* params defining source-location and detector */
{ /* </lalVerbatim> */

  LIGOTimeGPS SSBofguess;
  LIGOTimeGPS GPSguess;
  INT4 iterations, E9=1000000000;
  INT8 delta, guess;

  INITSTATUS( stat, "ConvertSSB2GPS", GENERATEPULSARSIGNALC );
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
      TRY ( LALConvertGPS2SSB (stat->statusPtr, &SSBofguess, GPSguess, params), stat);

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
    ABORT ( stat, GENERATEPULSARSIGNALH_ESSBCONVERT, GENERATEPULSARSIGNALH_MSGESSBCONVERT);
  }

  /* Now that we've found the GPS time that corresponds to the given SSB time */
  *GPSout = GPSguess;
  
  
  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* LALConvertSSB2GPS() */




/*----------------------------------------------------------------------
 * create one SFT-struct
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="GeneratePulsarSignalCP"> */
void
LALCreateSFTtype (LALStatus *stat, 
		  SFTtype **output, 	/* output: allocated SFT-struct */
		  UINT4 SFTlen)		/* number of frequency-bins */
{ /* </lalVerbatim> */
  SFTtype *sft = NULL;

  INITSTATUS( stat, "LALCreateSFTtype", GENERATEPULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );

  ASSERT (output != NULL, stat, GENERATEPULSARSIGNALH_ENULL,  GENERATEPULSARSIGNALH_MSGENULL);
  ASSERT (*output == NULL, stat, GENERATEPULSARSIGNALH_ENONULL,  GENERATEPULSARSIGNALH_MSGENONULL);

  sft = LALCalloc (1, sizeof(*sft) );
  if (sft == NULL) {
    ABORT (stat, GENERATEPULSARSIGNALH_EMEM, GENERATEPULSARSIGNALH_MSGEMEM);
  }
  LALCCreateVector (stat->statusPtr, &(sft->data), SFTlen);
  BEGINFAIL (stat) { 
    LALFree (sft);
  } ENDFAIL (stat);

  *output = sft;

  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* LALCreateSFTVector() */


/*----------------------------------------------------------------------
 * create a whole vector of <numSFT> SFTs with <SFTlen> frequency-bins 
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="GeneratePulsarSignalCP"> */
void
LALCreateSFTVector (LALStatus *stat, 
		    SFTVector **output, /* output: allocated SFT-vector */
		    UINT4 numSFTs, 	/* number of SFTs */
		    UINT4 SFTlen)	/* number of frequency-bins per SFT */
{ /* </lalVerbatim> */
  UINT4 iSFT, j;
  SFTVector *vect;	/* vector to be returned */
  COMPLEX8Vector *data = NULL;

  INITSTATUS( stat, "LALCreateSFTVector", GENERATEPULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );

  ASSERT (output != NULL, stat, GENERATEPULSARSIGNALH_ENULL,  GENERATEPULSARSIGNALH_MSGENULL);
  ASSERT (*output == NULL, stat, GENERATEPULSARSIGNALH_ENONULL,  GENERATEPULSARSIGNALH_MSGENONULL);

  vect = LALCalloc (1, sizeof(*vect) );
  if (vect == NULL) {
    ABORT (stat, GENERATEPULSARSIGNALH_EMEM, GENERATEPULSARSIGNALH_MSGEMEM);
  }

  vect->length = numSFTs;

  vect->data = LALCalloc (1, numSFTs * sizeof ( *vect->data ) );
  if (vect->data == NULL) {
    LALFree (vect);
    ABORT (stat, GENERATEPULSARSIGNALH_EMEM, GENERATEPULSARSIGNALH_MSGEMEM);
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

/*----------------------------------------------------------------------
 * destroy an SFT-struct
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="GeneratePulsarSignalCP"> */
void
LALDestroySFTtype (LALStatus *stat, 
		   SFTtype **sft)
{ /* </lalVerbatim> */

  INITSTATUS( stat, "LALDestroySFTtype", GENERATEPULSARSIGNALC);
  ATTATCHSTATUSPTR (stat);

  ASSERT (sft != NULL, stat, GENERATEPULSARSIGNALH_ENULL,  GENERATEPULSARSIGNALH_MSGENULL);
  ASSERT (*sft != NULL, stat, GENERATEPULSARSIGNALH_ENULL,  GENERATEPULSARSIGNALH_MSGENULL);
  

  LALCDestroyVector (stat->statusPtr, &((*sft)->data) );
  LALFree ( (*sft) );

  *sft = NULL;

  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* LALDestroySFTtype() */


/*----------------------------------------------------------------------
 * destroy an SFT-vector
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="GeneratePulsarSignalCP"> */
void
LALDestroySFTVector (LALStatus *stat, 
		     SFTVector **vect)	/* the SFT-vector to free */
{ /* </lalVerbatim> */
  UINT4 i;

  INITSTATUS( stat, "LALDestroySFTVector", GENERATEPULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );

  ASSERT (vect != NULL, stat, GENERATEPULSARSIGNALH_ENULL,  GENERATEPULSARSIGNALH_MSGENULL);
  ASSERT (*vect != NULL, stat, GENERATEPULSARSIGNALH_ENULL,  GENERATEPULSARSIGNALH_MSGENULL);

  
  for (i=0; i < (*vect)->length; i++)
    LALCDestroyVector (stat->statusPtr, &((*vect)->data[i].data) );

  LALFree ( (*vect)->data );
  LALFree ( *vect );

  *vect = NULL;

  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* LALDestroySFTVector() */



/*----------------------------------------------------------------------
 * allocate a LIGOTimeGPSVector
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="GeneratePulsarSignalCP"> */
void
LALCreateTimestampVector (LALStatus *stat, LIGOTimeGPSVector **vect, UINT4 length)
{ /* </lalVerbatim> */
  LIGOTimeGPSVector *out = NULL;

  INITSTATUS( stat, "LALCreateTimestampVector", GENERATEPULSARSIGNALC);

  ASSERT (vect != NULL, stat, GENERATEPULSARSIGNALH_ENULL,  GENERATEPULSARSIGNALH_MSGENULL);
  ASSERT (*vect == NULL, stat, GENERATEPULSARSIGNALH_ENONULL,  GENERATEPULSARSIGNALH_MSGENONULL);

  out = LALCalloc (1, sizeof(LIGOTimeGPSVector));
  if (out == NULL) {
    ABORT (stat,  GENERATEPULSARSIGNALH_EMEM,  GENERATEPULSARSIGNALH_MSGEMEM);
  }
  out->length = length;
  out->data = LALCalloc (1, length * sizeof(LIGOTimeGPS));
  if (out->data == NULL) {
    LALFree (out);
    ABORT (stat,  GENERATEPULSARSIGNALH_EMEM,  GENERATEPULSARSIGNALH_MSGEMEM);
  }

  *vect = out;

  RETURN (stat);
  
} /* LALCreateTimestampVector() */



/*----------------------------------------------------------------------
 * de-allocate a LIGOTimeGPSVector
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="GeneratePulsarSignalCP"> */
void
LALDestroyTimestampVector (LALStatus *stat, LIGOTimeGPSVector **vect)
{ /* </lalVerbatim> */
  INITSTATUS( stat, "LALDestroyTimestampVector", GENERATEPULSARSIGNALC);

  ASSERT (vect != NULL, stat, GENERATEPULSARSIGNALH_ENULL,  GENERATEPULSARSIGNALH_MSGENULL);
  ASSERT (*vect != NULL, stat, GENERATEPULSARSIGNALH_ENULL,  GENERATEPULSARSIGNALH_MSGENULL);

  LALFree ( (*vect)->data);
  LALFree ( *vect );
  
  *vect = NULL;

  RETURN (stat);
  
} /* LALDestroyTimestampVector() */



/************************************************************************
 * the following are INTERNAL FUNCTIONS not to be called outside of this 
 * module 
 ************************************************************************/

/*----------------------------------------------------------------------
 * given a start-time, duration and Tsft, returns a list of timestamps
 * covering this time-stretch.
 * returns NULL on out-of-memory
 *----------------------------------------------------------------------*/
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

/*----------------------------------------------------------------------
 *   check that all timestamps given lie within the range [t0, t1] 
 *    return: 0 if ok, -1 if not
 *----------------------------------------------------------------------*/
int
check_timestamp_bounds (const LIGOTimeGPSVector *timestamps, LIGOTimeGPS t0, LIGOTimeGPS t1)
{
  LALStatus status1, status2;
  REAL8 diff0, diff1;
  UINT4 i;

  status1 = status2 = emptyStatus;

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

} /* check_timestamp_bounds() */

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

  INITSTATUS( stat, "checkNoiseSFTs", GENERATEPULSARSIGNALC);

  for (i=0; i < sfts->length; i++)
    {
      thisSFT = &(sfts->data[i]);
      deltaFn = thisSFT->deltaF;
      fn0 = thisSFT->f0;
      fn1 = f0 + thisSFT->data->length * deltaFn;
      
      if (deltaFn != deltaF) {
	ABORT (stat,  GENERATEPULSARSIGNALH_ENOISEDELTAF,  GENERATEPULSARSIGNALH_MSGENOISEDELTAF);
      }

      if ( (f0 < fn0) || (f1 > fn1) ) {
	ABORT (stat, GENERATEPULSARSIGNALH_ENOISEBAND, GENERATEPULSARSIGNALH_MSGENOISEBAND);
      }
      
      shift = (fn0 - f0) / deltaF;
      /* frequency bins have to coincide! ==> check that shift is integer!  */
      nshift = (UINT4)(shift+0.5);
      if ( fabs( nshift - shift) > eps ) {
	ABORT (stat, GENERATEPULSARSIGNALH_ENOISEBINS, GENERATEPULSARSIGNALH_MSGENOISEBINS);
      }

    } /* for i < numSFTs */

  RETURN (stat);

} /* checkNoiseSFTs() */



/*----------------------------------------------------------------------
 * Yousuke's phase-correction function, taken from makefakedata_v2
 *----------------------------------------------------------------------*/
void
correct_phase (LALStatus* stat, SFTtype *sft, LIGOTimeGPS tHeterodyne) 
{
  UINT4 i;
  REAL8 cosx,sinx;
  COMPLEX8 fvec1;
  REAL8 deltaT;

  INITSTATUS( stat, "correct_phase", GENERATEPULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );

  TRY (LALDeltaFloatGPS(stat->statusPtr, &deltaT, &(sft->epoch), &tHeterodyne), stat);
  deltaT *= sft->f0; 

  /* check if we really need to do anything here? (i.e. is deltaT an integer?) */
  if ( fabs (deltaT - (INT4) deltaT ) > eps )
    {
      deltaT *= LAL_TWOPI;

      cosx = cos (deltaT);
      sinx = sin(deltaT);

      for (i = 0; i < sft->data->length; i++)
	{
	  fvec1 = sft->data->data[i];
	  sft->data->data[i].re = fvec1.re * cosx - fvec1.im * sinx;
	  sft->data->data[i].im = fvec1.im * cosx + fvec1.re * sinx;
	} /* for i < length */

    } /* if deltaT/2pi not integer */

  DETATCHSTATUSPTR( stat );
  RETURN (stat);
  
} /* correct_phase() */
