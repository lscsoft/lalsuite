/************************************ <lalVerbatim file="GeneratePulsarSignalCV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/* NOTES: */
/* 07/14/04 gam; add functions LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
/* 10/08/04 gam; fix indexing into trig lookup tables (LUTs) by having table go from -2*pi to 2*pi */
/* 10/12/04 gam; When computing fCross and fPlus need to use 2.0*psi. */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{GeneratePulsarSignal.c}}
\label{ss:GeneratePulsarSignal.c}

Module to generate simulated pulsar-signals.

\subsubsection*{Prototypes}
\idx{LALGeneratePulsarSignal()}
\idx{LALSignalToSFTs()}
\idx{LALComputeSkyAndZeroPsiAMResponse()}
\idx{LALFastGeneratePulsarSFTs()}
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
\end{itemize}

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

\subsubsection*{Use of LALFastGeneratePulsarSFTs()}

The functions \verb+LALComputeSkyAndZeroPsiAMResponse()+ and \verb+LALFastGeneratePulsarSFTs()+ 
use approximate analytic formulas to generate SFTs.  This should be significantly
faster than \verb+LALGeneratePulsarSignal()+ and \verb+LALSignalToSFTs()+, which generate the
time series data and then FFT it.  Simple tests performed by the code in
\verb+GeneratePulsarSignalTest.c+ indicate that the maximum modulus of the SFTs output
by the approximate and exact codes differs by less that 10\%.  Since the tests
are not exhaustive, the user should use caution and conduct their own test
to compare \verb+LALComputeSkyAndZeroPsiAMResponse()+ and \verb+LALFastGeneratePulsarSFTs()+ with
\verb+LALGeneratePulsarSignal()+ and \verb+LALSignalToSFTs()+.

The strain of a periodic signal at the detector is given by
$$
h(t) = F_+(t) A_+ {\rm cos}\Phi(t) + F_\times(t) A_\times {\rm sin}\Phi(t),
$$
where $F_+$ and $F_\times$ are the usual beam pattern response functions, 
$A_+$ and $A_\times$ are the amplitudes of the gravitational wave for the
plus and cross polarizations, and $\Phi$ is the phase.  The phase contains modulations
from doppler shifts due to the relative motion between the source and the
detector and the spin evolution of the source.  (The functions discussed here
support both isolated sources and those in binary systems. The binary case
has not been tested.)

If we Taylor expand the phase out to first order about the time at the midpoint of
an SFT and approximate $F_+$ and $F_\times$ as constants, for one SFT we can write
$$
\Phi(t) \approx \Phi_{1/2} + 2\pi f_{1/2}(t - t_{1/2}).
$$
The strain at discrete time $t_j$, measured from the start of the SFT, can 
thus be approximated as
$$
h_j \approx F_{+ 1/2} A_+ {\rm cos} [\Phi_{1/2} + 2\pi f_{1/2}(t_0 + t_j - t_{1/2})]
+ F_{\times 1/2} A_\times {\rm sin} [\Phi_{1/2} + 2\pi f_{1/2}(t_0 + t_j - t_{1/2})],
$$
where $t_0$ is the time as the start of the SFT, and $t_{1/2} - t_0 = T_{\rm sft}/2$,
where $T_{\rm sft}$ is the duration of one SFT.  This simplifies to
$$
h_j \approx F_{+ 1/2} A_+ {\rm cos} (\Phi_0 + 2\pi f_{1/2}t_j)
+ F_{\times 1/2} A_\times {\rm sin} (\Phi_0 + 2\pi f_{1/2}t_j),
$$
where $\Phi_0$ is the phase at the start of the SFT
(not the initial phase at the start of the observation), i.e.,
$$
\Phi_0 = \Phi_{1/2} - 2 \pi f_{1/2} (T_{\rm sft} / 2).
$$
Note that these are the same approximations used by LALDemod.

One can show that the Discrete Fourier Transform (DFT) of $h_j$ above is:
$$
\tilde{h}_k = e^{i\Phi_0}  { ( F_{+ 1/2} A_+ - i F_{\times 1/2} A_\times) \over 2 } 
{ 1 - e^{2\pi i (\kappa - k)} \over 1 - e^{2\pi i (\kappa - k)/N} } 
\\
+ e^{-i\Phi_0}  { ( F_{+ 1/2} A_+ + i F_{\times 1/2} A_\times) \over 2 }
{ 1 - e^{-2\pi i (\kappa + k)} \over 1 - e^{-2\pi i (\kappa + k)/N} }
$$
where $N$ is the number of time samples used to find the
DFT (i.e., the sample rate times $T_{\rm sft}$), and
$$
\kappa \equiv f_{1/2} T_{\rm sft},
$$
is usually not an integer.

Note that the factor $e^{\pm 2\pi i k}$ in the numerators of the equation for $\tilde{h}_k$
equals 1.  Furthermore, for $0 < \kappa < N/2$ and $|\kappa - k| << N$ the first term
dominates and can be Taylor expanded to give:
$$
\tilde{h}_k = N e^{i\Phi_0} { ( F_{+ 1/2} A_+ - i F_{\times 1/2} A_\times) \over 2 }
\left [ \, { {\rm sin} (2\pi\kappa) \over 2 \pi (\kappa - k) } \,
+ \, i { 1 - {\rm cos} (2\pi\kappa) \over 2 \pi (\kappa - k) } \, \right ]
$$
Note that the last factor in square brackets is $P_{\alpha k}^*$ and
$e^{i\Phi_0} = Q_{\alpha}^*$ used by LALDemod.

\subsubsection*{Example pseudocode}

The structs used by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs 
are given in previous sections, and make use of those used by
LALGeneratePulsarSignal and LALSignalToSFTs plus a small number of
additional parameters.  Thus it is fairly easy to change between the above
approximate routines the exact routines. See GeneratePulsarSignalTest.c for
an example implementation of the code.

Note that one needs to call LALComputeSkyAndZeroPsiAMResponse once per sky position,
and then call LALFastGeneratePulsarSFTs for each set of signal parameters at that
sky position.  Thus, one could perform a Monte Carlo simulation, as shown
by the pseudo code:

\begin{verbatim}

loop over sky positions {
   ...
   LALComputeSkyAndZeroPsiAMResponse();
   ...
   loop over spindown {
      ...
      loop over frequencies {
         ...
         LALFastGeneratePulsarSFTs();
         ...
      }
      ...  
   }
   ...
}

\end{verbatim}

\subsubsection*{Notes on LALFastGeneratePulsarSFTs}

1) If \verb+*outputSFTs+ sent to \verb+LALFastGeneratePulsarSFTs()+ is \verb+NULL+ then
\verb+LALFastGeneratePulsarSFTs()+ allocates memory for the output SFTs; otherwise it assumes
memory has already been allocated.  Thus, the user does not have to deallocate
memory for the SFTs until all calls to \verb+LALFastGeneratePulsarSFTs()+ are completed.

\noindent 2) \verb+fHeterodyne+ and 0.5*\verb+samplingRate+ set in the \verb+PulsarSignalParams+ struct
give the start frequency and frequency band of the SFTs output from \verb+LALFastGeneratePulsarSFTs+.

\noindent 3) If \verb+resTrig+ is set to zero in the \verb+SFTandSignalParams+ struct, then
the C math libary \verb+cos+ and \verb+sin+ functions are called, else lookup tables (LUTs) are used
for calls to trig functions.  There may be a slight speedup in using LUTs.

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
  /* set source position: make sure it's "normalized", i.e. [0<=alpha<2pi]x[-pi/2<=delta<=pi/2] */
  TRY( LALNormalizeSkyPosition(stat->statusPtr, &(sourceParams.position), &(params->pulsar.position)), stat);

  /* if pulsar is in binary-orbit, set binary parameters */
  if (params->orbit)
    {
      /*------------------------------------------------------------ */
      /* temporary fix for comparison with Chris' code */
      /*
	TRY (LALConvertGPS2SSB (stat->statusPtr, &tmpTime, params->orbit->orbitEpoch, params), stat);
	sourceParams.orbitEpoch = tmpTime;
      */
      sourceParams.orbitEpoch =  params->orbit->orbitEpoch;
      /* ------------------------------------------------------------*/
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
      TRY ( LALConvertGPS2SSB(stat->statusPtr, &tmpTime, params->startTimeGPS, params), stat);
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
  TRY (LALConvertGPS2SSB(stat->statusPtr, &t0, params->startTimeGPS, params), stat);
  t0.gpsSeconds -= (UINT4)sourceParams.deltaT; /* start one time-step earlier to be safe */

  /* end time in SSB */
  t1 = params->startTimeGPS;
  TRY ( LALAddFloatToGPS(stat->statusPtr, &t1, &t1, params->duration), stat);
  TRY (LALConvertGPS2SSB(stat->statusPtr, &t1, t1, params), stat);	 /* convert time to SSB */

  /* get duration of source-signal */
  TRY (LALDeltaFloatGPS(stat->statusPtr, &SSBduration, &t1, &t0), stat);
  SSBduration += 2.0 * sourceParams.deltaT; /* add two time-steps to be safe */

  sourceParams.epoch = t0; 
  sourceParams.length = (UINT4) ceil( SSBduration / sourceParams.deltaT );

  /* we use frequency-spindowns, but GenerateSpinOrbitCW wants it f0-normalized,
     so we have to do that here: */
  sourceParams.f = NULL;
  if (params->pulsar.spindown)
    {
      TRY ( LALDCreateVector(stat->statusPtr, &(sourceParams.f), params->pulsar.spindown->length), stat);
      for (i=0; i < sourceParams.f->length; i++)
	sourceParams.f->data[i] = params->pulsar.spindown->data[i] / params->pulsar.f0;
    }

  /* finally, call the function to generate the source waveform */
  TRY ( LALGenerateSpinOrbitCW(stat->statusPtr, &sourceSignal, &sourceParams), stat);

  /* free spindown-vector right away, so we don't forget */
  if (sourceParams.f) {
    TRY ( LALDDestroyVector(stat->statusPtr, &(sourceParams.f) ), stat);
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

  /* NOTE: a timeseries of length N*dT has no timestep at N*dT !! (convention) */
  LALCreateVector(stat->statusPtr, &(output->data), (UINT4) ceil( params->samplingRate * params->duration) );
  BEGINFAIL(stat) {
    LALFree (output);
  } ENDFAIL(stat);

  output->deltaT = 1.0 / params->samplingRate;
  output->f0 = params->fHeterodyne;
  output->epoch = params->startTimeGPS;

  
  TRY ( LALSimulateCoherentGW(stat->statusPtr, output, &sourceSignal, &detector ), stat );

			       
  /*----------------------------------------------------------------------*/
  /* Free all allocated memory that is not returned */
  TRY (LALSDestroyVectorSequence( stat->statusPtr, &(sourceSignal.a->data)), stat);
  LALFree( sourceSignal.a );
  TRY (LALSDestroyVector(stat->statusPtr, &(sourceSignal.f->data ) ), stat);
  LALFree(sourceSignal.f);
  TRY (LALDDestroyVector(stat->statusPtr, &(sourceSignal.phi->data )), stat);
  LALFree(sourceSignal.phi);

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
  UINT4 numSFTsamples;			/* number of time-samples in an Tsft */
  UINT4 iSFT;
  REAL8 Band, SFTsamples, f0, deltaF;
  LIGOTimeGPSVector *timestamps = NULL;
  REAL4Vector timeStretch = {0,0};
  LIGOTimeGPS tStart;			/* start time of input time-series */
  LIGOTimeGPS tLast;			/* start-time of last _possible_ SFT */
  LIGOTimeGPS tmpTime;
  REAL8 duration, delay;
  UINT4 SFTlen;				/* number of frequency-bins in an SFT */
  UINT4 index0n;			/* first frequency-bin to use from noise-SFT */
  SFTtype *thisSFT, *thisNoiseSFT;	/* SFT-pointers */
  REAL4 renorm;				/* renormalization-factor of taking only part of an SFT */
  SFTVector *sftvect = NULL;		/* return value. For better readability */
  UINT4 j;
  RealFFTPlan *pfwd = NULL;		/* FFTW plan */
  LIGOTimeGPS tPrev;			/* real timstamp of previous SFT */
  UINT4 totalIndex;			/* timestep-index to start next FFT from */
  INT4 relIndexShift;			/* relative index-shift from previous SFT (number of timesteps) */

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
    TRY (checkNoiseSFTs(stat->statusPtr, params->noiseSFTs, f0, f0 + Band, deltaF), stat);
  }
    
  /* make sure that number of timesamples/SFT is an integer (up to possible rounding errors */
  SFTsamples = params->Tsft / signal->deltaT;		/* this is a float!*/
  numSFTsamples = (UINT4) (SFTsamples + 0.5);		/* round to closest int */
  ASSERT ( fabs(SFTsamples - numSFTsamples)/SFTsamples < eps, stat, 
	   GENERATEPULSARSIGNALH_EINCONSBAND, GENERATEPULSARSIGNALH_MSGEINCONSBAND);
  
  /* Prepare FFT: compute plan for FFTW */
  TRY (LALCreateForwardRealFFTPlan(stat->statusPtr, &pfwd, numSFTsamples, 0), stat); 	

  /* get some info about time-series */
  tStart = signal->epoch;					/* start-time of time-series */

  /* get last possible start-time for an SFT */
  duration =  (UINT4) (1.0* signal->data->length * signal->deltaT +0.5); /* total duration rounded to seconds */
  TRY ( LALAddFloatToGPS(stat->statusPtr, &tLast, &tStart, duration - params->Tsft), stat);

  /* for simplicity we _always_ work with timestamps. 
   * Therefore, we have to generate them now if none have been provided by the user. */
  if (params->timestamps == NULL)
    {
      timestamps = make_timestamps(&tStart, duration, params->Tsft );
      if (timestamps == NULL) {
	ABORT (stat, GENERATEPULSARSIGNALH_EMEM, GENERATEPULSARSIGNALH_MSGEMEM);
      }
    }
  else	/* if given, use those, and check they are valid */
    {
      timestamps = params->timestamps;
      /* check that all timestamps lie within [tStart, tLast] */
      if ( check_timestamp_bounds(timestamps, tStart, tLast) != 0) {
	ABORT (stat, GENERATEPULSARSIGNALH_ETIMEBOUND, GENERATEPULSARSIGNALH_MSGETIMEBOUND);
      }
    }

  /* prepare SFT-vector for return */
  numSFTs = timestamps->length;			/* number of SFTs to produce */
  SFTlen = (UINT4)(numSFTsamples/2) + 1;	/* number of frequency-bins per SFT */

  LALCreateSFTVector (stat->statusPtr, &sftvect, numSFTs, SFTlen);
  BEGINFAIL (stat) {
    if (params->timestamps == NULL)
      LALDestroyTimestampVector(stat->statusPtr, &timestamps);
  } ENDFAIL (stat);

  tPrev = tStart;	/* initialize */
  totalIndex = 0;	/* start from first timestep by default */

  /* main loop: apply FFT the requested time-stretches */
  for (iSFT = 0; iSFT < numSFTs; iSFT++)
    {
      thisSFT = &(sftvect->data[iSFT]);	/* point to current SFT-slot */

      /* find the start-bin for this SFT in the time-series */
      TRY ( LALDeltaFloatGPS(stat->statusPtr, &delay, &(timestamps->data[iSFT]), &tPrev), stat);
      relIndexShift = (INT4) (delay / signal->deltaT + 0.5);	/* round properly: picks *closest* timestep (==> "nudging") !!  */
      totalIndex += relIndexShift;

      timeStretch.length = numSFTsamples;    
      timeStretch.data = signal->data->data + totalIndex;	/* point to the right sample-bin */

      /* fill the header of the i'th output SFT */
      TRY ( LALAddFloatToGPS(stat->statusPtr, &tmpTime, &tPrev, relIndexShift * signal->deltaT), stat);
      thisSFT->epoch = tmpTime;			/* set the ACTUAL timestamp! (can be different from requested one ==> "nudging") */
      thisSFT->f0 = signal->f0;			/* minimum frequency */
      thisSFT->deltaF = 1.0 / params->Tsft;	/* frequency-spacing */

      tPrev = tmpTime;				/* prepare next loop */

      /* ok, issue at least a warning if we have "nudged" an SFT-timestamp */
      if (lalDebugLevel > 0)
	{
	  REAL8 diff;
	  TRY ( LALDeltaFloatGPS(stat->statusPtr, &diff, &(timestamps->data[iSFT]), &tmpTime), stat);
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
      LALForwardRealFFT(stat->statusPtr, thisSFT->data, &timeStretch, pfwd);
      BEGINFAIL(stat) {
	LALDestroySFTVector(stat->statusPtr, &sftvect);
      } ENDFAIL(stat);


      /* correct heterodyning-phase, IF NECESSARY */
      if ( ( (INT4)signal->f0 != signal->f0  )
	   || (signal->epoch.gpsNanoSeconds != 0)
	   || (thisSFT->epoch.gpsNanoSeconds != 0) )
	{
	  correct_phase(stat->statusPtr, thisSFT, signal->epoch);	/* theterodyne = signal->epoch!*/
	  BEGINFAIL (stat) {
	    LALDestroySFTVector(stat->statusPtr, &sftvect);
	  } ENDFAIL (stat);
	} /* if phase-correction necessary */

      /* Now add the noise-SFTs if given */
      if (params->noiseSFTs)
	{
	  COMPLEX8 *data, *noise;
	  thisNoiseSFT = &(params->noiseSFTs->data[iSFT]);
	  index0n = (UINT4) ((thisSFT->f0 - thisNoiseSFT->f0) / thisSFT->deltaF);

	  /* The renormalization follows strictly makefakedata_v2. */
	  renorm = 1.0*SFTlen/(thisNoiseSFT->data->length);

	  data = &(thisSFT->data->data[0]);
	  noise = &(thisNoiseSFT->data->data[index0n]);
	  for (j=0; j < SFTlen; j++)
	    {
	      data->re += renorm * noise->re;
	      data->im += renorm * noise->im;
	      data++;
	      noise++;
	    } /* for j < SFTlen */
	
	} /* if noiseSFTs */

    } /* for iSFT < numSFTs */ 

  /* free stuff */
  LALDestroyRealFFTPlan(stat->statusPtr, &pfwd);

  /* did we get timestamps or did we make them? */
  if (params->timestamps == NULL)
    {
      LALFree(timestamps->data);
      LALFree(timestamps);
      timestamps = NULL;
    }

  *outputSFTs = sftvect;

  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* LALSignalToSFTs() */



/* 07/14/04 gam */
/*--------------------------------------------------------------------------
 * Wrapper for ComputeSky and  LALComputeDetAMResponse that finds the sky
 * constants and F_+ and F_x for use with LALFastGeneratePulsarSFTs. Uses
 * the output of LALComputeSkyAndZeroPsiAMResponse and the same inputs as
 * LALGeneratePulsarSignal and LALSignalToSFTs.
 * This function used ComputeSkyBinary if params->pSigParams->orbit is not
 * NULL, else it uses ComputeSky to find the skyConsts.
 * NOTE THAT THIS FUNCTION COMPUTES F_+ and F_x for ZERO Psi!!!
 * LALFastGeneratePulsarSFTs used these to find F_+ and F_x for NONZERO Psi.
 *--------------------------------------------------------------------------*/
/* <lalVerbatim file="GeneratePulsarSignalCP"> */
void
LALComputeSkyAndZeroPsiAMResponse (LALStatus *stat,
                                   SkyConstAndZeroPsiAMResponse *output,
                                   const SFTandSignalParams *params)
{ /* </lalVerbatim> */
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
  LALGPSandAcc   timeAndAcc;  /* input for LALComputeDetAMResponse */
  REAL8 halfTsft;             /* half the time of one SFT */
  LIGOTimeGPS midTS;          /* midpoint time for an SFT */

  INITSTATUS( stat, "LALComputeSkyAndZeroPsiAMResponse", GENERATEPULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );
  
  numSFTs = params->pSFTParams->timestamps->length; /* number of SFTs */
  halfTsft = 0.5*params->pSFTParams->Tsft;          /* half the time of one SFT */

  /* setup baryinput for ComputeSky */
  baryinput.site = *(params->pSigParams->site);
  /* account for a quirk in LALBarycenter(): -> see documentation of type BarycenterInput */
  baryinput.site.location[0] /= LAL_C_SI;
  baryinput.site.location[1] /= LAL_C_SI;
  baryinput.site.location[2] /= LAL_C_SI;
  if (params->pSigParams->pulsar.position.system != COORDINATESYSTEM_EQUATORIAL) {
      ABORT (stat, GENERATEPULSARSIGNALH_EBADCOORDS, GENERATEPULSARSIGNALH_MSGEBADCOORDS);
  }
  TRY( LALNormalizeSkyPosition (stat->statusPtr, &tmp, &(params->pSigParams->pulsar.position)), stat);
  baryinput.alpha = tmp.longitude;
  baryinput.delta = tmp.latitude;
  baryinput.dInv = 0.e0;      /* following makefakedata_v2 */
  
  if (params->pSigParams->orbit) {
    /* ComputeSkyBinary parameters */
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
    csbParams->OrbitalEccentricity = 1.0 - params->pSigParams->orbit->oneMinusEcc; /* Orbital eccentricy */
    csbParams->ArgPeriapse = params->pSigParams->orbit->omega;       /* argument of periapsis (radians) */
    csbParams->TperiapseSSB = params->pSigParams->orbit->orbitEpoch; /* time of periapsis passage (in SSB) */
    /* compute semi-major axis and orbital period */
    csbParams->SemiMajorAxis = params->pSigParams->orbit->rPeriNorm/params->pSigParams->orbit->oneMinusEcc; 
    csbParams->OrbitalPeriod = ( ((REAL8)LAL_TWOPI) / params->pSigParams->orbit->angularSpeed ) * sqrt( (1.0 + csbParams->OrbitalEccentricity) / pow(params->pSigParams->orbit->oneMinusEcc,3.0) );
    csbParams->baryinput=&baryinput;
    csbParams->emit = &emit;
    csbParams->earth = &earth;
    csbParams->edat=params->pSigParams->ephemerides;

    /* Call ComputeSkyBinary */
    TRY ( ComputeSkyBinary (stat->statusPtr, output->skyConst, 0, csbParams), stat);
    LALFree(csbParams->skyPos);
    LALFree(csbParams);
  } else {
    /* ComputeSky parameters */
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

    /* Call ComputeSky */
    TRY ( ComputeSky (stat->statusPtr, output->skyConst, 0, csParams), stat);
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
  timeAndAcc.accuracy=LALLEAPSEC_STRICT;
 
  /* loop that calls LALComputeDetAMResponse to find F_+ and F_x at the midpoint of each SFT for ZERO Psi */
  for(i=0; i<numSFTs; i++) {
      /* Find mid point from timestamp, half way through SFT. */
      TRY ( LALAddFloatToGPS (stat->statusPtr, &midTS, &(params->pSFTParams->timestamps->data[i]), halfTsft), stat);
      timeAndAcc.gps=midTS;
      TRY ( LALComputeDetAMResponse(stat->statusPtr, &response, das, &timeAndAcc), stat);
      output->fPlusZeroPsi[i] = response.plus;
      output->fCrossZeroPsi[i] = response.cross;
  }
  LALFree(das->pSource);
  LALFree(das);

  DETATCHSTATUSPTR( stat );
  RETURN (stat);
} /* LALComputeSkyAndZeroPsiAMResponse */

/* 07/14/04 gam */
/*--------------------------------------------------------------------------
 * Fast generation of Fake SFTs for a pure pulsar signal.
 * Uses the output of LALComputeSkyAndZeroPsiAMResponse and the same inputs
 * as LALGeneratePulsarSignal and LALSignalToSFTs.  The fake signal is
 * Taylor expanded to first order about the midpoint time of each SFT.
 * Analytic expressions are used to find each SFT
 *--------------------------------------------------------------------------*/
/* <lalVerbatim file="GeneratePulsarSignalCP"> */
void 
LALFastGeneratePulsarSFTs (LALStatus *stat,
                           SFTVector **outputSFTs,
                           const SkyConstAndZeroPsiAMResponse *input,
                           const SFTandSignalParams *params)
{ /* </lalVerbatim> */
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
  REAL8 smallX=0.000000001;
  /* Next are for LUT for trig calls */
  INT4 indexTrig;
  REAL8 halfResTrig = ((REAL8)params->resTrig)/2.0; /* 10/08/04 gam; fix indexing into trig lookup tables (LUTs) by having table go from -2*pi to 2*pi */
  REAL8 varTmp, dTmp, dTmp2, sinTmp, cosTmp;

  INITSTATUS( stat, "LALFastGeneratePulsarSFTs", GENERATEPULSARSIGNALC);
  ATTATCHSTATUSPTR( stat );
 
  /* fprintf(stdout,"\n Hello from LALFastGeneratePulsarSFTs \n");
  fflush(stdout); */
       
  ASSERT (outputSFTs != NULL, stat, GENERATEPULSARSIGNALH_ENULL, GENERATEPULSARSIGNALH_MSGENULL);
  /* ASSERT (*outputSFTs == NULL, stat,  GENERATEPULSARSIGNALH_ENONULL,  GENERATEPULSARSIGNALH_MSGENONULL); */
  ASSERT (params != NULL, stat, GENERATEPULSARSIGNALH_ENULL, GENERATEPULSARSIGNALH_MSGENULL);  
  ASSERT (input != NULL, stat, GENERATEPULSARSIGNALH_ENULL, GENERATEPULSARSIGNALH_MSGENULL);

  if ( params->pSFTParams->timestamps && params->pSFTParams->noiseSFTs) {
    ASSERT ( params->pSFTParams->timestamps->length == params->pSFTParams->noiseSFTs->length, stat,  
      GENERATEPULSARSIGNALH_ENUMSFTS,  GENERATEPULSARSIGNALH_MSGENUMSFTS);
  }
  
  /* 10/08/04 gam; fix indexing into trig lookup tables (LUTs) by having table go from -2*pi to 2*pi */
  if (params->resTrig > 0) {
     ASSERT ( fabs( ( params->trigArg[0] + ((REAL8)LAL_TWOPI) )/ ((REAL8)LAL_TWOPI) ) < ( 2.0e-6*((REAL8)LAL_TWOPI) / params->resTrig ), stat,
           GENERATEPULSARSIGNALH_ELUTS,  GENERATEPULSARSIGNALH_MSGELUTS);
     ASSERT ( fabs( ( params->trigArg[params->resTrig] - ((REAL8)LAL_TWOPI) ) / ((REAL8)LAL_TWOPI)  ) < ( 2.0e-6*((REAL8)LAL_TWOPI) / params->resTrig ), stat,
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
  
  /* prepare SFT-vector for return */
  if (*outputSFTs == NULL) {
    TRY (LALCreateSFTVector (stat->statusPtr, &sftvect, numSFTs, SFTlen), stat);
  } else {
    sftvect = *outputSFTs;  /* Assume memory already allocated for SFTs */
  }
    
  /* if noiseSFTs are given: check they are consistent with signal! */
  if (params->pSFTParams->noiseSFTs) {
    TRY (checkNoiseSFTs (stat->statusPtr, params->pSFTParams->noiseSFTs, f0, f0 + band, deltaF), stat);
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
      
      /* fill in the data */
      for (j=0; j<SFTlen; j++) {
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
          thisSFT->data->data[j].re = (REAL4)(realTmp*realA - imagTmp*imagA);
          thisSFT->data->data[j].im = (REAL4)(realTmp*imagA + imagTmp*realA);
      } /* for j < SFTlen */
            
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
           thisSFT->data->data[j].re += thisNoiseSFT->data->data[index0n + j].re;
           thisSFT->data->data[j].im += thisNoiseSFT->data->data[index0n + j].im;
        } /* for j < SFTlen */
      }
  } /* for iSFT < numSFTs */ 

  /* prepare SFT-vector for return */
  if (*outputSFTs == NULL) {
    *outputSFTs = sftvect;
  } /* else sftvect already points to same memory as *outputSFTs */
          
  DETATCHSTATUSPTR( stat );
  RETURN (stat);
} /* LALFastGeneratePulsarSFTs () */



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

  /* be flexible: if points to null, nothing to do here.. (like 'free()') */
  if (*sft == NULL) {
    DETATCHSTATUSPTR( stat );
    RETURN (stat);
  }
    
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
  LALStatus status = emptyStatus;
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

#define oneBillion 1000000000;
/*----------------------------------------------------------------------
 *   check that all timestamps given lie within the range [t0, t1] 
 *    return: 0 if ok, -1 if not
 *----------------------------------------------------------------------*/
int
check_timestamp_bounds (const LIGOTimeGPSVector *timestamps, LIGOTimeGPS t0, LIGOTimeGPS t1)
{
  LALStatus status1, status2;
  INT4 diff0_s, diff0_ns, diff1_s, diff1_ns;
  UINT4 i;
  LIGOTimeGPS *ti;

  status1 = status2 = emptyStatus;

  for (i = 0; i < timestamps->length; i ++)
    {
      ti = &(timestamps->data[i]);
      diff0_s = ti->gpsSeconds - t0.gpsSeconds;
      diff0_ns = ti->gpsNanoSeconds - t0.gpsNanoSeconds;
      if ( diff0_ns < 0 )
	{
	  diff0_ns += oneBillion;
	  diff0_s -= 1;
	}


      diff1_s = t1.gpsSeconds - ti->gpsSeconds;
      diff1_ns = t1.gpsNanoSeconds - ti->gpsNanoSeconds;
      if (diff1_ns < 0)
	{
	  diff1_ns += oneBillion;
	  diff0_s -= 1;
	}

      /* check timestamps-bounds */
      if ( (diff0_s < 0) || (diff1_s < 0) ) 
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
  REAL8 relError;
  volatile REAL8 bin1, bin2;	/* keep compiler from optimizing these away! */

  INITSTATUS( stat, "checkNoiseSFTs", GENERATEPULSARSIGNALC);

  for (i=0; i < sfts->length; i++)
    {
      thisSFT = &(sfts->data[i]);
      deltaFn = thisSFT->deltaF;
      fn0 = thisSFT->f0;
      fn1 = f0 + thisSFT->data->length * deltaFn;
      
      if (deltaFn != deltaF) {
	if (lalDebugLevel) 
	  LALPrintError ("\n\nTime-base of noise-SFTs Tsft_n=%f differs from signal-SFTs Tsft=%f\n", 1.0/deltaFn, 1.0/deltaF);
	ABORT (stat,  GENERATEPULSARSIGNALH_ENOISEDELTAF,  GENERATEPULSARSIGNALH_MSGENOISEDELTAF);
      }

      if ( (f0 < fn0) || (f1 > fn1) ) {
	if (lalDebugLevel) 
	  LALPrintError ("\n\nSignal frequency-band [%f,%f] is not contained in noise SFTs [%f,%f]\n", f0, f1, fn0, fn1);
	ABORT (stat, GENERATEPULSARSIGNALH_ENOISEBAND, GENERATEPULSARSIGNALH_MSGENOISEBAND);
      }
      
      bin1 = f0 / deltaF;	/* exact division if f is an integer frequency-bin! */
      bin2 = fn0 / deltaF;
      shift = bin1 - bin2;
      /* frequency bins have to coincide! ==> check that shift is integer!  */
      nshift = (UINT4)(shift+0.5);
      relError = fabs( nshift - shift);
      if ( relError > eps ) {
	if (lalDebugLevel) 
	  LALPrintError ("\n\nNoise frequency-bins don't coincide with signal-bins. Relative deviation=%g\n", relError);
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
      if (lalDebugLevel)
	LALPrintError ("Warning: we need to apply  heterodyning phase-correction now: correct_phase()\n");

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
