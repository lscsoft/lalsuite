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


#include "GeneratePulsarSignal.h"

/*----------------------------------------------------------------------*/
/* prototypes for internal functions */
static void
ConvertGPS2SSB (LALStatus* stat, LIGOTimeGPS *SSBout, LIGOTimeGPS GPSin, PulsarSignalParams *params);

/*----------------------------------------------------------------------*/

NRCSID( GENERATEPULSARSIGNAL, "$Id$");

extern INT4 lalDebugLevel;

void
LALGeneratePulsarSignal (LALStatus *stat, REAL4TimeSeries **signal, PulsarSignalParams *params)
{
  static SpinOrbitCWParamStruc spinorbit;
  static CoherentGW sourceSignal;
  LIGOTimeGPS time = {0,0};

  INITSTATUS( stat, "LALGeneratePulsarSignal", GENERATEPULSARSIGNAL );
  ATTATCHSTATUSPTR (stat);

  ASSERT (signal != NULL, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT (*signal == NULL, stat,  DOPPLERSCANH_ENONULL,  DOPPLERSCANH_MSGENONULL);

  /*----------------------------------------------------------------------*/
  /* Prepare call to GenerateSpinOrbitCW() */
  /*----------------------------------------------------------------------*/
  spinorbit.position.system = COORDINATESYSTEM_EQUATORIAL;
  spinorbit.position.longitude = params->pulsar.Alpha;
  spinorbit.position.latitude =  params->pulsar.Delta;
  spinorbit.psi = params->pulsar.psi;
  spinorbit.aPlus = params->pulsar.aPlus;
  spinorbit.aCross = params->pulsar.aCross;
  spinorbit.phi0 = params->pulsar.phi0;
  spinorbit.f0 = params->pulsar.f0;
  spinorbit.f = params->pulsar.f;

  /* if pulsar is in binary-orbit, set binary parameters */
  if (params->orbit)
    {
      spinorbit.orbitEpoch = params->orbit->orbitEpoch;
      spinorbit.omega = params->orbit->omega;
      spinorbit.rPeriNorm = params->orbit->rPeriNorm;
      spinorbit.oneMinusEcc = params->orbit->oneMinusEcc;
      spinorbit.angularSpeed = params->orbit->angularSpeed;
    }
  else
    spinorbit.rPeriNorm = 0.0;		/* this defines an isolated pulsar */

  /* start-time in GPS time */
  spinorbit.epoch = params->startTimeGPS;
  /* pulsar reference-time in SSB frame !*/
  if (params->pulsar.TRefSSB.gpsSeconds != 0)
    spinorbit.spinEpoch = params->pulsar.TRefSSB;
  else
    {
      TRY (ConvertGPS2SSB (stat->statusPtr, &time, params->startTimeGPS, params), stat);
      spinorbit.spinEpoch = time;
    }
  /* sampling-timestep and length for source-parameters */
  spinorbit.deltaT = 60;	/* in seconds; hardcoded default from makefakedata_v2 */
  spinorbit.length = (UINT4)( 1.0* params->duration / spinorbit.deltaT );

  /*
   * finally, call the function to generate the source waveform 
   */
  TRY ( LALGenerateSpinOrbitCW (stat->statusPtr, &sourceSignal, &spinorbit), stat);
  /* check that sampling interval was short enough */
  if ( spinorbit.dfdt > 2.0 )  /* taken from makefakedata_v2 */
    {
      LALPrintError ("GenerateSpinOrbitCW() returned df*dt = %f > 2.0", spinorbit.dfdt);
      ABORT (stat, DOPPLERSCANH_ESAMPLING, DOPPLERSCANH_MSGESAMPLING);
    }
  /*----------------------------------------------------------------------*/


  
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


  INITSTATUS( stat, "ConvertGPS2SSB", GENERATEPULSARSIGNAL );
  ATTATCHSTATUSPTR (stat);

  ASSERT (SSBout != NULL, stat,  DOPPLERSCANH_ENULL,  DOPPLERSCANH_MSGENULL);
  ASSERT (params != NULL, stat,  DOPPLERSCANH_ENULL,  DOPPLERSCANH_MSGENULL);

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

  printf ("\nDEBUG: checking  ConvertGPS2SSB(): difference = %d s, %d ns\n", 
	  ssb.gpsSeconds - emit.te.gpsSeconds, ssb.gpsNanoSeconds - emit.te.gpsNanoSeconds );

  *SSBout = emit.te;

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* ConvertGPS2SSB() */
