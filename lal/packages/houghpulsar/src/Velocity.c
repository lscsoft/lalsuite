/*****************************************************************
 *
 *  File Name: Velocity.c
 *
 *  $Id$
 *
 *  Authors: Badri Krishnan, Alicia Sintes
 *           May 2003
 *
 *
 ******************************************************************/

/************************************<lalVerbatim file="VelocityCV">
Authors: Krishnan, B., Sintes, A.M.
$Id$
*************************************</lalVerbatim> */

/* <lalLaTeX>  **********************************************

\subsection{Module \texttt{Velocity.c}}
\label{ss:Velocity.c}
Computation of instant and averaged velocities of a given detector and the
like.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{VelocityD}
\index{\verb&LALDetectorVel()&}
\index{\verb&LALAvgDetectorVel()&}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}
The function \verb@LALDetectorVel@ finds the velocity of a given detector at a
given time. It is basically a wrapper for LALBarycenter.
The output is of the form \verb@REAL8  v[3]@, and the input is a
time \verb@LIGOTimeGPS  *time@, the detector \verb@LALDetector detector, 
and the ephemeris \verb@EphemerisData *edat from LALInitBarycenter

The function \verb@LALAvgDetectorVel@ outputs the
   average velocity \verb@REAL8 v[3]@ of the detector during a time interval by using
   the trapeziodal rule. The input structure is of type \verb@VelocityPar *in@
   containing all the required parmaters. \verb@vTol@ is the fractional accuracy 
   required for velocity. Recomended value $0.01$.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
\begin{verbatim}
LALFloatToGPS()
LALMalloc() 
LALBarycenterEarth() 
LALBarycenter() 
LALFree()
\end{verbatim}

%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%
\vfill{\footnotesize\input{VelocityCV}}
********************************************** </lalLaTeX> */


#include <lal/Velocity.h> 
/* #include "./Velocity.h" */

NRCSID (VELOCITYC, "$Id$");

/*
 * The functions that make up the guts of this module
 */

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* Given a detector and a time interval, this function outputs the 
   average velocity of the detector during the interval by using 
   the trapeziodal rule */
   
/* *******************************  <lalVerbatim file="VelocityD"> */   
void LALAvgDetectorVel( LALStatus *status, 
		        REAL8 v[3], 
		        VelocityPar *in)
{ /*-------------------------------------------------</lalVerbatim> */

  REAL8           trapSum[3], average[3], tempVel[3]; 
  REAL8           tBase, vTol, rTol, oldVeloMod, newVeloMod;
  LIGOTimeGPS     t0gps, tgps;
  REAL8           t0; 
  REAL8           t, ts, tn; 
  INT4            n, k, points; 
  LALDetector     detector;
  EphemerisData   *edat;

  /*---------------------------------------------------------------*/

  INITSTATUS (status, "LALAvgDetectorVel", VELOCITYC);
  ATTATCHSTATUSPTR (status);

  /* check that arguments are not null */
  ASSERT(in, status, VELOCITYH_ENULL, VELOCITYH_MSGENULL);

  /* copy input values */
  tBase = in->tBase;
  vTol = in->vTol;
  t0gps = in->startTime;
  detector = in->detector;
  edat = in->edat;

  /* basic checks that input values are sane */
  ASSERT(tBase > 0.0, status, VELOCITYH_EVAL, VELOCITYH_MSGEVAL);
  ASSERT(vTol > 0.0, status, VELOCITYH_EVAL, VELOCITYH_MSGEVAL);
  ASSERT(vTol < 1.0, status, VELOCITYH_EVAL, VELOCITYH_MSGEVAL);

  /* calculate starting time as float */
  ts = (REAL8)(t0gps.gpsSeconds) * 1.00;
  tn = (REAL8)(t0gps.gpsNanoSeconds) * 1.00E-9;
  t0 = ts + tn;

  /* calculate finish time */
  t = t0 + tBase;
  TRY( LALFloatToGPS( status->statusPtr, &tgps, &t), status);

  /* The first approximation: (b-a)^-1 * int(f(x),x=a..b) = 0.5*(f(a)+f(b)) */
  /* calculate velocity at starting time */
  TRY( LALDetectorVel( status->statusPtr, tempVel, &t0gps, detector, edat), status);
  for (n=0; n<3; n++) trapSum[n] = 0.5 * tempVel[n]; 

  /*calculate velocity at finish time */
  
  TRY( LALDetectorVel( status->statusPtr, tempVel, &tgps, detector, edat), status);
  
  /* first approximation to average */  
  for (n=0; n<3; n++) 
    {
      trapSum[n] += 0.5 * tempVel[n];
      average[n] = trapSum[n];
    }
  
  points = 1;
  /* now add more points and stop when desired accuracy is reached*/ 
  do {
    points *= 2;
    for (k=1; k<points; k+=2)
      {
	t = t0 + 1.0 * k * tBase / (1.0 * points);
        TRY( LALFloatToGPS( status->statusPtr, &tgps, &t), status);
	TRY( LALDetectorVel( status->statusPtr, tempVel, &tgps, detector, edat), status);
	for (n=0; n<3; n++) trapSum[n] += tempVel[n];
      }
    oldVeloMod = newVeloMod = 0.0;
    for (n=0; n<3; n++) 
      {
	oldVeloMod += average[n]*average[n];
	average[n] = trapSum[n] / (1.0*points);
        newVeloMod += average[n]*average[n];
      }
    /* now calculate the fractional change in magnitude of average */
    /* is it sufficient to require the magnitude to converge or should
       we look at the individual components? */
    rTol = fabs((sqrt(oldVeloMod) - sqrt(newVeloMod))) / (sqrt(oldVeloMod));      
  } while (rTol > vTol);
    
  /* copy the result to the output structure */
  for (n=0; n<3; n++) v[n] = average[n]; 

  DETATCHSTATUSPTR (status);

  RETURN (status);
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* Given a detector and a time interval, this function outputs the 
   average position of the detector during the interval by using 
   the trapeziodal rule */
   
/* *******************************  <lalVerbatim file="VelocityD"> */   
void LALAvgDetectorPos( LALStatus *status, 
		        REAL8 x[3], 
		        VelocityPar *in)
{ /*-------------------------------------------------</lalVerbatim> */

  REAL8           trapSum[3], average[3], tempVel[3]; 
  REAL8           tBase, vTol, rTol, oldVeloMod, newVeloMod;
  LIGOTimeGPS     t0gps, tgps;
  REAL8           t0; 
  REAL8           t, ts, tn; 
  INT4            n, k, points; 
  LALDetector     detector;
  EphemerisData   *edat;



  /*---------------------------------------------------------------*/

  INITSTATUS (status, "LALAvgDetectorPos", VELOCITYC);
  ATTATCHSTATUSPTR (status);

  /* check that arguments are not null */
  ASSERT(in, status, VELOCITYH_ENULL, VELOCITYH_MSGENULL);

  /* copy input values */
  tBase = in->tBase;
  vTol = in->vTol;
  t0gps = in->startTime;
  detector = in->detector;
  edat = in->edat;

  /* check that input values make sense */
  ASSERT(tBase > 0.0, status, VELOCITYH_EVAL, VELOCITYH_MSGEVAL);
  ASSERT(vTol > 0.0, status, VELOCITYH_EVAL, VELOCITYH_MSGEVAL);
  ASSERT(vTol < 1.0, status, VELOCITYH_EVAL, VELOCITYH_MSGEVAL);

  /* calculate starting time as float */
  ts = (REAL8)(t0gps.gpsSeconds) * 1.00;
  tn = (REAL8)(t0gps.gpsNanoSeconds) * 1.00E-9;
  t0 = ts + tn;

  /* calculate finish time */
  t = t0 + tBase;
  TRY( LALFloatToGPS( status->statusPtr, &tgps, &t), status);

  /* The first approximation: (b-a)^-1 * int(f(x),x=a..b) = 0.5*(f(a)+f(b)) */
  /* calculate velocity at starting time */
  TRY( LALDetectorVel( status->statusPtr, tempVel, &t0gps, detector, edat), status);
  for (n=0; n<3; n++) trapSum[n] = 0.5 * tempVel[n]; 

  /*calculate velocity at finish time */
  
  TRY( LALDetectorVel( status->statusPtr, tempVel, &tgps, detector, edat), status);
  
  /* first approximation to average */  
  for (n=0; n<3; n++) 
    {
      trapSum[n] += 0.5 * tempVel[n];
      average[n] = trapSum[n];
    }
  
  points = 1;
  /* now add more points and stop when desired accuracy is reached*/ 
  do {
    points *= 2;
    for (k=1; k<points; k+=2)
      {
	t = t0 + 1.0 * k * tBase / (1.0 * points);
        TRY( LALFloatToGPS( status->statusPtr, &tgps, &t), status);
	TRY( LALDetectorPos( status->statusPtr, tempVel, &tgps, detector, edat), status);
	for (n=0; n<3; n++) trapSum[n] += tempVel[n];
      }
    oldVeloMod = newVeloMod = 0.0;
    for (n=0; n<3; n++) 
      {
	oldVeloMod += average[n]*average[n];
	average[n] = trapSum[n] / (1.0*points);
        newVeloMod += average[n]*average[n];
      }
    /* now calculate the fractional change in magnitude of average */
    /* is it sufficient to require the magnitude to converge or should
       we look at the individual components? */
    rTol = fabs((sqrt(oldVeloMod) - sqrt(newVeloMod))) / (sqrt(oldVeloMod));      
  } while (rTol > vTol);
    
  /* copy the result to the output structure */
  for (n=0; n<3; n++) x[n] = average[n]; 

  DETATCHSTATUSPTR (status);

  RETURN (status);
}



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* This finds velocity of a given detector at a given time */
/* *******************************  <lalVerbatim file="VelocityD"> */
void LALDetectorVel(LALStatus    *status, 
		    REAL8        v[3], 
		    LIGOTimeGPS  *time0, 
		    LALDetector  detector,
		    EphemerisData *edat)
{ /*-------------------------------------------------</lalVerbatim> */

  INT4  i;
  EmissionTime     *emit=NULL;
  EarthState       *earth=NULL;
  BarycenterInput  baryinput;

  /*----------------------------------------------------------------*/
 
  INITSTATUS ( status, "LALDetectorVel", VELOCITYC);
  ATTATCHSTATUSPTR (status);

  ASSERT (time0, status, VELOCITYH_ENULL, VELOCITYH_MSGENULL);
  ASSERT (edat, status, VELOCITYH_ENULL, VELOCITYH_MSGENULL);

  emit = (EmissionTime *)LALMalloc(sizeof(EmissionTime));
  earth = (EarthState *)LALMalloc(sizeof(EarthState));

  /* detector info */
  baryinput.site.location[0] = detector.location[0]/LAL_C_SI;
  baryinput.site.location[1] = detector.location[1]/LAL_C_SI;
  baryinput.site.location[2] = detector.location[2]/LAL_C_SI;
  
  /* set other barycentering info */
  baryinput.tgps = *time0;

  /* for the purposes of calculating the velocity of the earth */
  /* at some given time, the position of the source in the sky should not matter. */
  /* So, for this function, we set alpha and delta to zero as inputs to the */
  /* barycenter routines  */  
  baryinput.alpha = 0.0;
  baryinput.delta = 0.0; 
  
  /* call barycentering routines to calculate velocities */
  TRY( LALBarycenterEarth( status->statusPtr, earth, time0, edat), status);
  TRY( LALBarycenter( status->statusPtr, emit, &baryinput, earth), status);
  
  /* set values of velocity for all the SFT's */
  for (i=0; i < 3; i++) { v[i] = emit->vDetector[i]; }

  LALFree(emit);
  LALFree(earth); 

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* This finds velocity of a given detector at a given time */
/* *******************************  <lalVerbatim file="VelocityD"> */
void LALDetectorPos(LALStatus    *status, 
		    REAL8        x[3], 
		    LIGOTimeGPS  *time0, 
		    LALDetector  detector,
		    EphemerisData *edat)
{ /*-------------------------------------------------</lalVerbatim> */

  INT4  i;
  EmissionTime     *emit=NULL;
  EarthState       *earth=NULL;
  BarycenterInput  baryinput;

  /*----------------------------------------------------------------*/
 
  INITSTATUS ( status, "LALDetectorPos", VELOCITYC);
  ATTATCHSTATUSPTR (status);

  ASSERT (time0, status, VELOCITYH_ENULL, VELOCITYH_MSGENULL);
  ASSERT (edat, status, VELOCITYH_ENULL, VELOCITYH_MSGENULL);

  emit = (EmissionTime *)LALMalloc(sizeof(EmissionTime));
  earth = (EarthState *)LALMalloc(sizeof(EarthState));

  /* detector info */
  baryinput.site.location[0] = detector.location[0]/LAL_C_SI;
  baryinput.site.location[1] = detector.location[1]/LAL_C_SI;
  baryinput.site.location[2] = detector.location[2]/LAL_C_SI;
  
  /* set other barycentering info */
  baryinput.tgps = *time0;

  /* for the purposes of calculating the velocity of the earth */
  /* at some given time, the position of the source in the sky should not matter. */
  /* So, for this function, we set alpha and delta to zero as inputs to the */
  /* barycenter routines  */  
  baryinput.alpha = 0.0;
  baryinput.delta = 0.0; 
  
  /* call barycentering routines to calculate velocities */
  TRY( LALBarycenterEarth( status->statusPtr, earth, time0, edat), status);
  TRY( LALBarycenter( status->statusPtr, emit, &baryinput, earth), status);
  
  /* set values of velocity for all the SFT's */
  for (i=0; i < 3; i++) { x[i] = emit->rDetector[i]; }

  LALFree(emit);
  LALFree(earth); 

  DETATCHSTATUSPTR (status);
  RETURN (status);
}









