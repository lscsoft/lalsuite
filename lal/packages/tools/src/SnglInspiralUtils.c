/*----------------------------------------------------------------------- 
 * 
 * File Name: SnglInspiralUtils.c
 *
 * Author: Brady, P. R., Brown, D. A., and Fairhurst, S
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="SnglInspiralUtilsCV">
Author: Brown, D. A.
$Id$
</lalVerbatim> 
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>

NRCSID( SNGLINSPIRALUTILSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{SnglInspiralUtils.c}}

\noindent Blah.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SnglInspiralUtilsCP}
\idx{LALSortSnglInspiral()}
\idx{LALCompareSnglInspiralByMass()}
\idx{LALCompareSnglInspiralByPsi()}
\idx{LALCompareSnglInspiralByTime()}

\subsubsection*{Description}

\noindent Blah.

\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent None.

\subsubsection*{Notes}
%% Any relevant notes.

\vfill{\footnotesize\input{SnglInspiralUtilsCV}}

</lalLaTeX>
#endif

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALSortSnglInspiral (
    LALStatus          *status,
    SnglInspiralTable **eventHead,
    int(*comparfunc)    (const void *, const void *)
    )
/* </lalVerbatim> */
{
  INT4                  i;
  INT4                  numEvents = 0;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable   **eventHandle = NULL;

  INITSTATUS( status, "LALSortSnglInspiral", SNGLINSPIRALUTILSC );

  /* count the number of events in the linked list */
  for ( thisEvent = *eventHead; thisEvent; thisEvent = thisEvent->next )
  {
    ++numEvents;
  }
  if ( ! numEvents )
  {
    LALWarning( status, "No events in list to sort" );
    RETURN( status );
  }

  /* allocate memory for an array of pts to sort and populate array */
  eventHandle = (SnglInspiralTable **) 
    LALCalloc( numEvents, sizeof(SnglInspiralTable *) );
  for ( i = 0, thisEvent = *eventHead; i < numEvents; 
      ++i, thisEvent = thisEvent->next )
  {
    eventHandle[i] = thisEvent;
  }

  /* qsort the array using the specified function */
  qsort( eventHandle, numEvents, sizeof(eventHandle[0]), comparfunc );

  /* re-link the linked list in the right order */
  thisEvent = *eventHead = eventHandle[0];
  for ( i = 1; i < numEvents; ++i )
  {
    thisEvent = thisEvent->next = eventHandle[i];
  }
  thisEvent->next = NULL;

  /* free the internal memory */
  LALFree( eventHandle );

  RETURN( status );
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
LALCompareSnglInspiralByMass (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  SnglInspiralTable *aPtr = *((SnglInspiralTable **)a);
  SnglInspiralTable *bPtr = *((SnglInspiralTable **)b);

  if ( aPtr->mass1 > bPtr->mass1 )
  {
    return 1;
  }
  else if ( aPtr->mass1 < bPtr->mass1 )
  {
    return -1;
  }
  else if ( aPtr->mass2 > bPtr->mass2 )
  {
    return 1;
  }
  else if ( aPtr->mass2 < bPtr->mass2 )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
LALCompareSnglInspiralByPsi (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  SnglInspiralTable *aPtr = *((SnglInspiralTable **)a);
  SnglInspiralTable *bPtr = *((SnglInspiralTable **)b);

  if ( aPtr->psi0 > bPtr->psi0 )
  {
    return 1;
  }
  else if ( aPtr->psi0 < bPtr->psi0 )
  {
    return -1;
  }
  else if ( aPtr->psi3 > bPtr->psi3 )
  {
    return 1;
  }
  else if ( aPtr->psi3 < bPtr->psi3 )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}



/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
LALCompareSnglInspiralByTime (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  LALStatus     status;
  SnglInspiralTable *aPtr = *((SnglInspiralTable **)a);
  SnglInspiralTable *bPtr = *((SnglInspiralTable **)b);
  INT8 ta, tb;

  memset( &status, 0, sizeof(LALStatus) );
  LALGPStoINT8( &status, &ta, &(aPtr->end_time) );
  LALGPStoINT8( &status, &tb, &(bPtr->end_time) );

  if ( ta > tb )
  {
    return 1;
  }
  else if ( ta < tb )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALCompareSnglInspiral (
    LALStatus                *status,
    SnglInspiralTable        *aPtr,
    SnglInspiralTable        *bPtr,
    SnglInspiralAccuracy     *params
    )
/* </lalVerbatim> */
{
  INT8 ta, tb;
  REAL4 dm1, dm2;
  REAL4 dmchirp, deta;
  REAL4 dpsi0, dpsi3;
  REAL4 sigmaRatio;

  INITSTATUS( status, "LALCompareSnglInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  params->match = 0;

  LALGPStoINT8( status->statusPtr, &ta, &(aPtr->end_time) );
  LALGPStoINT8( status->statusPtr, &tb, &(bPtr->end_time) );

  /* compare on trigger time coincidence */
  if ( labs( ta - tb ) < params->dt )
  {
    if ( params->test == psi0_and_psi3 )
    {
      dpsi0 = fabs( aPtr->psi0 - bPtr->psi0 );
      dpsi3 = fabs( aPtr->psi3 - bPtr->psi3 );

      if ( dpsi0 <= params->dpsi0 && dpsi3 <= params->dpsi3 )
      {
	params->match = 1;
	LALInfo( status, "Triggers are coincident in psi0 and psi3" );
      }
      else
      {
	LALInfo( status, "Triggers are not coincident in psi0 and psi3" );
      }
    }
    else if ( params->test == m1_and_m2 )
    {  
      dm1 = fabs( aPtr->mass1 - bPtr->mass1 );
      dm2 = fabs( aPtr->mass2 - bPtr->mass2 );

      /* compare mass1 and mass2 parameters */
      if ( dm1 <= params->dm && dm2 <= params->dm )
      {
	if ( fabs( 
	      (aPtr->eff_distance - bPtr->eff_distance) / aPtr->eff_distance
	      ) < params->epsilon / bPtr->snr + params->kappa )
	{
	  params->match = 1;
	  LALInfo( status, "Triggers are coincident in eff_distance" );
	}
	else
	{
	  LALInfo( status, "Triggers fail eff_distance coincidence test" );
	}
      }
      else
      {
	LALInfo( status, "Triggers fail mass coincidence test" );
      }
    }
    else if ( params->test == mchirp_and_eta )
    {  
      dmchirp = fabs( aPtr->mchirp - bPtr->mchirp );
      deta = fabs( aPtr->eta - bPtr->eta );

      /* compare mchirp and eta parameters */
      if ( dmchirp <= params->dmchirp && deta <= params->deta )
      {
	params->match = 1;
	LALInfo( status, "Triggers are coincident in mchirp and eta" );
      }
      else
      {
	LALInfo( status, "Triggers fail mchirp, eta coincidence test" );
      }
    }
    else
    {
      LALInfo( status, "error: unknown test\n" );
    }
  }
  else
  {
    LALInfo( status, "Triggers fail time coincidence test" );
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALClusterSnglInspiralTable (
    LALStatus                  *status,
    SnglInspiralTable          *inspiralEvent,
    INT8                        dtimeNS,
    SnglInspiralClusterChoice   clusterchoice
    )
/* </lalVerbatim> */
{
  SnglInspiralTable     *thisEvent=NULL;
  SnglInspiralTable     *prevEvent=NULL;

  INITSTATUS( status, "LALClusterSnglInspiralTable", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  thisEvent = inspiralEvent->next;
  prevEvent = inspiralEvent;

  while ( thisEvent )
  {
    INT8 currTime;
    INT8 prevTime;

    /* compute the time in nanosec for each event trigger */
    LALGPStoINT8(status->statusPtr, &currTime, &(thisEvent->end_time));
    CHECKSTATUSPTR(status);

    LALGPStoINT8(status->statusPtr, &prevTime, &(prevEvent->end_time));
    CHECKSTATUSPTR(status);

    /* find events within the cluster window */
    if ( (currTime - prevTime) < dtimeNS )
    {
      /* displace previous event in cluster */
      if ( 
	  (
	   (clusterchoice == snr_and_chisq) && 
	   (thisEvent->snr > prevEvent->snr) && 
	   (thisEvent->chisq < prevEvent->chisq)
	  ) || (
	    (clusterchoice == snrsq_over_chisq) &&
	    (thisEvent->snr)*(thisEvent->snr)/(thisEvent->chisq) > 
	    (prevEvent->snr)*(prevEvent->snr)/(prevEvent->chisq)
	    ) || (
	      (clusterchoice == snr) && (thisEvent->snr > prevEvent->snr)
	      )
	 )
      {
	memcpy( prevEvent, thisEvent, sizeof(SnglInspiralTable) );
      }

      /* otherwise just dump this event from cluster */
      prevEvent->next = thisEvent->next;
      LALFree( thisEvent );
      thisEvent = prevEvent->next;
    }
    else 
    {
      /* otherwise we keep this unique event trigger */
      prevEvent = thisEvent;
      thisEvent = thisEvent->next;
    }
  }

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALGalacticInspiralParamsToSimInspiralTable(
    LALStatus                  *status,
    SimInspiralTable           *output,
    GalacticInspiralParamStruc *input,
    RandomParams               *params
    )
/* </lalVerbatim> */
{
  PPNParamStruc         ppnParams;
  LALMSTUnitsAndAcc     gmstUnits = { MST_HRS, LALLEAPSEC_STRICT };
  LALGPSandAcc          gpsAndAcc;
  SkyPosition	        skyPos;
  LALSource             source;
  LALPlaceAndGPS        placeAndGPS;
  DetTimeAndASource     detTimeAndSource;
  LALDetector           lho = lalCachedDetectors[LALDetectorIndexLHODIFF];
  LALDetector           llo = lalCachedDetectors[LALDetectorIndexLLODIFF];
  LALDetAndSource       detAndSource;
  LALDetAMResponse      resp;
  REAL8			time_diff_ns;
  REAL4                 splus, scross, cosiota;

  INITSTATUS( status, "LALGalacticParamsToSimInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( output, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( input, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( params, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );


  /*
   *
   * compute sky position and inspiral params
   *
   */


  /* generate the ppn inspiral params */
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  LALGetInspiralParams( status->statusPtr, &ppnParams, input, params );
  CHECKSTATUSPTR( status );

  if ( ppnParams.position.system != COORDINATESYSTEM_EQUATORIAL )
  {
    ABORT( status, LIGOMETADATAUTILSH_ECOOR, LIGOMETADATAUTILSH_MSGECOOR );
  }

  /* copy the inspiral data into sim_inspiral table */
  output->mass1        = input->m1;
  output->mass2        = input->m2;
  output->eta          = ppnParams.eta;
  output->distance     = ppnParams.d / (1.0e6 * LAL_PC_SI); /* Mpc */
  output->longitude    = ppnParams.position.longitude;
  output->latitude     = ppnParams.position.latitude;
  output->inclination  = ppnParams.inc;
  output->coa_phase    = ppnParams.phi;
  output->polarization = ppnParams.psi;

  /* populate geocentric end time */
  output->geocent_end_time = input->geocentEndTime;

  /* populate gmst field */
  LALGPStoGMST1( status->statusPtr, &(output->end_time_gmst), 
      &(output->geocent_end_time), &gmstUnits );
  CHECKSTATUSPTR( status );

  /* set up params for the site end times and detector response */
  memset( &skyPos, 0, sizeof(SkyPosition) );
  memset( &source, 0, sizeof(LALSource) );
  memset( &placeAndGPS, 0, sizeof(LALPlaceAndGPS) );
  memset( &detTimeAndSource, 0, sizeof(DetTimeAndASource) );
  memset( &detAndSource, 0, sizeof(LALDetAndSource) );

  skyPos.longitude = output->longitude;
  skyPos.latitude  = output->latitude;
  skyPos.system    = COORDINATESYSTEM_EQUATORIAL;

  source.equatorialCoords = skyPos;
  source.orientation      = output->polarization;

  placeAndGPS.p_gps = &(output->geocent_end_time);

  detTimeAndSource.p_det_and_time = &placeAndGPS;
  detTimeAndSource.p_source = &skyPos;

  detAndSource.pSource = &source;

  gpsAndAcc.accuracy = LALLEAPSEC_STRICT;
  gpsAndAcc.gps = output->geocent_end_time;


  /*
   *
   * compute site end times
   *
   */


  /* initialize end times with geocentric value */
  output->h_end_time = output->l_end_time = input->geocentEndTime;

  /* ligo hanford observatory */
  placeAndGPS.p_detector = &lho;
  LALTimeDelayFromEarthCenter( status->statusPtr, &time_diff_ns, 
      &detTimeAndSource );
  CHECKSTATUSPTR( status );
  LALAddFloatToGPS( status->statusPtr, &(output->h_end_time),
      &(output->h_end_time), time_diff_ns );
  CHECKSTATUSPTR( status );

  /* ligo livingston observatory */
  placeAndGPS.p_detector = &llo;
  LALTimeDelayFromEarthCenter( status->statusPtr, &time_diff_ns, 
      &detTimeAndSource );
  CHECKSTATUSPTR( status );
  LALAddFloatToGPS( status->statusPtr, &(output->l_end_time), 
      &(output->l_end_time), time_diff_ns );
  CHECKSTATUSPTR( status );


  /*
   *
   * compute the effective distance of the inspiral
   *
   */


  /* initialize distances with real distance and compute splus and scross */
  output->eff_dist_h = output->eff_dist_l = 2.0 * output->distance;
  cosiota = cos( output->inclination );
  splus = -( 1.0 + cosiota * cosiota );
  scross = -2.0 * cosiota;

  /* compute the response of the LHO detectors */
  detAndSource.pDetector = &lho;
  LALComputeDetAMResponse( status->statusPtr, &resp, &detAndSource,
      &gpsAndAcc );
  CHECKSTATUSPTR( status );

  /* compute the effective distance for LHO */
  output->eff_dist_h /= sqrt( 
      splus*splus*resp.plus*resp.plus + scross*scross*resp.cross*resp.cross );

  /* compute the response of the LLO detector */
  detAndSource.pDetector = &llo;
  LALComputeDetAMResponse( status->statusPtr, &resp, &detAndSource,
      &gpsAndAcc );
  CHECKSTATUSPTR( status );

  /* compute the effective distance for LLO */
  output->eff_dist_l /= sqrt( 
      splus*splus*resp.plus*resp.plus + scross*scross*resp.cross*resp.cross );


  /* 
   *
   * normal exit
   *
   */


  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/* function to compute the site end time and effective distance of an event */
static int site_time_and_dist( 
    LALStatus         *status,
    LALDetector       *detector,
    LIGOTimeGPS       *end_time,
    REAL4             *eff_dist,
    SkyPosition       *skyPos,
    SimInspiralTable  *sim_inspiral)
{
  LALMSTUnitsAndAcc     gmstUnits = { MST_HRS, LALLEAPSEC_STRICT };
  LALGPSandAcc          gpsAndAcc;
  LALSource             source;
  LALPlaceAndGPS        placeAndGPS;
  DetTimeAndASource     detTimeAndSource;
  LALDetAndSource       detAndSource;
  LALDetAMResponse      resp;
  REAL8			time_diff_ns;
  REAL4                 splus, scross, cosiota;

  memset( &source, 0, sizeof(LALSource) );
  memset( &placeAndGPS, 0, sizeof(LALPlaceAndGPS) );
  memset( &detTimeAndSource, 0, sizeof(DetTimeAndASource) );
  memset( &detAndSource, 0, sizeof(LALDetAndSource) );


  source.equatorialCoords = *skyPos;
  source.orientation      = sim_inspiral->polarization;

  placeAndGPS.p_gps = &(sim_inspiral->geocent_end_time);

  detTimeAndSource.p_det_and_time = &placeAndGPS;
  detTimeAndSource.p_source = skyPos;
  detTimeAndSource.p_det_and_time->p_detector = detector;

  detAndSource.pSource = &source;
  detAndSource.pDetector = detector;

  gpsAndAcc.accuracy = LALLEAPSEC_STRICT;
  gpsAndAcc.gps = sim_inspiral->geocent_end_time;

  /* initialize end time with geocentric value */
  *end_time = sim_inspiral->geocent_end_time;

  /* calculate the detector end time */
  LALTimeDelayFromEarthCenter( status->statusPtr, &time_diff_ns, 
      &detTimeAndSource );
  CHECKSTATUSPTR( status );
  LALAddFloatToGPS( status->statusPtr, end_time,
      end_time, time_diff_ns );
  CHECKSTATUSPTR( status );

  /* initialize distance with real distance and compute splus and scross */
  *eff_dist = 2.0 * sim_inspiral->distance;
  cosiota = cos( sim_inspiral->inclination );
  splus = -( 1.0 + cosiota * cosiota );
  scross = -2.0 * cosiota;

  /* compute the response of the detector */
  LALComputeDetAMResponse( status->statusPtr, &resp, &detAndSource,
      &gpsAndAcc );
  CHECKSTATUSPTR( status );

  /* compute the effective distance */
  *eff_dist /= sqrt( 
      splus*splus*resp.plus*resp.plus + scross*scross*resp.cross*resp.cross );

  return(0);
}



/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALPopulateSimInspiralSiteInfo(
    LALStatus                  *status,
    SimInspiralTable           *output
    )
/* </lalVerbatim> */
{
  SkyPosition	        skyPos;
  LALDetector           detector; 
  REAL4                *eff_dist;
  LIGOTimeGPS          *end_time;


  INITSTATUS( status, "LALPopulateSimInspiralSiteInfo", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( output, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  /* set up params for the geocent end time and source location */
  memset( &skyPos, 0, sizeof(SkyPosition) );

  skyPos.longitude = output->longitude;
  skyPos.latitude  = output->latitude;
  skyPos.system    = COORDINATESYSTEM_EQUATORIAL;

  /* LIGO Hanford observatory*/
  detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  end_time = &(output->h_end_time);
  eff_dist = &(output->eff_dist_h);
  site_time_and_dist(status, &detector, end_time, eff_dist, 
      &skyPos, output);

  /* LIGO Livingston observatory*/
  detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  end_time = &(output->l_end_time);
  eff_dist = &(output->eff_dist_l);
  site_time_and_dist(status, &detector, end_time, eff_dist, 
      &skyPos, output);

  /* GEO observatory*/
  detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  end_time = &(output->g_end_time);
  eff_dist = &(output->eff_dist_g);
  site_time_and_dist(status, &detector, end_time, eff_dist, 
      &skyPos, output);

  /* TAMA observatory*/
  detector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  end_time = &(output->t_end_time);
  eff_dist = &(output->eff_dist_t);
  site_time_and_dist(status, &detector, end_time, eff_dist, 
      &skyPos, output);

  /* Virgo observatory*/
  detector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  end_time = &(output->v_end_time);
  eff_dist = &(output->eff_dist_v);
  site_time_and_dist(status, &detector, end_time, eff_dist, 
      &skyPos, output);

  /* 
   *
   * normal exit
   *
   */


  DETATCHSTATUSPTR (status);
  RETURN (status);
}
