/*
*  Copyright (C) 2007 Drew Keppel, George Birthisel, Patrick Brady, Peter Shawhan, Stephen Fairhurst
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

/*-----------------------------------------------------------------------
 *
 * File Name: SimInspiralUtils.c
 *
 * Author: Brady, P. R., Brown, D. A., and Fairhurst, S
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="SimInspiralUtilsCV">
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

NRCSID( SIMINSPIRALUTILSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{SimInspiralUtils.c}}

Provides a set of utilities for manipulating \texttt{simInspiralTable}s.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SimInspiralUtilsCP}
\idx{LALGalacticInspiralParamsToSimInspiralTable()}
\idx{LALInspiralSiteTimeAndDist()}
\idx{LALPopulateSimInspiralSiteInfo()}

\subsubsection*{Description}

The function \texttt{LALInspiralSiteTimeAndDist()} calculates detector end
time (\texttt{endTime}) and effective distance (\texttt{effDist}) for an
inspiral signal from a specific location in the sky (\texttt{skyPos}) assumed
to be given in equatorial coordinates.  The detector end time is obtained by
using \texttt{LALTimeDelayFromEarthCenter()}, while the effective distance
requires calculation of the detector response, calculated using
\texttt{LALComputeDetAMResponse()}.

The function \texttt{LALPopulateSimInspiralSiteInfo()} populates the end time
and effective distance for each of the interferometer sites.  The sky location
(in equatorial coordinates) is assumed to be already contained in the input
\texttt{SimInspiralTable}.  The end time and effective distance for each site
is calculated by calling \texttt{LALInspiralSiteTimeAndDist()} once for each
of the detectors, and setting the \texttt{detector} appropriately.



\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent LALGetInspiralParams, LALGPStoGMST1, LALTimeDelayFromEarthCenter,
  LALAddFloatToGPS, LALComputeDetAMResponse.

  \subsubsection*{Notes}
  %% Any relevant notes.

  \vfill{\footnotesize\input{SimInspiralUtilsCV}}

  </lalLaTeX>
#endif

  /* a few useful static functions */
static INT8 geocent_end_time(const SimInspiralTable *x)
{
  return(XLALGPSToINT8NS(&x->geocent_end_time));
}


/* <lalVerbatim file="SimInspiralUtilsCP"> */
void
XLALPlayTestSimInspiral(
    SimInspiralTable         **eventHead,
    LALPlaygroundDataMask      *dataType
    )
/* </lalVerbatim> */
{
  SimInspiralTable    *inspiralEventList = NULL;
  SimInspiralTable    *thisEvent = NULL;
  SimInspiralTable    *prevEvent = NULL;

  INT8 triggerTime = 0;
  INT4 isPlay = 0;
  INT4 numTriggers;

  /* Remove all the triggers which are not of the desired type */

  numTriggers = 0;
  thisEvent = *eventHead;

  if ( (*dataType == playground_only) || (*dataType == exclude_play) )
  {
    while ( thisEvent )
    {
      SimInspiralTable *tmpEvent = thisEvent;
      thisEvent = thisEvent->next;

      triggerTime = XLALGPSToINT8NS( &(tmpEvent->geocent_end_time) );
      isPlay = XLALINT8NanoSecIsPlayground( &triggerTime );

      if ( ( (*dataType == playground_only)  && isPlay ) ||
          ( (*dataType == exclude_play) && ! isPlay) )
      {
        /* keep this trigger */
        if ( ! inspiralEventList  )
        {
          inspiralEventList = tmpEvent;
        }
        else
        {
          prevEvent->next = tmpEvent;
        }
        tmpEvent->next = NULL;
        prevEvent = tmpEvent;
        ++numTriggers;
      }
      else
      {
        /* discard this template */
        XLALFreeSimInspiral ( &tmpEvent );
      }
    }
    *eventHead = inspiralEventList;
    if ( *dataType == playground_only )
    {
      XLALPrintInfo( "Kept %d playground triggers \n", numTriggers );
    }
    else if ( *dataType == exclude_play )
    {
      XLALPrintInfo( "Kept %d non-playground triggers \n", numTriggers );
    }
  }
  else if ( *dataType == all_data )
  {
    XLALPrintInfo(
        "XLALPlayTestSimInspiral: Keeping all triggers\n" );
  }
  else
  {
    XLALPrintInfo(
        "XLALPlayTestSimInspiral: Unknown data type, returning no triggers\n"
        );
    *eventHead = NULL;
  }

}


/* <lalVerbatim file="SimInspiralUtilsCP"> */
int
XLALSimInspiralInSearchedData(
    SimInspiralTable         **eventHead,
    SearchSummaryTable       **summList
    )
/* </lalVerbatim> */
{
  SearchSummaryTable   *thisSearchSumm = NULL;

  SimInspiralTable *eventList = NULL;
  SimInspiralTable *prevEvent = NULL;
  SimInspiralTable *thisEvent = NULL;

  int numInj = 0;

  XLALTimeSortSearchSummary( summList, LALCompareSearchSummaryByOutTime );
  XLALSortSimInspiral( eventHead, XLALCompareSimInspiralByGeocentEndTime );

  thisEvent = *eventHead;
  thisSearchSumm = *summList;


  while ( thisEvent )
  {
    SimInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    while ( thisSearchSumm )
    {
      if ( geocent_end_time(tmpEvent) <
          XLALGPSToINT8NS( &(thisSearchSumm->out_start_time) ))
      {
        XLALPrintInfo(
            "XLALSimInspiralInSearchedData: Discarding injection\n" );
        LALFree( tmpEvent );
        break;
      }
      else
      {
        if ( geocent_end_time(tmpEvent) <
            XLALGPSToINT8NS( &(thisSearchSumm->out_end_time) ))
        {
          XLALPrintInfo(
              "XLALSimInspiralInSearchedData: Keeping injection\n" );
          numInj++;
          if ( ! eventList )
          {
            eventList = tmpEvent;
          }
          else
          {
            prevEvent->next = tmpEvent;
          }
          tmpEvent->next = NULL;
          prevEvent = tmpEvent;
          break;
        }
      }

      thisSearchSumm = thisSearchSumm->next;
    }

    if ( !thisSearchSumm )
    {
      XLALPrintInfo(
          "XLALSimInspiralInSearchedData: Discarding injection\n" );
      LALFree( tmpEvent );
    }
  }

  *eventHead = eventList;

  return(numInj);
}


/* <lalVerbatim file="SimInspiralUtilsCP"> */
int
XLALSimInspiralChirpMassCut(
    SimInspiralTable   **eventHead,
    REAL4                minChirpMass,
    REAL4                maxChirpMass
    )
/* </lalVerbatim> */
{
  SimInspiralTable    *inspiralEventList = NULL;
  SimInspiralTable    *thisEvent = NULL;
  SimInspiralTable    *prevEvent = NULL;

  INT4 numInj;

  /* Remove all the triggers which are not of the desired type */

  numInj = 0;
  thisEvent = *eventHead;

  while ( thisEvent )
  {
    SimInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( ( tmpEvent->mchirp >= minChirpMass ) &&
         ( tmpEvent->mchirp < maxChirpMass ) )
    {
      /* keep this trigger */
      if ( ! inspiralEventList  )
      {
        inspiralEventList = tmpEvent;
      }
      else
      {
        prevEvent->next = tmpEvent;
      }
      tmpEvent->next = NULL;
      prevEvent = tmpEvent;
      ++numInj;
    }
    else
    {
      /* discard this template */
      XLALFreeSimInspiral ( &tmpEvent );
    }
  }

  *eventHead = inspiralEventList;

  return(numInj);
}


/* <lalVerbatim file="SimInspiralUtilsCP"> */
int
XLALSimInspiralCompMassCut(
    SimInspiralTable   **eventHead,
    REAL4                minCompMass,
    REAL4                maxCompMass,
    REAL4                minCompMass2,
    REAL4                maxCompMass2
    )
/* </lalVerbatim> */
{
  SimInspiralTable    *inspiralEventList = NULL;
  SimInspiralTable    *thisEvent = NULL;
  SimInspiralTable    *prevEvent = NULL;

  INT4 numInj;
  REAL4 mass1, mass2;


  /* Remove all the triggers which are not of the desired type */

  numInj = 0;
  thisEvent = *eventHead;

  while ( thisEvent )
  {
    SimInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( tmpEvent->mass1 >= tmpEvent->mass2 )
    {
      mass1 = tmpEvent->mass1;
      mass2 = tmpEvent->mass2;
    }
    else
    {
      mass1 = tmpEvent->mass2;
      mass2 = tmpEvent->mass1;
    }

    if ( ( mass1 >= minCompMass ) &&
         ( mass1 < maxCompMass ) &&
         ( mass2 >= minCompMass2 ) &&
         ( mass2 < maxCompMass2 ) )
    {
      /* keep this trigger */
      if ( ! inspiralEventList  )
      {
        inspiralEventList = tmpEvent;
      }
      else
      {
        prevEvent->next = tmpEvent;
      }
      tmpEvent->next = NULL;
      prevEvent = tmpEvent;
      ++numInj;
    }
    else
    {
      /* discard this template */
      XLALFreeSimInspiral ( &tmpEvent );
    }
  }

  *eventHead = inspiralEventList;

  return(numInj);
}


/* <lalVerbatim file="SimInspiralUtilsCP"> */
int
XLALSimInspiralTotalMassCut(
    SimInspiralTable   **eventHead,
    REAL4                minTotalMass,
    REAL4                maxTotalMass
    )
/* </lalVerbatim> */
{
  SimInspiralTable    *inspiralEventList = NULL;
  SimInspiralTable    *thisEvent = NULL;
  SimInspiralTable    *prevEvent = NULL;

  INT4 numInj;

  /* Remove all the triggers which are not of the desired type */

  numInj = 0;
  thisEvent = *eventHead;

  while ( thisEvent )
  {
    SimInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( ( (tmpEvent->mass1 + tmpEvent->mass2) >= minTotalMass ) &&
         ( (tmpEvent->mass1 + tmpEvent->mass2) < maxTotalMass ) )
    {
      /* keep this trigger */
      if ( ! inspiralEventList  )
      {
        inspiralEventList = tmpEvent;
      }
      else
      {
        prevEvent->next = tmpEvent;
      }
      tmpEvent->next = NULL;
      prevEvent = tmpEvent;
      ++numInj;
    }
    else
    {
      /* discard this template */
      XLALFreeSimInspiral ( &tmpEvent );
    }
  }

  *eventHead = inspiralEventList;

  return(numInj);
}


/* <lalVerbatim file="SimInspiralUtilsCP"> */
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
  SkyPosition           skyPos;
  LALSource             source;
  LALPlaceAndGPS        placeAndGPS;
  DetTimeAndASource     detTimeAndSource;
  LALDetector           lho = lalCachedDetectors[LALDetectorIndexLHODIFF];
  LALDetector           llo = lalCachedDetectors[LALDetectorIndexLLODIFF];
  LALDetAndSource       detAndSource;
  LALDetAMResponse      resp;
  REAL8     time_diff_ns;
  REAL4                 splus, scross, cosiota;

  INITSTATUS( status, "LALGalacticParamsToSimInspiral", SIMINSPIRALUTILSC );
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

/* <lalVerbatim file="SimInspiralUtilsCP"> */
void
LALInspiralSiteTimeAndDist(
    LALStatus         *status,
    SimInspiralTable  *output,
    LALDetector       *detector,
    LIGOTimeGPS       *endTime,
    REAL4             *effDist,
    SkyPosition       *skyPos
    )
/* </lalVerbatim> */
{
  LALGPSandAcc          gpsAndAcc;
  LALSource             source;
  LALPlaceAndGPS        placeAndGPS;
  DetTimeAndASource     detTimeAndSource;
  LALDetAndSource       detAndSource;
  LALDetAMResponse      resp;
  REAL8                 time_diff_ns;
  REAL4                 splus, scross, cosiota;

  INITSTATUS( status, "LALInspiralSiteTimeAndDist", SIMINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  /* check that the arguments are not null */
  ASSERT( output, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( detector, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( endTime, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( effDist, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( skyPos, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  memset( &source, 0, sizeof(LALSource) );
  memset( &placeAndGPS, 0, sizeof(LALPlaceAndGPS) );
  memset( &detTimeAndSource, 0, sizeof(DetTimeAndASource) );
  memset( &detAndSource, 0, sizeof(LALDetAndSource) );


  source.equatorialCoords = *skyPos;
  source.orientation      = output->polarization;

  placeAndGPS.p_gps = &(output->geocent_end_time);

  detTimeAndSource.p_det_and_time = &placeAndGPS;
  detTimeAndSource.p_source = skyPos;
  detTimeAndSource.p_det_and_time->p_detector = detector;

  detAndSource.pSource = &source;
  detAndSource.pDetector = detector;

  gpsAndAcc.accuracy = LALLEAPSEC_STRICT;
  gpsAndAcc.gps = output->geocent_end_time;

  /* initialize end time with geocentric value */
  *endTime = output->geocent_end_time;

  /* calculate the detector end time */
  LALTimeDelayFromEarthCenter( status->statusPtr, &time_diff_ns,
      &detTimeAndSource );
  CHECKSTATUSPTR( status );
  LALAddFloatToGPS( status->statusPtr, endTime,
      endTime, time_diff_ns );
  CHECKSTATUSPTR( status );

  /* initialize distance with real distance and compute splus and scross */
  *effDist = 2.0 * output->distance;
  cosiota = cos( output->inclination );
  splus = -( 1.0 + cosiota * cosiota );
  scross = -2.0 * cosiota;

  /* compute the response of the detector */
  LALComputeDetAMResponse( status->statusPtr, &resp, &detAndSource,
      &gpsAndAcc );
  CHECKSTATUSPTR( status );

  /* compute the effective distance */
  *effDist /= sqrt(
      splus*splus*resp.plus*resp.plus + scross*scross*resp.cross*resp.cross );

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/* <lalVerbatim file="SimInspiralUtilsCP"> */
void
LALPopulateSimInspiralSiteInfo(
    LALStatus                  *status,
    SimInspiralTable           *output
    )
/* </lalVerbatim> */
{
  SkyPosition           skyPos;
  LALDetector           detector;
  REAL4                *eff_dist;
  LIGOTimeGPS          *end_time;


  INITSTATUS( status, "LALPopulateSimInspiralSiteInfo", SIMINSPIRALUTILSC );
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
  LALInspiralSiteTimeAndDist(status->statusPtr, output, &detector, end_time,
      eff_dist, &skyPos);

  /* LIGO Livingston observatory*/
  detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  end_time = &(output->l_end_time);
  eff_dist = &(output->eff_dist_l);
  LALInspiralSiteTimeAndDist(status->statusPtr, output, &detector, end_time,
      eff_dist, &skyPos);

  /* GEO observatory*/
  detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  end_time = &(output->g_end_time);
  eff_dist = &(output->eff_dist_g);
  LALInspiralSiteTimeAndDist(status->statusPtr, output, &detector, end_time,
      eff_dist, &skyPos);

  /* TAMA observatory*/
  detector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  end_time = &(output->t_end_time);
  eff_dist = &(output->eff_dist_t);
  LALInspiralSiteTimeAndDist(status->statusPtr, output, &detector, end_time,
      eff_dist, &skyPos);

  /* Virgo observatory*/
  detector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  end_time = &(output->v_end_time);
  eff_dist = &(output->eff_dist_v);
  LALInspiralSiteTimeAndDist(status->statusPtr, output, &detector, end_time,
      eff_dist, &skyPos);

  /*
   *
   * normal exit
   *
   */


  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/* <lalVerbatim file="SimInspiralUtilsCP"> */
void
XLALSortSimInspiral(
    SimInspiralTable **head,
    int (*comparefunc)(const SimInspiralTable * const *,
      const SimInspiralTable * const *)
    )
/* </lalVerbatim> */
{
  INT4 i;
  INT4 length;
  SimInspiralTable *event;
  SimInspiralTable **array;

  /* empty list --> no-op */
  if(!head || !*head)
    return;

  /* count the number of events in the list */
  for(length = 0, event = *head; event; event = event->next)
    length++;

  /* construct an array of pointers into the list */
  array = LALCalloc(length, sizeof(*array));
  for(i = 0, event = *head; event; event = event->next)
    array[i++] = event;

  /* sort the array using the specified function */
  qsort(array, length, sizeof(*array),
      (int(*)(const void *, const void *)) comparefunc);

  /* re-link the list according to the sorted array */
  for(i = 0; i < length; i++, head = &(*head)->next)
    *head = array[i];
  *head = NULL;

  /* free the array */
  LALFree(array);
}

/* <lalVerbatim file="SimInspiralUtilsCP"> */
int
XLALCompareSimInspiralByGeocentEndTime(
    const SimInspiralTable * const *a,
    const SimInspiralTable * const *b
    )
/* </lalVerbatim> */
{
  INT8 ta, tb;
  INT8 epsilon = 10;	/* nanoseconds */

  ta = geocent_end_time(*a);
  tb = geocent_end_time(*b);

  if(ta > tb + epsilon)
    return(1);
  if(ta < tb - epsilon)
    return(-1);
  return(0);
}


/* <lalVerbatim file="SimInspiralUtilsCP"> */
int
XLALFreeSimInspiral (
    SimInspiralTable **eventHead
    )
/* </lalVerbatim> */
{
  EventIDColumn        *eventId;

  while ( (*eventHead)->event_id )
  {
    /* free any associated event_id's */
    eventId = (*eventHead)->event_id;
    (*eventHead)->event_id = (*eventHead)->event_id->next;
    LALFree( eventId );
  }
  LALFree( *eventHead );

  return (0);
}


/* <lalVerbatim file="SimInspiralUtilsCP"> */
INT8
XLALReturnSimInspiralEndTime (
    SimInspiralTable *event,
    CHAR             *ifo
    )
/* </lalVerbatim> */
{
  static const char *func = "ReturnSimInspiralEndTime";
  if ( ! strcmp( "L1", ifo ) )
  {
    return( XLALGPSToINT8NS(&(event->l_end_time) ) );
  }
  else if ( ! strcmp( "H1", ifo ) ||
      ! strcmp( "H2", ifo ) )
  {
    return( XLALGPSToINT8NS(&(event->h_end_time) ) );
  }
  else if ( ! strcmp( "G1", ifo ) )
  {
    return(XLALGPSToINT8NS(&(event->g_end_time) ) );
  }
  else if ( ! strcmp( "T1", ifo ) )
  {
    return( XLALGPSToINT8NS(&(event->t_end_time) ) );
  }
  else if ( ! strcmp( "V1", ifo ) )
  {
    return( XLALGPSToINT8NS(&(event->v_end_time) ) );
  }
  else
  {
    XLAL_ERROR(func,XLAL_EIO);
  }

}

/* <lalVerbatim file="SimInspiralUtilsCP"> */
int
XLALSnglSimInspiralTest (
    SimInspiralTable  **simHead,
    SnglInspiralTable **eventHead,
    SimInspiralTable  **missedSimHead,
    SnglInspiralTable **missedSnglHead,
    INT8                injectWindowNS
    )
/* </lalVerbatim> */
{

  /* Note: we are assuming that both the inspiral and */
  /* injection events are time sorted                 */
  SimInspiralTable *thisSimEvent = *simHead;
  SimInspiralTable *thisMissedSim= NULL;
  SimInspiralTable *prevSimEvent = NULL;
  SnglInspiralTable *thisEvent   = *eventHead;
  SnglInspiralTable *prevEvent   = NULL;
  SnglInspiralTable *thisMissed  = NULL;
  EventIDColumn     *thisId      = NULL;

  int numSimFound  = 0;
  int coincidence = 0;

  INT8 simGeocentTime, simSiteTime, inspiralTime;
  INT8 earthRadiusNS = (INT8) ( 1e9 * 2 * LAL_REARTH_SI / LAL_C_SI );

  *simHead     = NULL;
  *eventHead   = NULL;


  if ( ! thisEvent )
  {
    XLALPrintInfo( "No triggers in input data, all injections missed\n" );

    *missedSimHead = thisSimEvent;
    return(0);
  }
  else
  {

    /* begin loop over the sim_inspiral events */
    while ( thisSimEvent )
    {
      coincidence = 0;
      /* find the end time of the SimEvent */
      simGeocentTime = geocent_end_time( thisSimEvent );

      /* find the first inspiral event after the current sim event */
      while ( thisEvent )
      {
        /* compute the time in nanosec for thisEvent */
        inspiralTime = XLALGPSToINT8NS( &(thisEvent->end_time) );

        if( inspiralTime < (simGeocentTime - earthRadiusNS - injectWindowNS ) )
        {
          /* discard this event and move on to the next one */
          if ( ! *missedSnglHead )
          {
            *missedSnglHead = thisMissed = thisEvent;
          }
          else
          {
            thisMissed = thisMissed->next = thisEvent;
          }
          if ( prevEvent ) prevEvent->next = thisEvent->next;
          thisEvent = thisEvent->next;
          thisMissed->next = NULL;
          XLALPrintInfo( "-" );
        }
        else
        {
          /* we have reached the negative coincincidence window */
          break;
        }
      }

      while ( thisEvent )
      {
        /* compute the time in nanosec for thisEvent */
        inspiralTime = XLALGPSToINT8NS( &(thisEvent->end_time) );

        if( inspiralTime < (simGeocentTime + earthRadiusNS + injectWindowNS ) )
        {
          /* this event may be in coincidence window, need to check site
           * end time */
          simSiteTime = XLALReturnSimInspiralEndTime( thisSimEvent,
              thisEvent->ifo );


          if ( (inspiralTime > (simSiteTime - injectWindowNS)) &&
              (inspiralTime < (simSiteTime + injectWindowNS)) )
          {
            /* this event is within the coincidence window  */

            /* store the sim inspiral in the event_id's for this sngl */
            thisId = thisEvent->event_id;
            while ( thisId )
            {
              thisId->simInspiralTable = thisSimEvent;
              thisId = thisId->next;
            }

            /* store this event and move on to the next one */
            if ( ! *eventHead ) *eventHead = thisEvent;
            prevEvent = thisEvent;
            thisEvent = thisEvent->next;
            coincidence = 1;
            XLALPrintInfo( "+" );
          }
          else
          {
            /* discard this event and move on to the next one */
            if ( ! *missedSnglHead )
            {
              *missedSnglHead = thisMissed = thisEvent;
            }
            else
            {
              thisMissed = thisMissed->next = thisEvent;
            }

            if ( prevEvent ) prevEvent->next = thisEvent->next;
            thisEvent = thisEvent->next;
            thisMissed->next = NULL;
            XLALPrintInfo( "-" );
          }
        }
        else
        {
          /* we have reached the end of the positive coincincidence window */
          break;
        }
      }

      if ( coincidence )
      {
        /* keep this sim event in the list and move to the next sim event */
        if ( ! *simHead ) *simHead = thisSimEvent;
        prevSimEvent = thisSimEvent;
        ++numSimFound;
        thisSimEvent = thisSimEvent->next;
        XLALPrintInfo( "F" );
      }
      else
      {
        /* save this sim event in the list of missed events... */
        if ( ! *missedSimHead )
        {
          *missedSimHead = thisMissedSim = thisSimEvent;
        }
        else
        {
          thisMissedSim = thisMissedSim->next = thisSimEvent;
        }

        /* ...and remove it from the list of found events */
        if ( prevSimEvent ) prevSimEvent->next = thisSimEvent->next;
        XLALPrintInfo( "M" );

        /* move to the next sim in the list */
        thisSimEvent = thisSimEvent->next;

        /* make sure the missed sim list is terminated */
        thisMissedSim->next = NULL;
      }

      if ( ! thisEvent )
      {
        /* these are no more events to process so all the rest of the */
        /* injections must be put in the missed injections list       */
        if ( ! *missedSimHead )
        {
          /* this and any subsequent events are in the missed sim list */
          if ( thisSimEvent ) thisMissedSim = *missedSimHead = thisSimEvent;
        }
        else
        {
          if ( thisSimEvent )
          {
            /* append the rest of the list to the list of missed injections */
            thisMissedSim = thisMissedSim->next = thisSimEvent;
          }
          else
          {
            /* there are no injections after this one */
            thisMissedSim = thisMissedSim->next = NULL;
          }
        }

        /* terminate the list of found injections correctly */
        if ( prevSimEvent ) prevSimEvent->next = NULL;

        while ( thisMissedSim )
        {
          /* count the number of injections just stuck in the missed list */
          XLALPrintInfo( "M" );
          thisMissedSim = thisMissedSim->next;
        }
        thisSimEvent = NULL;
        break;
      }
    }

    if ( thisEvent )
    {
      while( thisEvent )
      {
        /* discard this event and move on to the next one */
        if ( ! *missedSnglHead )
        {
          *missedSnglHead = thisMissed = thisEvent;
        }
        else
        {
          thisMissed = thisMissed->next = thisEvent;
        }
        if ( prevEvent ) prevEvent->next = thisEvent->next;
        thisEvent = thisEvent->next;
        thisMissed->next = NULL;
        XLALPrintInfo( "-" );
      }
    }
  }
  XLALPrintInfo( "\n" );
  return( numSimFound );
}


/* <lalVerbatim file="SimInspiralUtilsCP"> */
int
XLALCoincSimInspiralTest (
    SimInspiralTable   **simHead,
    CoincInspiralTable **coincHead,
    SimInspiralTable   **missedSimHead,
    CoincInspiralTable **missedCoincHead
    )
/* </lalVerbatim> */
{
  CoincInspiralTable    *thisCoinc       = *coincHead;
  CoincInspiralTable    *prevCoinc       = NULL;
  CoincInspiralTable    *thisMissedCoinc = NULL;
  SimInspiralTable      *thisSim         = NULL;
  SimInspiralTable      *prevSim         = NULL;
  SimInspiralTable      *thisMissedSim   = NULL;
  SnglInspiralTable     *thisSngl        = NULL;
  EventIDColumn         *thisId          = NULL;

  InterferometerNumber   ifoInCoinc = LAL_UNKNOWN_IFO;
  int                    numSimFound = 0;

  if ( !*coincHead )
  {
    XLALPrintInfo(
        "XLALCoincSimInspiral: Empty coincInspiral passed as input" );
    *missedSimHead = *simHead;
    *simHead = NULL;
    return( 0 );
  }

  *coincHead = NULL;

  while( thisCoinc )
  {
    thisSim = NULL;
    /* loop over the interferometers to get the event_id*/

    for ( ifoInCoinc = 0; ifoInCoinc < LAL_NUM_IFO; ifoInCoinc++)
    {
      if ( (thisSngl = thisCoinc->snglInspiral[ifoInCoinc]) )
      {
        thisSim = thisSngl->event_id->simInspiralTable;
        break;
      }
    }

    for ( ; ifoInCoinc < LAL_NUM_IFO; ifoInCoinc++)
    {
      if ( (thisSngl = thisCoinc->snglInspiral[ifoInCoinc]) &&
          (thisSim != thisSngl->event_id->simInspiralTable) )
      {
        thisSim = NULL;
        break;
      }
    }

    if ( thisSim )
    {
      /* thisCoinc is coincident with a thisSim */
      thisCoinc->simInspiral = thisSim;

      /* set the event_id's */
      if ( !thisSim->event_id )
      {
        thisId = thisSim->event_id = LALCalloc( 1, sizeof(EventIDColumn) );
      }
      else
      {
        for ( thisId = thisSim->event_id; thisId->next; thisId = thisId->next);
        thisId = thisId->next = LALCalloc( 1, sizeof(EventIDColumn) );
      }
      thisId->simInspiralTable = thisSim;
      thisId->coincInspiralTable = thisCoinc;

      if ( ! *coincHead )
      {
        *coincHead = thisCoinc;
      }

      XLALPrintInfo( "+" );
      /* move on to next coinc */
      prevCoinc = thisCoinc;
      thisCoinc = thisCoinc->next;
    }
    else
    {
      /* discard this event and move on to the next one */
      if ( ! *missedCoincHead )
      {
        *missedCoincHead = thisMissedCoinc = thisCoinc;
      }
      else
      {
        thisMissedCoinc = thisMissedCoinc->next = thisCoinc;
      }

      if ( prevCoinc ) prevCoinc->next = thisCoinc->next;
      thisCoinc = thisCoinc->next;
      XLALPrintInfo( "-" );

      /* terminate the missed list */
      thisMissedCoinc->next = NULL;
    }
  }

  /* run through simInspirals, keeping only those in coincs */

  thisSim = *simHead;
  *simHead = NULL;

  while( thisSim )
  {
    if( thisSim->event_id )
    {
      /* keep this event in the list and move to the next sim event */
      if ( ! *simHead ) *simHead = thisSim;
      prevSim = thisSim;
      ++numSimFound;
      thisSim = thisSim->next;
      XLALPrintInfo( "F" );
    }
    else
    {
      /* save this sim event in the list of missed events... */
      if ( ! *missedSimHead )
      {
        *missedSimHead = thisMissedSim = thisSim;
      }
      else
      {
        thisMissedSim = thisMissedSim->next = thisSim;
      }

      /* ...and remove it from the list of found events */
      if ( prevSim ) prevSim->next = thisSim->next;
      XLALPrintInfo( "M" );

      /* move to the next sim in the list */
      thisSim = thisSim->next;

      /* make sure the missed sim list is terminated */
      thisMissedSim->next = NULL;
    }
  }
  XLALPrintInfo( "\n" );
  return( numSimFound );
}

