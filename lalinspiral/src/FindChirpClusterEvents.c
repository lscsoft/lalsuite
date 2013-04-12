/*
*  Copyright (C) 2007 Stas Babak, Chad Hanna, Darren Woods, Duncan Brown, Patrick Brady
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
 * File Name: FindChirpClusterEvents.c
 *
 * Author: Brown, D. A.
 *
 *
 *-----------------------------------------------------------------------
 */

/**

\author Brown D. A., Woods D.
\file

*/

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>

void
LALFindChirpClusterEvents (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,
    FindChirpBankVetoData      *bankVetoData,
    UINT4                       subBankIndex,
    int                         writeCData,
    InspiralTemplate           *bankCurrent
    )

{

  int                   xlalRetCode = 0;
  INT8                  bvTimeNS = 0;
  UINT4                 numPoints = 0;
  UINT4                 ignoreIndex = 0;
  UINT4                 eventStartIdx = 0;
  UINT4  		deltaEventIndex = 0;
  UINT4                 lowerIndex = 0;
  UINT4                 upperIndex = 0;
  UINT4                 buffer = 0;
  UINT4                 bvTimeIndex = 0;
  UINT4                 j, kmax;
  UINT4 		doChisqFlag = 1;
  REAL4                 norm = 0;
  REAL4                 deltaT;
  REAL8                 deltaF;
  REAL4                 modqsqThresh = 0;
  REAL4                 chisqThreshFac = 0;
  UINT4                 numChisqBins = 0;
  COMPLEX8             *q = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  CHAR                  searchName[LIGOMETA_SEARCH_MAX];
  UINT4			bvDOF = 0;
  REAL4			bvChisq = 0;
  UINT4                 ccDOF = 0;
  REAL4                 ccChisq = 0;
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* check the allowed approximants */
  xlalRetCode = XLALInspiralGetApproximantString( searchName, LIGOMETA_SEARCH_MAX,
                 params->approximant, params->order );

  if ( xlalRetCode == XLAL_FAILURE )
  {
    ABORTXLAL( status );
  }

  /* make sure the approximant in the tmplt and segment agree */
  if ( params->approximant != input->fcTmplt->tmplt.approximant ||
      params->approximant != input->segment->approximant )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }


  /*
   *
   * set up the variables needed to cluster
   *
   */

  q = params->qVec->data;
  numPoints = params->qVec->length;
  ignoreIndex = params->ignoreIndex;
  lowerIndex = ignoreIndex;
  upperIndex = numPoints - ignoreIndex;
  deltaT = params->deltaT;
  deltaF = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints );
  kmax = input->fcTmplt->tmplt.fFinal / deltaF < numPoints/2 ?
    input->fcTmplt->tmplt.fFinal / deltaF : numPoints/2;

  /* normalisation */
  norm = input->fcTmplt->norm;

  /* normalised snr threhold */
  modqsqThresh = params->rhosqThresh / norm;

  /* the length of the chisq bin vec is the number of bin boundaries so the */
  /* number of chisq bins is (length - 1) or 0 if there are no boundaries   */
  numChisqBins = input->segment->chisqBinVec->length ?
    input->segment->chisqBinVec->length - 1 : 0;


  /* we threshold on the "modified" chisq threshold computed from       */
  /*   chisqThreshFac = chisqDelta * norm / p                           */
  /*                                                                    */
  /*   rho^2 = norm * modqsq                                            */
  /*                                                                    */
  /* So we actually threshold on                                        */
  /*                                                                    */
  /*    r^2 < chisqThresh * ( 1 + modqsq * chisqThreshFac )             */
  /*                                                                    */
  /* which is the same as thresholding on                               */
  /*    r^2 < chisqThresh * ( 1 + rho^2 * chisqDelta / p )              */
  /* and since                                                          */
  /*    chisq = p r^2                                                   */
  /* this is equivalent to thresholding on                              */
  /*    chisq < chisqThresh * ( p + rho^2 chisqDelta )                  */
  /*                                                                    */
  /* The raw chisq is stored in the database. this quantity is chisq    */
  /* distributed with 2p-2 degrees of freedom.                          */
  chisqThreshFac = norm * params->chisqDelta / (REAL4) numChisqBins;


  /*
   *
   * Apply one of the clustering algorithms
   *
   */


   /* set deltaEventIndex depending on clustering method used */
   if ( params->clusterMethod == FindChirpClustering_tmplt )
   {
     deltaEventIndex =
       (UINT4) rint( (input->fcTmplt->tmplt.tC / deltaT) + 1.0 );
   }
   else if ( params->clusterMethod == FindChirpClustering_window )
   {
     deltaEventIndex =
       (UINT4) rint( (params->clusterWindow / deltaT) + 1.0 );
   }
   else if ( params->clusterMethod == FindChirpClustering_tmpltwindow )
   {
     if ( input->fcTmplt->tmplt.tC > params->clusterWindow )
     {
       deltaEventIndex =
        (UINT4) rint( (input->fcTmplt->tmplt.tC / deltaT) + 1.0 );
     }
     else
     {
       deltaEventIndex =
        (UINT4) rint( (params->clusterWindow / deltaT) + 1.0 );
     }
   }

  /* In coherent stage cluster around trigbank-coherent triggers */
  if (writeCData) {
    /* set the event LIGO GPS time of the bankVetoTrigger */
    bvTimeNS = 1000000000L * (INT8) (bankCurrent->end_time.gpsSeconds);
    bvTimeNS += (INT8) (bankCurrent->end_time.gpsNanoSeconds);
    bvTimeNS -= XLALGPSToINT8NS( &(input->segment->data->epoch) );
    bvTimeIndex = (UINT4) rint( ((REAL8) bvTimeNS)/ (deltaT*1.0e9) );

    buffer = (UINT4) rint( 64.0 / deltaT );
    if ( (bvTimeIndex < buffer ) || (bvTimeIndex > (numPoints-buffer) ) ) {
      upperIndex = 0;
      lowerIndex = 0;
    }
    else {
      lowerIndex = ( ( bvTimeIndex > deltaEventIndex ) ? (bvTimeIndex) : 0 );
      upperIndex = bvTimeIndex + deltaEventIndex;
    }
  }

  /* look for an events in the filter output */
  for ( j = lowerIndex; j < upperIndex; ++j )
  {
    REAL4 modqsq = crealf(q[j]) * crealf(q[j]) + q[j].im * q[j].im;

    /* if snrsq exceeds threshold at any point */
    if ( modqsq > modqsqThresh )
    {
      /* If it crosses the threshold see if we need to do a chisq test
         since this is no longer computed in FindChirpFilterSegment */
      if ( input->segment->chisqBinVec->length && doChisqFlag)
      {
        /* compute the chisq vector for this segment */
        memset( params->chisqVec->data, 0,
          params->chisqVec->length * sizeof(REAL4) );

        /* pointers to chisq input */
        params->chisqInput->qtildeVec = params->qtildeVec;
        params->chisqInput->qVec      = params->qVec;

        /* pointer to the chisq bin vector in the segment */
        params->chisqParams->chisqBinVec = input->segment->chisqBinVec;
        params->chisqParams->norm        = norm;

        /* compute the chisq bin boundaries for this template */
        if (params->chisqParams->chisqBinVec->data)
        {
          LALFree(params->chisqParams->chisqBinVec->data);
          params->chisqParams->chisqBinVec->data = NULL;
        }

        LALFindChirpComputeChisqBins( status->statusPtr,
            params->chisqParams->chisqBinVec, input->segment, kmax );
            CHECKSTATUSPTR( status );

        /* compute the chisq threshold: this is slow! */
        LALFindChirpChisqVeto( status->statusPtr, params->chisqVec,
          params->chisqInput, params->chisqParams );
        CHECKSTATUSPTR (status);
        doChisqFlag = 0;
      }

      /* if we have don't have a chisq or the chisq drops below the       */
      /* modified chisq threshold, start processing events                */
      if ( ! input->segment->chisqBinVec->length ||
          params->chisqVec->data[j] <
          (params->chisqThresh * ( 1.0 + modqsq * chisqThreshFac )) )
      {
        if (1) /* eventually check bank veto ! */
        {

          /*
           *
           * find the local maximum of the event
           *
           */


          if ( ! *eventList )
          {
            /* store the start of the crossing */
            eventStartIdx = j;

            /* if this is the first event, start the list */
            thisEvent = *eventList = (SnglInspiralTable *)
              LALCalloc( 1, sizeof(SnglInspiralTable) );
            if ( ! thisEvent )
            {
              ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
            }

            /* record the data that we need for the clustering algorithm */
            thisEvent->end_time.gpsSeconds = j;
            thisEvent->snr = modqsq;
          }
          else if (  ! params->clusterMethod == FindChirpClustering_none  &&
              j <= thisEvent->end_time.gpsSeconds + deltaEventIndex &&
             modqsq > thisEvent->snr )
          {
            /* if this is the same event, update the maximum */
            thisEvent->end_time.gpsSeconds = j;
            thisEvent->snr = modqsq;
          }
          else if ( j > thisEvent->end_time.gpsSeconds + deltaEventIndex ||
                params->clusterMethod == FindChirpClustering_none  )
          {
            /* clean up this event */
            SnglInspiralTable *lastEvent;
            if ( bankVetoData->length > 1 )
            {
              bvChisq = XLALComputeBankVeto( bankVetoData, subBankIndex,
					     thisEvent->end_time.gpsSeconds, deltaT, &bvDOF);
            }

	    if ( !writeCData ) {
	      ccChisq = XLALComputeFullChisq(bankVetoData,input,params,q,
                subBankIndex, thisEvent->end_time.gpsSeconds, &ccDOF, norm);
	    }

            LALFindChirpStoreEvent(status->statusPtr, input, params,
                thisEvent, q, kmax, norm, eventStartIdx, numChisqBins,
                searchName );
            CHECKSTATUSPTR( status );

            /* Set bank_chisq and cont_chisq */
            thisEvent->bank_chisq_dof = bvDOF;
            thisEvent->bank_chisq = bvChisq;
            thisEvent->cont_chisq_dof = ccDOF;
            thisEvent->cont_chisq = ccChisq;

            if ( writeCData ) {
              if ( !thisEvent->event_id )
                thisEvent->event_id = (EventIDColumn *) LALCalloc(1, sizeof(EventIDColumn) );
              thisEvent->event_id->id = bankVetoData->fcInputArray[subBankIndex]->fcTmplt->tmplt.event_id->id;
            }

            /* store the start of the crossing */
            eventStartIdx = j;

            /* allocate memory for the newEvent */
            lastEvent = thisEvent;

            lastEvent->next = thisEvent = (SnglInspiralTable *)
              LALCalloc( 1, sizeof(SnglInspiralTable) );
            if ( ! lastEvent->next )
            {
              ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
            }

            /* stick minimal data into the event */
            thisEvent->end_time.gpsSeconds = j;
            thisEvent->snr = modqsq;
          }
        } /* end if bank veto */
      } /* end if chisq */
    }
  }

  /*
   *
   * clean up last event
   *
   */

  if ( thisEvent )
  {
    if ( bankVetoData->length > 1 )
    {
      bvChisq = XLALComputeBankVeto( bankVetoData, subBankIndex,
				     thisEvent->end_time.gpsSeconds, deltaT, &bvDOF);
    }

    if ( !writeCData ) {
      ccChisq = XLALComputeFullChisq(bankVetoData, input,params,q,
            subBankIndex, thisEvent->end_time.gpsSeconds, &ccDOF, norm);
    }

    LALFindChirpStoreEvent(status->statusPtr, input, params,
         thisEvent, q, kmax, norm, eventStartIdx, numChisqBins,
         searchName );

    thisEvent->bank_chisq_dof = bvDOF;
    thisEvent->bank_chisq = bvChisq;
    thisEvent->cont_chisq_dof = ccDOF;
    thisEvent->cont_chisq = ccChisq;

    if ( writeCData ) {
      if ( !thisEvent->event_id )
	thisEvent->event_id = (EventIDColumn *) LALCalloc(1, sizeof(EventIDColumn) );
      
      thisEvent->event_id->id = bankVetoData->fcInputArray[subBankIndex]->fcTmplt->tmplt.event_id->id;
    }

    CHECKSTATUSPTR( status );
  }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
