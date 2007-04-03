/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpClusterEvents.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpClusterEventsCV">
Author: Brown D. A., Woods D.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpClusterEvents.c}}
\label{ss:FindChirpClusterEvents.c}

</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>

#define rint(x) (floor((x)+0.5))

NRCSID (FINDCHIRPCLUSTEREVENTSC, "$Id$"); 

/* <lalVerbatim file="FindChirpClusterEventsCP"> */
void
LALFindChirpClusterEvents (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,
    FindChirpBankVetoData      *bankVetoData,
    UINT4                       subBankIndex
    )
/* </lalVerbatim> */
{
  UINT4                 numPoints = 0;
  UINT4                 ignoreIndex = 0;
  UINT4                 eventStartIdx = 0;
  UINT4  		deltaEventIndex = 0;
  UINT4                 j, kmax;
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
  INITSTATUS( status, "LALFindChirpClusterEvents", FINDCHIRPCLUSTEREVENTSC );
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
  switch ( params->approximant )
  {
    case TaylorT1:
      LALSnprintf( searchName, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "TaylorT1twoPN" );
      break;

    case TaylorT2:
      LALSnprintf( searchName, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "TaylorT2twoPN" );
      break;

    case TaylorT3:
      LALSnprintf( searchName, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "TaylorT3twoPN" );
      break;

    case TaylorF2:
      LALSnprintf( searchName, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "TaylorF2twoPN" );
      break;

    case GeneratePPN:
      LALSnprintf( searchName, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "GeneratePPNtwoPN" );
      break;

    case PadeT1:
      LALSnprintf( searchName, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "PadeT1twoPN" );
      break;

    case EOB:
      LALSnprintf( searchName, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "EOBtwoPN" );
      break;

    case FindChirpSP:
      LALSnprintf( searchName, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "FindChirpSPtwoPN" );
      break;

    default:
      ABORT( status, FINDCHIRPH_EUAPX, FINDCHIRPH_MSGEUAPX );
      break;
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
  deltaT = params->deltaT;
  deltaF = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints );
  kmax = input->fcTmplt->tmplt.fFinal / deltaF < numPoints/2 ?
    input->fcTmplt->tmplt.fFinal / deltaF : numPoints/2;

  /* normalisation */
  norm = input->fcTmplt->norm; 
  /* 4.0 * (deltaT / (REAL4)numPoints) / 
         input->segment->segNorm->data[kmax]; */


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
   if ( params->clusterMethod == tmplt )
   {
     deltaEventIndex = 
       (UINT4) rint( (input->fcTmplt->tmplt.tC / deltaT) + 1.0 );
   }
   else if ( params->clusterMethod == window )
   {
     deltaEventIndex = 
       (UINT4) rint( (params->clusterWindow / deltaT) + 1.0 );
   }


  /* look for an events in the filter output */
  for ( j = ignoreIndex; j < numPoints - ignoreIndex; ++j )
  {
    REAL4 modqsq = q[j].re * q[j].re + q[j].im * q[j].im;

    /* if snrsq exceeds threshold at any point */
    if ( modqsq > modqsqThresh )
    {

      /* if we have don't have a chisq or the chisq drops below the       */
      /* modified chisq threshold, start processing events                */
      if ( ! input->segment->chisqBinVec->length ||
          params->chisqVec->data[j] < 
          (params->chisqThresh * ( 1.0 + modqsq * chisqThreshFac )) )
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
        else if (  ! params->clusterMethod == noClustering  &&
            j <= thisEvent->end_time.gpsSeconds + deltaEventIndex &&
            modqsq > thisEvent->snr )
        {
          /* if this is the same event, update the maximum */
          thisEvent->end_time.gpsSeconds = j;
          thisEvent->snr = modqsq;
        }
        else if ( j > thisEvent->end_time.gpsSeconds + deltaEventIndex ||
              params->clusterMethod == noClustering  )
        {
          /* clean up this event */
          SnglInspiralTable *lastEvent;
          if ( bankVetoData->length > 1 ) 
          {
            bvChisq = XLALComputeBankVeto( bankVetoData, subBankIndex,
                             thisEvent->end_time.gpsSeconds, &bvDOF);
          }

          LALFindChirpStoreEvent(status->statusPtr, input, params,
              thisEvent, q, kmax, norm, eventStartIdx, numChisqBins, 
              searchName );
          CHECKSTATUSPTR( status );

          if ( bankVetoData->length > 1 )
          {
            thisEvent->chisq_dof = bvDOF;
            thisEvent->chisq = bvChisq;
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
      }
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
                             thisEvent->end_time.gpsSeconds, &bvDOF);
    }

    LALFindChirpStoreEvent(status->statusPtr, input, params,
         thisEvent, q, kmax, norm, eventStartIdx, numChisqBins, 
         searchName );

    if ( bankVetoData->length > 1 )
    {
      thisEvent->chisq_dof = bvDOF;
      thisEvent->chisq = bvChisq;
    }
    CHECKSTATUSPTR( status );
  }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
