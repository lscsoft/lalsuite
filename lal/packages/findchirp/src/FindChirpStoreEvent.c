/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpStoreEvent.c
 *
 * Author: Brown, D. A.,  Woods D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpStoreEventCV">
Author: Brown D. A., Woods D.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpStoreEvent.c}}
\label{ss:FindChirpStoreEvent.c}

</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>

NRCSID (FINDCHIRPSTOREEVENTC, "$Id$"); 

/* <lalVerbatim file="FindChirpStoreEventCP"> */
void
LALFindChirpStoreEvent (
    LALStatus                  *status,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,
    SnglInspiralTable          *thisEvent,
    COMPLEX8                   *q,
    UINT4                       kmax,
    REAL4                       norm,
    UINT4                       eventStartIdx,
    UINT4                       numChisqBins,
    CHAR                        searchName[LIGOMETA_SEARCH_MAX]
)
/* </lalVerbatim> */

{
  INT8                       timeNS;
  INT4                       timeIndex;
  REAL4                      deltaT;
  LALMSTUnitsAndAcc          gmstUnits;

  INITSTATUS( status, "LALFindChirpStoreEvent", FINDCHIRPSTOREEVENTC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure the output handle exists, but points to a null pointer */
  /*ASSERT( eventList, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );*/

  /* make sure the event handle exists, but points to a null pointer */
  /*ASSERT( thisEvent, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*thisEvent, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );*/

  /* make sure that the input structure contains some input */
  /*ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );*/

  /* check that the workspace vector contains some input */
  /*ASSERT( q, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );*/


  /*
   *
   * store the event
   *
   */


  /* point local pointer to event and params pointers */
  /* note: we expect the gps seconds to be set before calling this routine */
  timeIndex = thisEvent->end_time.gpsSeconds;
  deltaT = params->deltaT;

  /* set the gmst units and strictness */
  gmstUnits.units = MST_HRS;
  gmstUnits.accuracy = LALLEAPSEC_STRICT;
 
  /* set the event LIGO GPS time of the event */
  timeNS = 1000000000L * 
    (INT8) (input->segment->data->epoch.gpsSeconds);
  timeNS += (INT8) (input->segment->data->epoch.gpsNanoSeconds);
  timeNS += (INT8) (1e9 * timeIndex * deltaT);
  thisEvent->end_time.gpsSeconds = (INT4) (timeNS/1000000000L);
  thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS%1000000000L);
  LALGPStoGMST1( status->statusPtr, &(thisEvent->end_time_gmst),
      &(thisEvent->end_time), &gmstUnits );
  CHECKSTATUSPTR( status );

  /* set the impulse time for the event */
  thisEvent->template_duration = (REAL8) input->fcTmplt->tmplt.tC;

  /* record the ifo name for the event */
  strncpy( thisEvent->ifo, input->segment->data->name, 
      2 * sizeof(CHAR) );
  strncpy( thisEvent->channel, input->segment->data->name + 3,
      (LALNameLength - 3) * sizeof(CHAR) );
  thisEvent->impulse_time = thisEvent->end_time;

  /* record the coalescence phase of the chirp */
  thisEvent->coa_phase = (REAL4)
      atan2( q[timeIndex].im, q[timeIndex].re );

  /* copy the template into the event */
  thisEvent->mass1   = (REAL4) input->fcTmplt->tmplt.mass1;
  thisEvent->mass2   = (REAL4) input->fcTmplt->tmplt.mass2;
  thisEvent->mchirp  = (REAL4) input->fcTmplt->tmplt.chirpMass;
  thisEvent->eta     = (REAL4) input->fcTmplt->tmplt.eta;
  thisEvent->tau0    = (REAL4) input->fcTmplt->tmplt.t0;
  thisEvent->tau2    = (REAL4) input->fcTmplt->tmplt.t2;
  thisEvent->tau3    = (REAL4) input->fcTmplt->tmplt.t3;
  thisEvent->tau4    = (REAL4) input->fcTmplt->tmplt.t4;
  thisEvent->tau5    = (REAL4) input->fcTmplt->tmplt.t5;
  thisEvent->ttotal  = (REAL4) input->fcTmplt->tmplt.tC;
  thisEvent->f_final = (REAL4) input->fcTmplt->tmplt.fFinal;

  /* set the type of the template used in the analysis */
  memcpy( thisEvent->search, searchName, 
      LIGOMETA_SEARCH_MAX * sizeof(CHAR) );

  /* set snrsq, chisq, sigma and effDist for this event */
  if ( input->segment->chisqBinVec->length )
  {
    /* we store chisq distributed with 2p - 2 degrees of freedom */
    /* in the database. params->chisqVec->data = r^2 = chisq / p */
    /* so we multiply r^2 by p here to get chisq                 */
    thisEvent->chisq = 
      params->chisqVec->data[timeIndex] * (REAL4) numChisqBins;
    thisEvent->chisq_dof = numChisqBins;
  }
  else
  {
    thisEvent->chisq     = 0;
    thisEvent->chisq_dof = 0;
  }
  thisEvent->sigmasq = norm * input->segment->segNorm->data[kmax] * 
    input->segment->segNorm->data[kmax] * input->fcTmplt->tmpltNorm;
  thisEvent->eff_distance = 
    (input->fcTmplt->tmpltNorm * input->segment->segNorm->data[kmax] * 
     input->segment->segNorm->data[kmax]) / thisEvent->snr;
  thisEvent->eff_distance = sqrt( thisEvent->eff_distance );
  thisEvent->snr *= norm;
  thisEvent->snr = sqrt( thisEvent->snr );

  /* compute the time since the snr crossing */
  thisEvent->event_duration = (REAL8) timeIndex - (REAL8) eventStartIdx;
  thisEvent->event_duration *= (REAL8) deltaT;
}
