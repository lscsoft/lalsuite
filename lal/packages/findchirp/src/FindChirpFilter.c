/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpFilter.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpFilterCV">
Author: Brown D. A.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpFilter.c}}
\label{ss:FindChirpFilter.c}

This module provides the core of the matched filter for binary inspiral
chirps.


\subsection{Matched Filtering Using Post-Newtonian Templates}

The gravitational wave strain induced in an interferometer by a binary 
inspiral may be written as
\begin{equation}
h(t) = \frac{A(t)}{\mathcal{D}} \cos\left( 2 \phi(t) - \theta \right),
\label{eq:rootwaveform}
\end{equation}
where
\begin{equation}
A(t) = - \frac{2G\mu}{c^4} \left[ \pi GM f(t) \right]^\frac{2}{3}
\end{equation}
and $\mathcal{D}$ is the \emph{effective distance}, given by
\begin{equation}
\mathcal{D} = \frac{r}{\sqrt{F_+^2 (1 + \cos^2 \iota)^2 + F_\times^2 4 \cos^2 \iota}}.
\end{equation}
The phase angle $\theta$ is
\begin{equation}
\tan \theta = \frac{F_\times 2\cos \iota}{F_+(1 + \cos^2 \iota)}
\end{equation}
and $\phi(t)$ is the phase evolution of the inspiral waveform.



\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpFilterCP}
\idx{LALFindChirpFilterSegment()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpFilterCV}}
</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>
#include <lal/FindChirpFilterOutputVeto.h>

double rint(double x);

NRCSID (FINDCHIRPFILTERC, "$Id$");


/* <lalVerbatim file="FindChirpFilterCP"> */
void
LALFindChirpFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )
/* </lalVerbatim> */
{
  UINT4                 j, k, kmax;
  UINT4                 numPoints;
  UINT4                 deltaEventIndex;
  UINT4                 ignoreIndex;
  REAL4                 deltaT;
  REAL8                 deltaF;
  REAL4                 norm;
  REAL4                 modqsqThresh;
  REAL4                 mismatch;
  REAL4                 chisqThreshFac;
  REAL4                 modChisqThresh;
  UINT4                 numChisqBins;
  UINT4                 eventStartIdx = 0;
  BOOLEAN               haveChisq     = 0;
  COMPLEX8             *qtilde        = NULL;
  COMPLEX8             *q             = NULL;
  COMPLEX8             *inputData     = NULL;
  COMPLEX8             *tmpltSignal   = NULL;
  SnglInspiralTable    *thisEvent     = NULL;
  LALMSTUnitsAndAcc     gmstUnits;
  CHAR                  searchName[LIGOMETA_SEARCH_MAX];

  INITSTATUS( status, "LALFindChirpFilter", FINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the filter parameters are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPH_EDTZO, FINDCHIRPH_MSGEDTZO );
  ASSERT( params->rhosqThresh >= 0, status,
      FINDCHIRPH_ERHOT, FINDCHIRPH_MSGERHOT );
  ASSERT( params->chisqThresh >= 0, status,
      FINDCHIRPH_ECHIT, FINDCHIRPH_MSGECHIT );

  /* check that the fft plan exists */
  ASSERT( params->invPlan, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the workspace vectors exist */
  ASSERT( params->qVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->qVec->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->qtildeVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->qtildeVec->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the chisq parameter and input structures exist */
  ASSERT( params->chisqParams, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInput, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec ) 
  {
    ASSERT( params->rhosqVec->data->data, status, 
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->rhosqVec->data, status, 
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  if ( params->cVec )
  {
    ASSERT( params->cVec->data->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL  );
    ASSERT( params->cVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* if a chisqVec vector has been created, check we can store data in it */
  if ( params->chisqVec ) 
  {
    ASSERT( params->chisqVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* make sure that the input structure exists */
  ASSERT( input, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the input structure contains some input */
  ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

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
   * point local pointers to input and output pointers
   *
   */


  /* workspace vectors */
  q = params->qVec->data;
  qtilde = params->qtildeVec->data;

  /* number of points in a segment */
  numPoints = params->qVec->length;

  /* template and data */
  inputData = input->segment->data->data->data;
  tmpltSignal = input->fcTmplt->data->data;
  deltaT = params->deltaT;
  deltaF = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints );
  kmax = input->fcTmplt->tmplt.fFinal / deltaF < numPoints/2 ? 
    input->fcTmplt->tmplt.fFinal / deltaF : numPoints/2;

  /* the length of the chisq bin vec is the number of bin boundaries so the */
  /* number of chisq bins is (length - 1) or 0 if there are no boundaries   */
  numChisqBins = input->segment->chisqBinVec->length ? 
    input->segment->chisqBinVec->length - 1 : 0;

  /* set the gmst units and strictness */
  gmstUnits.units = MST_HRS;
  gmstUnits.accuracy = LALLEAPSEC_STRICT;


  /*
   *
   * compute viable search regions in the snrsq vector
   *
   */


  if ( input->fcTmplt->tmplt.tC <= 0 )
  {
    ABORT( status, FINDCHIRPH_ECHTZ, FINDCHIRPH_MSGECHTZ );
  }

  deltaEventIndex = (UINT4) rint( (input->fcTmplt->tmplt.tC / deltaT) + 1.0 );

  /* ignore corrupted data at start and end */
  ignoreIndex = ( input->segment->invSpecTrunc / 2 ) + deltaEventIndex;

  if ( lalDebugLevel & LALINFO )
  {
    CHAR infomsg[256];

    LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
        "m1 = %e m2 = %e => %e seconds => %d points\n"
        "invSpecTrunc = %d => ignoreIndex = %d\n", 
        input->fcTmplt->tmplt.mass1, input->fcTmplt->tmplt.mass2, 
        input->fcTmplt->tmplt.tC , deltaEventIndex, 
        input->segment->invSpecTrunc, ignoreIndex );
    LALInfo( status, infomsg );
  }

  /* XXX check that we are not filtering corrupted data XXX */
  /* XXX this is hardwired to 1/4 segment length        XXX */
  if ( ignoreIndex > numPoints / 4 )
  {
    ABORT( status, FINDCHIRPH_ECRUP, FINDCHIRPH_MSGECRUP );
  }
  /* XXX reset ignoreIndex to one quarter of a segment XXX */
  ignoreIndex = numPoints / 4;

  if ( lalDebugLevel & LALINFO )
  {
    CHAR infomsg[256];

    LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg), 
        "filtering from %d to %d\n",
        ignoreIndex, numPoints - ignoreIndex );
    LALInfo( status, infomsg );
  }

  /*
   *
   * compute qtilde and q
   *
   */


  memset( qtilde, 0, numPoints * sizeof(COMPLEX8) );

  /* qtilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
  {
    REAL4 r = inputData[k].re;
    REAL4 s = inputData[k].im;
    REAL4 x = tmpltSignal[k].re;
    REAL4 y = 0 - tmpltSignal[k].im;       /* note complex conjugate */

    qtilde[k].re = r*x - s*y;
    qtilde[k].im = r*y + s*x;
  }

  /* inverse fft to get q */
  LALCOMPLEX8VectorFFT( status->statusPtr, params->qVec, params->qtildeVec, 
      params->invPlan );
  CHECKSTATUSPTR( status );


  /* 
   *
   * calculate signal to noise squared 
   *
   */


  /* if full snrsq vector is required, set it to zero */
  if ( params->rhosqVec )
    memset( params->rhosqVec->data->data, 0, numPoints * sizeof( REAL4 ) );

  if (params->cVec )
    memset( params->cVec->data->data, 0, numPoints * sizeof( COMPLEX8 ) ); 

  /* normalisation */
  params->norm = norm = 
    4.0 * (deltaT / (REAL4)numPoints) / input->segment->segNorm->data[kmax];

  /* normalised snr threhold */
  modqsqThresh = params->rhosqThresh / norm;

  /* we threshold on the "modified" chisq threshold computed from       */
  /*   chisqThreshFac = delta^2 * norm / p                              */
  /*   rho^2 = norm * modqsq                                            */
  /*                                                                    */
  /* So we actually threshold on                                        */
  /*                                                                    */
  /*    r^2 < chisqThresh * ( 1 + modqsq * chisqThreshFac )             */
  /*                                                                    */
  /* which is the same as thresholding on                               */
  /*    r^2 < chisqThresh * ( 1 + rho^2 * delta^2 / p )                 */
  /* and since                                                          */
  /*    chisq = p r^2                                                   */
  /* this is equivalent to thresholding on                              */
  /*    chisq < chisqThresh * ( p + rho^2 delta^2 )                     */
  /*                                                                    */
  /* The raw chisq is stored in the database. this quantity is chisq    */
  /* distributed with 2p-2 degrees of freedom.                          */
  mismatch = 1.0 - input->fcTmplt->tmplt.minMatch;
  chisqThreshFac = norm * mismatch * mismatch / (REAL4) numChisqBins;

  /* if full snrsq vector is required, store the snrsq */
  if ( params->rhosqVec ) 
  {
    memcpy( params->rhosqVec->name, input->segment->data->name,
        LALNameLength * sizeof(CHAR) );
    memcpy( &(params->rhosqVec->epoch), &(input->segment->data->epoch), 
        sizeof(LIGOTimeGPS) );
    params->rhosqVec->deltaT = input->segment->deltaT;

    for ( j = 0; j < numPoints; ++j )
    {
      REAL4 modqsq = q[j].re * q[j].re + q[j].im * q[j].im;
      params->rhosqVec->data->data[j] = norm * modqsq;
    }
  }

 if ( params->cVec ) 
 {
   memcpy( params->cVec->name, input->segment->data->name,
       LALNameLength * sizeof(CHAR) );
   memcpy( &(params->cVec->epoch), &(input->segment->data->epoch), 
       sizeof(LIGOTimeGPS) );
   params->cVec->deltaT = input->segment->deltaT;

   for ( j = 0; j < numPoints; ++j )
   {
     params->cVec->data->data[j].re = sqrt(norm) * q[j].re;
     params->cVec->data->data[j].im = sqrt(norm) * q[j].im;
   }
 }

  /* look for an events in the filter output */
  for ( j = ignoreIndex; j < numPoints - ignoreIndex; ++j )
  {
    REAL4 modqsq = q[j].re * q[j].re + q[j].im * q[j].im;

    /* if snrsq exceeds threshold at any point */
    if ( modqsq > modqsqThresh )
    {
      /* compute chisq vector if it does not exist and we want it */
      if ( ! haveChisq  && input->segment->chisqBinVec->length )
      {
        memset( params->chisqVec->data, 0, 
            params->chisqVec->length * sizeof(REAL4) );

        /* pointers to chisq input */
        params->chisqInput->qtildeVec = params->qtildeVec;
        params->chisqInput->qVec      = params->qVec;

        /* pointer to the chisq bin vector in the segment */
        params->chisqParams->chisqBinVec = input->segment->chisqBinVec;
        params->chisqParams->norm        = norm;

        /* compute the chisq bin boundaries for this template */
        if ( ! params->chisqParams->chisqBinVec->data )
        {
          LALFindChirpComputeChisqBins( status->statusPtr, 
              params->chisqParams->chisqBinVec, input->segment, kmax );
          CHECKSTATUSPTR( status );
        }

        /* compute the chisq threshold: this is slow! */
        LALFindChirpChisqVeto( status->statusPtr, params->chisqVec, 
            params->chisqInput, params->chisqParams );
        CHECKSTATUSPTR (status); 

        haveChisq = 1;
      }

      /* if we have don't have a chisq or the chisq drops below the       */
      /* modified chisq threshold, start processing events                */
      if ( ! input->segment->chisqBinVec->length ||
          params->chisqVec->data[j] < 
          (params->chisqThresh * ( 1.0 + modqsq * chisqThreshFac )) )
      {
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
        else if ( params->maximiseOverChirp &&
            j <= thisEvent->end_time.gpsSeconds + deltaEventIndex &&
            modqsq > thisEvent->snr )
        {
          /* if this is the same event, update the maximum */
          thisEvent->end_time.gpsSeconds = j;
          thisEvent->snr = modqsq;
        }
        else if ( j > thisEvent->end_time.gpsSeconds + deltaEventIndex ||
            ! params->maximiseOverChirp )
        {
          /* clean up this event */
          SnglInspiralTable *lastEvent;
          INT8               timeNS;
          INT4               timeIndex = thisEvent->end_time.gpsSeconds;

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

          /* set the impuse time for the event */
          thisEvent->template_duration = (REAL8) input->fcTmplt->tmplt.tC;

          /* record the ifo and channel name for the event */
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
   * clean up the last event if there is one
   *
   */


  if ( thisEvent )
  {
    INT8           timeNS;
    INT4           timeIndex = thisEvent->end_time.gpsSeconds;

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

    /* set the impuse time for the event */
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


  /*
   *
   * check the events pass the filter output veto
   *
   */

  
  if ( params->filterOutputVetoParams )
  {
    LALFindChirpFilterOutputVeto( status->statusPtr, eventList, params->qVec, 
        norm, params->filterOutputVetoParams );
    CHECKSTATUSPTR( status );
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
