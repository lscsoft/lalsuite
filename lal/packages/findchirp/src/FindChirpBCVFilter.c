/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpBCVFilter.c
 *
 * Author: Brown, D. A., BCV-Modifications by Messaritaki E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpBCVFilterCV">
Author: Brown D. A., BCV-Modifications by Messaritaki E.
$Id$
</lalVerbatim>

<lalLaTeX>
\input{FindChirpBCVFilterCDoc}

\vfill{\footnotesize\input{FindChirpBCVFilterCV}}
</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpBCV.h>

double rint(double x);

NRCSID (FINDCHIRPBCVFILTERC, "$Id$");


/* <lalVerbatim file="FindChirpBCVFilterCP"> */
void
LALFindChirpBCVFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )
/* </lalVerbatim> */
{
  UINT4                 j, k;
  UINT4                 numPoints;
  UINT4                 deltaEventIndex;
  UINT4                 ignoreIndex;
  REAL4                 deltaT;
  REAL4                 norm;
  REAL4                 modqsqThresh;
  REAL4                 rhosqThresh;
  REAL4                 mismatch;
  REAL4                 chisqThreshFac;
  REAL4                 modChisqThresh;
  UINT4                 numChisqBins;
  UINT4                 eventStartIdx = 0;
  REAL4                 chirpTime     = 0;
  BOOLEAN               haveChisq     = 0;
  COMPLEX8             *qtilde        = NULL; 
  COMPLEX8             *qtildeBCV     = NULL; 
  COMPLEX8             *q             = NULL; 
  COMPLEX8             *qBCV          = NULL;
  COMPLEX8             *inputData     = NULL;
  COMPLEX8             *inputDataBCV  = NULL;
  COMPLEX8             *tmpltSignal   = NULL;
  SnglInspiralTable    *thisEvent     = NULL;
  LALMSTUnitsAndAcc     gmstUnits;
  REAL4                 a1;
  REAL4                 b1;                  
  REAL4                 b2;                  
  REAL4                 templateNorm;
  REAL4                 modqsq;
  REAL4                 omega;
  REAL4                 Num1;
  REAL4                 Num2;
  REAL4                 Den1;
  REAL4                 Den2;
  REAL4                 InvTan1;
  REAL4                 InvTan2;
  REAL4                 m;

#if 0
     FindChirpChisqInput  *chisqInput;
     FindChirpChisqInput  *chisqInputBCV;
#endif

  INITSTATUS( status, "LALFindChirpBCVFilter", FINDCHIRPBCVFILTERC );
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
  ASSERT(params->qVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVec->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec->data,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL);
  ASSERT(params->qVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVecBCV->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV->data,status, FINDCHIRPH_ENULL, 
      FINDCHIRPH_MSGENULL);

  /* check that the chisq parameter and input structures exist */
  ASSERT( params->chisqParams, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInput,   status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInputBCV,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec )
  {
    ASSERT( params->rhosqVec->data->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->rhosqVec->data, status,
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
  ASSERT( input->tmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure the filter has been initialized for the correct approximant */
  if ( params->approximant != BCV )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }

  /* make sure that the template and the segment are both BCV */
  ASSERT( input->fcTmplt->approximant == BCV, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  ASSERT( input->segment->approximant == BCV, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );


  /*
   *
   * point local pointers to input and output pointers
   *
   */


  /* workspace vectors */
  q    = params->qVec->data;
  qBCV = params->qVecBCV->data;
  qtilde    = params->qtildeVec->data;
  qtildeBCV = params->qtildeVecBCV->data;


  /* template and data */
  inputData     = input->segment->data->data->data;
  inputDataBCV  = input->segment->dataBCV->data->data;
  tmpltSignal   = input->fcTmplt->data->data;
  templateNorm  = input->fcTmplt->tmpltNorm;
  deltaT        = params->deltaT;

  /* the length of the chisq bin vec is the number of bin   */
  /* _boundaries_ so the number of chisq bins is length - 1 */
  numChisqBins = input->segment->chisqBinVec->length - 1;

  /* number of points in a segment */
  numPoints = params->qVec->length;

  /* set the gmst units and strictness */
  gmstUnits.units = MST_HRS;
  gmstUnits.accuracy = LALLEAPSEC_STRICT;


  /*
   *
   * compute viable search regions in the snrsq vector
   *
   */


  {
    /* length of the chirp:                                            */
    /* calculated according to chirplen, using the values of M and eta */
    /* for the BCV tempaltes, as calculated using psi0 and psi3        */ 
    REAL4 psi0 = input->tmplt->psi0;
    REAL4 psi3 = input->tmplt->psi3;
    REAL4 fFinal = input->tmplt->fFinal;

    /* REAL4 eta = input->tmplt->eta; */

    /* m1 and m2 are currently equal to 0 in the BCV template bank */
    REAL4 m1 = input->tmplt->mass1;
    REAL4 m2 = input->tmplt->mass2;
    REAL4 fmin = input->segment->fLow;
    /* total mass and eta, for use in chirpTime calculation */
    REAL4 m =  fabs(psi3) / (16.0 * LAL_MTSUN_SI * LAL_PI * LAL_PI * psi0) ;
    REAL4 eta = 3.0 / (128.0*psi0 * pow( (m*LAL_MTSUN_SI*LAL_PI), (5.0/3.0)) );
    REAL4 c0 = 5*m*LAL_MTSUN_SI/(256*eta);
    REAL4 c2 = 743.0/252.0 + eta*11.0/3.0;
    REAL4 c3 = -32*LAL_PI/3;
    REAL4 c4 = 3058673.0/508032.0 + eta*(5429.0/504.0 + eta*617.0/72.0);
    REAL4 x  = pow(LAL_PI*m*LAL_MTSUN_SI*fmin, 1.0/3.0);
    REAL4 x2 = x*x;
    REAL4 x3 = x*x2;
    REAL4 x4 = x2*x2;
    REAL4 x8 = x4*x4;

    /* k that corresponds to fFinal, currently not used      */
    /* UINT4 kFinal = floor( numPoints * deltaT * fFinal ); */  
    /* BCV normalization parameters */

    chirpTime = fabs(c0*(1 + c2*x2 + c3*x3 + c4*x4)/x8);
    deltaEventIndex = (UINT4) rint( (chirpTime / deltaT) + 1.0 );

    /* ignore corrupted data at start and end */
    ignoreIndex = ( input->segment->invSpecTrunc / 2 ) + deltaEventIndex;

    if ( lalDebugLevel & LALINFO )
    {
      CHAR infomsg[256];

      LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
          "m = %e eta = %e => %e seconds => %d points\n"
          "invSpecTrunc = %d => ignoreIndex = %d\n",
          m, eta, chirpTime, deltaEventIndex,
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
  }

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
   * compute qtilde, qtildeBCV, and q, qBCV 
   *
   */


  memset( qtilde,    0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCV, 0, numPoints * sizeof(COMPLEX8) );

  /* qtilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
  {
    REAL4 r    = inputData[k].re;
    REAL4 s    = inputData[k].im;    
    REAL4 rBCV = inputDataBCV[k].re;
    REAL4 sBCV = inputDataBCV[k].im; 
    REAL4 x = tmpltSignal[k].re;
    REAL4 y = 0.0 - tmpltSignal[k].im; /* note complex conjugate */     

    qtilde[k].re = r * x - s * y ;
    qtilde[k].im = r * y + s * x ;
    qtildeBCV[k].re = rBCV * x - sBCV * y ;
    qtildeBCV[k].im = rBCV * y + sBCV * x ;
  }

  /* qtilde negative frequency only: not DC or nyquist */
  if ( params->computeNegFreq )
  {
    for ( k = numPoints/2 + 2; k < numPoints - 1; ++k )
    {
      REAL4 r = inputData[k].re;
      REAL4 s = inputData[k].im;    
      REAL4 rBCV = inputDataBCV[k].re;
      REAL4 sBCV = inputDataBCV[k].im; 
      REAL4 x = tmpltSignal[k].re;
      REAL4 y = 0.0 - tmpltSignal[k].im; /* note complex conjugate */

      qtilde[k].re = r * x - s * y ;
      qtilde[k].im = r * y + s * x ;
      qtildeBCV[k].re = rBCV * x - sBCV * y ;
      qtildeBCV[k].im = rBCV * y + sBCV * x ;
    }
  }

  /* inverse fft to get q, and qBCV */
  LALCOMPLEX8VectorFFT( status->statusPtr, params->qVec, 
      params->qtildeVec, params->invPlan );
  CHECKSTATUSPTR( status );
  LALCOMPLEX8VectorFFT( status->statusPtr, params->qVecBCV, 
      params->qtildeVecBCV, params->invPlan );
  CHECKSTATUSPTR( status );




  /* 
   *
   * calculate signal to noise squared 
   *
   */



  /* if full snrsq vector is required, set it to zero */
  if ( params->rhosqVec )
    memset( params->rhosqVec->data->data, 0, numPoints * sizeof( REAL4 ) );

  /* normalisation */
  /*  For the BCV templates, templateNorm is equal to                      */
  /*  (5*eta/96) * (M/pi^2)^(2/3) (Tsun/Dt)^(-1/3) (2*Msun(L)*d/1Mpc)^2,   */
  /* the square of the normalization factor that multiplies the template   */

  rhosqThresh = params->rhosqThresh;

  params->norm = norm = deltaT / ((REAL4) numPoints) ;


  /* normalized snr threhold */
  modqsqThresh = rhosqThresh / norm ;

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

  mismatch = 1.0 - input->tmplt->minMatch;
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
      REAL4 modqsqSP  = q[j].re * q[j].re + q[j].im * q[j].im ;
      REAL4 modqsqBCV = qBCV[j].re * qBCV[j].re + qBCV[j].im * qBCV[j].im ;
      REAL4 ImProd = 2.0 * ( - q[j].re * qBCV[j].im + qBCV[j].re * q[j].im ) ;

      REAL4 modqsq = ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
          0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) *
        ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
          0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) ;

      params->rhosqVec->data->data[j] = norm * modqsq;   
    }
  }


  /* look for an events in the filter output */
  for ( j = ignoreIndex; j < numPoints - ignoreIndex; ++j )
  {
    REAL4 modqsqSP  = q[j].re * q[j].re + q[j].im * q[j].im ;
    REAL4 modqsqBCV = qBCV[j].re * qBCV[j].re + qBCV[j].im * qBCV[j].im ;
    REAL4 ImProd = 2.0 * ( - q[j].re * qBCV[j].im + qBCV[j].re * q[j].im ) ;

    REAL4 modqsq = ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
        0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) *
      ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
        0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) ;


    /* if snrsq exceeds threshold at any point */
    if ( modqsq > modqsqThresh )                  
    {
      /* compute chisq vector if it does not exist and we want it */
      if ( ! haveChisq  && input->segment->chisqBinVec->length )
      {
        memset( params->chisqVec->data, 0,
            params->chisqVec->length * sizeof(REAL4) );

        /* pointers to chisq input */
        params->chisqInput->qtildeVec    = params->qtildeVec;
        params->chisqInput->qVec         = params->qVec;
        params->chisqInputBCV->qtildeVec = params->qtildeVecBCV;
        params->chisqInputBCV->qVec      = params->qVecBCV;

        /* pointer to the chisq bin vector in the segment */
        params->chisqParams->chisqBinVec    = input->segment->chisqBinVec;
        params->chisqParams->chisqBinVecBCV = input->segment->chisqBinVecBCV;
        params->chisqParams->norm           = norm;
        params->chisqParams->a1             = input->segment->a1 ;
        params->chisqParams->b1             = input->segment->b1 ;
        params->chisqParams->b2             = input->segment->b2 ;
#if 0
        params->chisqParams->bankMatch   = input->tmplt->minMatch;
#endif

        /* compute the chisq threshold: this is slow! */
        LALFindChirpBCVChisqVeto( status->statusPtr, params->chisqVec,
            params->chisqInput, params->chisqInputBCV, params->chisqParams );
        CHECKSTATUSPTR (status);

        haveChisq = 1;
      }


      /* if we don't have a chisq or the chisq drops below the            */
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
        else if (j > thisEvent->end_time.gpsSeconds + deltaEventIndex ||
              ! params->maximiseOverChirp)
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
          thisEvent->template_duration = (REAL8) chirpTime;

          /* record the ifo and channel name for the event */
          strncpy( thisEvent->ifo, input->segment->data->name,
              2 * sizeof(CHAR) );
          strncpy( thisEvent->channel, input->segment->data->name + 3,
              (LALNameLength - 3) * sizeof(CHAR) );
          thisEvent->impulse_time = thisEvent->end_time;

          /* record coalescence phase and alpha */

          /* calculate the numerators and the denominators */
          Num1 = qBCV[timeIndex].re + q[timeIndex].im ;
          Num2 = qBCV[timeIndex].re - q[timeIndex].re ;
          Den1 = q[timeIndex].re - qBCV[timeIndex].im ;
          Den2 = q[timeIndex].re + qBCV[timeIndex].im ;

          InvTan1 = (REAL4) atan2(Num1, Den1);
          InvTan2 = (REAL4) atan2(Num2, Den2);

          thisEvent->coa_phase = - 0.5 * InvTan1 + 0.5 * InvTan2 ;
          omega = 0.5 * InvTan1 + 0.5 * InvTan2 ;
          thisEvent->alpha = - input->segment->b2 * tan(omega) / 
            ( input->segment->a1 + input->segment->b1*tan(omega) );
          /* actually record alpha * fcut^(2/3) which must be b/t 0 and 1 */
          thisEvent->alpha *= pow((input->tmplt->fFinal) , 2.0/3.0);   

          /* copy the template into the event */
          thisEvent->psi0   = (REAL4) input->tmplt->psi0; 
          thisEvent->psi3   = (REAL4) input->tmplt->psi3;
          /* chirp mass in units of M_sun */
          thisEvent->mchirp = (1.0 / LAL_MTSUN_SI) * LAL_1_PI *
            pow( 3.0 / 128.0 / input->tmplt->psi0 , 3.0/5.0 );
          m =  fabs(thisEvent->psi3) / 
            (16.0 * LAL_MTSUN_SI * LAL_PI * LAL_PI * thisEvent->psi0) ;
          thisEvent->eta = 3.0 / (128.0*thisEvent->psi0 * 
              pow( (m*LAL_MTSUN_SI*LAL_PI), (5.0/3.0)) );
          thisEvent->f_final  = (REAL4) input->tmplt->fFinal ;

          /* set the type of the template used in the analysis */
          LALSnprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
              "FindChirpBCV" );

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
#if 0
          thisEvent->sigmasq = norm * input->segment->segNorm *
            input->segment->segNorm * input->fcTmplt->tmpltNorm;
          thisEvent->eff_distance =
            (input->fcTmplt->tmpltNorm * input->segment->segNorm *
             input->segment->segNorm) / thisEvent->snr;
          thisEvent->eff_distance = sqrt( thisEvent->eff_distance );
#endif

          thisEvent->snr *= norm;      
          thisEvent->snr = sqrt( thisEvent->snr );

          /* compute the time since the snr crossing */
          thisEvent->event_duration = 
            (REAL8) timeIndex - (REAL8) eventStartIdx;
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
    thisEvent->template_duration = (REAL8) chirpTime; 

    /* record the ifo name for the event */
    strncpy( thisEvent->ifo, input->segment->data->name,
        2 * sizeof(CHAR) );
    strncpy( thisEvent->channel, input->segment->data->name + 3,
        (LALNameLength - 3) * sizeof(CHAR) );
    thisEvent->impulse_time = thisEvent->end_time;

    /* record coalescence phase and alpha */

    /* calculate the numerators and the denominators */
    Num1 = qBCV[timeIndex].re + q[timeIndex].im ;
    Num2 = qBCV[timeIndex].re - q[timeIndex].re ;
    Den1 = q[timeIndex].re - qBCV[timeIndex].im ;
    Den2 = q[timeIndex].re + qBCV[timeIndex].im ;

    InvTan1 = (REAL4) atan2(Num1, Den1);
    InvTan2 = (REAL4) atan2(Num2, Den2 );


    thisEvent->coa_phase = - 0.5 * InvTan1 + 0.5 * InvTan2 ;
    omega = 0.5 * InvTan1 + 0.5 * InvTan2 ;
    thisEvent->alpha = - input->segment->b2 * tan(omega) /
      ( input->segment->a1 + input->segment->b1*tan(omega) );
    /* actually record alpha * fcut^(2/3) which must be b/t 0 and 1 */
    thisEvent->alpha *= pow( (input->tmplt->fFinal), 2.0/3.0);   


    /* copy the template into the event */
    thisEvent->psi0   = (REAL4) input->tmplt->psi0;   
    thisEvent->psi3   = (REAL4) input->tmplt->psi3;  
    /* chirp mass in units of M_sun */
    thisEvent->mchirp = (1.0 / LAL_MTSUN_SI) * LAL_1_PI *
      pow( 3.0 / 128.0 / input->tmplt->psi0, 3.0/5.0 );
    thisEvent->f_final  = (REAL4) input->tmplt->fFinal;
    m =  fabs(thisEvent->psi3) /
          (16.0 * LAL_MTSUN_SI * LAL_PI * LAL_PI * thisEvent->psi0) ;
    thisEvent->eta = 3.0 / (128.0*thisEvent->psi0 *
          pow( (m*LAL_MTSUN_SI*LAL_PI), (5.0/3.0)) );



    /* set the type of the template used in the analysis */
    LALSnprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
        "FindChirpBCV" );

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
#if 0    
    thisEvent->sigmasq = norm * input->segment->segNorm *
      input->segment->segNorm * input->fcTmplt->tmpltNorm;
    thisEvent->eff_distance =
      (input->fcTmplt->tmpltNorm * input->segment->segNorm *
       input->segment->segNorm) / thisEvent->snr;
    thisEvent->eff_distance = sqrt( thisEvent->eff_distance ); 
#endif      

    thisEvent->snr *=  norm ;   
    thisEvent->snr = sqrt( thisEvent->snr );

    /* compute the time since the snr crossing */
    thisEvent->event_duration = (REAL8) timeIndex - (REAL8) eventStartIdx;
    thisEvent->event_duration *= (REAL8) deltaT;
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
