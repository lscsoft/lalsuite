/*----------------------------------------------------------------------------
 *
 * File Name: CoincInspiralEllipsoid.c
 *
 * Author: Craig Robinson
 *
 * $Id$
 *
 *---------------------------------------------------------------------------*/

#if 0
<lalVerbatim file="CoincInspiralEllipsoidCV">
Author: Craig Robinson
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
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/XLALError.h>
#include <lal/EllipsoidOverlapTools.h>
#include <lal/CoincInspiralEllipsoid.h>

NRCSID( COINCINSPIRALELLIPSOIDC, "$Id$" );


static REAL8 getTimeError(const SnglInspiralTable *table, REAL8 eMatch);

/* <lalVerbatim file="CoincInspiralEllipsoidCP"> */
void
LALCreateTwoIFOCoincListEllipsoid(
    LALStatus                  *status,
    CoincInspiralTable        **coincOutput,
    SnglInspiralTable          *snglInput,
    InspiralAccuracyList       *accuracyParams
    )
/* </lalVerbatim> */
{
  SnglInspiralTable            *currentTrigger = NULL;
  INT8                          currentTriggerNS[2];
  CoincInspiralTable           *coincHead = NULL;
  CoincInspiralTable           *thisCoinc = NULL;
  INT4                          numEvents = 0;
  INT8                          maxTimeDiff = 0;
  TriggerErrorList             *errorListHead = NULL;
  TriggerErrorList             *thisErrorList = NULL;
  TriggerErrorList             *currentError[2];
  fContactWorkSpace            *workSpace;
  REAL8                        timeError = 0.0;


  INITSTATUS( status, "LALCreateTwoIFOCoincList", COINCINSPIRALELLIPSOIDC );
  ATTATCHSTATUSPTR( status );

  ASSERT( snglInput, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( coincOutput, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( ! *coincOutput, status, 
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );

  memset( currentTriggerNS, 0, 2 * sizeof(INT8) );

  printf("e-Match = %e\n", accuracyParams->eMatch);
/*  XLALClearErrno();*/
  /* Loop through triggers and assign each of them an error ellipsoid */
  for (currentTrigger = snglInput; currentTrigger;
      currentTrigger = currentTrigger->next)
  {
    REAL8 thisTimeError;

    if (!errorListHead)
    {
      errorListHead = (TriggerErrorList *) LALCalloc(1, sizeof(TriggerErrorList));
      thisErrorList = errorListHead;
    }
    else
    {
      thisErrorList->next = (TriggerErrorList *) LALCalloc(1, sizeof(TriggerErrorList));
      thisErrorList = thisErrorList->next;
    }
    thisErrorList->trigger    = currentTrigger;
    thisErrorList->err_matrix = XLALGetErrorMatrixFromSnglInspiral( currentTrigger,
                                  accuracyParams->eMatch );
    if (!thisErrorList->err_matrix)
    {
      XLALClearErrno();
      ABORTXLAL( status );
    }
    thisErrorList->position   = XLALGetPositionFromSnglInspiral( currentTrigger );
    thisTimeError = getTimeError(currentTrigger, accuracyParams->eMatch);
    if (thisTimeError > timeError)
      timeError = thisTimeError;
  }


  /* Initialise the workspace for ellipsoid overlaps */
  workSpace = XLALInitFContactWorkSpace( 3, NULL, NULL, gsl_min_fminimizer_brent, 1.0e-2 );
  if (!workSpace)
  {
    XLALClearErrno();
    ABORTXLAL( status );
  }

  /* calculate the maximum time delay 
   * set it equal to 2 * worst IFO timing accuracy plus
   * light travel time for earth's diameter 
   * (detectors can't be further apart than this) */
  
  maxTimeDiff = (INT8) (1e9 * 2.0 * timeError);    
  maxTimeDiff += (INT8) ( 1e9 * 2 * LAL_REARTH_SI / LAL_C_SI );
  
  for ( currentError[0] = errorListHead; currentError[0]->next;
      currentError[0] = currentError[0]->next)
  {

    /* calculate the time of the trigger */
    currentTriggerNS[0] = XLALGPStoINT8(
                 &(currentError[0]->trigger->end_time) );

    /* set next trigger for comparison */
    currentError[1] = currentError[0]->next;
    currentTriggerNS[1] = XLALGPStoINT8(
                 &(currentError[1]->trigger->end_time) );

    while ( (currentTriggerNS[1] - currentTriggerNS[0]) < maxTimeDiff )
    {

      INT2 match;

      /* test whether we have coincidence */
      match = XLALCompareInspiralsEllipsoid( currentError[0],
                 currentError[1], workSpace, accuracyParams );
      if ( match < 0 )
      {
        /* Error in the comparison function */
        thisErrorList = errorListHead;
        while (thisErrorList)
        {
          errorListHead = thisErrorList->next;
          gsl_matrix_free( thisErrorList->err_matrix );
          gsl_vector_free( thisErrorList->position );
          LALFree( thisErrorList );
          thisErrorList = errorListHead;
        }

        XLALFreeFContactWorkSpace( workSpace );
        ABORTXLAL( status );
      }

      /* Check whether the event was coincident */
      if ( match )
      {
        /* create a 2 IFO coinc and store */
        if ( ! coincHead  )
        {
          coincHead = thisCoinc = (CoincInspiralTable *) 
            LALCalloc( 1, sizeof(CoincInspiralTable) );
        }
        else
        {
          thisCoinc = thisCoinc->next = (CoincInspiralTable *) 
            LALCalloc( 1, sizeof(CoincInspiralTable) );
        }

        /* Add the two triggers to the coinc */
        LALAddSnglInspiralToCoinc( status->statusPtr, &thisCoinc,
            currentError[0]->trigger );
        LALAddSnglInspiralToCoinc( status->statusPtr, &thisCoinc,
            currentError[1]->trigger );

        ++numEvents;
         
      }

      /* scroll on to the next sngl inspiral */
      
      if ( (currentError[1] = currentError[1]->next) )
      {
        LALGPStoINT8( status->statusPtr, &currentTriggerNS[1], 
            &(currentError[1]->trigger->end_time) );
      }
      else
      {
        LALInfo(status, "Second trigger has reached end of list");
        break;
      }
    }
  }

  *coincOutput = coincHead;

  /* Free all the memory allocated for the ellipsoid overlap */
  thisErrorList = errorListHead;
  while (thisErrorList)
  {
    errorListHead = thisErrorList->next;
    gsl_matrix_free( thisErrorList->err_matrix );
    gsl_vector_free( thisErrorList->position );
    LALFree( thisErrorList );
    thisErrorList = errorListHead;
  }

  XLALFreeFContactWorkSpace( workSpace );

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="CoincInspiralEllipsoidCP"> */
void
LALCreateNIFOCoincListEllipsoid(
    LALStatus                  *status,
    CoincInspiralTable        **coincHead,
    InspiralAccuracyList       *accuracyParams,
    INT4                        N
    )
/* </lalVerbatim> */
{
  INT4                          numEvents  = 0;
  InterferometerNumber          ifoNumber  = LAL_UNKNOWN_IFO;
  InterferometerNumber          ifoNum     = LAL_UNKNOWN_IFO;
  InterferometerNumber          firstEntry = LAL_UNKNOWN_IFO;
  CoincInspiralTable           *thisCoinc     = NULL;
  CoincInspiralTable           *lastCoinc     = NULL;
  CoincInspiralTable           *otherCoinc    = NULL;
  CoincInspiralTable           *nIfoCoincHead = NULL;
  CoincInspiralTable           *thisNIfoCoinc = NULL;
  EventIDColumn                *eventIDHead   = NULL;
  EventIDColumn                *thisID        = NULL;


  INITSTATUS( status, "LALCreateNIFOCoincList", COINCINSPIRALELLIPSOIDC );
  ATTATCHSTATUSPTR( status );

  /* loop over all the coincidences in the list */
  for( thisCoinc = *coincHead; thisCoinc; thisCoinc = thisCoinc->next)
  {
    lastCoinc = thisCoinc;

    /* check that this is an (N-1) coinc */
    if ( thisCoinc->numIfos == N - 1 )
    {
      /* look up the first single inspiral */
      for ( firstEntry = 0; firstEntry < LAL_NUM_IFO; firstEntry++)
      {
        if ( thisCoinc->snglInspiral[firstEntry] )
        {
          LALInfo( status, "Found the first entry in the coinc" );
          break;
        }
      }

      /* get the list of event IDs for this first entry */
      eventIDHead = thisCoinc->snglInspiral[firstEntry]->event_id;

      /* loop over the (N-1) ifo coincs that first entry is a member of
       * and try to find an N ifo coinc */
      for( thisID = eventIDHead; thisID; thisID = thisID->next )
      {
        otherCoinc = thisID->coincInspiralTable;

        if( otherCoinc->numIfos == N - 1 )
        {
          /* loop over all singles which are alphabetically before the
           * first one in thisCoinc */
          for( ifoNumber = 0; ifoNumber < firstEntry; ifoNumber++ )
          {
            /* test whether we have an N ifo coincidence */
            accuracyParams->match = 0;

            if ( otherCoinc->snglInspiral[ifoNumber] )
            {
              XLALSnglInspiralCoincTestEllipsoid( thisCoinc,
                  otherCoinc->snglInspiral[ifoNumber], accuracyParams );
            }

            if ( accuracyParams->match )
            {
              LALInfo( status, "We have found an N ifo coinc, storing");
              ++numEvents;

              /* create a N IFO coinc and store */
              if ( ! nIfoCoincHead  )
              {
                nIfoCoincHead = thisNIfoCoinc = (CoincInspiralTable *)
                  LALCalloc( 1, sizeof(CoincInspiralTable) );
              }
              else
              {
                thisNIfoCoinc = thisNIfoCoinc->next = (CoincInspiralTable *)
                  LALCalloc( 1, sizeof(CoincInspiralTable) );
              }

              /* add the single to the new N coinc */
              LALAddSnglInspiralToCoinc( status->statusPtr, &thisNIfoCoinc,
                  otherCoinc->snglInspiral[ifoNumber] );

              /* add the triggers from the (N-1) coinc to the new N coinc */
              for( ifoNum = 0; ifoNum < LAL_NUM_IFO; ifoNum++ )
              {
                if( thisCoinc->snglInspiral[ifoNum] )
                {
                  LALAddSnglInspiralToCoinc( status->statusPtr, &thisNIfoCoinc,
                      thisCoinc->snglInspiral[ifoNum] );
                }
              }
            } /* closes: if ( accuracyParams->match ) */
          }
        }
      } /* closes: for( thisID = eventIDHead; thisID; thisID->next ) */
    }
  } /* closes: for( thisCoinc = coincHead; thisCoinc;
     *              thisCoinc = thisCoinc->next) */

  /* append the N ifo coincs to the end of the linked list */
  if ( lastCoinc )
  {
    lastCoinc->next = nIfoCoincHead;
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);

}

/* <lalVerbatim file="CoincInspiralEllipsoidCP"> */
void
XLALSnglInspiralCoincTestEllipsoid(
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral,
    InspiralAccuracyList       *accuracyParams
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *thisCoincEntry;
  INT4                  match = 1;
  INT4                  ifoNumber = 0;

  static const char *func = "XLALSnglInspiralCoincTest";


  printf("Calling the ellipsoid NIfo coinc function");

  /* Loop over sngl_inspirals contained in coinc_inspiral */
  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    thisCoincEntry = coincInspiral->snglInspiral[ifoNumber];

    if ( thisCoincEntry )
    {
      /* snglInspiral entry exists for this IFO, perform coincidence test */
      if ( ifoNumber == XLALIFONumber(snglInspiral->ifo) )
      {
        XLALPrintInfo( "We already have a coinc from this IFO" );
        accuracyParams->match = 0;
      }

      else
      {
        EventIDColumn  *eventIDHead = thisCoincEntry->event_id;
        EventIDColumn  *thisID      = NULL;
        CoincInspiralTable *thisCoinc;

        for ( thisID = eventIDHead; thisID; thisID = thisID->next )
        {
          thisCoinc = thisID->coincInspiralTable;
          if ( thisCoinc->snglInspiral[XLALIFONumber(snglInspiral->ifo)] == snglInspiral )
          {
            accuracyParams->match = 1;
            break;
          }
        }  
      }
      /* set match to zero if no match.  Keep same if match */
      match *= accuracyParams->match;
    }
  }
  /* returm errorParams->match to be 1 if we match, zero otherwise */
  accuracyParams->match = match;
  if ( accuracyParams->match == 0 )
    XLALPrintInfo( "Coincidence test failed" );
  if ( accuracyParams->match == 1 )
    XLALPrintInfo( "Coincidence test passed" );

  return;
}

/* <lalVerbatim file="CoincInspiralEllipsoidCP"> */
INT2 XLALCompareInspiralsEllipsoid(
      TriggerErrorList              *aPtr,
      TriggerErrorList              *bPtr,
      fContactWorkSpace             *workSpace,
      InspiralAccuracyList          *params
      )
/* </lalVerbatim> */
{
 
  static const char *func = "XLALCompareInspiralsEllipsoid";

  INT2  isCoinc   = 0;

  REAL8 overlap;

  if (strcmp(aPtr->trigger->ifo, bPtr->trigger->ifo))
  {
    INT8 curTimeANS, curTimeBNS;
    REAL8 originalTimeA;
    REAL8 originalTimeB;
    REAL8 travelTime;
    REAL8 timeShift = 0.0;
    InterferometerNumber ifoaNum, ifobNum;

    ifoaNum = XLALIFONumber( aPtr->trigger->ifo );
    ifobNum = XLALIFONumber( bPtr->trigger->ifo );

    originalTimeA = gsl_vector_get( aPtr->position, 0 );
    originalTimeB = gsl_vector_get( bPtr->position, 0 );

    travelTime = params->lightTravelTime[ifoaNum][ifobNum] * 1.0e-9;

    curTimeANS = XLALGPStoINT8( &(aPtr->trigger->end_time) );
    curTimeBNS = XLALGPStoINT8( &(bPtr->trigger->end_time) );

    /* Reset the times to avoid any precision problems */
    XLALSetTimeInPositionVector( aPtr->position, 0.0 );
    XLALSetTimeInPositionVector( bPtr->position,
          (curTimeBNS - curTimeANS) * 1.0e-9 );

    /* Loop over the time shift to sweep the light travel time */
    while (timeShift <= travelTime)
    {
      /* check for the intersection of the ellipsoids */
      workSpace->invQ1 = aPtr->err_matrix;
      workSpace->invQ2 = bPtr->err_matrix;
      overlap = XLALCheckOverlapOfEllipsoids( aPtr->position,
                   bPtr->position, workSpace );
      if (XLAL_IS_REAL8_FAIL_NAN(overlap))
      {
        XLAL_ERROR( func, XLAL_EFUNC );
      }

      /* test whether we have coincidence */
      if ( overlap <= 1.0 )
      {
        isCoinc = 1;
        break;
      }

      timeShift += 1.0 / 4096.0;
      XLALSetTimeInPositionVector( aPtr->position, timeShift );
    }

    /* Set the times back to their correct values */
    XLALSetTimeInPositionVector( aPtr->position, originalTimeA );
    XLALSetTimeInPositionVector( bPtr->position, originalTimeB );
  }
  params->match = isCoinc;
  return isCoinc;
}

/* 
 * This function returns the largest time error associated with a
 * particular error ellipsoid.
 */
static REAL8 getTimeError(const SnglInspiralTable *table, REAL8 eMatch)
{
  REAL8 a11;
  REAL8 a23; 
  REAL8 a22;
  REAL8 a33;
  REAL8 a12;
  REAL8 a13;
  REAL8 x;
  REAL8 denom;
  
  a11 = table->Gamma[0] / (1.0 - eMatch);
  a12 = table->Gamma[1] / (1.0 - eMatch);
  a13 = table->Gamma[2] / (1.0 - eMatch);
  a22 = table->Gamma[3] / (1.0 - eMatch);
  a23 = table->Gamma[4] / (1.0 - eMatch);
  a33 = table->Gamma[5] / (1.0 - eMatch);

  x = (a23 * a23 - a22 * a33) * a22;
  denom = (a12*a23 - a22*a13) * (a12*a23 - a22*a13)
              - (a23*a23 - a22*a33) * (a12*a12 - a22*a11);

  return ( sqrt( x / denom ));
}
