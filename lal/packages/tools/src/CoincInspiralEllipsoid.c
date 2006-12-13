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


/* <lalVerbatim file="CoincInspiralEllipsoidCP"> */
REAL8 XLALCalculateEThincaParameter( 
          SnglInspiralTable *table1,
          SnglInspiralTable *table2
             )
/* </lalVerbatim> */
{

   static const char *func = "XLALCalculateEThincaParameter";

   TriggerErrorList   errorList[2];

   REAL8 loMatch, hiMatch, midMatch;
   
   INT4 ifoNumber, ifoTwo, i;
   INT2  isOverlap;
   LALDetector aDet, bDet;
   InspiralAccuracyList accuracyParams;
   fContactWorkSpace    *workSpace;
   

   loMatch = 0.0;
   hiMatch = 0.9999;

  if ( !table1 || !table2 )
    XLAL_ERROR_REAL8( func, XLAL_EFAULT );

  workSpace = XLALInitFContactWorkSpace( 3, NULL, NULL, gsl_min_fminimizer_brent, 1.0e-2 );
  if (!workSpace)
  {
    XLAL_ERROR_REAL8( func, XLAL_EFUNC );
  }

  memset( &accuracyParams, 0, sizeof( InspiralAccuracyList ));

  /* Populate the lightTravel matrix */
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    XLALReturnDetector( &aDet, ifoNumber );

    for ( ifoTwo = 0; ifoTwo < LAL_NUM_IFO; ifoTwo++)
    {
      XLALReturnDetector( &bDet, ifoTwo );

      accuracyParams.lightTravelTime[ ifoNumber][ ifoTwo ] =
      XLALLightTravelTime( &aDet, &bDet );
    }
  }

  /* Set up the trigger lists */
  if ( XLALGPStoINT8( &(table1->end_time) ) < XLALGPStoINT8( &(table2->end_time )) )
  {
     errorList[0].trigger = table1;
     errorList[1].trigger = table2;
  }
  else
  {
     errorList[0].trigger = table2;
     errorList[1].trigger = table1;
  }

  for (i = 0; i < 2; i++ )
  {
    errorList[i].position = XLALGetPositionFromSnglInspiral( errorList[i].trigger );

    errorList[i].err_matrix = XLALGetErrorMatrixFromSnglInspiral( errorList[i].trigger,
                                 hiMatch );
  }

  /* Check to see if they overlap for the highest match possible. If so there is
   * nothing more to be done.
   */
   isOverlap = XLALCompareInspiralsEllipsoid( &errorList[0], &errorList[1], workSpace,
                                        &accuracyParams );

   for (i = 0; i < 2; i++ )
   {
     gsl_matrix_free( errorList[i].err_matrix );
   }

   if (isOverlap)
   {
     for (i = 0; i < 2; i++ )
     {
       gsl_vector_free( errorList[i].position );
     }
     XLALFreeFContactWorkSpace( workSpace );
     return hiMatch;
   }

  /* Now test for the lowest match */
  for (i = 0; i < 2; i++ )
  {

    errorList[i].err_matrix = XLALGetErrorMatrixFromSnglInspiral( errorList[i].trigger,
                                 loMatch );
  }
   
   isOverlap = XLALCompareInspiralsEllipsoid( &errorList[0], &errorList[1], workSpace,
                                        &accuracyParams );

   for (i = 0; i < 2; i++ )
   {
     gsl_matrix_free( errorList[i].err_matrix );
   }
   if (!isOverlap)
   {
     for (i = 0; i < 2; i++ )
     {
       gsl_vector_free( errorList[i].position );
     }
     XLALFreeFContactWorkSpace( workSpace );
     XLALPrintError("The two triggers provided are NOT coincident!!");
     XLAL_ERROR_REAL8( func, XLAL_EINVAL );
   }

   /* Now onto the algorithm proper */
   while (fabs( hiMatch - loMatch ) > 0.001 )
   {
     midMatch = (hiMatch + loMatch) / 2.0;


     for (i = 0; i < 2; i++ )
     {

      errorList[i].err_matrix = XLALGetErrorMatrixFromSnglInspiral( errorList[i].trigger,
                                 midMatch );
     }
     isOverlap = XLALCompareInspiralsEllipsoid( &errorList[0], &errorList[1], workSpace,
                                        &accuracyParams );

     for (i = 0; i < 2; i++ )
     {
       gsl_matrix_free( errorList[i].err_matrix );
     }
     if (isOverlap)
     {
       loMatch = midMatch;
     }
     else
     {
       hiMatch = midMatch;
     }
   }

   /* Free all accocated memory */
   for (i = 0; i < 2; i++ )
   {
     gsl_vector_free( errorList[i].position );
   }
   XLALFreeFContactWorkSpace( workSpace );

   return midMatch;

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
