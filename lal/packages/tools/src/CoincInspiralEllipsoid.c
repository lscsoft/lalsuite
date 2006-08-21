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
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/XLALError.h>
#include <lal/EllipsoidOverlapTools.h>

NRCSID( COINCINSPIRALUTILSC, "$Id$" );

typedef struct tagTriggerErrorList
{
  SnglInspiralTable          *trigger;
  gsl_matrix                 *err_matrix;
  gsl_vector                 *position;
  struct tagTriggerErrorList *next;
}
TriggerErrorList;

/* <lalVerbatim file="CoincInspiralEllipsoidCP"> */
void
LALCreateTwoIFOCoincListEllipsoid(
    LALStatus                  *status,
    CoincInspiralTable        **coincOutput,
    SnglInspiralTable          *snglInput
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

  INITSTATUS( status, "LALCreateTwoIFOCoincList", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( snglInput, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( coincOutput, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( ! *coincOutput, status, 
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );

  memset( currentTriggerNS, 0, 2 * sizeof(INT8) );

  XLALClearErrno();
  /* Loop through triggers and assign each of them an error ellipsoid */
  for (currentTrigger = snglInput; currentTrigger;
      currentTrigger = currentTrigger->next)
  {
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
    thisErrorList->err_matrix = XLALGetErrorMatrixFromSnglInspiral( currentTrigger );
    if (!thisErrorList->err_matrix)
    {
      XLALClearErrno();
      ABORTXLAL( status );
    }
    thisErrorList->position   = XLALGetPositionFromSnglInspiral( currentTrigger );
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
  
  maxTimeDiff = (INT8) 2e7;    
  maxTimeDiff += (INT8) ( 1e9 * 2 * LAL_REARTH_SI / LAL_C_SI );
  
  for ( currentError[0] = errorListHead; currentError[0]->next;
      currentError[0] = currentError[0]->next)
  {

    /* calculate the time of the trigger */
    LALGPStoINT8( status->statusPtr, &currentTriggerNS[0], 
        &(currentError[0]->trigger->end_time) );

    /* set next trigger for comparison */
    currentError[1] = currentError[0]->next;
    LALGPStoINT8( status->statusPtr, &currentTriggerNS[1], 
          &(currentError[1]->trigger->end_time) );

    while ( (currentTriggerNS[1] - currentTriggerNS[0]) < maxTimeDiff )
    {
      REAL8 overlap;

      if (strcmp(currentError[0]->trigger->ifo, currentError[1]->trigger->ifo))
      {

        /* check for the intersection of the ellipsoids */
        workSpace->invQ1 = currentError[0]->err_matrix;
        workSpace->invQ2 = currentError[1]->err_matrix;
        overlap = XLALCheckOverlapOfEllipsoids( currentError[0]->position, 
                     currentError[1]->position, workSpace );
        if (XLAL_IS_REAL8_FAIL_NAN(overlap))
        {
          XLALClearErrno();
          ABORTXLAL( status ); 
        }

        /* test whether we have coincidence */
        if ( overlap <= 1.0 )
        {
          LALInfo( status, "Found double coincident trigger,");
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
