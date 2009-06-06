/*
*  Copyright (C) 2007 Alexander Dietz, Drew Keppel, Nickolas Fotopoulos, Reinhard Prix, Craig Robinson
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

#if 0
<lalLaTeX>
\subsection{Module \texttt{CoincInspiralEllipsoid.c}}

\noindent Blah.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{CoincInspiralEllipsoidCP}
\idx{LALCreateTwoIFOCoincListEllipsoid()}
\idx{XLALSnglInspiralCoincTestEllipsoid()}
\idx{XLALCompareInspiralsEllipsoid()}
\idx{XLALCalculateEThincaParameter()}
\idx{XLALCalculateEThincaParameterForInjection()}

\subsubsection*{Description}

\texttt{LALCreateTwoIFOCoincListEllipsoid()} takes in a linked list of
single inspiral tables and returns a list of two instrument coincidences.
To determine coincidence, the triggers are modelled as ellipsoids in the
parameter space. Triggers are deemed to be coincident if these ellipsoids
are found to overlap.The ellipsoid scaling factor is given within the
\texttt{accuracyParams} structure. When single inspirals from two different
instruments are found to be coincident, the code creates a new
\texttt{coincInspiralTable} and uses \texttt{LALAddSnglInspiralToCoinc()}
to add the single inspirals to the coinc. The function returns
\texttt{coincOutput} which is a pointer to the head of a linked list of
\texttt{CoincInspiralTable}s.

\texttt{XLALSnglInspiralCoincTestEllipsoid()} is used in the creation of
multiple IFO coincident events. It is called by \texttt{LALCreateNIFOCoincList()}
when the coincidence test is set to be ellipsoid. Unlike in other coincidence
tests, coincidence here is determined by the use of event ids as opposed to
calling the comparison function. This is because the test for ellipsoid overlap
uses matrix inversions and function maximizations, which are potentially costly
operations. If all members of the coinc are found to be
coincident with the single, then \texttt{accuracyParams.match} is set to 1,
otherwise to 0.

\texttt{XLALCompareInspiralsEllipsoid()} checks for the overlap of ellipsoids
associated with two single inspiral tables. The ellipsoid scaling factor is
provided by \texttt{accuracyParams}. If the ellipsoids are found to overlap,
1 is returned; otherwise 0 is returned.

\texttt{XLALCalculateEThincaParameter()} calculates the maximum e-thinca
parameter between two single inspiral tables. It does this using a bisection
method, and uses \texttt{XLALCompareInspiralsEllipsoid()} to check for overlap.
The maximum value for the e-thinca parameter is returned. If the two triggers
do not overlap for an e-thinca parameter of 2.0, the triggers are not
coincident, and an error is thrown.

\texttt{XLALCalculateEThincaParameterForInjection()} takes in a
\texttt{SnglInspiralTable} and a \texttt{SimInspiralTable}, and returns the
e-thinca parameter between the trigger and the injection. This amounts to
calculating the square of the metric distance between the two points in
$(t_C, \tau_0, \tau_3)$ space.

\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent

\subsubsection*{Notes}
%% Any relevant notes.

\vfill{\footnotesize\input{CoincInspiralEllipsoidCV}}

</lalLaTeX>
#endif

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
   * light travel time for earths diameter
   * (detectors cant be further apart than this) */

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

        accuracyParams->match = 0;
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

    /* if analyzing a grb the position and time-delay is KNOWN */
    if (params->exttrig)
    {
      timeShift=travelTime;
    }
    else
    {
      timeShift = -travelTime;
    }

    /* Reset the times to avoid any precision problems */
    XLALSetTimeInPositionVector( aPtr->position, timeShift );
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
          SnglInspiralTable *table2,
          InspiralAccuracyList* accuracyParams
             )
/* </lalVerbatim> */
{

   static const char *func = "XLALCalculateEThincaParameter";

   TriggerErrorList   errorList[2];

   REAL8 loMatch, hiMatch, midMatch;
   INT4 i;
   INT2  isOverlap;
   fContactWorkSpace    *workSpace;

   loMatch = 0.0001;
   hiMatch = 2.0;

  if ( !table1 || !table2 || !accuracyParams)
    XLAL_ERROR_REAL8( func, XLAL_EFAULT );

  workSpace = XLALInitFContactWorkSpace( 3, NULL, NULL, gsl_min_fminimizer_brent, 1.0e-2 );
  if (!workSpace)
  {
    XLAL_ERROR_REAL8( func, XLAL_EFUNC );
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
                                 loMatch );
  }

  /* Check to see if they overlap for the lowest difference possible. If so there is
   * nothing more to be done.
   */
   isOverlap = XLALCompareInspiralsEllipsoid( &errorList[0], &errorList[1], workSpace,
                                        accuracyParams );

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
     return loMatch;
   }

  /* Now test for the largest error */
  for (i = 0; i < 2; i++ )
  {
    errorList[i].err_matrix = XLALGetErrorMatrixFromSnglInspiral( errorList[i].trigger,
                                 hiMatch );
  }

   isOverlap = XLALCompareInspiralsEllipsoid( &errorList[0], &errorList[1], workSpace,
                                        accuracyParams );

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
                                        accuracyParams );

     for (i = 0; i < 2; i++ )
     {
       gsl_matrix_free( errorList[i].err_matrix );
     }
     if (isOverlap)
     {
       hiMatch = midMatch;
     }
     else
     {
       loMatch = midMatch;
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


/* This function returns the e-thinca parameter between a trigger and an injection */
/* <lalVerbatim file="CoincInspiralEllipsoidCP"> */
REAL8 XLALEThincaParameterForInjection(
                    SimInspiralTable  *injection,
                    SnglInspiralTable *trigger
                    )
/* </lalVerbatim> */
{

  static const char *func = "XLALEThincaParameterForInjection";

  /* Trigger parameters */
  REAL8 fLower;
  REAL8 mTotal;
  REAL8 tau0;
  REAL8 eta;

  /* Inj parameters */
  InspiralTemplate injTmplt;   /* Used to calculate parameters */
  INT8 injEndTime;

  /* Parameter differences */
  REAL8 dtC, dt0, dt3;

  REAL8 eMatch;

  LALStatus status;

  if ( !injection || !trigger )
    XLAL_ERROR_REAL8( func, XLAL_EFAULT );

  memset( &status, 0, sizeof(LALStatus));

  /* We need the tau0 and tau3 for the injection */
  /* To do this, we will need to calculate the lower frequency used in the
   * case of the inspiral trigger */

  mTotal = trigger->mass1 + trigger->mass2;
  eta    = trigger->eta;
  tau0   = trigger->tau0;

  mTotal = mTotal * LAL_MTSUN_SI;

  fLower = 5.0 / (256.0 * eta * pow(mTotal, 5.0/3.0) * tau0 );
  fLower = pow(fLower, 3.0/8.0) / LAL_PI;

  XLALPrintInfo("%s: fLower found to be %e\n", func, fLower );

  /* Now populate the inspiral template with relevant parameters */
  injTmplt.mass1      = injection->mass1;
  injTmplt.mass2      = injection->mass2;
  injTmplt.massChoice = m1Andm2;
  injTmplt.fLower     = fLower;
  injTmplt.order      = LAL_PNORDER_THREE_POINT_FIVE;

  LALInspiralParameterCalc( &status, &injTmplt );

  /* Get the GPS time from the injection*/
  injEndTime = XLALReturnSimInspiralEndTime( injection, trigger->ifo );

  dtC = ( injEndTime - XLALGPStoINT8( &(trigger->end_time) ) ) * 1.0e-9;
  dt0 = injTmplt.t0 - trigger->tau0;
  dt3 = injTmplt.t3 - trigger->tau3;

  eMatch = trigger->Gamma[0] * dtC * dtC
         + 2.0 * trigger->Gamma[1] * dtC * dt0
         + 2.0 * trigger->Gamma[2] * dtC * dt3
         + trigger->Gamma[3] * dt0 * dt0
         + 2.0 * trigger->Gamma[4] * dt0 * dt3
         + trigger->Gamma[5] * dt3 * dt3;

  return eMatch;
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

  a11 = table->Gamma[0] / eMatch;
  a12 = table->Gamma[1] / eMatch;
  a13 = table->Gamma[2] / eMatch;
  a22 = table->Gamma[3] / eMatch;
  a23 = table->Gamma[4] / eMatch;
  a33 = table->Gamma[5] / eMatch;

  x = (a23 * a23 - a22 * a33) * a22;
  denom = (a12*a23 - a22*a13) * (a12*a23 - a22*a13)
              - (a23*a23 - a22*a33) * (a12*a12 - a22*a11);

  return ( sqrt( x / denom ));
}
