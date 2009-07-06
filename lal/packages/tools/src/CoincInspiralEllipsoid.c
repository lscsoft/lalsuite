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
 * $Id: CoincInspiralEllipsoid.c,v 1.19 2007/06/08 14:41:57 bema Exp $
 *
 *---------------------------------------------------------------------------*/

#if 0
<lalVerbatim file="CoincInspiralEllipsoidCV">
Author: Craig Robinson
$Id: CoincInspiralEllipsoid.c,v 1.19 2007/06/08 14:41:57 bema Exp $
</lalVerbatim>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALErrno.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/XLALError.h>
#include <lal/EllipsoidOverlapTools.h>
#include <lal/TrigScanEThincaCommon.h>
#include <lal/CoincInspiralEllipsoid.h>

NRCSID( COINCINSPIRALELLIPSOIDC, "$Id: CoincInspiralEllipsoid.c,v 1.19 2007/06/08 14:41:57 bema Exp $" );

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

typedef struct tagEThincaMinimizer
{
   fContactWorkSpace  *workSpace;
   TriggerErrorList   *aPtr;
   TriggerErrorList   *bPtr;
}
EThincaMinimizer;

static REAL8 minimizeEThincaParameterOverTimeDiff( REAL8 timeShift,
                                                   void *minimizer
                                                 );

REAL8 XLALMinimizeEThincaParameterOverTravelTime( REAL8 travelTime,
                                                  EThincaMinimizer *minimizer,
                                                  INT4   exttrig
                                                );


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

  /* Loop through triggers and assign each of them an error ellipsoid */
  errorListHead = XLALCreateTriggerErrorList( snglInput, accuracyParams->eMatch, &timeError );
  if ( !errorListHead )
  {
    ABORTXLAL( status );
  }

  /* Initialise the workspace for ellipsoid overlaps */
  workSpace = XLALInitFContactWorkSpace( 3, NULL, NULL, gsl_min_fminimizer_brent, 1.0e-2 );
  if (!workSpace)
  {
    XLALDestroyTriggerErrorList( errorListHead );
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
    currentTriggerNS[0] = XLALGPSToINT8NS(
                 &(currentError[0]->trigger->end_time) );

    /* set next trigger for comparison */
    currentError[1] = currentError[0]->next;
    currentTriggerNS[1] = XLALGPSToINT8NS(
                 &(currentError[1]->trigger->end_time) );

    while ( (currentTriggerNS[1] - currentTriggerNS[0]) < maxTimeDiff )
    {

      INT2 match;

      /* test whether we have coincidence */
      match = XLALCompareInspiralsEllipsoid( currentError[0],
                 currentError[1], workSpace, accuracyParams );
      if ( match == XLAL_FAILURE )
      {
        /* Error in the comparison function */
        XLALDestroyTriggerErrorList( errorListHead );
        XLALFreeFContactWorkSpace( workSpace );
        ABORTXLAL( status );
      }

      /* Check whether the event was coincident */
      if ( match )
      {
        REAL8 etp = XLALCalculateEThincaParameter( currentError[0]->trigger,  currentError[1]->trigger, accuracyParams );
        printf( "%e\n", etp );
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
        if ( !thisCoinc )
        {
          /* Error allocating memory */
          thisCoinc = coincHead;
          while ( thisCoinc )
          {
            coincHead = thisCoinc->next;
            LALFree( thisCoinc );
            thisCoinc = coincHead;
          }
          XLALDestroyTriggerErrorList( errorListHead );
          XLALFreeFContactWorkSpace( workSpace );
          ABORT( status, LAL_NOMEM_ERR, LAL_NOMEM_MSG );
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
        currentTriggerNS[1] = XLALGPSToINT8NS( &(currentError[1]->trigger->end_time) );
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
  XLALDestroyTriggerErrorList( errorListHead );
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

  /*static const char *func = "XLALSnglInspiralCoincTest";*/


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

#ifndef LAL_NDEBUG
  if (!aPtr || !aPtr->trigger || !aPtr->position ||
      !bPtr || !bPtr->trigger || !bPtr->position ||
      !workSpace || !params)
    XLAL_ERROR( func, XLAL_EFAULT);
#endif

  if (strcmp(aPtr->trigger->ifo, bPtr->trigger->ifo))
  {
    REAL8 travelTime;
    InterferometerNumber ifoaNum, ifobNum;

    EThincaMinimizer minimizer;

    memset( &minimizer, 0, sizeof(EThincaMinimizer) );

    ifoaNum = XLALIFONumber( aPtr->trigger->ifo );
    ifobNum = XLALIFONumber( bPtr->trigger->ifo );

    travelTime = params->lightTravelTime[ifoaNum][ifobNum] * 1.0e-9;

    minimizer.workSpace = workSpace;
    minimizer.aPtr      = aPtr;
    minimizer.bPtr      = bPtr;

    overlap = XLALMinimizeEThincaParameterOverTravelTime( travelTime, &minimizer, params->exttrig );
    if ( XLAL_IS_REAL8_FAIL_NAN( overlap ) )
    {
      XLAL_ERROR( func, XLAL_EFUNC );
    }

    /* test whether we have coincidence */
    if ( overlap <= 1.0 )
    {
      isCoinc = 1;
    }
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

   TriggerErrorList * errorList[2];

   REAL8 travelTime;
   INT4 i;
   INT4 ifoANum, ifoBNum;
   fContactWorkSpace    *workSpace;

   REAL8 ethinca;
   EThincaMinimizer *minimizer = NULL;

#ifndef LAL_NDEBUG
  if ( !table1 || !table2 || !accuracyParams)
    XLAL_ERROR_REAL8( func, XLAL_EFAULT );
#endif

  memset( errorList, 0, 2 * sizeof(TriggerErrorList *) );

  /* Before we do anything, check we have triggers from different ifos */
  if ( !strcmp(table1->ifo, table2->ifo) )
  {
    XLALPrintError("%s: Triggers provided are from the same ifo!\n", func );
    XLAL_ERROR_REAL8( func, XLAL_EINVAL );
  }

  /* Allocate memory */
  errorList[0] = LALCalloc( 1, sizeof( TriggerErrorList ));
  if (errorList[0])
  {
    errorList[1] = errorList[0]->next = LALCalloc( 1, sizeof( TriggerErrorList ));
  }
  if (!errorList[0] || !errorList[1] )
  {
    if (errorList[0])
      XLALDestroyTriggerErrorList( errorList[0] );
    XLAL_ERROR_REAL8( func, XLAL_ENOMEM );
  }

  workSpace = XLALInitFContactWorkSpace( 3, NULL, NULL, gsl_min_fminimizer_brent, 1.0e-5 );
  if (!workSpace)
  {
    XLALDestroyTriggerErrorList( errorList[0] );
    XLAL_ERROR_REAL8( func, XLAL_EFUNC | XLALClearErrno() );
  }

  /* Set up the trigger lists */
  if ( XLALGPSToINT8NS( &(table1->end_time) ) < XLALGPSToINT8NS( &(table2->end_time )) )
  {
     errorList[0]->trigger = table1;
     errorList[1]->trigger = table2;
  }
  else
  {
     errorList[0]->trigger = table2;
     errorList[1]->trigger = table1;
  }

  for (i = 0; i < 2; i++ )
  {
    errorList[i]->position = XLALGetPositionFromSnglInspiral( errorList[i]->trigger );

    errorList[i]->err_matrix = XLALGetErrorMatrixFromSnglInspiral( errorList[i]->trigger,
                                 1.0 );
    if ( !errorList[i]->position || !errorList[i]->err_matrix )
    {
      XLALDestroyTriggerErrorList( errorList[0] );
      XLALFreeFContactWorkSpace( workSpace );
      XLAL_ERROR_REAL8( func, XLAL_ENOMEM );
    }
  }


  /* Set the travel time */
  ifoANum = XLALIFONumber( errorList[0]->trigger->ifo );
  ifoBNum = XLALIFONumber( errorList[1]->trigger->ifo );

  travelTime = accuracyParams->lightTravelTime[ifoANum][ifoBNum] * 1.0e-9;

  /* Create the e-thinca minimizer */
  minimizer = LALCalloc( 1, sizeof(EThincaMinimizer) );
  if ( !minimizer )
  {
    XLALDestroyTriggerErrorList( errorList[0] );
    XLALFreeFContactWorkSpace( workSpace );
    XLAL_ERROR_REAL8( func, XLAL_ENOMEM );
  }
  minimizer->workSpace = workSpace;
  minimizer->aPtr = errorList[0];
  minimizer->bPtr = errorList[1];

  ethinca = XLALMinimizeEThincaParameterOverTravelTime( travelTime, minimizer, accuracyParams->exttrig );
  if ( XLAL_IS_REAL8_FAIL_NAN( ethinca ) )
  {
    LALFree( minimizer );
    XLALDestroyTriggerErrorList( errorList[0] );
    XLALFreeFContactWorkSpace( workSpace );
    XLAL_ERROR_REAL8( func, XLAL_EFUNC );
  }

  LALFree ( minimizer );
  XLALDestroyTriggerErrorList( errorList[0] );
  XLALFreeFContactWorkSpace( workSpace );
  XLALPrintInfo( " I leave here and return an ethinca of %e.\n", ethinca );
  return ethinca;

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

#ifndef LAL_NDEBUG
  if ( !injection || !trigger )
    XLAL_ERROR_REAL8( func, XLAL_EFAULT );
#endif

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

  dtC = ( injEndTime - XLALGPSToINT8NS( &(trigger->end_time) ) ) * 1.0e-9;
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

REAL8 XLALMinimizeEThincaParameterOverTravelTime( REAL8 travelTime,
                                                  EThincaMinimizer *minimizer,
                                                  INT4   exttrig
                                                )
{

  static const char func[] = "XLALMinimizeEThincaParameterOverTravelTime";

  REAL8 ethinca;


  /* If colocated detectors or known sky position, just return the e-thinca parameter */
  if (travelTime == 0.0 || exttrig )
  {
    ethinca = minimizeEThincaParameterOverTimeDiff( travelTime, minimizer );
    if ( XLAL_IS_REAL8_FAIL_NAN(ethinca) )
    {
      XLAL_ERROR_REAL8( func, XLAL_EFUNC );
    }
    return ethinca;
  }
  else
  {
    gsl_function        F;
    INT4                min_status;
    INT4                iter = 0;
    const INT4          max_iter = 100;
    REAL8               epsilon = 1.0 / 16384.0;
    REAL8               m = 0.0;
    REAL8               a = - travelTime, b = travelTime; /* Upper and lower bounds */
    REAL8               minEThinca, maxEThinca;
    REAL8               midEThinca;
    gsl_min_fminimizer  *s = gsl_min_fminimizer_alloc( minimizer->workSpace->T );

    if ( !s )
    {
      XLAL_ERROR_REAL8( func, XLAL_ENOMEM );
    }


    F.function = &minimizeEThincaParameterOverTimeDiff;
    F.params   = minimizer;

    /* Calculate e-thinca parameter at start, end and mid points */
    minEThinca = minimizeEThincaParameterOverTimeDiff( a, minimizer );
    maxEThinca = minimizeEThincaParameterOverTimeDiff( b, minimizer );
    midEThinca = minimizeEThincaParameterOverTimeDiff( m, minimizer );
    if ( XLAL_IS_REAL8_FAIL_NAN(minEThinca) || XLAL_IS_REAL8_FAIL_NAN(maxEThinca)
         || XLAL_IS_REAL8_FAIL_NAN(midEThinca) )
    {
      gsl_min_fminimizer_free( s );
      XLAL_ERROR_REAL8( func, XLAL_EFUNC );
    }

    /* Check we have contained a minimum. Otherwise take appropriate action */
    if ( midEThinca >= minEThinca || midEThinca >= maxEThinca )
    {
      REAL8 testEThinca; /* To contain the lowest end-point */
      if ( minEThinca < maxEThinca )
      {
        testEThinca = minEThinca;
        m           = a + 2.0 * epsilon;
      }
      else
      {
        testEThinca = maxEThinca;
        m = b - 2.0 * epsilon;
      }
      midEThinca = minimizeEThincaParameterOverTimeDiff( m, minimizer );
      if ( XLAL_IS_REAL8_FAIL_NAN(midEThinca) )
      {
        gsl_min_fminimizer_free( s );
        XLAL_ERROR_REAL8( func, XLAL_EFUNC );
      }

      /* If we still don't have the minimum return the lowest end-point */
      if ( midEThinca >= testEThinca )
      {
        gsl_min_fminimizer_free( s );
        return testEThinca;
      }
    }

    /* Set up the GSL minimizer */
    XLAL_CALLGSL( min_status = gsl_min_fminimizer_set_with_values(s, &F,
                       m, midEThinca, a, minEThinca, b, maxEThinca) );
    if ( min_status != GSL_SUCCESS )
    {
      gsl_min_fminimizer_free( s );
      XLAL_ERROR_REAL8( func, XLAL_EFUNC );
    }

    /* Loop to perform the minimization */
    do
    {
        iter++;
        XLAL_CALLGSL( min_status = gsl_min_fminimizer_iterate (s) );
        if (min_status != GSL_SUCCESS )
        {
            gsl_min_fminimizer_free( s );
            XLAL_ERROR_REAL8( func, XLAL_EFUNC );
        }

        m = gsl_min_fminimizer_x_minimum (s);
        a = gsl_min_fminimizer_x_lower (s);
        b = gsl_min_fminimizer_x_upper (s);

        XLAL_CALLGSL( min_status = gsl_min_test_interval (a, b, epsilon, 0.0) );
        if (min_status != GSL_CONTINUE && min_status != GSL_SUCCESS )
        {
          gsl_min_fminimizer_free( s );
          XLAL_ERROR_REAL8( func, XLAL_EFUNC );
        }
    }
    while ( min_status == GSL_CONTINUE && iter < max_iter );
    /* End of minimization routine */

    /* Throw an error if max iterations would have been exceeded */
    if ( iter == max_iter && min_status == GSL_CONTINUE )
    {
      gsl_min_fminimizer_free( s );
      XLAL_ERROR_REAL8( func, XLAL_EMAXITER );
    }

    /* Get the minimum e-thinca param, and free memory for minimizer */
    ethinca = gsl_min_fminimizer_f_minimum( s );
    gsl_min_fminimizer_free( s );
    XLALPrintInfo( "%s: Number of iterations = %d\n", func, iter);
  }

  /* Return the required e-thinca value */
  return ethinca;
}

/* The following function would be called by the GSL minimizer */
static REAL8 minimizeEThincaParameterOverTimeDiff( REAL8 timeShift,
                                                   void *minimizer
                                                 )
{

  EThincaMinimizer *params = (EThincaMinimizer *) minimizer;


  INT8 curTimeANS, curTimeBNS;
  REAL8 originalTimeA;
  REAL8 originalTimeB;
  REAL8 overlap;

  originalTimeA = gsl_vector_get( params->aPtr->position, 0 );
  originalTimeB = gsl_vector_get( params->bPtr->position, 0 );

  curTimeANS = XLALGPSToINT8NS( &(params->aPtr->trigger->end_time) );
  curTimeBNS = XLALGPSToINT8NS( &(params->bPtr->trigger->end_time) );

  /* Reset the times to avoid any precision problems */
  XLALSetTimeInPositionVector( params->aPtr->position, timeShift );
  XLALSetTimeInPositionVector( params->bPtr->position,
          (curTimeBNS - curTimeANS) * 1.0e-9 );

  /* check for the intersection of the ellipsoids */
  params->workSpace->invQ1 = params->aPtr->err_matrix;
  params->workSpace->invQ2 = params->bPtr->err_matrix;
  overlap = XLALCheckOverlapOfEllipsoids( params->aPtr->position,
               params->bPtr->position, params->workSpace );
  if (XLAL_IS_REAL8_FAIL_NAN(overlap))
  {
     /* Set the times back to their correct values */
     XLALSetTimeInPositionVector( params->aPtr->position, originalTimeA );
     XLALSetTimeInPositionVector( params->bPtr->position, originalTimeB );
     XLAL_ERROR_REAL8( "minimizeEThincaParameterOverTimeDiff", XLAL_EFUNC );
  }

  /* Set the times back to their correct values */
  XLALSetTimeInPositionVector( params->aPtr->position, originalTimeA );
  XLALSetTimeInPositionVector( params->bPtr->position, originalTimeB );

  return overlap;
}
