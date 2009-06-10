/*
*  Copyright (C) 2007 Anand Sengupta, Craig Robinson
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
 * File Name: LALTrigScanCluster.c
 *
 * Author: Sengupta, Anand. S. and Gupchup, Jayant A.
 *
 * $Id: LALTrigScanCluster.c,v 1.13 2007/10/02 13:48:53 spxcar Exp $
 *
 *---------------------------------------------------------------------------*/

#if 0
<lalVerbatim file="LALTrigScanClusterCV">
Author: Sengupta, Anand. S. and Gupchup, Jayant A.
$Id: LALTrigScanCluster.c,v 1.13 2007/10/02 13:48:53 spxcar Exp $
</lalVerbatim>
#endif

#include <lal/TrigScanEThincaCommon.h>
#include <lal/LALTrigScanCluster.h>

NRCSID (LALTRIGSCANCLUSTERC,
        "$Id: LALTrigScanCluster.c,v 1.13 2007/10/02 13:48:53 spxcar Exp $");

#if 0
<lalLaTeX>
\subsection{Module \texttt{LALTrigScanCluster.c}}

\noindent Functions for trigScan clustering

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALTrigScanClusterCP}
\idx{XLALTrigScanClusterTriggers()}
\idx{XLALTrigScanCreateCluster()}
\idx{XLALTrigScanRemoveStragglers()}
\idx{XLALTrigScanKeepLoudestTrigger()}
\idx{XLALTrigScanReLinkLists()}
\idx{XLALTrigScanDestroyCluster()}
\subsubsection*{Description}

\texttt{XLALTrigScanClusterTriggers()} is the main function for invoking 
trigScan clustering. It takes in the \texttt{SnglInspiralTable} to be clustered,
the method to be applied, the metric scaling factor, and a flag to say whether
to append stragglers (i.e. clusters of only 1 trigger). Upon success, the
return value will be XLAL_SUCCESS, with the texttt{SnglInspiralTable} having
been clustered. At present, the only clustering method implemented is T0T3Tc.

\texttt{XLALTrigScanCreateCluster()} takes in a \texttt{TriggerErrorList}
containing the triggers, their position vectors and ellipsoid matrices. It
creates the cluster by agglomerating triggers by checking for overlap of 
ellipsoids. Upon ellipsoid overlap, the trigger is added to the cluster, and
removed from the unclustered list.

\texttt{XLALTrigScanRemoveStragglers()} takes in a linked list of
\texttt{TrigScanCluster}s. It removes all clusters of 1 element from the list.

\texttt{XLALTrigScanKeepLoudestTrigger()} performs the actual clustering,
freeing all triggers in the cluster except the loudest. Currently, loudest
means the trigger with the highest SNR, but this could be extended to other
statistics.

\texttt{XLALTrigScanReLinkLists()} re-links all the inter-cluster lists after
the clustering has been performed, in preparation for returning the clustered
\texttt{SnglInspiralTable} to the program.

\texttt{XLALTrigScanDestroyCluster()} frees memory associated with a cluster.
It has two modes of operation, specified by the \texttt{TrigScanStatus}. If this
is TRIGSCAN_FAILURE, the SnglInspiralTable will also be freed. If it is
TRIGSCAN_SUCCESS, the SnglInspiralTable will be kept for returning to the 
calling program.
</lalLaTeX>
#endif 

static int CompareErrorListsBySNR( const void *errorA, const void *errorB);

/* <lalVerbatim file="LALTrigScanClusterCP"> */
int XLALTrigScanClusterTriggers( SnglInspiralTable **table,
                                 trigScanType      method,
                                 REAL8             scaleFactor,
                                 INT4              appendStragglers )
/* </lalVerbatim> */
{

  static const char func[] = "XLALTrigScanClusterTriggers";

  SnglInspiralTable *tableHead     = NULL;
  SnglInspiralTable *thisTable     = NULL;
  TriggerErrorList  *errorList     = NULL;
  TrigScanCluster   *clusterHead   = NULL;
  TrigScanCluster   *thisCluster   = NULL;

  /* The maximum time difference associated with an ellipsoid */
  REAL8 tcMax;

#ifndef LAL_NDEBUG
  if ( !table )
  {
    XLAL_ERROR( func, XLAL_EFAULT );
  }

  if ( (UINT4) method >= (UINT4) NUM_TRIGSCAN_TYPE )
  {
    XLAL_ERROR( func, XLAL_EINVAL );
  }
#endif

  tableHead = *table;

  if ( !tableHead )
  {
    XLALPrintWarning( "No triggers to cluster.\n" );
    return XLAL_SUCCESS;
  }

  if ( method == trigScanNone )
  {
    XLALPrintWarning( "No clustering requested.\n" );
    return XLAL_SUCCESS;
  }

  /* TrigScan only currently implemented for tau0/tau3 */
  if ( method != T0T3Tc )
  {
    XLALPrintError( "TrigScan only currently implemented for tau0/tau3!\n" );
    XLAL_ERROR( func, XLAL_EINVAL );
  }

  if ( scaleFactor <= 0.0 )
  {
    XLALPrintError( "TrigScan metric scaling must be > 0: %e given.\n", scaleFactor );
    XLAL_ERROR( func, XLAL_EINVAL );
  }

  /* TrigScan requires triggers to be time-ordered. Make sure this is the case */
  /* and if not, sort the triggers */
  for ( thisTable = tableHead; thisTable->next; thisTable = thisTable->next )
  {
    if ( XLALGPStoINT8( &(thisTable->end_time) ) 
           > XLALGPStoINT8( &(thisTable->next->end_time) ) )
    {
      *table = tableHead = XLALSortSnglInspiral( tableHead, LALCompareSnglInspiralByTime );
      break;
    }
  }

  /* Firstly, create the matrices, etc required for the clustering */
  errorList = XLALCreateTriggerErrorList( tableHead, scaleFactor, &tcMax );
  if ( !errorList )
  {
    XLAL_ERROR( func, XLAL_EFUNC );
  }

  /* Now create the list of clusters. Keep going until errorlist is exhausted */
  while ( errorList )
  {
    TrigScanCluster *newCluster = NULL;
    
    newCluster = XLALTrigScanCreateCluster( &errorList, tcMax );
    /* The next line is to keep track of memory in case of failure */
    if ( errorList )
      *table = errorList->trigger;
    else
      *table = NULL;

    if ( !newCluster )
    {

      thisCluster = clusterHead;

      while ( thisCluster )
      {
         TrigScanCluster *tmpCluster = thisCluster;
         thisCluster = thisCluster->next;
         XLALTrigScanDestroyCluster( tmpCluster, TRIGSCAN_ERROR );
      }
      if ( errorList ) XLALDestroyTriggerErrorList( errorList );
      
      XLAL_ERROR( func, XLAL_EFUNC );
    }
    /* Add the cluster to the list */
    if ( !clusterHead )
    {
      clusterHead = thisCluster = newCluster;
    }
    else
    {
      thisCluster = thisCluster->next = newCluster;
    }    
  }

  /* Remove stragglers if necessary */
  if ( !appendStragglers )
  {
    if ( XLALTrigScanRemoveStragglers( &clusterHead ) == XLAL_FAILURE )
    {
      thisCluster = clusterHead;
 
      while ( thisCluster )
      { 
         TrigScanCluster *tmpCluster = thisCluster;
         thisCluster = thisCluster->next;
         XLALTrigScanDestroyCluster( tmpCluster, TRIGSCAN_ERROR );
      } 
       
      XLAL_ERROR( func, XLAL_EFUNC );
    }

    if ( !clusterHead )
    {
      XLALPrintWarning( "All triggers were stragglers! All have been removed.\n" );
      return XLAL_SUCCESS;
    }
  }

  /* Keep the loudest trigger in each cluster */
  for ( thisCluster = clusterHead; thisCluster; thisCluster = thisCluster->next )
  {
    if ( XLALTrigScanKeepLoudestTrigger( thisCluster ) == XLAL_FAILURE )
    {
      thisCluster = clusterHead;
  
      while ( thisCluster ) 
      {  
         TrigScanCluster *tmpCluster = thisCluster;
         thisCluster = thisCluster->next;
         XLALTrigScanDestroyCluster( tmpCluster, TRIGSCAN_ERROR );
      }  
        
      XLAL_ERROR( func, XLAL_EFUNC );
    }
  }

  /* Re-link the lists ready for returning */
  if ( XLALTrigScanReLinkLists( clusterHead ) == XLAL_FAILURE )
  {
    thisCluster = clusterHead;

    while ( thisCluster )
    {
       TrigScanCluster *tmpCluster = thisCluster;
       thisCluster = thisCluster->next;
       XLALTrigScanDestroyCluster( tmpCluster, TRIGSCAN_ERROR );
    }
       
    XLAL_ERROR( func, XLAL_EFUNC );
  }

  *table = clusterHead->element->trigger;

  /* Since trigScan can have multiple clusters at similar times */
  /* We sort the list to ensure time-ordering */
  *table = XLALSortSnglInspiral( *table, LALCompareSnglInspiralByTime );
  XLALPrintInfo( "Returning %d clustered triggers.\n", XLALCountSnglInspiral( *table ) );

  /* Free the memory */
  thisCluster = clusterHead;
  while ( thisCluster )
  {
     TrigScanCluster *tmpCluster = thisCluster;
     thisCluster = thisCluster->next;
     XLALTrigScanDestroyCluster( tmpCluster, TRIGSCAN_SUCCESS );
  }

  return XLAL_SUCCESS;
}

/* <lalVerbatim file="LALTrigScanClusterCP"> */
TrigScanCluster * XLALTrigScanCreateCluster( TriggerErrorList **errorListHead,
                                             REAL8            tcMax )
/* </lalVerbatim> */
{
  static const char func[] = "XLALTrigScanCreateCluster";

  TrigScanCluster *cluster = NULL;

  /* Pointers to the main trigger error list */
  TriggerErrorList *thisErrorList     = NULL;
  TriggerErrorList *previousErrorList = NULL;

  /* Pointers to the trigger error list within the cluster */
  TriggerErrorList *thisClusterList   = NULL;
  TriggerErrorList *endClusterList    = NULL;

  INT8              maxTimeDiff;

  /* Stuff for checking ellipsoid overlap */
  fContactWorkSpace *workSpace        = NULL;
  REAL8             fContactValue;

#ifndef LAL_NDEBUG
  if ( !errorListHead )
  {
    XLAL_ERROR_NULL( func, XLAL_EFAULT );
  }

  if ( !(*errorListHead) )
  {
    XLAL_ERROR_NULL( func, XLAL_EFAULT );
  }

  if ( tcMax <= 0 )
  {
    XLAL_ERROR_NULL( func, XLAL_EINVAL );
  }
#endif

  /* Allocate memory for the cluster */
  cluster = LALCalloc( 1, sizeof( TrigScanCluster ) );
  if ( !cluster )
  {
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  }

  /* Create the workspace for checking ellipsoid overlap */
  workSpace = XLALInitFContactWorkSpace( 3, NULL, NULL, gsl_min_fminimizer_brent, 1.0e-2 );
  if ( !workSpace )
  {
    LALFree( cluster );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }

  /* Set the first trigger in the list to be part of the cluster */
  cluster->element                = *errorListHead;
  endClusterList                  = cluster->element;
  cluster->nelements              = 1;
  *errorListHead                  = cluster->element->next;
  cluster->element->next          = NULL;
  cluster->element->trigger->next = NULL;

  maxTimeDiff = (INT8)( (2.0 * tcMax + 1.0e-5) * 1.0e9 );

  thisClusterList = cluster->element;
  /* Now we go through the agglomeration procedure */
  while ( thisClusterList )
  {
    /* Timing info of the trigger in the cluster */
    INT8  endTimeA;
    REAL8 originalTimeA;

    /* Timing info of the trigger we wish to compare */
    INT8  endTimeB;
    REAL8 originalTimeB;

    endTimeA = XLALGPStoINT8( &(thisClusterList->trigger->end_time) );
    XLAL_CALLGSL( originalTimeA = gsl_vector_get( thisClusterList->position, 0 ) );

    /* Reset the time to avoid precision problems */
    XLALSetTimeInPositionVector( thisClusterList->position, 0 ); 

    /* Loop through the list of triggers */
    thisErrorList = *errorListHead;
    previousErrorList = NULL;
    while ( thisErrorList )
    {

      endTimeB = XLALGPStoINT8( &(thisErrorList->trigger->end_time) );
      XLAL_CALLGSL( originalTimeB = gsl_vector_get( thisErrorList->position, 0 ) );

      /* If the triggers are more than twice the max time error apart, no need to proceed */
      if ( endTimeB - endTimeA > maxTimeDiff )
      {
        break;
      }

      XLALSetTimeInPositionVector( thisErrorList->position, 
                (REAL8) ( ( endTimeB - endTimeA ) * 1.0e-9 ) );


      /* check for the intersection of the ellipsoids */
      workSpace->invQ1 = thisClusterList->err_matrix;
      workSpace->invQ2 = thisErrorList->err_matrix;
      fContactValue = XLALCheckOverlapOfEllipsoids( thisClusterList->position,
                   thisErrorList->position, workSpace );
      if (XLAL_IS_REAL8_FAIL_NAN(fContactValue))
      {
        /* The triggers in the cluster have been removed from the main list */
        /* so they must be freed here */
        thisClusterList = cluster->element;
        while ( thisClusterList )
        {
          TriggerErrorList *tmpClusterList = thisClusterList;
          thisClusterList = thisClusterList->next;
          XLALFreeSnglInspiral( &(tmpClusterList->trigger) );
          XLAL_CALLGSL( gsl_matrix_free( tmpClusterList->err_matrix ) );
          XLAL_CALLGSL( gsl_vector_free( tmpClusterList->position ) );
          LALFree( tmpClusterList );
        }
        XLALFreeFContactWorkSpace( workSpace );
        XLAL_ERROR_NULL( func, XLAL_EFUNC );
      }
      /* Reset the time to its original value */
      XLALSetTimeInPositionVector( thisErrorList->position, originalTimeB );

      /* test whether we have coincidence */
      if ( fContactValue <= 1.0 )
      {
        /* Add the trigger to the cluster, and pull it off the main list */
        if ( previousErrorList )
        {
          if ( thisErrorList->next )
          {
            previousErrorList->trigger->next = thisErrorList->next->trigger;
          }
          else
          {
            previousErrorList->trigger->next = NULL;
          }
          previousErrorList->next = thisErrorList->next;
        }
        else
        {
          *errorListHead = thisErrorList->next;
        }
        endClusterList->trigger->next = thisErrorList->trigger;
        endClusterList = endClusterList->next = thisErrorList;
        thisErrorList = thisErrorList->next;
        endClusterList->next = NULL;
        endClusterList->trigger->next = NULL;
        cluster->nelements++;
      }
      else
      {
        /* No coincidence, so we go on */
        previousErrorList = thisErrorList;
        thisErrorList = thisErrorList->next;
      }
    }
    /* Reset the time in the cluster list trigger */
    XLALSetTimeInPositionVector( thisClusterList->position, originalTimeA );
    thisClusterList = thisClusterList->next;
  }

  XLALFreeFContactWorkSpace( workSpace );

  /* We have now clustered the triggers - return the result */
  XLALPrintInfo( "Returning a cluster containing %d triggers.\n", cluster->nelements );
  return cluster;
}

/* <lalVerbatim file="LALTrigScanClusterCP"> */      
int XLALTrigScanRemoveStragglers( TrigScanCluster **clusters )
/* </lalVerbatim> */
{

  static const char func[] = "XLALTrigScanRemoveStragglers";

  TrigScanCluster *previous    = NULL; /* Keeping track of the previous element */
  TrigScanCluster *thisCluster = NULL;

#ifndef LAL_NDEBUG
  if ( !clusters )
  {
    XLAL_ERROR( func, XLAL_EFAULT );
  }
#endif

  /* Loop through the list and remove all clusters containing 1 trigger */
  while ( thisCluster )
  {
    if ( thisCluster->nelements == 1 )
    {
      TrigScanCluster *tmpCluster = thisCluster;

      thisCluster = thisCluster->next;

      if ( !previous )
      {
        *clusters = tmpCluster->next;
      }
      else
      {
        previous->next = tmpCluster->next;
      }
      XLALFreeSnglInspiral( &(tmpCluster->element->trigger) );
      XLAL_CALLGSL( gsl_matrix_free( tmpCluster->element->err_matrix ) );
      XLAL_CALLGSL( gsl_vector_free( tmpCluster->element->position ) );
      LALFree( tmpCluster );
    }
    else
    {
      previous = thisCluster;
      thisCluster = thisCluster->next;
    }
  }

  return XLAL_SUCCESS;
}

/* <lalVerbatim file="LALTrigScanClusterCP"> */
int XLALTrigScanKeepLoudestTrigger( TrigScanCluster *cluster )
{
  static const char func[] = "TrigScanKeepLoudestTrigger";

  TriggerErrorList *triggerToKeep;
  TriggerErrorList *thisTrigger;

#ifndef LAL_NDEBUG
  if ( !cluster )
  {
    XLAL_ERROR( func, XLAL_EFAULT );
  }

  if ( cluster->nelements < 1 )
  {
    XLALPrintError( "Invalid number of triggers in cluster: %d\n", cluster->nelements );
    XLAL_ERROR( func, XLAL_EINVAL );
  }
#endif

  if ( cluster->nelements == 1 )
  {
    /* No need to do anything */
    return XLAL_SUCCESS;
  }

  triggerToKeep = cluster->element;

  /* Find the loudest trigger */
  for ( thisTrigger = cluster->element->next; thisTrigger; thisTrigger = thisTrigger->next )
  {
    if ( thisTrigger->trigger->snr > triggerToKeep->trigger->snr )
    {
      triggerToKeep = thisTrigger;
    }
  }

  /* Keep only the loudest trigger */
  thisTrigger = cluster->element;
  while ( thisTrigger )
  {
    TriggerErrorList *tmpTrigger = thisTrigger;
    thisTrigger = thisTrigger->next;

    if ( tmpTrigger != triggerToKeep )
    {
      XLALFreeSnglInspiral( &(tmpTrigger->trigger ) );
      XLAL_CALLGSL( gsl_matrix_free( tmpTrigger->err_matrix ) );
      XLAL_CALLGSL( gsl_vector_free( tmpTrigger->position ) );
      LALFree( tmpTrigger );
    }
  }

  cluster->element = triggerToKeep;
  cluster->element->trigger->next = NULL;
  cluster->element->next = NULL;

  return XLAL_SUCCESS;
}

/* <lalVerbatim file="LALTrigScanClusterCP"> */
int XLALTrigScanReLinkLists( TrigScanCluster *clusterHead )
/* </lalVerbatim> */
{

  static const char func[] = "XLALTrigScanReLinkLists";
  
  TrigScanCluster *thisCluster;

  if ( ! clusterHead )
    XLAL_ERROR( func, XLAL_EFAULT );

  for ( thisCluster = clusterHead; thisCluster->next; thisCluster = thisCluster->next )
  {
    thisCluster->element->trigger->next = thisCluster->next->element->trigger;
  }
  /* Set the last one to null */
  thisCluster->element->trigger->next = NULL;

  return XLAL_SUCCESS;
}

/* <lalVerbatim file="LALTrigScanClusterCP"> */
void XLALTrigScanDestroyCluster( TrigScanCluster *cluster,
                                TrigScanStatus   status
                              )
/* </lalVerbatim> */
{

  static const char func[] = "XLALTrigScanDestroyCluster";

  TriggerErrorList *thisList;

#ifndef LAL_NDEBUG
  if ( !cluster )
    XLAL_ERROR_VOID( func, XLAL_EFAULT );

  if ( (UINT4) status >= (UINT4) TRIGSCAN_NUM_STATUS )
    XLAL_ERROR_VOID( func, XLAL_EINVAL );
#endif

  /* If something has failed, we need to free the SnglInspirals */
  if ( status == TRIGSCAN_ERROR )
  {
    for ( thisList = cluster->element; thisList; thisList = thisList->next )
    {
      XLALFreeSnglInspiral( &(thisList->trigger ) );
    }
  }

  XLALDestroyTriggerErrorList( cluster->element );
  LALFree( cluster );

  return;
}
