/*
*  Copyright (C) 2007 Lisa M. Goggin
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

/*-----------------------------------------------------------------------
 *
 * File Name: CoincRingdownUtils.c
 *
 * Author: Goggin, L. M. based on CoincInspiralUtils.c by  Brady, P. R., Brown, D. A.,
 * and Fairhurst, S
 *
 *
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GenerateRing.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/XLALError.h>

/**
 * \defgroup CoincRingdownUtils_c Module CoincRingdownUtils.c
 * \ingroup pkg_CBC_NEW
 * \author Fairhurst, S.
 * \brief Blah.
 */
/*@{*/

/**
 * Takes in a linked list of single inspiral
 * tables and returns a list of two instrument coincidences.  The coincidence
 * requirements are given by the \c accuracyParams.  When single inspirals
 * from two different instruments are found to be coincident, the code creates a
 * new \c coincInspiralTable and uses <tt>LALAddSnglInspiralToCoinc()</tt>
 * to add the single inspirals to the coinc.  The function returns
 * \c qcoincOutput which is a pointer to the head of a linked list of
 * \c CoincInspiralTables.
 */
int
XLALCoincRingdownIfosDiscard(
    CoincRingdownTable **coincHead,
    char                *ifos
    )

{
  CoincRingdownTable    *prevCoinc = NULL;
  CoincRingdownTable    *thisCoinc = NULL;
  int                    numCoinc = 0;

  thisCoinc = *coincHead;
  *coincHead = NULL;

  while ( thisCoinc )   {
    CoincRingdownTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    if ( XLALCoincRingdownIfos( tmpCoinc, ifos ) )
    {
      /* ifos match so discard tmpCoinc */
      XLALFreeCoincRingdown( &tmpCoinc );
    }
    else
    {
      /* keep tmpCoinc */
      if ( ! *coincHead  )
      {
        *coincHead = tmpCoinc;
      }
      else
      {
        prevCoinc->next = tmpCoinc;
      }
      tmpCoinc->next = NULL;
      prevCoinc = tmpCoinc;
      ++numCoinc;
    }
  }

  return( numCoinc );
}



void
LALCreateTwoIFORingdownCoincList(
    LALStatus                  *status,
    CoincRingdownTable        **coincOutput,
    SnglRingdownTable          *snglInput,
    RingdownAccuracyList       *accuracyParams
    )

{
  SnglRingdownTable            *currentTrigger[2];
  INT8                          currentTriggerNS[2];
  CoincRingdownTable           *coincHead = NULL;
  CoincRingdownTable           *thisCoinc = NULL;
  INT4                          numEvents = 0;
  INT4                          ifoNumber;
  INT8                          maxTimeDiff = 0;
  REAL8                         ds2;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( coincOutput, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( ! *coincOutput, status,
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );

  memset( currentTriggerNS, 0, 2 * sizeof(INT8) );
  memset( currentTrigger, 0, 2 * sizeof(SnglRingdownTable *) );


  /* calculate the maximum time delay
   * set it equal to 2 * worst IFO timing accuracy plus
   * light travel time for earth's diameter
   * (detectors can't be further apart than this) */

  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    maxTimeDiff = (maxTimeDiff > accuracyParams->ifoAccuracy[ifoNumber].dt) ?
      maxTimeDiff : accuracyParams->ifoAccuracy[ifoNumber].dt;
  }
  maxTimeDiff *= 2;
  maxTimeDiff += (INT8) ( 1e9 * 2 * LAL_REARTH_SI / LAL_C_SI );

  for ( currentTrigger[0] = snglInput; currentTrigger[0]->next;
      currentTrigger[0] = currentTrigger[0]->next)
  {

    /* calculate the time of the trigger */
    currentTriggerNS[0] = XLALGPSToINT8NS( &(currentTrigger[0]->start_time) );

    /* set next trigger for comparison */
    currentTrigger[1] = currentTrigger[0]->next;
    currentTriggerNS[1] = XLALGPSToINT8NS( &(currentTrigger[1]->start_time) );

    while ( (currentTriggerNS[1] - currentTriggerNS[0]) < maxTimeDiff )
    {
      /* check that triggers pass coincidence test */
      ds2 = XLALCompareRingdowns( currentTrigger[0],
          currentTrigger[1], accuracyParams );

      /* test whether we have coincidence */
      if ( accuracyParams->match )
      {
        LALInfo( status, "Found double coincident trigger,");
        /* create a 2 IFO coinc and store */
        if ( ! coincHead  )
        {
          coincHead = thisCoinc = (CoincRingdownTable *)
            LALCalloc( 1, sizeof(CoincRingdownTable) );
        }
        else
        {
          thisCoinc = thisCoinc->next = (CoincRingdownTable *)
            LALCalloc( 1, sizeof(CoincRingdownTable) );
        }

        /* Store the appropriate ds2 columns in the coinc table */
        if ( (strcmp(currentTrigger[0]->ifo,"H1")==0 && strcmp(currentTrigger[1]->ifo,"H2")==0)
             ||(strcmp(currentTrigger[0]->ifo,"H2")==0 && strcmp(currentTrigger[1]->ifo,"H1")==0) )
        {
          thisCoinc->ds2_H1H2 = ds2;
        }

        else if ( (strcmp(currentTrigger[0]->ifo,"H1")==0 && strcmp(currentTrigger[1]->ifo,"L1")==0)
             ||(strcmp(currentTrigger[0]->ifo,"L1")==0 && strcmp(currentTrigger[1]->ifo,"H1")==0) )
        {
          thisCoinc->ds2_H1L1 = ds2;
        }

        else if ( (strcmp(currentTrigger[0]->ifo,"H1")==0 && strcmp(currentTrigger[1]->ifo,"V1")==0)
             ||(strcmp(currentTrigger[0]->ifo,"V1")==0 && strcmp(currentTrigger[1]->ifo,"H1")==0) )
        {
          thisCoinc->ds2_H1V1 = ds2;
        }

        else if ( (strcmp(currentTrigger[0]->ifo,"H2")==0 && strcmp(currentTrigger[1]->ifo,"L1")==0)
             ||(strcmp(currentTrigger[0]->ifo,"L1")==0 && strcmp(currentTrigger[1]->ifo,"H2")==0) )
        {
          thisCoinc->ds2_H2L1 = ds2;
        }

        else if ( (strcmp(currentTrigger[0]->ifo,"H2")==0 && strcmp(currentTrigger[1]->ifo,"V1")==0)
             ||(strcmp(currentTrigger[0]->ifo,"V1")==0 && strcmp(currentTrigger[1]->ifo,"H2")==0) )
        {
          thisCoinc->ds2_H2V1 = ds2;
        }

        else if ( (strcmp(currentTrigger[0]->ifo,"L1")==0 && strcmp(currentTrigger[1]->ifo,"V1")==0)
             ||(strcmp(currentTrigger[0]->ifo,"V1")==0 && strcmp(currentTrigger[1]->ifo,"L1")==0) )
        {
          thisCoinc->ds2_L1V1 = ds2;
        }

        else
        {
          LALInfo( status, "Unknown pair of ifo's");
        }

        /* Add the two triggers to the coinc */
        LALAddSnglRingdownToCoinc( status->statusPtr, &thisCoinc,
            currentTrigger[0] );
        LALAddSnglRingdownToCoinc( status->statusPtr, &thisCoinc,
            currentTrigger[1] );

        ++numEvents;

      }

      /* scroll on to the next sngl ringdown */

      if ( (currentTrigger[1] = currentTrigger[1]->next) )
      {
        currentTriggerNS[1] = XLALGPSToINT8NS( &(currentTrigger[1]->start_time) );
      }
      else
      {
        LALInfo(status, "Second trigger has reached end of list");
        break;
      }
    }
  }

  *coincOutput = coincHead;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/**
 * Takes linked list of
 * \c CoincInspiralTables, assumed to contain (N-1) ifo coincidences and
 * creates all N ifo coincidences.  Both the input and output list of
 * \c CoincInspiralTables are passed as \c coincHead.
 */
void
LALCreateNIFORingdownCoincList(
    LALStatus                  *status,
    CoincRingdownTable        **coincHead,
    RingdownAccuracyList       *accuracyParams,
    INT4                        N
    )

{
  REAL8                         ds2        = 0;
  INT4                          numEvents  = 0;
  InterferometerNumber          ifoNumber  = LAL_UNKNOWN_IFO;
  InterferometerNumber          ifoNum     = LAL_UNKNOWN_IFO;
  InterferometerNumber          firstEntry = LAL_UNKNOWN_IFO;
  CoincRingdownTable           *thisCoinc     = NULL;
  CoincRingdownTable           *lastCoinc     = NULL;
  CoincRingdownTable           *otherCoinc    = NULL;
  CoincRingdownTable           *nIfoCoincHead = NULL;
  CoincRingdownTable           *thisNIfoCoinc = NULL;
  EventIDColumn                *eventIDHead   = NULL;
  EventIDColumn                *thisID        = NULL;


  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* loop over all the coincidences in the list */
  for( thisCoinc = *coincHead; thisCoinc; thisCoinc = thisCoinc->next)
  {
    lastCoinc = thisCoinc;

    /* check that this is an (N-1) coinc */
    if ( thisCoinc->numIfos == N - 1 )
    {
      /* look up the first single ringdown */
      for ( firstEntry = (InterferometerNumber) 0; firstEntry < LAL_NUM_IFO; firstEntry++)
      {
        if ( thisCoinc->snglRingdown[firstEntry] )
        {
          LALInfo( status, "Found the first entry in the coinc" );
          break;
        }
      }

      /* get the list of event IDs for this first entry */
      eventIDHead = thisCoinc->snglRingdown[firstEntry]->event_id;

      /* loop over the (N-1) ifo coincs that first entry is a member of
       * and try to find an N ifo coinc */
      for( thisID = eventIDHead; thisID; thisID = thisID->next )
      {
        otherCoinc = thisID->coincRingdownTable;

        if( otherCoinc->numIfos == N - 1 )
        {
          /* loop over all singles which are alphabetically before the
           * first one in thisCoinc */
          for( ifoNumber = (InterferometerNumber) 0; ifoNumber < firstEntry; ifoNumber++ )
          {
            /* test whether we have an N ifo coincidence */
            accuracyParams->match = 0;

            if ( otherCoinc->snglRingdown[ifoNumber] )
            {
              ds2 = XLALSnglRingdownCoincTest( thisCoinc, 
                otherCoinc->snglRingdown[ifoNumber], accuracyParams );
            }

            if ( accuracyParams->match )
            {
              LALInfo( status, "We have found an N ifo coinc, storing");
              ++numEvents;

              /* create a N IFO coinc and store */
              if ( ! nIfoCoincHead  )
              {
                nIfoCoincHead = thisNIfoCoinc = (CoincRingdownTable *)
                  LALCalloc( 1, sizeof(CoincRingdownTable) );
              }
              else
              {
                thisNIfoCoinc = thisNIfoCoinc->next = (CoincRingdownTable *)
                  LALCalloc( 1, sizeof(CoincRingdownTable) );
              }

              /* Store the appropriate ds2 columns in the coinc table */
              if ( (strcmp(otherCoinc->snglRingdown[ifoNumber]->ifo,"H1")==0 && strcmp(thisCoinc->snglRingdown[firstEntry]->ifo,"H2")==0)
                ||(strcmp(otherCoinc->snglRingdown[ifoNumber]->ifo,"H2")==0 && strcmp(thisCoinc->snglRingdown[firstEntry]->ifo,"H1")==0) )
              {
                thisNIfoCoinc->ds2_H1H2 = ds2;
              }

              else if ( (strcmp(otherCoinc->snglRingdown[ifoNumber]->ifo,"H1")==0 && strcmp(thisCoinc->snglRingdown[firstEntry]->ifo,"L1")==0)
                ||(strcmp(otherCoinc->snglRingdown[ifoNumber]->ifo,"L1")==0 && strcmp(thisCoinc->snglRingdown[firstEntry]->ifo,"H1")==0) )
              {
                thisNIfoCoinc->ds2_H1L1 = ds2;
              }

              else if ( (strcmp(otherCoinc->snglRingdown[ifoNumber]->ifo,"H1")==0 && strcmp(thisCoinc->snglRingdown[firstEntry]->ifo,"V1")==0)
                ||(strcmp(otherCoinc->snglRingdown[ifoNumber]->ifo,"V1")==0 && strcmp(thisCoinc->snglRingdown[firstEntry]->ifo,"H1")==0) )
              {
                thisNIfoCoinc->ds2_H1V1 = ds2;
              }

              else if ( (strcmp(otherCoinc->snglRingdown[ifoNumber]->ifo,"H2")==0 && strcmp(thisCoinc->snglRingdown[firstEntry]->ifo,"L1")==0)
                ||(strcmp(otherCoinc->snglRingdown[ifoNumber]->ifo,"L1")==0 && strcmp(thisCoinc->snglRingdown[firstEntry]->ifo,"H2")==0) )
              {
                thisNIfoCoinc->ds2_H2L1 = ds2;
              }

              else if ( (strcmp(otherCoinc->snglRingdown[ifoNumber]->ifo,"H2")==0 && strcmp(thisCoinc->snglRingdown[firstEntry]->ifo,"V1")==0)
                ||(strcmp(otherCoinc->snglRingdown[ifoNumber]->ifo,"V1")==0 && strcmp(thisCoinc->snglRingdown[firstEntry]->ifo,"H2")==0) )
              {
                thisNIfoCoinc->ds2_H2V1 = ds2;
              }

              else if ( (strcmp(otherCoinc->snglRingdown[ifoNumber]->ifo,"L1")==0 && strcmp(thisCoinc->snglRingdown[firstEntry]->ifo,"V1")==0)
                ||(strcmp(otherCoinc->snglRingdown[ifoNumber]->ifo,"V1")==0 && strcmp(thisCoinc->snglRingdown[firstEntry]->ifo,"L1")==0) )
              {
                thisNIfoCoinc->ds2_L1V1 = ds2;
              }

              else
              {
                LALInfo( status, "Unknown pair of ifo's");
              }

              /* add the single to the new N coinc */
              LALAddSnglRingdownToCoinc( status->statusPtr, &thisNIfoCoinc,
                  otherCoinc->snglRingdown[ifoNumber] );

              /* add the triggers from the (N-1) coinc to the new N coinc */
              for( ifoNum = (InterferometerNumber) 0; ifoNum < LAL_NUM_IFO; ifoNum++ )
              {
                if( thisCoinc->snglRingdown[ifoNum] )
                {
                  LALAddSnglRingdownToCoinc( status->statusPtr, &thisNIfoCoinc,
                      thisCoinc->snglRingdown[ifoNum] );
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

/**
 * Will remove any lower order coincidences
 * if they are contained in a higher order coincidence.  For example, if an H1-L1
 * double coincident trigger is also part of an H1-H2-L1 triple coincident
 * trigger, the double coincident trigger will be removed.  The head of the list
 * of coincident triggers is passed and returned as \c coincHead.
 */
void
LALRemoveRepeatedRingdownCoincs(
    LALStatus                  *status,
    CoincRingdownTable        **coincHead
    )
{
  INT4                          removeThisCoinc  = 0;
  InterferometerNumber          ifoNumber  = LAL_UNKNOWN_IFO;
  InterferometerNumber          firstEntry = LAL_UNKNOWN_IFO;
  CoincRingdownTable           *thisCoinc     = NULL;
  CoincRingdownTable           *prevCoinc     = NULL;
  CoincRingdownTable           *otherCoinc    = NULL;
  EventIDColumn                *eventIDHead   = NULL;
  EventIDColumn                *thisID        = NULL;


  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* loop over all the coincidences in the list */
  thisCoinc = *coincHead;

  while( thisCoinc )
  {
    /* look up the first single ringdown */
    for ( firstEntry = (InterferometerNumber) 0; firstEntry < LAL_NUM_IFO; firstEntry++)
    {
      if ( thisCoinc->snglRingdown[firstEntry] )
      {
        LALInfo( status, "Found the first entry in the coinc" );
        break;
      }
    }

    /* get the list of event IDs for this first entry */
    eventIDHead = thisCoinc->snglRingdown[firstEntry]->event_id;

    /* loop over the coincs that firstEntry is a member of and see if
     * thisCoinc is a subset of a higher order coinc */

    removeThisCoinc = 0;

    for( thisID = eventIDHead; thisID; thisID = thisID->next )
    {
      otherCoinc = thisID->coincRingdownTable;

      if( otherCoinc->numIfos >= thisCoinc->numIfos &&
          otherCoinc != thisCoinc )
      {
        /* we have a higher (or equal) coinc, thisCoinc could be a subset
         * test whether all sngls in thisCoinc are also in otherCoinc */

        for( ifoNumber = (InterferometerNumber) ( firstEntry + 1 ); ifoNumber < LAL_NUM_IFO;
            ifoNumber++ )
        {
          if ( thisCoinc->snglRingdown[ifoNumber] &&
              !(thisCoinc->snglRingdown[ifoNumber] ==
                otherCoinc->snglRingdown[ifoNumber]) )
          {
            LALInfo( status, "No Match");
            break;
          }
        }

        if ( ifoNumber == LAL_NUM_IFO )
        {
          LALInfo( status, "Removing lower order coinc");
          removeThisCoinc = 1;
          break;
        }
      }
    }

    if ( removeThisCoinc )
    {
      if ( !prevCoinc )
      {
        *coincHead = thisCoinc->next;
        LALFreeCoincRingdown( status->statusPtr, &thisCoinc );
        thisCoinc = *coincHead;
      }
      else
      {
        prevCoinc->next = thisCoinc->next;
        LALFreeCoincRingdown( status->statusPtr, &thisCoinc );
        thisCoinc = prevCoinc->next;
      }
    }
    else
    {
      prevCoinc = thisCoinc;
      thisCoinc = thisCoinc->next;
    }
  }


  DETATCHSTATUSPTR (status);
  RETURN (status);

}

void
LALFreeCoincRingdown(
    LALStatus                  *status,
    CoincRingdownTable        **coincPtr
    )
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  XLALFreeCoincRingdown( coincPtr );


  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/**
 * Free the
 * memory associated to the \c CoincInspiralTable pointed to by
 * \c coincPtr.  This entails freeing the \c CoincInspiralTable as
 * well as any \c eventIds which point to the coinc.
 */
void
XLALFreeCoincRingdown(
    CoincRingdownTable        **coincPtr
    )
{
  InterferometerNumber          ifoNumber  = LAL_UNKNOWN_IFO;
  EventIDColumn                *prevID     = NULL;
  EventIDColumn                *thisID     = NULL;
  SnglRingdownTable            *thisSngl   = NULL;

  for ( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( (thisSngl = (*coincPtr)->snglRingdown[ifoNumber]) )
    {
      /* loop over the list of eventID's until we get to the one that
       * points to thisCoinc */
      thisID = thisSngl->event_id;
      prevID = NULL;

      while ( thisID )
      {
        /* test if thisID points to our coinc */
        if ( thisID->coincRingdownTable == *coincPtr )
        {
          if ( !prevID )
          {
            thisSngl->event_id = thisID->next;
          }
          else
          {
            prevID->next = thisID->next;
          }
          LALFree(thisID);
          break;
        }
        else
        {
          prevID = thisID;
          thisID = thisID->next;
        }
      }
    }
  }
  LALFree(*coincPtr);
}

/**
 * Adds a pointer to a single inspiral table
 * to a coinc inspiral table.  Upon entry, if \c coincPtr points to a
 * \c NULL coinc inspiral table, the table is created before a pointer to
 * the single inspiral table is added.  Additionally, an \c eventId table is
 * created for the single inspiral table.  This points to both the single and
 * coinc inspirals.  If an \c eventId already exists for the single
 * inspiral, another eventId table is added to the linked list.  The linked list
 * of \c eventIds associated to a single inspiral table allow us to easily
 * determine which coincident events each single is a part of.
 */
void
LALAddSnglRingdownToCoinc(
    LALStatus                  *status,
    CoincRingdownTable        **coincPtr,
    SnglRingdownTable          *snglRingdown
    )

{
  /*
  CoincRingdownTable  *coincRingdown = NULL;
  EventIDColumn       *eventId = NULL;
  */

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( coincPtr, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( snglRingdown, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  *coincPtr = XLALAddSnglRingdownToCoinc(*coincPtr, snglRingdown);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}



CoincRingdownTable *
XLALAddSnglRingdownToCoinc(
    CoincRingdownTable         *coincRingdown,
    SnglRingdownTable          *snglRingdown
    )

{
  EventIDColumn     *eventId = NULL;

  /* allocate memory for new coinc if it doesn't exist */
  if (! coincRingdown )
  {
    coincRingdown = (CoincRingdownTable *)
      LALCalloc( 1, sizeof(CoincRingdownTable) );
    if ( !coincRingdown )
    {
      LALFree( coincRingdown );
      XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
  }

  switch ( (snglRingdown->ifo)[0] )
  {
    case 'H':
      if ( !strcmp( snglRingdown->ifo, "H1" ) )
      {
        coincRingdown->snglRingdown[LAL_IFO_H1] = snglRingdown;
      }
      else if (!strcmp( snglRingdown->ifo, "H2" ) )
      {
        coincRingdown->snglRingdown[LAL_IFO_H2] = snglRingdown;
      }
      else
      {
        /* Invalid Hanford Detector */
        XLALPrintError( "Invalid ifo in input snglInspiral" );
        XLAL_ERROR_NULL(XLAL_EIO);
      }
      break;

    case 'L':
      coincRingdown->snglRingdown[LAL_IFO_L1] = snglRingdown;
      break;

    case 'V':
      coincRingdown->snglRingdown[LAL_IFO_V1] = snglRingdown;
      break;

    default:
      /* Invalid Detector Site */
      XLALPrintError( "Invalid ifo in input snglInspiral" );
      XLAL_ERROR_NULL(XLAL_EIO);
  }

  ++(coincRingdown->numIfos);

  /* create an eventId for the single, populate it with the single and coinc */
  if ( ! snglRingdown->event_id )
  {
    eventId = (EventIDColumn *) LALCalloc( 1, sizeof(EventIDColumn) );
    if ( !eventId )
    {
      LALFree(eventId);
      XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
    snglRingdown->event_id = eventId;
  }
  else
  {
     for( eventId = snglRingdown->event_id; eventId->next;
         eventId = eventId->next);
     eventId = eventId->next = (EventIDColumn *)
         LALCalloc( 1, sizeof(EventIDColumn) );
    if ( !eventId )
    {
      LALFree(eventId);
      XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
  }
  eventId->snglRingdownTable = snglRingdown;
  eventId->coincRingdownTable = coincRingdown;

  return coincRingdown;
}

/**
 * Tests for coincidence between a single
 * inspiral and a coinc inspiral.  It works by testing for coincidence between
 * each non-null entry in the coinc inspiral and the single.  This is done using
 * <tt>LALCompareSnglInspiral()</tt>.  If all members of the coinc are found to be
 * coincident with the single, the <tt>accuracyParams.match</tt> is set to 1,
 * otherwise to 0.
 */
REAL8
XLALSnglRingdownCoincTest(
    CoincRingdownTable         *coincRingdown,
    SnglRingdownTable          *snglRingdown,
    RingdownAccuracyList       *accuracyParams
    )

{
  SnglRingdownTable    *thisCoincEntry;
  INT4                  match = 1;
  INT4                  ifoNumber = 0;
  REAL8                 ds2 = 0;


  /* Loop over sngl_ringdowns contained in coinc_ringdown */
  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    thisCoincEntry = coincRingdown->snglRingdown[ifoNumber];

    if ( thisCoincEntry )
    {
      /* snglRingdown entry exists for this IFO, perform coincidence test */
      if ( ifoNumber == XLALIFONumber(snglRingdown->ifo) )
      {
        XLALPrintInfo( "We already have a coinc from this IFO" );
        accuracyParams->match = 0;
      }

      else
      {
        ds2 = XLALCompareRingdowns ( snglRingdown,
            thisCoincEntry, accuracyParams );
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


  return ds2;
}


/**
 * Extracts the information from a
 * linked list of \c coincInspiralTables and returns it as a linked list of
 * \c snglInspiralTables.  Thus, the output \c snglPtr is a pointer to
 * a linked list of single inspiral tables.  That list contains only single
 * inspirals which are found in coincidence.  In order to preserve the coincidence
 * information, we assign to each coincident event an integer value.  This is
 * stored in the <tt>UINT8 id</tt> field of the \c eventIDColumn of each
 * single inspiral which forms part of the coincidence.  The \c id is set
 * equal to \f$10^{9} \times\f$ \c gpsStartTime \f$+ 10^{5} \times\f$
 * \c slideNum \f$+\f$ event number. We do not assign multiple \c id
 * values to a given single inspiral table, but instead make multiple copies of
 * the table, each with a unique \c id.
 */
void
LALExtractSnglRingdownFromCoinc(
    LALStatus                  *status,
    SnglRingdownTable         **snglPtr,
    CoincRingdownTable         *coincRingdown,
    LIGOTimeGPS                *gpsStartTime,
    INT4                        slideNum
    )

{
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  *snglPtr = XLALExtractSnglRingdownFromCoinc( coincRingdown, gpsStartTime,
      slideNum );


  DETATCHSTATUSPTR (status);
  RETURN (status);
}


SnglRingdownTable *
XLALExtractSnglRingdownFromCoinc(
    CoincRingdownTable         *coincRingdown,
    LIGOTimeGPS                *gpsStartTime,
    INT4                        slideNum
    )

{
  SnglRingdownTable  *snglHead = NULL;
  SnglRingdownTable  *thisSngl = NULL;
  SnglRingdownTable  *thisCoincEntry = NULL;
  CoincRingdownTable *thisCoinc = NULL;
  EventIDColumn      *eventId = NULL;
  UINT4               eventNum = 1;
  INT4                j;

  if ( !coincRingdown )
  {
    XLALPrintInfo(
        "XLALExtractSnglRingdownFromCoinc: Empty coincRingdown passed as input"
        );
    return( NULL );
  }

    /* loop over the linked list of coinc ringdown */
  for( thisCoinc = coincRingdown; thisCoinc; thisCoinc = thisCoinc->next,
      ++eventNum)
  {
    /* loop over the interferometers */
    for ( j = 0; j < LAL_NUM_IFO; j++)
    {
      thisCoincEntry = thisCoinc->snglRingdown[j];

      if ( thisCoincEntry )
      {
        /* allocate memory for a new sngl ringdown */
        if ( !snglHead )
        {
          thisSngl = snglHead = (SnglRingdownTable *)
            LALCalloc( 1, sizeof(SnglRingdownTable) );
        }
        else
        {
          thisSngl = thisSngl->next = (SnglRingdownTable *)
            LALCalloc( 1, sizeof(SnglRingdownTable) );
        }

        /* copy thisCoincEntry into our list */
        memcpy( thisSngl, thisCoincEntry, sizeof(SnglRingdownTable) );
        thisSngl->next = NULL;

        /* create an eventId and populate the id */
        eventId = (EventIDColumn *) LALCalloc( 1, sizeof(EventIDColumn) );
        if ( thisCoincEntry->event_id->id )
        {
          /* event id number exists, use it */
          eventId->id = thisCoincEntry->event_id->id;
        }
        else if ( gpsStartTime )
        {
          eventId->id = LAL_INT8_C(1000000000) *
            (INT8) gpsStartTime->gpsSeconds + (INT8) eventNum;
        }
        else
        {
          XLALPrintError(
              "Event does not have id and no GPS start time given" );
          while ( snglHead )
          {
            thisSngl = snglHead;
            snglHead = snglHead->next;
            XLALFreeSnglRingdown( &thisSngl );
          }
          XLAL_ERROR_NULL(XLAL_EIO);
        }

        if ( slideNum < 0 )
        {
          eventId->id += LAL_INT8_C(100000)* (-1 *slideNum + 5000);
        }
        else
        {
          eventId->id += LAL_INT8_C(100000) * slideNum;
        }
        thisSngl->event_id = eventId;
        eventId->snglRingdownTable = thisSngl;

        /* copy the ds2 values from the coinc_table */
        thisSngl->ds2_H1H2 = thisCoinc->ds2_H1H2;
        thisSngl->ds2_H1L1 = thisCoinc->ds2_H1L1;
        thisSngl->ds2_H1V1 = thisCoinc->ds2_H1V1;
        thisSngl->ds2_H2L1 = thisCoinc->ds2_H2L1;
        thisSngl->ds2_H2V1 = thisCoinc->ds2_H2V1;
        thisSngl->ds2_L1V1 = thisCoinc->ds2_L1V1;
      }
    }
  }

  return( snglHead );

}



int
XLALCoincRingdownIfos (
    CoincRingdownTable  *coincRingdown,
    char                *ifos
    )

{
  InterferometerNumber  ifoNumber  = LAL_UNKNOWN_IFO;
  int                   ifosMatch  = 1;
  CHAR                  ifo[LIGOMETA_IFO_MAX];

  if ( !coincRingdown )
  {
    return ( 0 );
  }

  for( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    XLALReturnIFO( ifo, ifoNumber);

    /* check that the coinc is of the correct type */
    if ( (coincRingdown->snglRingdown[ifoNumber] &&  !strstr(ifos,ifo)) ||
        (!coincRingdown->snglRingdown[ifoNumber] &&  strstr(ifos,ifo)) )
    {
      ifosMatch = 0;
      break;
    }
  }
  return( ifosMatch );
}


int
XLALCoincRingdownIfosCut(
    CoincRingdownTable **coincHead,
    char                *ifos
    )
{
  CoincRingdownTable    *prevCoinc = NULL;
  CoincRingdownTable    *thisCoinc = NULL;
  int                    numCoinc = 0;

  thisCoinc = *coincHead;
  *coincHead = NULL;

  while ( thisCoinc )
  {
    CoincRingdownTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    if ( XLALCoincRingdownIfos( tmpCoinc, ifos ) )
    {
      /* ifos match so keep tmpCoinc */
      if ( ! *coincHead  )
      {
        *coincHead = tmpCoinc;
      }
      else
      {
        prevCoinc->next = tmpCoinc;
      }
      tmpCoinc->next = NULL;
      prevCoinc = tmpCoinc;
      ++numCoinc;
    }
    else
    {
      /* discard tmpCoinc */
      XLALFreeCoincRingdown( &tmpCoinc );
    }
  }

  return( numCoinc );
}



UINT8
XLALCoincRingdownIdNumber (
    CoincRingdownTable  *coincRingdown
    )

{
  SnglRingdownTable    *thisSngl = NULL;
  InterferometerNumber  ifoNumber  = LAL_UNKNOWN_IFO;

  if ( !coincRingdown )
  {
    XLAL_ERROR(XLAL_EIO);
  }

  for( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    EventIDColumn *thisID = NULL;
    if ( (thisSngl = coincRingdown->snglRingdown[ifoNumber]) )
    {
      /* loop over the list of eventID's until we get to the one that
       * points to thisCoinc */
      thisID = thisSngl->event_id;

      while ( thisID )
      {
        /* test if thisID points to our coinc */
        if ( thisID->coincRingdownTable == coincRingdown )
        {
          return( thisID->id );
          break;
        }
      }
    }
  }
  /* should never get here */
  XLALPrintError( "Unable to find id associated to this event" );
  XLAL_ERROR(XLAL_EIO);
}


CoincRingdownTable *
XLALCoincRingdownSlideCut(
    CoincRingdownTable **coincHead,
    int                  slideNum
    )
{
  CoincRingdownTable    *prevCoinc      = NULL;
  CoincRingdownTable    *thisCoinc      = NULL;
  CoincRingdownTable    *slideHead      = NULL;
  CoincRingdownTable    *thisSlideCoinc = NULL;

  UINT8 idNumber = 0;

  if( slideNum < 0 )
  {
    slideNum = 5000 - slideNum;
  }

  thisCoinc = *coincHead;
  *coincHead = NULL;

  while ( thisCoinc )
  {
    idNumber = XLALCoincRingdownIdNumber( thisCoinc );

    if ( (int) ((idNumber % 1000000000) / 100000) == slideNum )
    {
      /* add thisCoinc to the slideCoinc list */
      if ( slideHead )
      {
        thisSlideCoinc = thisSlideCoinc->next = thisCoinc;
      }
      else
      {
        slideHead = thisSlideCoinc = thisCoinc;
      }

      /* remove from coincHead list */
      if ( prevCoinc )
      {
        prevCoinc->next = thisCoinc->next;
      }

      thisCoinc = thisCoinc->next;
      thisSlideCoinc->next = NULL;
    }
    else
    {
      /* move along the list */
      if( ! *coincHead )
      {
        *coincHead = thisCoinc;
      }

      prevCoinc = thisCoinc;
      thisCoinc = thisCoinc->next;
    }
  }
  return( slideHead );
}






INT4 XLALCountCoincRingdown( CoincRingdownTable *head )

{
  INT4 length;
  CoincRingdownTable *event;

  if ( !head )
  {
    return( 0 );
  }

  /* count the number of events in the list */
  for(length = 0, event = head; event; event = event->next)
    length++;

  return length;
}



CoincRingdownTable *
XLALStatCutCoincRingdown (
    CoincRingdownTable         *eventHead,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    REAL4                       statCut
    )

{
  CoincRingdownTable    *thisEvent = NULL;
  CoincRingdownTable    *prevEvent = NULL;


  thisEvent = eventHead;
  eventHead = NULL;

  while ( thisEvent )
  {
    CoincRingdownTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( XLALCoincRingdownStat(tmpEvent,coincStat,bittenLParams) >= statCut )
    {
      /* keep this template */
      if ( ! eventHead  )
      {
        eventHead = tmpEvent;
      }
      else
      {
        prevEvent->next = tmpEvent;
      }
      tmpEvent->next = NULL;
      prevEvent = tmpEvent;
    }
    else
    {
      /* discard this template */
      XLALFreeCoincRingdown ( &tmpEvent );
    }
  }
  return( eventHead );
}




SnglRingdownTable *
XLALCompleteCoincRingdown (
    CoincRingdownTable         *eventHead,
    int                         ifoList[LAL_NUM_IFO]
    )

{
  CoincRingdownTable    *thisCoinc = NULL;
  SnglRingdownTable     *snglHead  = NULL;
  SnglRingdownTable     *thisSngl   = NULL;
  InterferometerNumber   ifoNumber  = LAL_UNKNOWN_IFO;
  InterferometerNumber   ifoNum  = LAL_UNKNOWN_IFO;

  for ( thisCoinc = eventHead; thisCoinc; thisCoinc = thisCoinc->next )
  {
    for ( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
    {
      if ( ifoList[ifoNumber] && !thisCoinc->snglRingdown[ifoNumber] )
      {
        /* we need to add a trigger for this ifo with zero snr,
         * but correct end time */
        if ( !snglHead )
        {
          snglHead = thisSngl = (SnglRingdownTable *)
              LALCalloc( 1, sizeof(SnglRingdownTable) );
        }
        else
        {
          thisSngl = thisSngl->next = (SnglRingdownTable *)
              LALCalloc( 1, sizeof(SnglRingdownTable) );
        }
        /* check that the sngl was allocated successfully */
        if ( !thisSngl )
        {
          while ( snglHead )
          {
            thisSngl = snglHead;
            snglHead = snglHead->next;
            LALFree(thisSngl);
          }
          XLAL_ERROR_NULL(XLAL_ENOMEM);
        }

        /* populate the ifo field */
        XLALReturnIFO(thisSngl->ifo,ifoNumber);
        XLALPrintInfo( "Appending a zero snr trigger for %s\n", thisSngl->ifo);

        /* obtain the end time */
        ifoNum = (InterferometerNumber) 0;
        while (!thisCoinc->snglRingdown[ifoNum]) ifoNum++;
        thisSngl->start_time = thisCoinc->snglRingdown[ifoNum]->start_time;

        /* add sngl to coinc */
        thisCoinc = XLALAddSnglRingdownToCoinc( thisCoinc, thisSngl );
      }
    }
  }
  return( snglHead );
}




CoincRingdownTable *
XLALPlayTestCoincRingdown(
    CoincRingdownTable         *eventHead,
    LALPlaygroundDataMask      *dataType
    )

{
  CoincRingdownTable    *coincEventList = NULL;
  CoincRingdownTable    *thisEvent = NULL;
  CoincRingdownTable    *prevEvent = NULL;

  INT8 triggerTime = 0;
  INT4 isPlay = 0;
  INT4 numTriggers;

  /* Remove all the triggers which are not of the desired type */

  numTriggers = 0;
  thisEvent = eventHead;

  if ( (*dataType == playground_only) || (*dataType == exclude_play) )
  {
    while ( thisEvent )
    {
      CoincRingdownTable *tmpEvent = thisEvent;
      thisEvent = thisEvent->next;

      triggerTime = XLALCoincRingdownTimeNS( tmpEvent );
      isPlay = XLALINT8NanoSecIsPlayground( triggerTime );

      if ( ( (*dataType == playground_only)  && isPlay ) ||
          ( (*dataType == exclude_play) && ! isPlay) )
      {
        /* keep this trigger */
        if ( ! coincEventList  )
        {
          coincEventList = tmpEvent;
        }
        else
        {
          prevEvent->next = tmpEvent;
        }
        tmpEvent->next = NULL;
        prevEvent = tmpEvent;
        ++numTriggers;
      }
      else
      {
        /* discard this template */
        XLALFreeCoincRingdown ( &tmpEvent );
      }
    }
    eventHead = coincEventList;
    if ( *dataType == playground_only )
    {
      XLALPrintInfo( "Kept %d playground triggers \n", numTriggers );
    }
    else if ( *dataType == exclude_play )
    {
      XLALPrintInfo( "Kept %d non-playground triggers \n", numTriggers );
    }
  }
  else if ( *dataType == all_data )
  {
    XLALPrintInfo( "Keeping all triggers since all_data specified\n" );
  }
  else
  {
    XLALPrintInfo( "Unknown data type, returning no triggers\n" );
    eventHead = NULL;
  }

  return(eventHead);
}




/**
 * Is used to recreate a list of coinc
 * inspirals from a list of \c snglInspiralTables with populated
 * \c eventIDColumn.  The code searches for entries in
 * \c snglInspiral which have the same numerical value of the \c id
 * field in the \c eventIDColumn.
 */
int
XLALRecreateRingdownCoincFromSngls(
    CoincRingdownTable        **coincPtr,
    SnglRingdownTable          *snglRingdown
    )

{
  SnglRingdownTable    *thisSngl  = NULL;
  CoincRingdownTable   *thisCoinc = NULL;
  CoincRingdownTable   *prevCoinc = NULL;
  CoincRingdownTable   *coincHead = NULL;
  UINT8                 eventId = 0;
  INT4                  numCoincs = 0;
  InterferometerNumber  ifoNumber = LAL_UNKNOWN_IFO;
  InterferometerNumber  ifoInCoinc = LAL_UNKNOWN_IFO;


  if ( !snglRingdown )
  {
    XLALPrintInfo(
      "XLALRecreateCoincFromSngls: Empty snglRingdown passed as input" );
    return( 0 );
  }

  /* loop over the linked list of sngl ringdowns */
  for( thisSngl = snglRingdown; thisSngl; thisSngl = thisSngl->next )
  {
    ifoNumber = (InterferometerNumber) XLALIFONumber( thisSngl->ifo );
    thisCoinc = coincHead;
    while ( thisCoinc )
    {
      /* loop over the interferometers to get the event_id*/
      for ( ifoInCoinc = (InterferometerNumber) 0; ifoInCoinc < LAL_NUM_IFO; ifoInCoinc++)
      {
        if ( thisCoinc->snglRingdown[ifoInCoinc] )
        {
          eventId = thisCoinc->snglRingdown[ifoInCoinc]->event_id->id;
          break;
        }
      }

      if ( thisSngl->event_id->id == eventId )
      {
        /* thisSngl is part of the coinc, so add it */
        if ( thisCoinc->snglRingdown[ifoNumber] )
        {
          /* already have an event for this ifo */
          XLALPrintError(
              "Already have a single from this ifo with event id %lld",
              eventId);
          /* free memory */
          while ( coincHead )
          {
            thisCoinc = coincHead;
            coincHead = coincHead->next;
            LALFree(thisCoinc);
          }
          XLAL_ERROR( XLAL_EDATA);
        }
        else
        {
          thisCoinc->snglRingdown[ifoNumber] = thisSngl;
          thisCoinc->numIfos += 1;
          thisSngl->event_id->coincRingdownTable = thisCoinc;
          break;
        }
      }

      /* proceed to the next coinc */
      prevCoinc = thisCoinc;
      thisCoinc = thisCoinc->next;
    }

    if ( thisSngl->event_id->id != eventId )
    {
      /* need to start a new coinc */
      if ( coincHead )
      {
        thisCoinc = prevCoinc->next =
          LALCalloc( 1, sizeof(CoincRingdownTable) );
      }
      else
      {
        thisCoinc = coincHead = LALCalloc( 1, sizeof(CoincRingdownTable) );
      }
      if ( !thisCoinc )
      {
        /* out of memory: free memory + exit*/
        while ( coincHead )
        {
          thisCoinc = coincHead;
          coincHead = coincHead->next;
          LALFree( thisCoinc );
        }
        XLAL_ERROR(XLAL_ENOMEM);
      }

      thisCoinc->snglRingdown[ifoNumber] = thisSngl;
      thisCoinc->numIfos = 1;
      thisSngl->event_id->coincRingdownTable = thisCoinc;
      numCoincs +=1;
    }
  }

  *coincPtr = coincHead;

  return( numCoincs );
}

#if 0

/**
 * Is used to generate a coherent bank from
 * a list of \c coincInspiralTables.  The coherent bank has the same mass
 * parameters for each ifo.  These are currently chosen as the mass parameters
 * of the trigger in the coinc with the highest \c snr.  If the
 * \c ifos field is not \c NULL, then a template is generated for
 * every ifo in \c ifos.  If it is \c NULL then templates are only
 * generated for those ifos which have triggers in the coinc.
 */
int
XLALGenerateCoherentBank(
    SnglRingdownTable         **coherentBank,
    CoincRingdownTable         *coincInput,
    CHAR                       *ifos
    )

{
  InterferometerNumber  ifoInCoinc = LAL_UNKNOWN_IFO;
  InterferometerNumber  ifoNumber  = LAL_UNKNOWN_IFO;
  InterferometerNumber  ifoMax  = LAL_UNKNOWN_IFO;
  SnglRingdownTable    *bankHead = NULL;
  SnglRingdownTable    *currentTrigger = NULL;
  CoincRingdownTable   *thisCoinc = NULL;
  INT4                  numTmplts = 0;

  if ( !coincInput )
  {
    XLALPrintInfo(
      "XLALGenerateCoherentBank: Empty coincInput passed as input" );
    return( 0 );
  }

  for ( thisCoinc = coincInput; thisCoinc; thisCoinc = thisCoinc->next )
  {
    REAL4 max_snr = 0;

    /* loop over the interferometers to get the highest snr*/
    for ( ifoInCoinc = 0; ifoInCoinc < LAL_NUM_IFO; ifoInCoinc++)
    {
      if (( thisCoinc->snglRingdown[ifoInCoinc] ) &&
        (thisCoinc->snglRingdown[ifoInCoinc]->snr > max_snr) )
      {
        max_snr = thisCoinc->snglRingdown[ifoInCoinc]->snr;
        ifoMax = ifoInCoinc;
      }
    }

    for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ++ifoNumber )
    {

      CHAR ifo[LIGOMETA_IFO_MAX];

      XLALReturnIFO( ifo, ifoNumber);

      /* decide whether we want a template for this ifo */
      if ( (thisCoinc->snglRingdown[ifoNumber] &&  !ifos) ||
           ( ifos && strstr(ifos,ifo)) )
      {
        numTmplts++;

        if( bankHead )
        {
          currentTrigger = currentTrigger->next =
            LALCalloc( 1, sizeof(SnglRingdownTable) );
        }
        else
        {
          bankHead = currentTrigger =
              LALCalloc( 1, sizeof(SnglRingdownTable) );
        }
        if ( !currentTrigger )
        {
          goto error;
        }
        /* copy the info from the loudest trigger */
        memcpy(currentTrigger, thisCoinc->snglRingdown[ifoMax],
            sizeof(SnglRingdownTable));
        /* terminate the list */
        currentTrigger->next = NULL;
        currentTrigger->event_id = NULL;
        /* set the ifo */
        snprintf( currentTrigger->ifo, LIGOMETA_IFO_MAX, ifo );
        /* set the event id */
        currentTrigger->event_id = LALCalloc( 1, sizeof(EventIDColumn) );
        if ( !(currentTrigger->event_id) )
        {
		      goto error;
        }
        currentTrigger->event_id->id =
          thisCoinc->snglRingdown[ifoMax]->event_id->id;
        currentTrigger->event_id->snglRingdownTable = currentTrigger;
      }
    }
  }

  *coherentBank = bankHead;
  return( numTmplts );

  error:
  while ( bankHead )
  {
    currentTrigger = bankHead;
    bankHead = bankHead->next;
    XLALFreeSnglRingdown( &currentTrigger );
  }
  XLAL_ERROR(XLAL_ENOMEM);

}
#endif


/**
 * Is used to perform a distance cut between
 * the triggers in a coincidence.
 */
CoincRingdownTable *
XLALRingdownDistanceCut(
    CoincRingdownTable        **coincRingdown,
    RingdownAccuracyList       *accuracyParams
    )

{
  InterferometerNumber  ifoA = LAL_UNKNOWN_IFO;
  InterferometerNumber  ifoB = LAL_UNKNOWN_IFO;
  CoincRingdownTable   *thisCoinc = NULL;
  CoincRingdownTable   *prevCoinc = NULL;
  CoincRingdownTable   *coincHead = NULL;

  thisCoinc = *coincRingdown;
  coincHead = NULL;

  while( thisCoinc )
  {
    INT4  discardTrigger = 0;
    REAL4 kappaA = 0, kappaB = 0;
    REAL4 distA = 0, distB = 0;
    REAL8 sigmasqA = 0, sigmasqB = 0;

    CoincRingdownTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    for ( ifoA = (InterferometerNumber) 0; ifoA < LAL_NUM_IFO; ifoA++ )
    {
      kappaA = accuracyParams->ifoAccuracy[ifoA].kappa;
      for ( ifoB = (InterferometerNumber) ( ifoA + 1 ); ifoB < LAL_NUM_IFO; ifoB++ )
      {
        kappaB = accuracyParams->ifoAccuracy[ifoB].kappa;

        if( tmpCoinc->snglRingdown[ifoA] && kappaA
            && tmpCoinc->snglRingdown[ifoB] && kappaB  )
        {
          /* perform the distance consistency test */
          sigmasqA = tmpCoinc->snglRingdown[ifoA]->sigma_sq;
          sigmasqB = tmpCoinc->snglRingdown[ifoB]->sigma_sq;
          distA = tmpCoinc->snglRingdown[ifoA]->eff_dist;
          distB = tmpCoinc->snglRingdown[ifoB]->eff_dist;

          if( ( sigmasqA > sigmasqB && distA/distB > kappaA ) ||
              ( sigmasqB > sigmasqA && distB/distA > kappaB ) )
          {
            discardTrigger = 1;
            break;
          }
        }
      }

      if ( discardTrigger )
      {
        break;
      }
    }


    if( discardTrigger )
    {
      XLALFreeCoincRingdown( &tmpCoinc );
    }
    else
    {
      if ( ! coincHead )
      {
        coincHead = tmpCoinc;
      }
      else
      {
        prevCoinc->next = tmpCoinc;
      }
      tmpCoinc->next = NULL;
      prevCoinc = tmpCoinc;
    }
  }
  *coincRingdown = coincHead;
  return( coincHead );
}


#if 0

/**
 * Extracts all single inspirals from a
 * specific ifo which are in coinc inspirals.  The output \c snglPtr is a
 * pointer to a linked list of single inspiral tables.  That list contains only
 * single inspirals from the specified \c ifo which are found in
 * coincidence.
 */
void
LALCoincCutSnglInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead
    )

{
  SnglInspiralTable  *eventList = NULL;
  SnglInspiralTable  *prevEvent = NULL;
  SnglInspiralTable  *thisEvent = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* check that eventHead is non-null */
  ASSERT( eventHead, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );


  /* Scan through a linked list of sngl_inspiral tables and return a
     pointer to the head of a linked list of tables which are in coincs */

  thisEvent = *eventHead;

  while ( thisEvent )
  {
    INT4           keepEvent = 0;
    EventIDColumn *eventID   = NULL;

    SnglInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( tmpEvent->event_id )
    {
      for( eventID = tmpEvent->event_id; tmpEvent; tmpEvent = tmpEvent->next )
      {
        /* there is an eventID, check if there's a coinc inspiral */
        eventID = tmpEvent->event_id;

        if ( eventID->coincInspiralTable )
        {
          keepEvent = 1;
          break;
        }

        eventID = eventID->next;
      }
    }

    if ( keepEvent )
    {
      /* event is in a coinc so keep it */
      if ( ! eventList  )
      {
        eventList = tmpEvent;
      }
      else
      {
        prevEvent->next = tmpEvent;
      }
      tmpEvent->next = NULL;
      prevEvent = tmpEvent;
    }
    else
    {
      /* discard this trigger since it's not in a coinc */
      LALFreeSnglInspiral ( status->statusPtr, &tmpEvent );
    }
  }

  *eventHead = eventList;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
#endif

INT8
XLALCoincRingdownTimeNS (
    const CoincRingdownTable         *coincRingdown
    )
{
  InterferometerNumber  ifoNumber;
  INT8 startTime = 0;

  for( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    if ( coincRingdown->snglRingdown[ifoNumber] )
    {
      startTime = XLALGPSToINT8NS(
          &(coincRingdown->snglRingdown[ifoNumber]->start_time) );
      return(startTime);
    }
  }
  XLAL_ERROR(XLAL_EIO);
}


REAL4
XLALCoincRingdownStat(
    CoincRingdownTable         *coincRingdown,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams
    )
{
  InterferometerNumber  ifoNumber;
  SnglRingdownTable    *snglRingdown;
  REAL4                 statValues[LAL_NUM_IFO];
  REAL4 statValue = 0;
  INT4  i;
  INT4  ifoCounter = 0;

  if( coincStat == no_stat )
  {
    return(0);
  }

  /* for bittenL only*/
  if( coincStat == bitten_l )
  {
    for ( i = 0; i < LAL_NUM_IFO ; i++)
    {
      statValues[i] = 1e9; /* sufficiently high values */
    }
  }


  for( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    if ( (snglRingdown = coincRingdown->snglRingdown[ifoNumber]) )
    {
      /* count the number of IFOs for this coincidence */
      ifoCounter++;

      if ( coincStat == snrsq )
      {
        statValue += snglRingdown->snr * snglRingdown->snr;
      }

      else if ( coincStat == bitten_l )
      {
        statValues[ifoNumber] = bittenLParams->param_a[ifoNumber]
                * snglRingdown->snr
                + bittenLParams->param_b[ifoNumber];
        statValue += snglRingdown->snr ;
      }

    }
  }

  /*    for the bitten L case only , we need to compare different
        values and keep the minimum one */
  if ( coincStat == bitten_l )
  {
    if (coincStat == bitten_l || ifoCounter<3) {
      for( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( (snglRingdown = coincRingdown->snglRingdown[ifoNumber]) )
        {
          if (statValues[ifoNumber] < statValue)
          {
           statValue = statValues[ifoNumber];
          }
        }
      }
    }
  }


  return( statValue );
}



int
XLALClusterCoincRingdownTable (
    CoincRingdownTable        **coincList,
    INT8                        dtimeNS,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams
    )

{
  CoincRingdownTable     *thisCoinc = NULL;
  CoincRingdownTable     *prevCoinc = NULL;
  CoincRingdownTable     *nextCoinc = NULL;
  int                     numCoincClust = 0;
  REAL4 thisStat = 0;
  REAL4 nextStat = 0;
  InterferometerNumber  ifoNumber;
  SnglRingdownTable    *snglRingdown;

  CoincInspiralStatistic  tripleStat = snrsq;

  if ( !coincList )
  {
    XLAL_ERROR(XLAL_EIO);
  }

  if ( ! *coincList )
  {
    XLALPrintInfo(
      "XLALClusterCoincRingdownTable: Empty coincList passed as input" );
  return( 0 );
  }

  thisCoinc = (*coincList);
  nextCoinc = (*coincList)->next;
  *coincList = NULL;

  while ( nextCoinc )
  {
    INT8 thisTime = XLALCoincRingdownTimeNS( thisCoinc );
    INT8 nextTime = XLALCoincRingdownTimeNS( nextCoinc );

    /* find events within the cluster window */
    if ( (nextTime - thisTime) < dtimeNS )
    {
      if ( thisCoinc->numIfos > 2)
      {
         thisStat = XLALCoincRingdownStat( thisCoinc, tripleStat, bittenLParams );
      }
      else
      {
        thisStat = XLALCoincRingdownStat( thisCoinc, coincStat, bittenLParams );
      }

      if ( nextCoinc->numIfos > 2)
      {
        nextStat = XLALCoincRingdownStat( nextCoinc, tripleStat, bittenLParams );
      }
      else
      {
        nextStat = XLALCoincRingdownStat( nextCoinc, coincStat, bittenLParams );
      }

      for( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( (snglRingdown = thisCoinc->snglRingdown[ifoNumber]) )
        {
          thisCoinc->snglRingdown[ifoNumber]->epsilon=thisStat;
        }
        if ( (snglRingdown = nextCoinc->snglRingdown[ifoNumber]) )
        {
          nextCoinc->snglRingdown[ifoNumber]->epsilon=nextStat;
        }
      }

      if ( nextStat > thisStat )
      {
        /* displace previous event in cluster */
        if( prevCoinc )
        {
          prevCoinc->next = nextCoinc;
        }
        XLALFreeCoincRingdown( &thisCoinc );
        thisCoinc = nextCoinc;
        nextCoinc = thisCoinc->next;
      }
      else
      {
        /* otherwise just dump next event from cluster */
        thisCoinc->next = nextCoinc->next;
        XLALFreeCoincRingdown ( &nextCoinc );
        nextCoinc = thisCoinc->next;
      }
    }
    else
    {
      /* otherwise we keep this unique event trigger */
      if ( ! (*coincList) )
      {
        *coincList = thisCoinc;
      }
      if ( thisCoinc->numIfos > 2)
      {
        thisStat = XLALCoincRingdownStat( thisCoinc, tripleStat, bittenLParams );
      }
      else
      {
        thisStat = XLALCoincRingdownStat( thisCoinc, coincStat, bittenLParams );
      }

      for( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( (snglRingdown = thisCoinc->snglRingdown[ifoNumber]) )
        {
          thisCoinc->snglRingdown[ifoNumber]->epsilon=thisStat;
        }
      }

      prevCoinc = thisCoinc;
      thisCoinc = thisCoinc->next;
      nextCoinc = thisCoinc->next;
      ++numCoincClust;
    }
  }

  /* store the last event */
  if ( ! (*coincList) )
  {
    *coincList = thisCoinc;
  }
  ++numCoincClust;

  return(numCoincClust);
}



int
XLALCompareCoincRingdownByTime (
    const void *a,
    const void *b
    )

{
  const CoincRingdownTable *aPtr = *((const CoincRingdownTable * const *)a);
  const CoincRingdownTable *bPtr = *((const CoincRingdownTable * const *)b);
  INT8 ta, tb;

  ta = XLALCoincRingdownTimeNS ( aPtr );
  tb = XLALCoincRingdownTimeNS ( bPtr );

  if ( ta > tb )
  {
    return 1;
  }
  else if ( ta < tb )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}



CoincRingdownTable *
XLALSortCoincRingdown (
    CoincRingdownTable  *eventHead,
    int(*comparfunc)    (const void *, const void *)
    )

{
  INT4                   i;
  INT4                   numEvents = 0;
  CoincRingdownTable    *thisEvent = NULL;
  CoincRingdownTable   **eventHandle = NULL;

  /* count the number of events in the linked list */
  for ( thisEvent = eventHead; thisEvent; thisEvent = thisEvent->next )
  {
    ++numEvents;
  }
  if ( ! numEvents )
  {
     XLALPrintInfo(
      "XLALSortCoincRingdown: Empty coincRingdown passed as input" );
    return( eventHead );
  }

  /* allocate memory for an array of pts to sort and populate array */
  eventHandle = (CoincRingdownTable **)
    LALCalloc( numEvents, sizeof(CoincRingdownTable *) );
  for ( i = 0, thisEvent = eventHead; i < numEvents;
      ++i, thisEvent = thisEvent->next )
  {
    eventHandle[i] = thisEvent;
  }

  /* qsort the array using the specified function */
  qsort( eventHandle, numEvents, sizeof(eventHandle[0]), comparfunc );

  /* re-link the linked list in the right order */
  thisEvent = eventHead = eventHandle[0];
  for ( i = 1; i < numEvents; ++i )
  {
    thisEvent = thisEvent->next = eventHandle[i];
  }
  thisEvent->next = NULL;

  /* free the internal memory */
  LALFree( eventHandle );

  return( eventHead );

}


void
LALRingdownH1H2Consistency(
    LALStatus                  *status,
    CoincRingdownTable        **coincRingdown,
    REAL4			H1snrCutThreshold,
    LALSegList                 *vetoSegsH1,
    LALSegList                 *vetoSegsH2
    )

{
  CoincRingdownTable   *thisCoinc = NULL;
  CoincRingdownTable   *prevCoinc = NULL;
  CoincRingdownTable   *coincHead = NULL;
  REAL4 snrH1;
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  thisCoinc = *coincRingdown;
  coincHead = NULL;

  while( thisCoinc )
  {
    INT4  discardTrigger = 0;
    CoincRingdownTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

/* Discard the H1 trigger if it's snr is above the threshold, */
/*   a H2 trigger does not exist and H2 is not vetoed.*/
    if( tmpCoinc->snglRingdown[LAL_IFO_H1] && !tmpCoinc->snglRingdown[LAL_IFO_H2])
    {
      snrH1 = tmpCoinc->snglRingdown[LAL_IFO_H1]->snr;
      if( snrH1 > H1snrCutThreshold )
      {
        if ( vetoSegsH2->initMagic == SEGMENTSH_INITMAGICVAL )
	{
	  if (!XLALSegListSearch( vetoSegsH2,
			      &(tmpCoinc->snglRingdown[LAL_IFO_H1]->start_time)))
          {
            discardTrigger =1;
          }
        }
        else
        {
          discardTrigger = 1;
        }
      }
    }

/* Discard the H2 trigger if a H1 trigger does not exist and H1 was not vetoed. */
    if( tmpCoinc->snglRingdown[LAL_IFO_H2] && !tmpCoinc->snglRingdown[LAL_IFO_H1])
    {
      if ( vetoSegsH1->initMagic == SEGMENTSH_INITMAGICVAL )
      {
        if (!XLALSegListSearch( vetoSegsH1,
                            &(tmpCoinc->snglRingdown[LAL_IFO_H2]->start_time)))
        {
          discardTrigger =1;
        }
      }
      else
      {
        discardTrigger = 1;
      }
    }

    if( discardTrigger )
    {
      XLALFreeCoincRingdown( &tmpCoinc );
    }
    else
    {
      if ( ! coincHead )
      {
        coincHead = tmpCoinc;
      }
      else
      {
        prevCoinc->next = tmpCoinc;
      }
      tmpCoinc->next = NULL;
      prevCoinc = tmpCoinc;
    }
  }
  *coincRingdown = coincHead;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
/*@}*/ /* end:CoincRingdownUtils_c */
