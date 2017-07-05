/*
*  Copyright (C) 2007 Stas Babak, Alexander Dietz, Drew Keppel, Gareth Jones, Jolien Creighton, Patrick Brady, Stephen Fairhurst, Craig Robinson , Thomas Cokelaer
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
 * File Name: CoincInspiralUtils.c
 *
 * Author: Brady, P. R., Brown, D. A., and Fairhurst, S
 *
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/XLALError.h>
#include <lal/CoincInspiralEllipsoid.h>

/**
 * \defgroup CoincInspiralUtils_c Module CoincInspiralUtils.c
 * \ingroup CoincInspiralEllipsoid_h
 * \author Fairhurst, S.
 *
 * \brief Blah.
 *
 * ### Description ###
 *
 * <tt>LALCreateNIFOCoincList()</tt> takes linked list of
 * \c CoincInspiralTables, assumed to contain (N-1) ifo coincidences and
 * creates all N ifo coincidences.  Both the input and output list of
 * \c CoincInspiralTables are passed as \c coincHead.
 *
 * <tt>LALRemoveRepeatedCoincs()</tt> will remove any lower order coincidences
 * if they are contained in a higher order coincidence.  For example, if an H1-L1
 * double coincident trigger is also part of an H1-H2-L1 triple coincident
 * trigger, the double coincident trigger will be removed.  The head of the list
 * of coincident triggers is passed and returned as \c coincHead.
 *
 * <tt>LALSnglInspiralCoincTest()</tt> tests for coincidence between a single
 * inspiral and a coinc inspiral.  It works by testing for coincidence between
 * each non-null entry in the coinc inspiral and the single.  This is done using
 * <tt>LALCompareSnglInspiral()</tt>.  If all members of the coinc are found to be
 * coincident with the single, the <tt>accuracyParams.match</tt> is set to 1,
 * otherwise to 0.
 *
 * <tt>XLALCountCoincInspiral()</tt> scans through a linked list of coincidence
 * inspiral table and counts the number of events. This count is returned
 * as \c numTrigs.
 */
/*@{*/


void
LALSnglInspiralCoincTest(
    LALStatus                  *status,
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral,
    InspiralAccuracyList       *accuracyParams
    )

{

  INITSTATUS(status);

  XLALSnglInspiralCoincTest( coincInspiral, snglInspiral, accuracyParams );

  RETURN( status );
}



void
XLALSnglInspiralCoincTest(
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral,
    InspiralAccuracyList       *accuracyParams
    )

{
  SnglInspiralTable    *thisCoincEntry;
  INT4                  match = 1;
  INT4                  ifoNumber = 0;


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
        XLALCompareInspirals ( snglInspiral,
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

  return;
}



int
XLALCreateCoincSlideTable(
    CoincInspiralSlideTable   **slideTableHead,
    INT4                        numSlides
    )

{
  CoincInspiralSlideTable  *thisSlideTable = NULL;
  INT4                      idx = 0;
  INT4                      slideNum = 0;

  for ( idx = 0; idx < 2*numSlides; idx++ )
  {
    if ( idx < numSlides )
    {
      slideNum = idx + 1;
    }
    else
    {
      slideNum = idx - 2*numSlides;
    }

    if ( *slideTableHead )
    {
      thisSlideTable->next = (CoincInspiralSlideTable*)
          LALCalloc( 1, sizeof(CoincInspiralSlideTable) );
      thisSlideTable = thisSlideTable->next;
    }
    else
    {
      *slideTableHead = thisSlideTable = (CoincInspiralSlideTable*)
          LALCalloc( 1, sizeof(CoincInspiralSlideTable) );
    }

    if ( !thisSlideTable )
    {
      /* out of memory: free memory + exit*/
      while ( *slideTableHead )
      {
        thisSlideTable = *slideTableHead;
        *slideTableHead = (*slideTableHead)->next;
        LALFree( thisSlideTable );
      }

      XLAL_ERROR(XLAL_ENOMEM );
    }

    thisSlideTable->coincInspiral = NULL;
    thisSlideTable->slideNum = slideNum;
    thisSlideTable->slideTimeAnalyzed = 0;
    thisSlideTable->currentRate = 0;
    thisSlideTable->next = NULL;
  }

  return( numSlides );
}


INT8
XLALCoincInspiralTimeNS (
    const CoincInspiralTable         *coincInspiral
    )
{
  InterferometerNumber  ifoNumber;
  INT8 endTime = 0;

  for( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    if ( coincInspiral->snglInspiral[ifoNumber] )
    {
      endTime = XLALGPSToINT8NS(
          &(coincInspiral->snglInspiral[ifoNumber]->end) );
      return(endTime);
    }
  }
  XLAL_ERROR(XLAL_EIO);
}

REAL4
XLALCoincInspiralStat(
    const CoincInspiralTable   *coincInspiral,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams
    )
{
  InterferometerNumber  ifoNumber;
  SnglInspiralTable    *snglInspiral;
  REAL4                 statValues[LAL_NUM_IFO];
  REAL4 statValue = 0;
  INT4  i;
  INT4  ifoCounter = 0;
  /* This replaces the 250 in the effective snr formula. */
  REAL4 eff_snr_denom_fac = bittenLParams->eff_snr_denom_fac;
  REAL4 chisq_index = bittenLParams->chisq_index;

  if( coincStat == no_stat )
  {
    return(0);
  }

  /* for bittenL only*/
  if( coincStat == bitten_l || coincStat == bitten_lsq)
  {
    for ( i = 0; i < LAL_NUM_IFO ; i++)
    {
      statValues[i] = 1e9; /* sufficiently high values */
    }
  }


  for( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    if ( (snglInspiral = coincInspiral->snglInspiral[ifoNumber]) )
    {
      /* count the number of IFOs for this coincidence */
      ifoCounter++;

      if ( coincStat == snrsq )
      {
        statValue += snglInspiral->snr * snglInspiral->snr;
      }
      else if ( coincStat == effective_snrsq )
      {
        REAL4 tmp_snr = snglInspiral->snr;
        REAL4 tmp_chisq = snglInspiral->chisq;
        /* XXX Assuming that chisq_dof contains the number of bins, not dof */
        REAL4 tmp_bins = snglInspiral->chisq_dof;

        statValue += tmp_snr * tmp_snr /
          sqrt ( tmp_chisq/(2*tmp_bins-2) * (1+tmp_snr*tmp_snr/eff_snr_denom_fac) ) ;
      }
      else if ( coincStat == new_snrsq )
      {
        REAL4 tmp_snr = snglInspiral->snr;
        REAL4 tmp_chisq = snglInspiral->chisq;
        /* XXX Assuming that chisq_dof contains the number of bins, not dof */
        REAL4 tmp_bins = snglInspiral->chisq_dof;
        REAL4 tmp_chisq_r = 1.0;

        tmp_chisq_r = tmp_chisq/(2*tmp_bins -2);

        if ( tmp_chisq_r > 1.0 ) {
          statValue += tmp_snr * tmp_snr /
            pow( 0.5*(1 + pow(tmp_chisq_r, 0.5*chisq_index)), 2.0/chisq_index);
        }
        else statValue += tmp_snr * tmp_snr;
      }
      else if ( coincStat == bitten_l || coincStat == bitten_lsq)
      {
        statValues[ifoNumber] = bittenLParams->param_a[ifoNumber]
                * snglInspiral->snr
                - bittenLParams->param_b[ifoNumber];
        statValue += snglInspiral->snr * snglInspiral->snr ;
      }
      else if ( coincStat == s3_snr_chi_stat )
      {
        REAL4 tmp_snr = snglInspiral->snr;
        REAL4 tmp_chisq = snglInspiral->chisq;

        statValue += tmp_snr * tmp_snr * tmp_snr * tmp_snr /
          ( tmp_chisq * ( 250 + tmp_snr * tmp_snr ) );
      }
      else if ( coincStat == ifar )
      {
        statValue = 1/snglInspiral->alpha;
      }

    }
  }

  /*    for the bitten L case only , we need to compare different
        values and keep the minimum one */
  if ( coincStat == bitten_l || coincStat == bitten_lsq)
  {
    statValue = sqrt(statValue);

    if (coincStat == bitten_l || ifoCounter<3) {
      for( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( (snglInspiral = coincInspiral->snglInspiral[ifoNumber]) )
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
XLALCompareCoincInspiralByTime (
    const void *a,
    const void *b
    )

{
  const CoincInspiralTable *aPtr = *((const CoincInspiralTable * const *)a);
  const CoincInspiralTable *bPtr = *((const CoincInspiralTable * const *)b);
  INT8 ta, tb;

  ta = XLALCoincInspiralTimeNS ( aPtr );
  tb = XLALCoincInspiralTimeNS ( bPtr );

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

int
XLALCompareCoincInspiralByEffectiveSnr (
    const void *a,
    const void *b
    )

{
  const CoincInspiralTable *aPtr = *((const CoincInspiralTable * const *)a);
  const CoincInspiralTable *bPtr = *((const CoincInspiralTable * const *)b);
  REAL4 ta, tb;

  CoincInspiralStatistic coincStat = effective_snrsq;
  CoincInspiralStatParams    bittenLParams;
  memset( &bittenLParams, 0, sizeof(CoincInspiralStatParams   ) );
  /* Default value of denom fac is 250 to preserve old functionality */
  bittenLParams.eff_snr_denom_fac = 250.0;
  ta = XLALCoincInspiralStat(aPtr,coincStat,&bittenLParams);
  tb = XLALCoincInspiralStat(bPtr,coincStat,&bittenLParams);

  if ( ta > tb )
  {
    return -1;
  }
  else if ( ta < tb )
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

int  XLALComputeAndStoreEffectiveSNR( CoincInspiralTable *head, CoincInspiralStatistic *stat, CoincInspiralStatParams *par)
  {
  while (head)
    {
    head->stat = XLALCoincInspiralStat(head, *stat, par);
    head = head->next;
    }
  return 0;
  }

int
XLALCompareCoincInspiralByStat (
    const void *a,
    const void *b
    )

{
  const CoincInspiralTable *aPtr = *((const CoincInspiralTable * const *)a);
  const CoincInspiralTable *bPtr = *((const CoincInspiralTable * const *)b);
  REAL4 ta, tb;
  ta = aPtr->stat;
  tb = bPtr->stat;

  if ( ta > tb )
  {
    return -1;
  }
  else if ( ta < tb )
  {
    return 1;
  }
  else
  {
    return 0;
  }
}


CoincInspiralTable *
XLALSortCoincInspiralByStat (
    CoincInspiralTable  *eventHead,
    int(*comparfunc)    (const void *, const void *),
    CoincInspiralStatParams *statParams,
    CoincInspiralStatistic *stat
    )
{

  INT4                   i;
  INT4                   numEvents = 0;
  CoincInspiralTable    *thisEvent = NULL;
  CoincInspiralTable   **eventHandle = NULL;

  /* count the number of events in the linked list */
  for ( thisEvent = eventHead; thisEvent; thisEvent = thisEvent->next )
  {
    ++numEvents;
  }

  if ( ! numEvents )
  {
     XLALPrintInfo(
      "XLALSortCoincInspiral: Empty coincInspiral passed as input" );
    return( eventHead );
  }

  XLALComputeAndStoreEffectiveSNR( eventHead, stat, statParams);

  /* allocate memory for an array of pts to sort and populate array */
  eventHandle = (CoincInspiralTable **)
    LALCalloc( numEvents, sizeof(CoincInspiralTable *) );
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



CoincInspiralTable *
XLALSortCoincInspiral (
    CoincInspiralTable  *eventHead,
    int(*comparfunc)    (const void *, const void *)
    )

{
  INT4                   i;
  INT4                   numEvents = 0;
  CoincInspiralTable    *thisEvent = NULL;
  CoincInspiralTable   **eventHandle = NULL;

  /* count the number of events in the linked list */
  for ( thisEvent = eventHead; thisEvent; thisEvent = thisEvent->next )
  {
    ++numEvents;
  }
  if ( ! numEvents )
  {
     XLALPrintInfo(
      "XLALSortCoincInspiral: Empty coincInspiral passed as input" );
    return( eventHead );
  }

  /* allocate memory for an array of pts to sort and populate array */
  eventHandle = (CoincInspiralTable **)
    LALCalloc( numEvents, sizeof(CoincInspiralTable *) );
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


int
XLALCoincInspiralIfos (
    CoincInspiralTable  *coincInspiral,
    const char          *ifos
    )

{
  InterferometerNumber  ifoNumber  = LAL_UNKNOWN_IFO;
  int                   ifosMatch  = 1;
  CHAR                  ifo[LIGOMETA_IFO_MAX];

  if ( !coincInspiral )
  {
    return ( 0 );
  }

  for( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    XLALReturnIFO( ifo, ifoNumber);

    /* check that the coinc is of the correct type */
    if ( (coincInspiral->snglInspiral[ifoNumber] &&  !strstr(ifos,ifo)) ||
        (!coincInspiral->snglInspiral[ifoNumber] &&  strstr(ifos,ifo)) )
    {
      ifosMatch = 0;
      break;
    }
  }
  return( ifosMatch );
}



INT4 XLALCountCoincInspiral( CoincInspiralTable *head )

{
  INT4 length;
  CoincInspiralTable *event;

  if ( !head )
  {
    return( 0 );
  }

  /* count the number of events in the list */
  for(length = 0, event = head; event; event = event->next)
    length++;

  return length;
}


int
XLALCalcExpFitNLoudestBackground (
    CoincInspiralTable         *coincSlideHead,
    int                         fitNum,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    REAL4                      *fitStat,
    REAL4                      *fitA,
    REAL4                      *fitB
    )

{
  CoincInspiralTable    *thisSlideEvent = coincSlideHead;
  int idx = 0;
  REAL4 Delta = 0;
  REAL4 X0 = 0;
  REAL4 X1 = 0;
  REAL4 X2 = 0;
  REAL4 Y0 = 0;
  REAL4 Y1 = 0;
  REAL4 thisStat = 0;

  for (idx = 1, thisSlideEvent = coincSlideHead; idx <= fitNum; idx++,
      thisSlideEvent = thisSlideEvent->next )
  {
    if ( ! thisSlideEvent )
    {
      /* should never get here */
      XLALPrintError( "Not enough Background Triggers: have %d, need %d",
          idx - 1, fitNum );
      XLAL_ERROR(XLAL_ERANGE);
    }

    thisStat = XLALCoincInspiralStat(thisSlideEvent, coincStat, bittenLParams);
    X0 += idx;
    X1 += thisStat * idx;
    X2 += pow(thisStat, 2.0) * idx;
    Y0 += idx * log(idx);
    Y1 += thisStat * idx * log(idx);
  }

  *fitStat = thisStat;

  Delta = X0 * X2 - X1 * X1;

  *fitA = X2 * Y0 - X1 * Y1;
  *fitA /= Delta;
  *fitA = exp(*fitA);

  *fitB = X0 * Y1 - X1 * Y0;
  *fitB /= Delta;

  return( fitNum );
}



REAL4
XLALRateCalcCoincInspiral (
    CoincInspiralTable         *coincZeroHead,
    CoincInspiralTable         *coincSlideHead,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    REAL4                       timeAnalyzed,
    REAL4                       fitStat,
    REAL4                       fitA,
    REAL4                       fitB
    )

{
  CoincInspiralTable    *thisSlideEvent = coincSlideHead;
  CoincInspiralTable    *thisEvent = coincZeroHead;

  REAL4  thisStat;
  REAL4  thisSlideStat;
  REAL4  thisRate = 0.;
  REAL4  loudestRate = -1;

  while ( thisEvent )
  {
    thisStat = XLALCoincInspiralStat(thisEvent, coincStat, bittenLParams);

    if ( (fitStat > 0) && (thisStat >= fitStat) )
    {
      /* put thisRate in sngl_inspiral->alpha for thisEvent[ifo]
       * for all ifos in thisEvent
       * FIXME in the future */
      InterferometerNumber ifoNumber = LAL_UNKNOWN_IFO;

      thisRate = fitA * exp(fitB * thisStat);

      for ( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( thisEvent->snglInspiral[ifoNumber] )
          thisEvent->snglInspiral[ifoNumber]->alpha = thisRate/timeAnalyzed;
      }
    }
    else
    {
      while ( thisSlideEvent )
      {
        thisSlideStat =
            XLALCoincInspiralStat(thisSlideEvent, coincStat, bittenLParams);

        if ( thisSlideStat >= thisStat )
        {
          /* count this slide coinc towards the rate */
          thisRate += 1.;
          thisSlideEvent = thisSlideEvent->next;
        }
        else
        {
          /* no more slide triggers with stat>=thisStat so break */
          break;
        }
      }
    {
      /* put thisRate in sngl_inspiral->alpha for thisEvent[ifo]
       * for all ifos in thisEvent
       * FIXME in the future */
      InterferometerNumber ifoNumber = LAL_UNKNOWN_IFO;

      for ( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( thisEvent->snglInspiral[ifoNumber] )
          thisEvent->snglInspiral[ifoNumber]->alpha = thisRate;
      }
    }
    }

    if ( loudestRate < 0 ) loudestRate = thisRate/timeAnalyzed;
    thisEvent = thisEvent->next;
  }

  return( loudestRate );
}



REAL4
XLALRateErrorCalcCoincInspiral (
    CoincInspiralTable         *coincZeroHead,
    CoincInspiralSlideTable    *slideHeads,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    int                         numSlides,
    REAL4                       timeAnalyzed,
    REAL4                       fitStat,
    REAL4                       fitA,
    REAL4                       fitB
    )

{
  CoincInspiralSlideTable    *thisSlideHead = NULL;
  CoincInspiralSlideTable    *thisHeadSlideHead = NULL;
  CoincInspiralSlideTable    *tmpSlideHead = NULL;
  CoincInspiralTable         *thisEvent = coincZeroHead;

  REAL4  thisStat;
  REAL4  thisSlideStat;
  REAL4  thisRate = 0.;
  REAL4  thisRateNum = 0.;
  REAL4  thisRateDenom = 0.;
  REAL4  thisRateError = 0.;
  REAL4  loudestRate = -1;

  XLALCreateCoincSlideTable( &thisHeadSlideHead, numSlides );
  thisSlideHead = thisHeadSlideHead;

  while ( thisSlideHead && slideHeads )
  {
    thisSlideHead->coincInspiral = slideHeads->coincInspiral;
    thisSlideHead->slideTimeAnalyzed = slideHeads->slideTimeAnalyzed;
    thisSlideHead->slideNum = slideHeads->slideNum;
    thisSlideHead = thisSlideHead->next;
    slideHeads = slideHeads->next;
  }

  while ( thisEvent )
  {
    thisStat = XLALCoincInspiralStat(thisEvent, coincStat, bittenLParams);

    if ( (fitStat > 0) && (thisStat >= fitStat) )
    {
      /* put thisRate in sngl_inspiral->alpha for thisEvent[ifo]
       * for all ifos in thisEvent
       * FIXME in the future */
      InterferometerNumber ifoNumber = LAL_UNKNOWN_IFO;

      thisRate = fitA * exp(fitB * thisStat);
      thisRateError = pow(thisRate, 0.5);
      for ( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( thisEvent->snglInspiral[ifoNumber] )
        {
          thisEvent->snglInspiral[ifoNumber]->alpha = thisRate;
          thisEvent->snglInspiral[ifoNumber]->alpha1 =
              thisRateError;
        }
      }
    }
    else
    {
      thisRate = 0.;
      thisRateNum = 0.;
      thisRateDenom = 0.;
      thisRateError = 0.;

      thisSlideHead = thisHeadSlideHead;
      while ( thisSlideHead )
      {
        while ( thisSlideHead->coincInspiral )
        {
          thisSlideStat = XLALCoincInspiralStat(thisSlideHead->coincInspiral,
              coincStat, bittenLParams);

          if ( thisSlideStat >= thisStat )
          {
            /* count this slide coinc towards the rate */
            thisSlideHead->currentRate += 1;
            thisSlideHead->coincInspiral = thisSlideHead->coincInspiral->next;
          }
          else
          {
            /* no more slide triggers in this slide number with
             * stat>=thisStat so break */
            break;
          }
        }
        /* add this slide's current rate to thisRate */
        thisRateNum += thisSlideHead->currentRate;
        thisRateDenom += thisSlideHead->slideTimeAnalyzed;

        /* move on to the next slide number */
        thisSlideHead = thisSlideHead->next;
      }

      thisRate = timeAnalyzed * thisRateNum / thisRateDenom;

      thisSlideHead = thisHeadSlideHead;
      while ( thisSlideHead )
      {
        /* calculate error on thisRate */
        thisRateError += pow( timeAnalyzed * thisSlideHead->currentRate / \
            thisSlideHead->slideTimeAnalyzed
            - thisRate, 2.0 ) / \
            (thisRateDenom / thisSlideHead->slideTimeAnalyzed);

        /* move on to the next slide number */
        thisSlideHead = thisSlideHead->next;
      }
      thisRateError = pow(thisRateError, 0.5);

      {
        /* put thisRate in sngl_inspiral->alpha for thisEvent[ifo]
         * for all ifos in thisEvent
         * FIXME in the future */
        InterferometerNumber ifoNumber = LAL_UNKNOWN_IFO;

        for ( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
        {
          if ( thisEvent->snglInspiral[ifoNumber] )
          {
            thisEvent->snglInspiral[ifoNumber]->alpha = thisRate;
            thisEvent->snglInspiral[ifoNumber]->alpha1 = thisRateError;
          }
        }
      }
    }

    if ( loudestRate < 0 ) loudestRate = thisRate;
    thisEvent = thisEvent->next;
  }

  if ( thisEvent )
  {
    /* should never get here */
    XLALPrintError( "Have events where FAR not calculated" );
    XLAL_ERROR(XLAL_EIO);
  }

  /* free the CoincInspiralSlideTable thisSlideHead */
  thisSlideHead = thisHeadSlideHead;
  while ( thisSlideHead )
  {
    tmpSlideHead = thisSlideHead;
    thisSlideHead = thisSlideHead->next;
    LALFree( tmpSlideHead );
  }

  return( loudestRate );
}


void
XLALPopulateAccuracyParams(
       InspiralAccuracyList  *accuracyParams
)

{
  INT4 ifoNumber, ifoTwo;
  LALDetector aDet, bDet;


  /* check that the accuracyParams structure is allocated */
  if ( accuracyParams == NULL )
  {
    XLAL_ERROR_VOID( XLAL_EFAULT );
  }

  /* Populate the lightTravel matrix */
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    XLALReturnDetector( &aDet, (InterferometerNumber) ifoNumber );

    for ( ifoTwo = 0; ifoTwo < LAL_NUM_IFO; ifoTwo++)
    {
      XLALReturnDetector( &bDet, (InterferometerNumber) ifoTwo );

      /* compute maximum light travel time */
      accuracyParams->lightTravelTime[ ifoNumber][ ifoTwo ] =
	XLALLightTravelTime( &aDet, &bDet );
    }
  }

  /* set the exttrig flag */
  accuracyParams->exttrig=0;
}




void
XLALPopulateAccuracyParamsExt(
       InspiralAccuracyList  *accuracyParams,
       const LIGOTimeGPS     *gpstime,
       const REAL8            ra_deg,
       const REAL8            dec_deg
)

{
  INT4 ifoNumber, ifoTwo;
  REAL8 timeDelay;
  REAL8 ra_radians, dec_radians;
  LALDetector aDet, bDet;

  /* check that the accuracyParams structure is allocated */
  if ( accuracyParams == NULL )
  {
    XLAL_ERROR_VOID( XLAL_EFAULT );
  }

  /* check the values given */
  if (ra_deg<0 || ra_deg > 360)
  {
    XLALPrintError("Right ascension value outside [0; 360]. Value given: %f\n", ra_deg);
    XLAL_ERROR_VOID( XLAL_EDATA );
  }
  if (dec_deg<-90 || dec_deg>90)
  {
    XLALPrintError("Declination value outside [-90; 90]. Value given: %f\n", dec_deg);
    XLAL_ERROR_VOID( XLAL_EDATA );
  }

  /* convert position */
  ra_radians  =  ra_deg * LAL_PI_180;
  dec_radians = dec_deg * LAL_PI_180;


  /* Populate the lightTravel matrix */
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    XLALReturnDetector( &aDet, (InterferometerNumber) ifoNumber );

    for ( ifoTwo = 0; ifoTwo < LAL_NUM_IFO; ifoTwo++)
    {
      XLALReturnDetector( &bDet, (InterferometerNumber) ifoTwo );

      /* compute signal travel time  */
      timeDelay=-XLALArrivalTimeDiff( aDet.location, bDet.location,
				      ra_radians, dec_radians, gpstime);

      accuracyParams->lightTravelTime[ ifoNumber][ ifoTwo ] =
	(INT8) 1e9*timeDelay;
    }
  }

  /* set the exttrig flag */
  accuracyParams->exttrig=1;

}

/*@}*/ /* end:CoincInspiralUtils_c */
