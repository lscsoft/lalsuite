/*
*  Copyright (C) 2007 Drew Keppel, Duncan Brown, Jolien Creighton, Patrick Brady, Stephen Fairhurst, Thomas Cokelaer
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
 * File Name: LIGOMetadataUtils.c
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown, D. A.
 * \file
 * \ingroup lalmetaio_general
 * \brief General routines for manipulating LIGO metadatabase tables.
 *
 * ### Description ###
 *
 * <tt>LALCheckOutTimeFromSearchSummary()</tt> verifies that all times
 * between the specified \c startTime and \c endTime have been
 * searched precisely once for the given \c ifo.
 *
 * ### Algorithm ###
 *
 * None.
 *
 * ### Uses ###
 *
 * LALCalloc, LALMalloc, LALFree.
 *
 * ### Notes ###
 *
 * %% Any relevant notes.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>

int XLALCountProcessTable(ProcessTable *head)

{
	int length;

	/* count the number of events in the list */
	for(length = 0; head; head = head->next)
		length++;

	return(length);
}


int XLALCountProcessParamsTable(ProcessParamsTable *head)

{
	int length;

	/* count the number of events in the list */
	for(length = 0; head; head = head->next)
		length++;

	return(length);
}



int
XLALIFONumber(
    const char *ifo
    )

{
  switch( ifo[0] )
  {
    case 'G':
      return LAL_IFO_G1;
      break;

    case 'H':
      if ( !strcmp( ifo, "H1" ) )
      {
        return LAL_IFO_H1;
      }
      else if (!strcmp( ifo, "H2" ) )
      {
        return LAL_IFO_H2;
      }
      else
      {
        /* Invalid Hanford Detector */
        return LAL_UNKNOWN_IFO ;
      }
      break;

    case 'L':
      return LAL_IFO_L1;
      break;

    case 'T':
      return LAL_IFO_T1;
      break;

    case 'V':
      return LAL_IFO_V1;
      break;

    default:
      /* Invalid Detector Site */
      return LAL_UNKNOWN_IFO ;
  }
}


void
XLALReturnIFO(
    char                *ifo,
    InterferometerNumber IFONumber
    )

{
  switch( IFONumber )
  {
    case LAL_IFO_G1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "G1");
      break;

    case LAL_IFO_H1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "H1");
      break;

    case LAL_IFO_H2:
      snprintf( ifo, LIGOMETA_IFO_MAX, "H2");
      break;

    case LAL_IFO_L1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "L1");
      break;

    case LAL_IFO_T1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "T1");
      break;

    case LAL_IFO_V1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "V1");
      break;

    default:
      /* Invalid Detector Site */
      snprintf( ifo, LIGOMETA_IFO_MAX, " ");
  }
}



void
XLALReturnDetector(
    LALDetector           *det,
    InterferometerNumber   IFONumber
    )

{
  switch( IFONumber )
  {
    case LAL_IFO_G1:
      *det = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
      break;

    case LAL_IFO_H1:
      *det = lalCachedDetectors[LALDetectorIndexLHODIFF];
      break;

    case LAL_IFO_H2:
      *det = lalCachedDetectors[LALDetectorIndexLHODIFF];
      break;

    case LAL_IFO_L1:
      *det = lalCachedDetectors[LALDetectorIndexLLODIFF];
      break;

    case LAL_IFO_T1:
      *det = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
      break;

    case LAL_IFO_V1:
      *det = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
      break;

    default:
      /* Invalid Detector Site */
      memset(det, 0, sizeof(LALDetector) );
  }
}


int
XLALCompareSearchSummaryByOutTime (
    const void *a,
    const void *b
    )

{
  const SearchSummaryTable *aPtr = *((const SearchSummaryTable * const *)a);
  const SearchSummaryTable *bPtr = *((const SearchSummaryTable * const *)b);

  INT8 ta = 0;
  INT8 tb = 0;

  /* determine the out start times */
  ta = XLALGPSToINT8NS( &(aPtr->out_start_time) );
  tb = XLALGPSToINT8NS( &(bPtr->out_start_time) );

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
    /* determine the out end times */
    ta = XLALGPSToINT8NS( &(aPtr->out_end_time) );
    tb = XLALGPSToINT8NS( &(bPtr->out_end_time) );

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
}


int
XLALTimeSortSearchSummary(
    SearchSummaryTable  **summHead,
    int(*comparfunc)    (const void *, const void *)
    )

{
  INT4                  i;
  INT4                  numSumms = 0;
  SearchSummaryTable    *thisSearchSumm = NULL;
  SearchSummaryTable   **summHandle = NULL;

  if ( !summHead )
  {
    XLAL_ERROR(XLAL_EIO);
  }


  /* count the number of summs in the linked list */
  for ( thisSearchSumm = *summHead; thisSearchSumm;
      thisSearchSumm = thisSearchSumm->next )
  {
    ++numSumms;
  }
  if ( ! numSumms )
  {
    return( 0 );
  }

  /* allocate memory for an array of ptrs to sort and populate array */
  summHandle = (SearchSummaryTable **)
    LALCalloc( numSumms, sizeof(SearchSummaryTable *) );
  for ( i = 0, thisSearchSumm = *summHead; i < numSumms;
      ++i, thisSearchSumm = thisSearchSumm->next )
  {
    summHandle[i] = thisSearchSumm;
  }

  /* qsort the array using the specified function */
  qsort( summHandle, numSumms, sizeof(summHandle[0]), comparfunc );

  /* re-link the linked list in the right order */
  thisSearchSumm = *summHead = summHandle[0];
  for ( i = 1; i < numSumms; ++i, thisSearchSumm = thisSearchSumm->next )
  {
    thisSearchSumm->next = summHandle[i];
  }
  thisSearchSumm->next = NULL;

  /* free the internal memory */
  LALFree( summHandle );

  return( 0 );
}


SearchSummaryTable *
XLALIfoScanSearchSummary(
    SearchSummaryTable         *input,
    CHAR                       *ifos
    )

{
  SearchSummaryTable    *output = NULL;
  SearchSummaryTable    *thisSearchSumm = NULL;
  SearchSummaryTable    *keptSumm = NULL;


  if ( !input )
  {
    XLAL_ERROR_NULL(XLAL_EIO);
  }

  /* Scan through a linked list of search_summary tables and return a
     pointer to the head of a linked list of tables for a specific IFO */

  for( thisSearchSumm = input; thisSearchSumm;
      thisSearchSumm = thisSearchSumm->next )
  {

    if ( !strcmp(thisSearchSumm->ifos, ifos) )
    {
      /* IFOs match so write this entry to the output table */
      if ( ! output  )
      {
        output = keptSumm = (SearchSummaryTable *)
          LALMalloc( sizeof(SearchSummaryTable) );
      }
      else
      {
        keptSumm = keptSumm->next = (SearchSummaryTable *)
          LALMalloc( sizeof(SearchSummaryTable) );
      }
      if ( !keptSumm )
      {
        while ( output )
        {
          thisSearchSumm = output;
          output = (output)->next;
          LALFree( thisSearchSumm );
        }
        XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      memcpy(keptSumm, thisSearchSumm, sizeof(SearchSummaryTable));
      keptSumm->next = NULL;
    }
  }
  return( output);
}


/**
 * Create a ProcessTable structure.
 */


ProcessTable *XLALCreateProcessTableRow(void)
{
  ProcessTable *new = XLALMalloc(sizeof(*new));

  if(!new)
    XLAL_ERROR_NULL(XLAL_EFUNC);

  new->next = NULL;
  memset(new->program, 0, sizeof(new->program));
  memset(new->version, 0, sizeof(new->version));
  memset(new->cvs_repository, 0, sizeof(new->cvs_repository));
  XLALGPSSet(&new->cvs_entry_time, 0, 0);
  memset(new->comment, 0, sizeof(new->comment));
  new->is_online = 0;
  memset(new->node, 0, sizeof(new->node));
  memset(new->username, 0, sizeof(new->username));
  XLALGPSSet(&new->start_time, 0, 0);
  XLALGPSSet(&new->end_time, 0, 0);
  new->jobid = 0;
  memset(new->domain, 0, sizeof(new->domain));
  new->unix_procid = 0;
  memset(new->ifos, 0, sizeof(new->ifos));
  new->process_id = -1;

  return new;
}


/**
 * Destroy a ProcessTable structure.
 */


void XLALDestroyProcessTableRow(ProcessTable *row)
{
  XLALFree(row);
}


/**
 * Destroy a ProcessTable linked list.
 */


void XLALDestroyProcessTable(ProcessTable *head)
{
  while(head)
  {
    ProcessTable *next = head->next;
    XLALDestroyProcessTableRow(head);
    head = next;
  }
}


/**
 * Return the next available process ID.
 */


long XLALProcessTableGetNextID(ProcessTable *head)
{
  long highest = -1;
  for(; head; head = head->next)
    if(head->process_id > highest)
      highest = head->process_id;
  return highest + 1;
}


/**
 * Create a ProcessParamsTable structure.
 */


ProcessParamsTable *XLALCreateProcessParamsTableRow(const ProcessTable *process)
{
  ProcessParamsTable *new = XLALMalloc(sizeof(*new));

  if(!new)
    XLAL_ERROR_NULL(XLAL_EFUNC);

  new->next = NULL;
  memset(new->program, 0, sizeof(new->program));
  if(process)
    new->process_id = process->process_id;
  else
    new->process_id = -1;
  memset(new->param, 0, sizeof(new->param));
  memset(new->type, 0, sizeof(new->type));
  memset(new->value, 0, sizeof(new->value));

  return new;
}


/**
 * Destroy a ProcessParamsTable structure.
 */


void XLALDestroyProcessParamsTableRow(ProcessParamsTable *row)
{
  XLALFree(row);
}


/**
 * Destroy a ProcessParamsTable linked list.
 */


void XLALDestroyProcessParamsTable(ProcessParamsTable *head)
{
  while(head)
  {
    ProcessParamsTable *next = head->next;
    XLALDestroyProcessParamsTableRow(head);
    head = next;
  }
}


/**
 * Create a TimeSlide structure.
 */


TimeSlide *XLALCreateTimeSlide(void)
{
  TimeSlide *new = XLALMalloc(sizeof(*new));

  if(!new)
    XLAL_ERROR_NULL(XLAL_EFUNC);

  new->next = NULL;
  new->process_id = -1;
  new->time_slide_id = -1;
  memset(new->instrument, 0, sizeof(new->instrument));
  new->offset = 0;

  return new;
}


/**
 * Destroy a TimeSlide structure.
 */


void XLALDestroyTimeSlide(TimeSlide *row)
{
  XLALFree(row);
}


/**
 * Destroy a TimeSlide linked list.
 */


void XLALDestroyTimeSlideTable(TimeSlide *head)
{
  while(head)
  {
    TimeSlide *next = head->next;
    XLALDestroyTimeSlide(head);
    head = next;
  }
}


/**
 * Find and return the address of the first element in the linked list of
 * TimeSlide objects whose time_slide_id and instrument name equal the
 * values given.  TimeSlide elements whose instrument pointer is NULL are
 * skipped.  Returns NULL if no matching row is found.  There are two
 * versions, one for cost * TimeSlide rows and one for non-const (neither
 * modifies the TimeSlide rows, the two versions are identical, they are
 * provided to allow the const'ness to be "passed" through the function).
 */


const TimeSlide *XLALTimeSlideConstGetByIDAndInstrument(const TimeSlide *time_slide, long time_slide_id, const char *instrument)
{
	for(; time_slide && (time_slide->time_slide_id != time_slide_id || strcmp(time_slide->instrument, instrument)); time_slide = time_slide->next);
	return time_slide;
}


TimeSlide *XLALTimeSlideGetByIDAndInstrument(TimeSlide *time_slide, long time_slide_id, const char *instrument)
{
	for(; time_slide && (time_slide->time_slide_id != time_slide_id || strcmp(time_slide->instrument, instrument)); time_slide = time_slide->next);
	return time_slide;
}


/**
 * Create a SearchSummaryTable structure.
 */


SearchSummaryTable *XLALCreateSearchSummaryTableRow(const ProcessTable *process)
{
  SearchSummaryTable *new = XLALMalloc(sizeof(*new));

  if(!new)
    XLAL_ERROR_NULL(XLAL_EFUNC);

  new->next = NULL;
  if(process)
    new->process_id = process->process_id;
  else
    new->process_id = -1;
  memset(new->comment, 0, sizeof(new->comment));
  XLALGPSSet(&new->in_start_time, 0, 0);
  XLALGPSSet(&new->in_end_time, 0, 0);
  XLALGPSSet(&new->out_start_time, 0, 0);
  XLALGPSSet(&new->out_end_time, 0, 0);
  new->nevents = -1;
  new->nnodes = -1;
  memset(new->ifos, 0, sizeof(new->ifos));

  return new;
}


/**
 * Destroy a SearchSummaryTable structure.
 */


void XLALDestroySearchSummaryTableRow(SearchSummaryTable *row)
{
  XLALFree(row);
}


/**
 * Destroy a SearchSummaryTable linked list.
 */


void XLALDestroySearchSummaryTable(SearchSummaryTable *head)
{
  while(head)
  {
    SearchSummaryTable *next = head->next;
    XLALDestroySearchSummaryTableRow(head);
    head = next;
  }
}

/*
 * Create a SegmentTable structure.
 */


SegmentTable *XLALCreateSegmentTableRow(const ProcessTable *process)
{
  SegmentTable *new = XLALMalloc(sizeof(*new));

  if(!new)
    XLAL_ERROR_NULL(XLAL_EFUNC);

  new->next = NULL;
  if(process)
    new->process_id = process->process_id;
  else
    new->process_id = -1;
  new->creator_db = -1;
  new->segment_id = -1;
  new->segment_def_id = -1;
  new->segment_def_cdb = -1;
  XLALGPSSet(&new->start_time, 0, 0);
  XLALGPSSet(&new->end_time, 0, 0);

  return new;
}


/*
 * Destroy a SegmentTable structure.
 */


void XLALDestroySegmentTableRow(SegmentTable *row)
{
  XLALFree(row);
}


/*
 * Destroy a SegmentTable linked list.
 */


void XLALDestroySegmentTable(SegmentTable *head)
{
  while(head)
  {
    SegmentTable *next = head->next;
    XLALDestroySegmentTableRow(head);
    head = next;
  }
}
