/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Lisa M. Goggin
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
 * File Name: LIGOLwXMLRingdownRead.c
 *
 * Author: Brown, D. A., and Goggin, L. M.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="LIGOLwXMLRingdownReadCV">
Author: Brown, D. A. and Goggin, L. M.
$Id$
</lalVerbatim> 
#endif

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/StringInput.h>

NRCSID( LIGOLWXMLRINGDOWNREADC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{LIGOLwXMLRingdownRead.c}}

Routines to read the various ringdown search XML data into LAL structures.

\subsubsection*{Prototypes}
\input{LIGOLwXMLRingdownReadCP}
\idx{XLALSnglRingdownTableFromLIGOLw()}
\idx{XLALSimRingdownTableFromLIGOLw()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

None.

\subsubsection*{Uses}
Functions in the Metaio library:
\begin{itemize}
\item \verb+MetaioFindColumn+
\item \verb+MetaioGetRow+
\item \verb+MetaioOpenTable+
\item \verb+MetaioClose+
\end{itemize}
\subsubsection*{Notes}

%% Any relevant notes.

\vfill{\footnotesize\input{LIGOLwXMLRingdownReadCV}}

</lalLaTeX>
#endif


#define XLAL_CLOBBER_EVENTS \
  while ( eventHead ); \
{ \
  thisEvent = eventHead; \
  eventHead = (eventHead)->next; \
  LALFree( thisEvent ); \
  thisEvent = NULL; \
}

/* <lalVerbatim file="LIGOLwXMLRingdownReadCP"> */
SnglRingdownTable* XLALSnglRingdownTableFromLIGOLw (
    CHAR               *fileName
    )
/* </lalVerbatim> */
{
  static const char *func = "XLALSnglRingdownTableFromLIGOLw";
  int                            i, j;
  int                            mioStatus = 0;
  INT4                           i4colData;
  INT8                           i8colData;
  REAL4                          r4colData;
  REAL8                          r8colData;
  SnglRingdownTable             *eventHead = NULL;
  SnglRingdownTable             *thisEvent = NULL;
  struct MetaioParseEnvironment  parseEnv;
  const  MetaioParseEnv          env = &parseEnv;
  MetaTableDirectory            *tableDir = NULL;

  /* open the sngl_ringdown XML file */
  mioStatus = MetaioOpenTable( env, fileName, "sngl_ringdown" );
  if ( mioStatus )
  {
    XLALPrintError( "XLAL Error - unable to open sngl_ringdown table: "
        "metaio error code %d\n", mioStatus );
    XLAL_ERROR_NULL( func, XLAL_EDATA ); 
   /* return 0;*/
  }

  /* create table directory to find columns in file */
  tableDir = XLALCreateMetaTableDir(env, sngl_ringdown_table);
  if ( ! tableDir )
  {
    XLALPrintError( "XLAL Error - "
        "unable to create sngl_ringdown table directory\n" );
    XLAL_ERROR_NULL( func, XLAL_EIO );
  }

  /* loop over the rows in the file */
  i = 0;
  while ( (mioStatus = MetaioGetRow(env)) == 1 ) 
  {
    /* count the rows in the file */
    i++;

    /* allocate memory for the ringdown structure we are about to read in */
    if ( ! eventHead )
    {
      thisEvent = eventHead = (SnglRingdownTable *) 
        LALCalloc( 1, sizeof(SnglRingdownTable) );
    }
    else
    {
      thisEvent = thisEvent->next = (SnglRingdownTable *) 
        LALCalloc( 1, sizeof(SnglRingdownTable) );
    }
    if ( ! thisEvent )
    {
      XLALPrintError( "XLAL Error - could not allocate sngl_ringdown table\n" );
      XLAL_CLOBBER_EVENTS;
      MetaioClose( env );
      XLAL_ERROR_NULL( func, XLAL_ENOMEM );
    }

    /* parse the contents of the row into the SnglRingdownTable structure */
    for ( j = 0; tableDir[j].name; ++j )
    {
      if ( tableDir[j].pos < 0 )
      {
        XLALPrintError( "XLAL Error - bad table directory for element %d\n", j );
        XLAL_CLOBBER_EVENTS;
        XLAL_ERROR_NULL( func, XLAL_EIO );
      }

      /* dereference the data stored in the table */
      r4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
      r8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_8; 
      i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;
      i8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_8s;

      if ( tableDir[j].idx == 0 )
      {
        LALSnprintf( thisEvent->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), 
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 1 )
      {
        LALSnprintf( thisEvent->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 2 )
      {
        thisEvent->start_time.gpsSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 3 )
      {
        thisEvent->start_time.gpsNanoSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 4 )
      {
        thisEvent->start_time_gmst = r8colData;
      }
      else if ( tableDir[j].idx == 5 )
      {
        thisEvent->frequency = r4colData;
      }
      else if ( tableDir[j].idx == 6 )
      {
        thisEvent->quality = r4colData;
      }
      else if ( tableDir[j].idx == 7 )
      {
        thisEvent->phase = r4colData;
      }
      else if ( tableDir[j].idx == 8 )
      {
        thisEvent->mass = r4colData;
      }
      else if ( tableDir[j].idx == 9 )
      {
        thisEvent->spin = r4colData;
      }
      else if ( tableDir[j].idx == 10 )
      {
        thisEvent->epsilon = r4colData;
      }
      else if ( tableDir[j].idx == 11 )
      {
        thisEvent->num_clust_trigs = i4colData;
      }
      else if ( tableDir[j].idx == 12 )
      {
        thisEvent->ds2_H1H2 = r4colData;
      }
      else if ( tableDir[j].idx == 13 )
      {
        thisEvent->ds2_H1L1 = r4colData;
      }
      else if ( tableDir[j].idx == 14 )
      {
        thisEvent->ds2_H2L1 = r4colData;
      }
      else if ( tableDir[j].idx == 15 )
      {
        thisEvent->amplitude = r4colData;
      }
      else if ( tableDir[j].idx == 16 )
      {
        thisEvent->snr = r4colData;
      }
      else if ( tableDir[j].idx == 17)
      {
        thisEvent->eff_dist = r4colData;
      }
      else if ( tableDir[j].idx == 18 )
      {
        thisEvent->sigma_sq = r8colData;
      }
      else if ( tableDir[j].idx == 19 )
      {
        /* JC: AVOID BUG IN METAIO -- BAD */
        union { const char *cs; const unsigned char *cus; } bad;
        thisEvent->event_id = (EventIDColumn *) 
          LALCalloc( 1, sizeof(EventIDColumn) );
        bad.cus = env->ligo_lw.table.elt[tableDir[j].pos].data.ilwd_char.data;
        sscanf( bad.cs, "sngl_ringdown:event_id:%" LAL_UINT8_FORMAT, &(thisEvent->event_id->id) );
        thisEvent->event_id->snglRingdownTable = thisEvent;
      }
      else
      {
        XLALPrintError( "XLAL Error - "
            "table directory index %d out of bounds\n", j );
        XLAL_CLOBBER_EVENTS;
        XLAL_ERROR_NULL( func, XLAL_EIO );
      }
    }
  }

  if ( mioStatus == -1 )
  {
    XLALPrintError( "XLAL Error - error parsing after row %d\n", i );
    XLAL_CLOBBER_EVENTS;
    MetaioClose( env );
    XLAL_ERROR_NULL( func, XLAL_EIO );
  }

  /* Normal exit */
  LALFree( tableDir );
  MetaioClose( env );

  return eventHead;
}


/* <lalVerbatim file="LIGOLwXMLRingdownReadCP"> */
SimRingdownTable* XLALSimRingdownTableFromLIGOLw (
    CHAR               *fileName,
    INT4                startTime,
    INT4                stopTime
    )
/* </lalVerbatim> */
{
  static const char *func = "XLALSimRingdownTableFromLIGOLw";
  int                            i, j;
  int                            mioStatus = 0;
  INT4                           i4colData;
  REAL4                          r4colData;
  REAL8                          r8colData;
  SimRingdownTable              *eventHead = NULL;
  SimRingdownTable              *thisEvent = NULL;
  struct MetaioParseEnvironment  parseEnv;
  const  MetaioParseEnv          env = &parseEnv;
  MetaTableDirectory            *tableDir = NULL;


  /* open the sim_ringdown XML file */
  mioStatus = MetaioOpenTable( env, fileName, "sim_ringdown" );
  if ( mioStatus )
  {
    XLALPrintError( "XLAL Error - unable to open sim_ringdown table: "
        "metaio error code %d\n", mioStatus );
    XLAL_ERROR_NULL( func, XLAL_EIO );
  }

  /* create table directory to find columns in file */
  tableDir = XLALCreateMetaTableDir(env, sim_ringdown_table);
  if ( ! tableDir )
  {
    XLALPrintError( "XLAL Error - "
        "unable to create sim_ringdown table directory\n" );
    XLAL_ERROR_NULL( func, XLAL_EIO );
  }

  /* loop over the rows in the file */
  i = 0;
  while ( (mioStatus = MetaioGetRow(env)) == 1 ) 
  {
    /* count the rows in the file */
    i++;

    /* check the injection time is withing the requested inteval */
    if ( tableDir[2].pos < 0 )
    {
      XLALPrintError( "XLAL Error - bad table directory for element %d\n", j );
      XLAL_CLOBBER_EVENTS;
      XLAL_ERROR_NULL( func, XLAL_EIO );
    }

    i4colData = env->ligo_lw.table.elt[tableDir[2].pos].data.int_4s;

    if ( ! stopTime || ( i4colData >= startTime && i4colData < stopTime ) )
    {
      /* allocate memory for the template we are about to read in */
      if ( ! eventHead )
      {
        thisEvent = eventHead = (SimRingdownTable *) 
          LALCalloc( 1, sizeof(SimRingdownTable) );
      }
      else
      {
        thisEvent = thisEvent->next = (SimRingdownTable *) 
          LALCalloc( 1, sizeof(SimRingdownTable) );
      }
      if ( ! thisEvent )
      {
        XLALPrintError( "XLAL Error - could not allocate sim_ringdown table\n" );
        XLAL_CLOBBER_EVENTS;
        MetaioClose( env );
        XLAL_ERROR_NULL( func, XLAL_ENOMEM );
      }

      /* parse the contents of the row into the SimRingdownTable structure */
      for ( j = 0; tableDir[j].name; ++j )
      {
        if ( tableDir[j].pos < 0 )
        {
          XLALPrintError( "XLAL Error - bad table directory for element %d\n", j );
          XLAL_CLOBBER_EVENTS;
          XLAL_ERROR_NULL( func, XLAL_EIO );
        }

        /* dereference the data stored in the table */
        i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;
        r4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
        r8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_8;

        if ( tableDir[j].idx == 0 )
        {
          LALSnprintf( thisEvent->waveform, 
              LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s", 
              env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 1 )
        {
          LALSnprintf( thisEvent->coordinates, 
              LIGOMETA_COORDINATES_MAX * sizeof(CHAR), 
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 2 )
        {
          thisEvent->geocent_start_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 3 )
        {
          thisEvent->geocent_start_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 4 )
        {
          thisEvent->h_start_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 5 )
        {
          thisEvent->h_start_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 6 )
        {
          thisEvent->l_start_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 7 )
        {
          thisEvent->l_start_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 8 )
        {
          thisEvent->start_time_gmst = r8colData;
        }
        else if ( tableDir[j].idx == 9 )
        {
          thisEvent->longitude = r4colData;
        }
        else if ( tableDir[j].idx == 10 )
        {
          thisEvent->latitude = r4colData;
        }
        else if ( tableDir[j].idx == 11 )
        {
          thisEvent->distance = r4colData;
        }
        else if ( tableDir[j].idx == 12 )
        {
          thisEvent->inclination = r4colData;
        }
        else if ( tableDir[j].idx == 13 )
        {
          thisEvent->polarization = r4colData;
        }
        else if ( tableDir[j].idx == 14 )
        {
          thisEvent->frequency = r4colData;
        }
        else if ( tableDir[j].idx == 15 )
        {
          thisEvent->quality = r4colData;
        }
        else if ( tableDir[j].idx == 16 )
        {
          thisEvent->phase = r4colData;
        }
        else if ( tableDir[j].idx == 17 )
        {
          thisEvent->mass = r4colData;
        }
        else if ( tableDir[j].idx == 18 )
        {
          thisEvent->spin = r4colData;
        }
        else if ( tableDir[j].idx == 19 )
        {
          thisEvent->epsilon= r4colData;
        }
        else if ( tableDir[j].idx == 20 )
        {
          thisEvent->amplitude = r4colData;
        }
        else if ( tableDir[j].idx == 21 )
        {
          thisEvent->eff_dist_h = r4colData;
        }
        else if ( tableDir[j].idx == 22 )
        {
          thisEvent->eff_dist_l = r4colData;
        }
        else if ( tableDir[j].idx == 23 )
        {
          thisEvent->hrss = r4colData;
        }
        else if ( tableDir[j].idx == 24 )
        {
          thisEvent->hrss_h = r4colData;
        }
        else if ( tableDir[j].idx == 25 )
        {
          thisEvent->hrss_l = r4colData;
        }
        else
        {
          XLALPrintError( "XLAL Error - "
              "table directory index %d out of bounds\n", j );
          XLAL_CLOBBER_EVENTS;
          XLAL_ERROR_NULL( func, XLAL_EIO );
        }
      }
    }
  }

  if ( mioStatus == -1 )
  {
    XLALPrintError( "XLAL Error - error parsing after row %d\n", i );
    XLAL_CLOBBER_EVENTS;
    MetaioClose( env );
    XLAL_ERROR_NULL( func, XLAL_EIO);
  }

  /* Normal exit */
  LALFree( tableDir );
  MetaioClose( env );

  return eventHead;
}

#undef XLAL_CLOBBER_EVENTS

/* <lalVerbatim file="LIGOLwXMLRingdownReadCP"> */
INT4 XLALReadRingdownTriggerFile (
    SnglRingdownTable    **ringdownEventList,
    SnglRingdownTable    **lastTrigger,
    SearchSummaryTable   **searchSummList,
    SearchSummvarsTable  **inputFileList,
    CHAR                  *fileName
    )
/* </lalVerbatim> */
{
  const char *func = "XLALReadRingdownTriggerFile";
  INT4                 numFileTriggers = 0;
  int 		       errnum;
  SnglRingdownTable   *inputData = NULL;
  SearchSummaryTable  *inputSummary = NULL;
  SearchSummaryTable  *thisSearchSumm = NULL;
  SearchSummvarsTable *thisInputFile = NULL;
  
  /* store the file name in search summvars */
  XLALPrintInfo(
      "XLALReadRingdownTriggerFile(): storing input file name %s\n"
      "in search summvars table\n", fileName );
  
  if ( ! *inputFileList )
  {
    *inputFileList = thisInputFile = (SearchSummvarsTable *)
      LALCalloc( 1, sizeof(SearchSummvarsTable) );
  }
  else
  {
    for ( thisInputFile = *inputFileList; thisInputFile->next;
        thisInputFile = thisInputFile->next );
    thisInputFile = thisInputFile->next = (SearchSummvarsTable *)
      LALCalloc( 1, sizeof(SearchSummvarsTable) );
  }
  if ( ! thisInputFile )
  {
    XLALPrintError( "XLAL Error - could not allocate search_summvars table\n" );
    XLAL_ERROR( func, XLAL_ENOMEM );
  }
  
  LALSnprintf( thisInputFile->name, LIGOMETA_NAME_MAX, "input_file" );
  LALSnprintf( thisInputFile->string, LIGOMETA_NAME_MAX, "%s", fileName );
  
  /* read in the search summary and store */
  XLALPrintInfo( "XLALReadRingdownTriggerFile(): "
      "Reading search_summary table\n");
  
  inputSummary = XLALSearchSummaryTableFromLIGOLw( fileName );
  if ( ! inputSummary )
  {
    LALFree( thisInputFile );
    XLALPrintError( "XLAL Error - error reading search_summary table from %s\n", 
        fileName );
    XLAL_ERROR( func, XLAL_EIO );
  }
  else
  {
    /* store the search summary table in searchSummList list */
    if ( ! *searchSummList )
    {
      *searchSummList = thisSearchSumm = inputSummary;
    }
    else
    {
      for ( thisSearchSumm = *searchSummList; thisSearchSumm->next;
          thisSearchSumm = thisSearchSumm->next);
      thisSearchSumm = thisSearchSumm->next = inputSummary;
    }
  }
  
  /* read in the triggers */
  XLAL_TRY( inputData = XLALSnglRingdownTableFromLIGOLw( fileName ), errnum);
  if ( ! inputData )
    switch ( errnum )
    {
      case XLAL_EDATA:
        XLALPrintError("Unable to read sngl_ringdown table from %s\n", fileName );
        /*LALFree(thisInputFile);*/
        XLALClearErrno();
        break;
      default:
        XLALSetErrno( errnum );
        XLAL_ERROR( func, XLAL_EFUNC );
  }
  else
  {
    /* store the triggers */
    if ( ! *ringdownEventList )
    {
      /* store the head of the linked list */
      *ringdownEventList = *lastTrigger = inputData;
    }
    else
    {
      /* append to the end of the linked list and set current    */
      /* trigger to the first trigger of the list being appended */
      *lastTrigger = (*lastTrigger)->next = inputData;
    }
  
    /* scroll to the end of the linked list of triggers */
    for ( ; (*lastTrigger)->next; *lastTrigger = (*lastTrigger)->next )
      numFileTriggers++;
  }
  
  return numFileTriggers;
}
