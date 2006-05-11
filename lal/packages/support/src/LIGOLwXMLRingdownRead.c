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
    XLAL_ERROR_NULL( func, XLAL_EIO );
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
        thisEvent->percent_mass_loss = r4colData;
      }
      else if ( tableDir[j].idx == 11 )
      {
        thisEvent->amplitude = r4colData;
      }
      else if ( tableDir[j].idx == 12 )
      {
        thisEvent->snr = r4colData;
      }
      else if ( tableDir[j].idx == 13)
      {
        thisEvent->eff_dist = r4colData;
      }
      else if ( tableDir[j].idx == 14 )
      {
        thisEvent->sigma_sq = r8colData;
      }
      else if ( tableDir[j].idx == 15 )
      {
        thisEvent->event_id = (EventIDColumn *) 
          LALCalloc( 1, sizeof(EventIDColumn) );
        thisEvent->event_id->id = i8colData;
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
        else if ( tableDir[j].idx == 10 )
        {
          thisEvent->longitude = r4colData;
        }
        else if ( tableDir[j].idx == 11 )
        {
          thisEvent->latitude = r4colData;
        }
        else if ( tableDir[j].idx == 13 )
        {
          thisEvent->distance = r4colData;
        }
        else if ( tableDir[j].idx == 14 )
        {
          thisEvent->inclination = r4colData;
        }
        else if ( tableDir[j].idx == 15 )
        {
          thisEvent->polarization = r4colData;
        }
        else if ( tableDir[j].idx == 18 )
        {
          thisEvent->frequency = r4colData;
        }
        else if ( tableDir[j].idx == 19 )
        {
          thisEvent->quality = r4colData;
        }
        else if ( tableDir[j].idx == 12 )
        {
          thisEvent->phase = r4colData;
        }
        else if ( tableDir[j].idx == 9 )
        {
          thisEvent->mass = r4colData;
        }
        else if ( tableDir[j].idx == 17 )
        {
          thisEvent->spin = r4colData;
        }
        else if ( tableDir[j].idx == 16 )
        {
          thisEvent->percent_mass_loss = r4colData;
        }
        else if ( tableDir[j].idx == 22 )
        {
          thisEvent->amplitude = r4colData;
        }
        else if ( tableDir[j].idx == 20 )
        {
          thisEvent->eff_dist_h = r4colData;
        }
        else if ( tableDir[j].idx == 21 )
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
  inputData = XLALSnglRingdownTableFromLIGOLw( fileName );
  if ( ! inputData )
  {
    XLALPrintError("Unable to read sngl_ringdown table from %s\n", fileName );
    LALFree(thisInputFile);
    XLAL_ERROR(func, XLAL_EIO);
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
