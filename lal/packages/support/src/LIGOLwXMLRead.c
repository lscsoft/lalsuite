/*----------------------------------------------------------------------- 
 * 
 * File Name: ligolwxmlread.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="LIGOLwXMLReadCV">
Author: Brown, D. A. and Fairhurst, S.
$Id$
</lalVerbatim> 
#endif

#include <lal/LALStdio.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXMLHeaders.h>
#include <lal/LIGOLwXMLRead.h>

NRCSID( LIGOLWXMLREADC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{LIGOLwXMLRead.c}}

Routines to write LIGO metadata database structures to LIGO lightweight XML
files.

\subsubsection*{Prototypes}
\input{LIGOLwXMLReadCP}
\idx{LALCreateMetaTableDir()}
\idx{LALSnglBurstTableFromLIGOLw()}
\idx{LALSimBurstTableFromLIGOLw()}
\idx{LALSnglInspiralTableFromLIGOLw()}
\idx{InspiralTmpltBankFromLIGOLw()}
\idx{SimInspiralTableFromLIGOLw()}
\idx{SearchSummaryTableFromLIGOLw()}
\idx{SummValueTableFromLIGOLw()}

    
\subsubsection*{Description}


\subsubsection*{Algorithm}

None.

\subsubsection*{Uses}

\subsubsection*{Notes}
 
%% Any relevant notes.
 
\vfill{\footnotesize\input{LIGOLwXMLReadCV}}

</lalLaTeX>
#endif

/* <lalVerbatim file="LIGOLwXMLReadCP"> */
void
LALCreateMetaTableDir(
    LALStatus              *status,
    MetaTableDirectory    **tableDir,
    const MetaioParseEnv    env,
    MetadataTableType       table
    )
/* </lalVerbatim> */
{
  INT4 i;

  INITSTATUS( status, "LALCreateMetaTableDir", LIGOLWXMLREADC );
  ATTATCHSTATUSPTR (status);

  /* check the inputs */
  ASSERT( !*tableDir, status, LIGOLWXMLREADH_ENNUL, LIGOLWXMLREADH_MSGENNUL );

  switch( table )
  {
    case no_table:
      ABORT( status, LIGOLWXMLREADH_ENTAB, LIGOLWXMLREADH_MSGENTAB );
      break;
    case process_table:
      break;
    case process_params_table:
      break;
    case search_summary_table:
      break;
    case search_summvars_table:
      break;
    case sngl_burst_table:
      {
        MetaTableDirectory tmpTableDir[] =
        {
          {"ifo",                     -1, 0},
          {"search",                  -1, 1},
          {"channel",                 -1, 2},
          {"start_time",              -1, 3},
          {"start_time_ns",           -1, 4},
          {"duration",                -1, 5},
          {"central_freq",            -1, 6},
          {"bandwidth",               -1, 7},
          {"amplitude",               -1, 8},
          {"snr",                     -1, 9},
          {"confidence",              -1, 10},
          {"peak_time",               -1, 11},
          {"peak_time_ns",            -1, 12},
          {NULL,                       0, 0}
        };
        for ( i=0 ; tmpTableDir[i].name; ++i )
        {
          if ( (tmpTableDir[i].pos = 
                MetaioFindColumn( env, tmpTableDir[i].name )) < 0 )
          {
            fprintf( stderr, "unable to find column %s\n", 
                tmpTableDir[i].name );
            ABORT(status,LIGOLWXMLREADH_ENCOL,LIGOLWXMLREADH_MSGENCOL);
          }
        }

        *tableDir = (MetaTableDirectory *) LALMalloc( (i+1) * 
            sizeof(MetaTableDirectory)) ;
        memcpy(*tableDir, tmpTableDir, (i+1)*sizeof(MetaTableDirectory) );
      }
      break;
    case sngl_inspiral_table:
      break;
    case multi_inspiral_table:
      break;
    case sim_inspiral_table:
      break;
    case sim_burst_table:
      {
        MetaTableDirectory tmpTableDir[] =
        {
          {"waveform",                     -1, 0},
          {"geocent_peak_time",            -1, 1},
          {"geocent_peak_time_ns",         -1, 2},
          {"h_peak_time",                  -1, 3},
          {"h_peak_time_ns",               -1, 4},
          {"l_peak_time",                  -1, 5},
          {"l_peak_time_ns",               -1, 6},
          {"peak_time_gmst",               -1, 7},
          {"dtplus",                       -1, 8},
          {"dtminus",                      -1, 9},
          {"longitude",                    -1, 10},
          {"latitude",                     -1, 11},
          {"coordinates",                  -1, 12},
          {"polarization",                 -1, 13},
          {"hrss",                         -1, 14},
          {"hpeak",                        -1, 15},
          {"freq",                         -1, 16},
          {"tau",                          -1, 17},
          {"zm_number",                    -1, 18},
          {NULL,                            0, 0}
        };
        for ( i=0 ; tmpTableDir[i].name; ++i )
        {
          if ( (tmpTableDir[i].pos = 
                MetaioFindColumn( env, tmpTableDir[i].name )) < 0 )
          {
            fprintf( stderr, "unable to find column %s\n", 
                tmpTableDir[i].name );
            ABORT(status,LIGOLWXMLREADH_ENCOL,LIGOLWXMLREADH_MSGENCOL);
          }
        }

        *tableDir = (MetaTableDirectory *) LALMalloc( (i+1) * 
            sizeof(MetaTableDirectory)) ;
        memcpy(*tableDir, tmpTableDir, (i+1)*sizeof(MetaTableDirectory) );
      }
      break;
    case summ_value_table:
      break;
    default:
      ABORT( status, LIGOLWXMLREADH_EUTAB, LIGOLWXMLREADH_MSGEUTAB );
  }

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN( status );
}

#define CLOBBER_EVENTS \
    while ( *eventHead ); \
    { \
      thisEvent = *eventHead; \
      *eventHead = (*eventHead)->next; \
      LALFree( thisEvent ); \
      thisEvent = NULL; \
    }

/* <lalVerbatim file="LIGOLwXMLReadCP"> */
void
LALSnglBurstTableFromLIGOLw (
    LALStatus          *status,
    SnglBurstTable    **eventHead,
    CHAR               *fileName
    )
/* </lalVerbatim> */
{
  int                                   i, j, nrows;
  int                                   mioStatus=0;
  SnglBurstTable                       *thisEvent = NULL;
  struct MetaioParseEnvironment         parseEnv;
  const  MetaioParseEnv                 env = &parseEnv;
  MetaTableDirectory                   *tableDir = NULL;

  INITSTATUS( status, "LALSnglBurstTableFromLIGOLw", LIGOLWXMLREADC );
  ATTATCHSTATUSPTR (status);

  /* check that the event handle and pointer are vaid */
  if ( ! eventHead )
  {
    ABORT(status, LIGOLWXMLREADH_ENULL, LIGOLWXMLREADH_MSGENULL);
  }
  if ( *eventHead )
  {
    ABORT(status, LIGOLWXMLREADH_ENNUL, LIGOLWXMLREADH_MSGENNUL);
  }

  /* open the sngl_burst XML file */
  mioStatus = MetaioOpenTable( env, fileName, "sngl_burst" );
  if ( mioStatus )
  {
    ABORT(status, LIGOLWXMLREADH_ENTAB, LIGOLWXMLREADH_MSGENTAB);
  }

  /* create table directory to find columns in file*/
  LALCreateMetaTableDir(status->statusPtr, &tableDir, env, sngl_burst_table);
  CHECKSTATUSPTR (status);

  /* loop over the rows in the file */
  i = nrows = 0;
  while ( (mioStatus = MetaioGetRow(env)) == 1 ) 
  {
    /* count the rows in the file */
    i++;

    /* allocate memory for the template we are about to read in */
    if ( ! *eventHead )
    {
      thisEvent = *eventHead = (SnglBurstTable *) 
        LALCalloc( 1, sizeof(SnglBurstTable) );
    }
    else
    {
      thisEvent = thisEvent->next = (SnglBurstTable *) 
        LALCalloc( 1, sizeof(SnglBurstTable) );
    }
    if ( ! thisEvent )
    {
      fprintf( stderr, "could not allocate burst event\n" );
      CLOBBER_EVENTS;
      MetaioClose( env );
      ABORT(status, LIGOLWXMLREADH_EALOC, LIGOLWXMLREADH_MSGEALOC);
    }

    /* parse the contents of the row into the InspiralTemplate structure */
    for ( j = 0; tableDir[j].name; ++j )
    {
      REAL4 r4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
      REAL8 r8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_8;
      INT4  i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;

      if ( tableDir[j].idx == 0 )
      {
        LALSnprintf( thisEvent->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), 
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 1 )
      {
        LALSnprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 2 )
      {
        LALSnprintf( thisEvent->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 3 )
      {
        thisEvent->start_time.gpsSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 4 )
      {
        thisEvent->start_time.gpsNanoSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 5 )
      {
        thisEvent->duration = r4colData;
      }
      else if ( tableDir[j].idx == 6 )
      {
        thisEvent->central_freq = r4colData;
      }
      else if ( tableDir[j].idx == 7 )
      {
        thisEvent->bandwidth = r4colData;
      }
      else if ( tableDir[j].idx == 8 )
      {
        thisEvent->amplitude = r4colData;
      }
      else if ( tableDir[j].idx == 9 )
      {
        thisEvent->snr = r4colData;
      }
      else if ( tableDir[j].idx == 10 )
      {
        thisEvent->confidence = r4colData;
      }
      else if ( tableDir[j].idx == 11 )
      {
        thisEvent->peak_time.gpsSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 12 )
      {
        thisEvent->peak_time.gpsNanoSeconds = i4colData;
      }
      else
      {
        CLOBBER_EVENTS;
        ABORT(status, LIGOLWXMLREADH_ENCOL, LIGOLWXMLREADH_MSGENCOL);
      }
    }
    /* count the number of triggers parsed */
    nrows++;
  }

  if ( mioStatus == -1 )
  {
    fprintf( stderr, "error parsing after row %d\n", i );
    CLOBBER_EVENTS;
    MetaioClose( env );
    ABORT(status, LIGOLWXMLREADH_EPARS, LIGOLWXMLREADH_MSGEPARS);
  }

  /* Normal exit */
  LALFree( tableDir );
  MetaioClose( env );
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/* <lalVerbatim file="LIGOLwXMLReadCP"> */
void
LALSimBurstTableFromLIGOLw (
    LALStatus          *status,
    SimBurstTable    **eventHead,
    CHAR               *fileName,
    INT4                startTime,
    INT4                stopTime
    )
/* </lalVerbatim> */
{
  int                                   i, j, nrows;
  int                                   mioStatus=0;
  SimBurstTable                        *thisEvent = NULL;
  struct MetaioParseEnvironment         parseEnv;
  const  MetaioParseEnv                 env = &parseEnv;
  MetaTableDirectory                   *tableDir = NULL;

  INITSTATUS( status, "LALSimBurstTableFromLIGOLw", LIGOLWXMLREADC );
  ATTATCHSTATUSPTR (status);

  /* check that the event handle and pointer are vaid */
  if ( ! eventHead )
  {
    ABORT(status, LIGOLWXMLREADH_ENULL, LIGOLWXMLREADH_MSGENULL);
  }
  if ( *eventHead )
  {
    ABORT(status, LIGOLWXMLREADH_ENNUL, LIGOLWXMLREADH_MSGENNUL);
  }

  /* open the sim_burst XML file */
  mioStatus = MetaioOpenTable( env, fileName, "sim_burst" );
  if ( mioStatus )
  {
    ABORT(status, LIGOLWXMLREADH_ENTAB, LIGOLWXMLREADH_MSGENTAB);
  }

  /* create table directory to find columns in file*/
  LALCreateMetaTableDir(status->statusPtr, &tableDir, env, sim_burst_table);
  CHECKSTATUSPTR (status);

  /* loop over the rows in the file */
  i = nrows = 0;
  while ( (mioStatus = MetaioGetRow(env)) == 1 ) 
  {
    INT4 geo_time = env->ligo_lw.table.elt[tableDir[1].pos].data.int_4s;

    /* count the rows in the file */
    i++;

    /* get the injetcion time and check that it is within the time window */
    if ( ! stopTime || geo_time > startTime && geo_time < stopTime )
    {

      /* allocate memory for the template we are about to read in */
      if ( ! *eventHead )
      {
        thisEvent = *eventHead = (SimBurstTable *) 
          LALCalloc( 1, sizeof(SimBurstTable) );
      }
      else
      {
        thisEvent = thisEvent->next = (SimBurstTable *) 
          LALCalloc( 1, sizeof(SimBurstTable) );
      }
      if ( ! thisEvent )
      {
        fprintf( stderr, "could not allocate burst event\n" );
        CLOBBER_EVENTS;
        MetaioClose( env );
        ABORT(status, LIGOLWXMLREADH_EALOC, LIGOLWXMLREADH_MSGEALOC);
      }

      /* parse the contents of the row into the InspiralTemplate structure */
      for ( j = 0; tableDir[j].name; ++j )
      {
        REAL4 r4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
        REAL8 r8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_8;
        INT4  i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;

        if ( tableDir[j].idx == 0 )
        {
          LALSnprintf( thisEvent->waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), 
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 1 )
        {
          thisEvent->geocent_peak_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 2 )
        {
          thisEvent->geocent_peak_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 3 )
        {
          thisEvent->h_peak_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 4 )
        {
          thisEvent->h_peak_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 5 )
        {
          thisEvent->l_peak_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 6 )
        {
          thisEvent->l_peak_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 7 )
        {
          thisEvent->peak_time_gmst = r8colData;
        }
        else if ( tableDir[j].idx == 8 )
        {
          thisEvent->dtminus = r4colData;
        }
        else if ( tableDir[j].idx == 9 )
        {
          thisEvent->dtplus = r4colData;
        }
        else if ( tableDir[j].idx == 10 )
        {
          thisEvent->longitude = r4colData;
        }
        else if ( tableDir[j].idx == 11 )
        {
          thisEvent->latitude = r4colData;
        }
        else if ( tableDir[j].idx == 12 )
        {
          LALSnprintf( thisEvent->coordinates, LIGOMETA_COORDINATES_MAX * sizeof(CHAR), 
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 13 )
        {
          thisEvent->polarization = r4colData;
        }
        else if ( tableDir[j].idx == 14 )
        {
          thisEvent->hrss = r4colData;
        }
        else if ( tableDir[j].idx == 15 )
        {
          thisEvent->hpeak = r4colData;
        }
        else if ( tableDir[j].idx == 16 )
        {
          thisEvent->freq = r4colData;
        }
        else if ( tableDir[j].idx == 17 )
        {
          thisEvent->tau = r4colData;
        }
        else if ( tableDir[j].idx == 18 )
        {
          thisEvent->zm_number = i4colData;
        }
        else
        {
          CLOBBER_EVENTS;
          ABORT(status, LIGOLWXMLREADH_ENCOL, LIGOLWXMLREADH_MSGENCOL);
        }
      }

      nrows++;
    }
  }

  if ( mioStatus == -1 )
  {
    fprintf( stderr, "error parsing after row %d\n", i );
    CLOBBER_EVENTS;
    MetaioClose( env );
    ABORT(status, LIGOLWXMLREADH_EPARS, LIGOLWXMLREADH_MSGEPARS);
  }

  /* Normal exit */
  LALFree( tableDir );
  MetaioClose( env );
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/* <lalVerbatim file="LIGOLwXMLReadCP"> */
int
LALSnglInspiralTableFromLIGOLw (
    SnglInspiralTable **eventHead,
    CHAR               *fileName,
    INT4                startEvent,
    INT4                stopEvent
    )
/* </lalVerbatim> */
{
  int                                   i, j, nrows;
  int                                   mioStatus;
  SnglInspiralTable                    *thisEvent = NULL;
  struct MetaioParseEnvironment         parseEnv;
  const  MetaioParseEnv                 env = &parseEnv;
  MetaTableDirectory tableDir[] =
  {
    {"ifo",                     -1, 0},
    {"search",                  -1, 1},
    {"channel",                 -1, 2},
    {"end_time",                -1, 3},
    {"end_time_ns",             -1, 4},
    {"end_time_gmst",           -1, 5},
    {"impulse_time",            -1, 6},
    {"impulse_time_ns",         -1, 7},
    {"template_duration",       -1, 8},
    {"event_duration",          -1, 9},
    {"amplitude",               -1, 10},
    {"eff_distance",            -1, 11},
    {"coa_phase",               -1, 12},
    {"mass1",                   -1, 13},
    {"mass2",                   -1, 14},
    {"mchirp",                  -1, 15},
    {"eta",                     -1, 16},
    {"tau0",                    -1, 17},
    {"tau2",                    -1, 18},
    {"tau3",                    -1, 19},
    {"tau4",                    -1, 20},
    {"tau5",                    -1, 21},
    {"ttotal",                  -1, 22},
    {"psi0",                    -1, 23},
    {"psi3",                    -1, 24},
    {"f_cut",                   -1, 25},
    {"snr",                     -1, 26},
    {"chisq",                   -1, 27},
    {"chisq_dof",               -1, 28},
    {"sigmasq",                 -1, 29},
    {"alpha",                   -1, 30},
    {NULL,                       0, 0}
  };

  /* check that the bank handle and pointer are vaid */
  if ( ! eventHead )
  {
    fprintf( stderr, "null pointer passed as handle to event list" );
    return -1;
  }
  if ( *eventHead )
  {
    fprintf( stderr, "non-null pointer passed as pointer to event list" );
    return -1;
  }

  /* open the sngl_inspiral table template bank file */
  mioStatus = MetaioOpenTable( env, fileName, "sngl_inspiral" );
  if ( mioStatus )
  {
    fprintf( stdout, "no sngl_inspiral table in file %s\n", fileName );
    return 0;
  }

  /* figure out the column positions of the template parameters */
  for ( i = 0; tableDir[i].name; ++i )
  {
    if ( (tableDir[i].pos = MetaioFindColumn( env, tableDir[i].name )) < 0 )
    {
      fprintf( stderr, "unable to find column %s\n", tableDir[i].name );
      MetaioClose(env);
      return -1;
    }
  }

  /* loop over the rows in the file */
  i = nrows = 0;
  while ( (mioStatus = MetaioGetRow(env)) == 1 ) 
  {
    /* count the rows in the file */
    i++;

    /* stop parsing if we have reach the last row requested */
    if ( stopEvent > -1 && i > stopEvent )
    {
      break;
    }

    /* if we have reached the first requested row, parse the row */
    if ( i > startEvent )
    {
      /* allocate memory for the template we are about to read in */
      if ( ! *eventHead )
      {
        thisEvent = *eventHead = (SnglInspiralTable *) 
          LALCalloc( 1, sizeof(SnglInspiralTable) );
      }
      else
      {
        thisEvent = thisEvent->next = (SnglInspiralTable *) 
          LALCalloc( 1, sizeof(SnglInspiralTable) );
      }
      if ( ! thisEvent )
      {
        fprintf( stderr, "could not allocate inspiral template\n" );
        CLOBBER_EVENTS;
        MetaioClose( env );
        return -1;
      }

      /* parse the contents of the row into the InspiralTemplate structure */
      for ( j = 0; tableDir[j].name; ++j )
      {
        REAL4 r4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
        REAL8 r8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_8;
        INT4  i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;

        if ( tableDir[j].idx == 0 )
        {
          LALSnprintf( thisEvent->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), 
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 1 )
        {
          LALSnprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 2 )
        {
          LALSnprintf( thisEvent->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 3 )
        {
          thisEvent->end_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 4 )
        {
          thisEvent->end_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 5 )
        {
          thisEvent->end_time_gmst = r8colData;
        }
        else if ( tableDir[j].idx == 6 )
        {
          thisEvent->impulse_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 7 )
        {
          thisEvent->impulse_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 8 )
        {
          thisEvent->template_duration = r8colData;
        }
        else if ( tableDir[j].idx == 9 )
        {
          thisEvent->event_duration = r8colData;
        }
        else if ( tableDir[j].idx == 10 )
        {
          thisEvent->amplitude = r4colData;
        }
        else if ( tableDir[j].idx == 11 )
        {
          thisEvent->eff_distance = r4colData;
        }
        else if ( tableDir[j].idx == 12 )
        {
          thisEvent->coa_phase = r4colData;
        }
        else if ( tableDir[j].idx == 13 )
        {
          thisEvent->mass1 = r4colData;
        }
        else if ( tableDir[j].idx == 14 )
        {
          thisEvent->mass2 = r4colData;
        }
        else if ( tableDir[j].idx == 15 )
        {
          thisEvent->mchirp = r4colData;
        }
        else if ( tableDir[j].idx == 16 )
        {
          thisEvent->eta = r4colData;
        }
        else if ( tableDir[j].idx == 17 )
        {
          thisEvent->tau0 = r4colData;
        }
        else if ( tableDir[j].idx == 18 )
        {
          thisEvent->tau2 = r4colData;
        }
        else if ( tableDir[j].idx == 19 )
        {
          thisEvent->tau3 = r4colData;
        }
        else if ( tableDir[j].idx == 20 )
        {
          thisEvent->tau4 = r4colData;
        }
        else if ( tableDir[j].idx == 21 )
        {
          thisEvent->tau5 = r4colData;
        }
        else if ( tableDir[j].idx == 22 )
        {
          thisEvent->ttotal = r4colData;
        }
        else if ( tableDir[j].idx == 23 )
        {
          thisEvent->psi0 = r4colData;
        }
        else if ( tableDir[j].idx == 24 )
        {
          thisEvent->psi3 = r4colData;
        }
        else if ( tableDir[j].idx == 25 )
        {
          thisEvent->f_cut = r4colData;
        }
        else if ( tableDir[j].idx == 26 )
        {
          thisEvent->snr = r4colData;
        }
        else if ( tableDir[j].idx == 27 )
        {
          thisEvent->chisq = r4colData;
        }
        else if ( tableDir[j].idx == 28 )
        {
          thisEvent->chisq_dof = i4colData;
        }
        else if ( tableDir[j].idx == 29 )
        {
          thisEvent->sigmasq = r8colData;
        }
        else if ( tableDir[j].idx == 30 )
        {
          thisEvent->alpha = r4colData;
        }
        else
        {
          CLOBBER_EVENTS;
          fprintf( stderr, "unknown column while parsing sngl_inspiral\n" );
          return -1;
        }
      }

      /* count the number of template parsed */
      nrows++;
    }
  }

  if ( mioStatus == -1 )
  {
    fprintf( stderr, "error parsing after row %d\n", i );
    CLOBBER_EVENTS;
    MetaioClose( env );
    return -1;
  }

  /* we have sucesfully parsed temples */
  MetaioClose( env );
  return nrows;  
}

#undef CLOBBER_EVENTS

#define CLOBBER_BANK \
    while ( *bankHead ); \
    { \
      thisTmplt = *bankHead; \
      *bankHead = (*bankHead)->next; \
      LALFree( thisTmplt ); \
      thisTmplt = NULL; \
    }

/* <lalVerbatim file="LIGOLwXMLReadCP"> */
int
InspiralTmpltBankFromLIGOLw (
    InspiralTemplate  **bankHead,
    CHAR               *fileName,
    INT4                startTmplt,
    INT4                stopTmplt
    )
/* </lalVerbatim> */
{
  int                                   i, j, nrows;
  int                                   mioStatus;
  InspiralTemplate                     *thisTmplt = NULL;
  struct MetaioParseEnvironment         parseEnv;
  const  MetaioParseEnv                 env = &parseEnv;
  int   pParParam;
  int   pParValue;
  REAL4 minMatch = 0;
  MetaTableDirectory tableDir[] =
  {
    {"mass1",   -1, 0},
    {"mass2",   -1, 1},
    {"mchirp",  -1, 2},
    {"eta",     -1, 3},
    {"tau0",    -1, 4},
    {"tau2",    -1, 5},
    {"tau3",    -1, 6},
    {"tau4",    -1, 7},
    {"tau5",    -1, 8},
    {"ttotal",  -1, 9},
    {"psi0",    -1, 10},
    {"psi3",    -1, 11},
    {"f_cut",   -1, 12},
    {NULL,      0, 0}
  };

  /* check that the bank handle and pointer are vaid */
  if ( ! bankHead )
  {
    fprintf( stderr, "null pointer passed as handle to template bank" );
    return -1;
  }
  if ( *bankHead )
  {
    fprintf( stderr, "non-null pointer passed as pointer to template bank" );
    return -1;
  }
  

  /* open the procress_params table from the bank file */
  mioStatus = MetaioOpenTable( env, fileName, "process_params" );
  if ( mioStatus )
  {
    fprintf( stderr, "error opening process_params table from file %s\n", 
        fileName );
    return -1;
  }

  /* figure out where the param and value columns are */
  if ( (pParParam = MetaioFindColumn( env, "param" )) < 0 )
  {
    fprintf( stderr, "unable to find column param in process_params\n" );
    MetaioClose(env);
    return -1;
  }
  if ( (pParValue = MetaioFindColumn( env, "value" )) < 0 )
  {
    fprintf( stderr, "unable to find column value in process_params\n" );
    MetaioClose(env);
    return -1;
  }

  /* get the minimal match of the bank from the process params */
  while ( (mioStatus = MetaioGetRow(env)) == 1 )
  {
    if ( ! strcmp( env->ligo_lw.table.elt[pParParam].data.lstring.data, 
          "--minimal-match" ) )
    {
      minMatch = (REAL4) 
        atof( env->ligo_lw.table.elt[pParValue].data.lstring.data );
    }
  }

  MetaioClose( env );
  
  /* open the sngl_inspiral table template bank file */
  mioStatus = MetaioOpenTable( env, fileName, "sngl_inspiral" );
  if ( mioStatus )
  {
    fprintf( stdout, "no sngl_inspiral table in file %s\n", fileName );
    return 0;
  }

  /* figure out the column positions of the template parameters */
  for ( i = 0; tableDir[i].name; ++i )
  {
    if ( (tableDir[i].pos = MetaioFindColumn( env, tableDir[i].name )) < 0 )
    {
      fprintf( stderr, "unable to find column %s\n", tableDir[i].name );
      MetaioClose(env);
      return -1;
    }
  }

  /* loop over the rows in the file */
  i = nrows = 0;
  while ( (mioStatus = MetaioGetRow(env)) == 1 ) 
  {
    /* count the rows in the file */
    i++;

    /* stop parsing if we have reach the last row requested */
    if ( stopTmplt > -1 && i > stopTmplt )
    {
      break;
    }

    /* if we have reached the first requested row, parse the row */
    if ( i > startTmplt )
    {
      /* allocate memory for the template we are about to read in */
      if ( ! *bankHead )
      {
        thisTmplt = *bankHead = (InspiralTemplate *) 
          LALCalloc( 1, sizeof(InspiralTemplate) );
      }
      else
      {
        thisTmplt = thisTmplt->next = (InspiralTemplate *) 
          LALCalloc( 1, sizeof(InspiralTemplate) );
      }
      if ( ! thisTmplt )
      {
        fprintf( stderr, "could not allocate inspiral template\n" );
        CLOBBER_BANK;
        MetaioClose( env );
        return -1;
      }
      
      /* parse the contents of the row into the InspiralTemplate structure */
      for ( j = 0; tableDir[j].name; ++j )
      {
        REAL4 colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
        if ( tableDir[j].idx == 0 )
        {
          thisTmplt->mass1 = colData;
        }
        else if ( tableDir[j].idx == 1 )
        {
          thisTmplt->mass2 = colData;
        }
        else if ( tableDir[j].idx == 2 )
        {
          thisTmplt->chirpMass = colData;
        }
        else if ( tableDir[j].idx == 3 )
        {
          thisTmplt->eta = colData;
        }
        else if ( tableDir[j].idx == 4 )
        {
          thisTmplt->t0 = colData;
        }
        else if ( tableDir[j].idx == 5 )
        {
          thisTmplt->t2 = colData;
        }
        else if ( tableDir[j].idx == 6 )
        {
          thisTmplt->t3 = colData;
        }
        else if ( tableDir[j].idx == 7 )
        {
          thisTmplt->t4 = colData;
        }
        else if ( tableDir[j].idx == 8 )
        {
          thisTmplt->t5 = colData;
        }
        else if ( tableDir[j].idx == 9 )
        {
          thisTmplt->tC = colData;
        }
        else if ( tableDir[j].idx == 10 )
        {
          thisTmplt->psi0 = colData;
        }
        else if ( tableDir[j].idx == 11 )
        {
          thisTmplt->psi3 = colData;
        }
        else if ( tableDir[j].idx == 12 )
        {
          thisTmplt->fFinal = colData;
        }
        else
        {
          CLOBBER_BANK;
          fprintf( stderr, "unknown column while parsing\n" );
          return -1;
        }
      }

      /* compute derived mass parameters */
      thisTmplt->totalMass = thisTmplt->mass1 + thisTmplt->mass2;
      if ( thisTmplt->totalMass > 0 )
      {
        thisTmplt->mu = thisTmplt->mass1 * thisTmplt->mass2 / 
          thisTmplt->totalMass;
      }
      
      /* set the match determined from the bank generation process params */
      thisTmplt->minMatch = minMatch;

      /* count the number of template parsed */
      thisTmplt->number = nrows++;
    }
  }

  if ( mioStatus == -1 )
  {
    fprintf( stderr, "error parsing after row %d\n", i );
    CLOBBER_BANK;
    MetaioClose( env );
    return -1;
  }

  /* we have sucesfully parsed temples */
  MetaioClose( env );
  return nrows;  
}

#define CLOBBER_SIM \
    while ( *simHead ); \
    { \
      thisSim = *simHead; \
      *simHead = (*simHead)->next; \
      LALFree( thisSim ); \
      thisSim = NULL; \
    }

/* <lalVerbatim file="LIGOLwXMLReadCP"> */
int
SimInspiralTableFromLIGOLw (
    SimInspiralTable   **simHead,
    CHAR                *fileName,
    INT4                 startTime,
    INT4                 endTime
    )
/* </lalVerbatim> */
{
  int                                   i, j, nrows;
  int                                   mioStatus;
  SimInspiralTable                     *thisSim = NULL;
  struct MetaioParseEnvironment         parseEnv;
  const  MetaioParseEnv                 env = &parseEnv;
  MetaTableDirectory tableDir[] =
  {
    {"waveform",	    -1, 0},
    {"geocent_end_time",    -1, 1},
    {"geocent_end_time_ns", -1, 2},
    {"h_end_time",	    -1, 3},
    {"h_end_time_ns",	    -1, 4},
    {"l_end_time",	    -1, 5},
    {"l_end_time_ns",	    -1, 6},
    {"end_time_gmst",	    -1,	7},
    {"source",		    -1, 8},
    {"mass1",               -1, 9},
    {"mass2",               -1, 10},
    {"eta",                 -1, 11},
    {"distance",            -1, 12},
    {"longitude",           -1, 13},
    {"latitude",            -1, 14},
    {"inclination",         -1, 15},
    {"coa_phase",           -1, 16},
    {"polarization",        -1, 17},
    {"eff_dist_h",	    -1, 18},
    {"eff_dist_l",	    -1, 19},
    {NULL,                   0, 0}
  };

  /* check that the bank handle and pointer are valid */
  if ( ! simHead )
  {
    fprintf( stderr, "null pointer passed as handle to simulation list" );
    return -1;
  }
  if ( *simHead )
  {
    fprintf( stderr, "non-null pointer passed as pointer to simulation list" );
    return -1;
  }
  
  /* open the sngl_inspiral table file */
  mioStatus = MetaioOpenTable( env, fileName, "sim_inspiral" );
  if ( mioStatus )
  {
    fprintf( stderr, "error opening sim_inspiral table from file %s\n", 
        fileName );
    return -1;
  }

  /* figure out the column positions of the template parameters */
  for ( i = 0; tableDir[i].name; ++i )
  {
    if ( (tableDir[i].pos = MetaioFindColumn( env, tableDir[i].name )) < 0 )
    {
      fprintf( stderr, "unable to find column %s\n", tableDir[i].name );
      MetaioClose(env);
      return -1;
    }
  }

  /* loop over the rows in the file */
  i = nrows = 0;
  while ( (mioStatus = MetaioGetRow(env)) == 1 ) 
  {
    INT4 geo_time = env->ligo_lw.table.elt[tableDir[1].pos].data.int_4s;

      /* get the injetcion time and check that it is within the time window */
      if ( ! endTime || geo_time > startTime && geo_time < endTime )
      {
        /* allocate memory for the template we are about to read in */
        if ( ! *simHead )
        {
          thisSim = *simHead = (SimInspiralTable *) 
            LALCalloc( 1, sizeof(SimInspiralTable) );
        }
        else
        {
          thisSim = thisSim->next = (SimInspiralTable *) 
            LALCalloc( 1, sizeof(SimInspiralTable) );
        }
        if ( ! thisSim )
        {
          fprintf( stderr, "could not allocate inspiral simulation\n" );
          CLOBBER_SIM;
          MetaioClose( env );
          return -1;
        }

        /* parse the row into the SimInspiralTable structure */
	for ( j = 0; tableDir[j].name; ++j )
	{
	  REAL4 r4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
	  REAL8 r8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_8;
	  INT4  i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;
	  if ( tableDir[j].idx == 0 )
	  {
	    LALSnprintf(thisSim->waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
	      "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
	  }    
	  else if ( tableDir[j].idx == 1 ) 	  
	  {
	    thisSim->geocent_end_time.gpsSeconds = i4colData;
	  }
	  else if ( tableDir[j].idx == 2 )
	  {
	    thisSim->geocent_end_time.gpsNanoSeconds = i4colData;
	  }
	  else if ( tableDir[j].idx == 3 ) 	  
	  {
	    thisSim->h_end_time.gpsSeconds = i4colData;
	  }
	  else if ( tableDir[j].idx == 4 )
	  {
	    thisSim->h_end_time.gpsNanoSeconds = i4colData;
	  }
	  else if ( tableDir[j].idx == 5 ) 	  
	  {
	    thisSim->l_end_time.gpsSeconds = i4colData;
	  }
	  else if ( tableDir[j].idx == 6 )
	  {
	    thisSim->l_end_time.gpsNanoSeconds = i4colData;
	  }
	  else if ( tableDir[j].idx == 7 )
          {
	    thisSim->end_time_gmst = r8colData;
	  }
	  else if ( tableDir[j].idx == 8 )
	  {
	    LALSnprintf(thisSim->source, LIGOMETA_SOURCE_MAX * sizeof(CHAR),
	      "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
	  }
	  else if ( tableDir[j].idx == 9 )
	  {
	    thisSim->mass1 = r4colData;
          }
	  else if ( tableDir[j].idx == 10 )
	  {
	    thisSim->mass2 = r4colData;
          }
          else if ( tableDir[j].idx == 11 )
          {
            thisSim->eta = r4colData;
          }
          else if ( tableDir[j].idx == 12 )
          {
	    thisSim->distance = r4colData;
          }
          else if ( tableDir[j].idx == 13 )
          {
            thisSim->longitude = r4colData;
          }
          else if ( tableDir[j].idx == 14 )
          {
            thisSim->latitude = r4colData;
          }
          else if ( tableDir[j].idx == 15 )
          {
            thisSim->inclination = r4colData;
          }
          else if ( tableDir[j].idx == 16 )
          {
            thisSim->coa_phase = r4colData;
          }
          else if ( tableDir[j].idx == 17 )
          {
            thisSim->polarization = r4colData;
          }
	  else if ( tableDir[j].idx == 18 )
          {
            thisSim->eff_dist_h = r4colData;
          }
	  else if ( tableDir[j].idx == 19 )
          {
            thisSim->eff_dist_l = r4colData;
          }
          else
          {
            CLOBBER_SIM;
            fprintf( stderr, "unknown column while parsing\n" );
            return -1;
          }
        }
        /* increase the count of rows parsed */
        ++nrows;       
      }
  }

  if ( mioStatus == -1 )
  {
    fprintf( stderr, "error parsing after row %d\n", i );
    CLOBBER_SIM;
    MetaioClose( env );
    return -1;
  }

  /* we have sucesfully parsed temples */
  MetaioClose( env );
  return nrows;  
}

#undef CLOBBER_SIM

/* <lalVerbatim file="LIGOLwXMLReadCP"> */
int
SearchSummaryTableFromLIGOLw (
    SearchSummaryTable **sumHead,
    CHAR                *fileName
    )
/* </lalVerbatim> */
{
  int                                   i, j, nrows;
  int                                   mioStatus;
  struct MetaioParseEnvironment         parseEnv;
  const  MetaioParseEnv                 env = &parseEnv;
  MetaTableDirectory tableDir[] =
  {
    {"comment",                 -1, 0},
    {"in_start_time",           -1, 1},
    {"in_start_time_ns",        -1, 2},
    {"in_end_time",             -1, 3},
    {"in_end_time_ns",          -1, 4},
    {"out_start_time",          -1, 5},
    {"out_start_time_ns",       -1, 6},
    {"out_end_time",            -1, 7},
    {"out_end_time_ns",         -1, 8},
    {"nevents",                 -1, 9},
    {"nnodes",                  -1, 10},
    {NULL,                       0, 0}
  };

  /* check that the bank handle and pointer are vaid */
  if ( ! sumHead )
  {
    fprintf( stderr, "null pointer passed as handle to search summary" );
    return -1;
  }
  if ( *sumHead )
  {
    fprintf( stderr, "non-null pointer passed as pointer to search summary" );
    return -1;
  }
  
  /* open the search_summary table in the file file */
  mioStatus = MetaioOpenTable( env, fileName, "search_summary" );
  if ( mioStatus )
  {
    fprintf( stderr, "error opening search_summary table from file %s\n", 
        fileName );
    return -1;
  }

  /* figure out the column positions */
  for ( i = 0; tableDir[i].name; ++i )
  {
    if ( (tableDir[i].pos = MetaioFindColumn( env, tableDir[i].name )) < 0 )
    {
      fprintf( stderr, "unable to find column %s\n", tableDir[i].name );
      MetaioClose(env);
      return -1;
    }
  }

  /* allocate memory for the table */
  *sumHead = (SearchSummaryTable *) LALCalloc( 1, sizeof(SearchSummaryTable) );

  /* loop over the rows in the file */
  i = nrows = 0;
  while ( (mioStatus = MetaioGetRow(env)) == 1 ) 
  {
    /* parse the rows into the SearhSummary structure */
    for ( j = 1; tableDir[j].name; ++j )
    {
      INT4 intData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;

      if ( tableDir[j].idx == 0 )
      {
        LALSnprintf( (*sumHead)->comment, LIGOMETA_COMMENT_MAX * sizeof(CHAR),
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 1 )
      {
        (*sumHead)->in_start_time.gpsSeconds = intData;
      }
      else if ( tableDir[j].idx == 2 )
      {
        (*sumHead)->in_start_time.gpsNanoSeconds = intData;
      }
      else if ( tableDir[j].idx == 3 )
      {
        (*sumHead)->in_end_time.gpsSeconds = intData;
      }
      else if ( tableDir[j].idx == 4 )
      {
        (*sumHead)->in_end_time.gpsNanoSeconds = intData;
      }
      else if ( tableDir[j].idx == 5 )
      {
        (*sumHead)->out_start_time.gpsSeconds = intData;
      }
      else if ( tableDir[j].idx == 6 )
      {
        (*sumHead)->out_start_time.gpsNanoSeconds = intData;
      }
      else if ( tableDir[j].idx == 7 )
      {
        (*sumHead)->out_end_time.gpsSeconds = intData;
      }
      else if ( tableDir[j].idx == 8 )
      {
        (*sumHead)->out_end_time.gpsNanoSeconds = intData;
      }
      else if ( tableDir[j].idx == 9 )
      {
        (*sumHead)->nevents = intData;
      }
      else if ( tableDir[j].idx == 10 )
      {
        (*sumHead)->nnodes = intData;
      }
      else
      {
        LALFree( *sumHead );
        fprintf( stderr, "unknown column while parsing\n" );
        return -1;
      }
    }

    /* increase the count of rows parsed */
    ++nrows;
  }

  if ( mioStatus == -1 )
  {
    fprintf( stderr, "error parsing after row %d\n", i );
    LALFree( *sumHead );
    MetaioClose( env );
    return -1;
  }

  /* we have sucesfully parsed table */
  MetaioClose( env );
  return nrows;  
}

/* <lalVerbatim file="LIGOLwXMLReadCP"> */
int
SummValueTableFromLIGOLw (
    SummValueTable **sumHead,
    CHAR           *fileName
    )
/* </lalVerbatim> */
{
  int                                   i, j, nrows;
  int                                   mioStatus;
  struct MetaioParseEnvironment         parseEnv;
  const  MetaioParseEnv                 env = &parseEnv;
  MetaTableDirectory tableDir[] =
  {
    {"program",		      -1, 0},
    {"process_id",	      -1, 1},
    {"start_time",	      -1, 2},
    {"start_time_ns",         -1, 3},
    {"end_time",	      -1, 4},
    {"end_time_ns",	      -1, 5},
    {"ifo",		      -1, 6},
    {"name",		      -1, 7},
    {"value",		      -1, 8},
    {"comment",               -1, 9},
    {NULL,                     0, 0}
  };

  /* check that the bank handle and pointer are vaid */
  if ( ! sumHead )
  {
    fprintf( stderr, "null pointer passed as handle to summ value" );
    return -1;
  }
  if ( *sumHead )
  {
    fprintf( stderr, "non-null pointer passed as pointer to summ value" );
    return -1;
  }
  
  /* open the summ_value table in the file file */
  mioStatus = MetaioOpenTable( env, fileName, "summ_value" );
  if ( mioStatus )
  {
    fprintf( stderr, "error opening summ_value table from file %s\n", 
        fileName );
    return -1;
  }

  /* figure out the column positions */
  for ( i = 0; tableDir[i].name; ++i )
  {
    if ( (tableDir[i].pos = MetaioFindColumn( env, tableDir[i].name )) < 0 )
    {
      fprintf( stderr, "unable to find column %s\n", tableDir[i].name );
      MetaioClose(env);
      return -1;
    }
  }

  /* allocate memory for the table */
  *sumHead = (SummValueTable *) LALCalloc( 1, sizeof(SummValueTable) );

  /* loop over the rows in the file */
  i = nrows = 0;
  while ( (mioStatus = MetaioGetRow(env)) == 1 ) 
  {
    /* parse the rows into the SearhSummary structure */
    for ( j = 1; tableDir[j].name; ++j )
    {
      REAL4 r4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
      INT4  i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;

      if ( tableDir[j].idx == 0 )
      {
        LALSnprintf( (*sumHead)->program, LIGOMETA_PROGRAM_MAX * sizeof(CHAR),
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 1 )
      {
        (*sumHead)->start_time.gpsSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 2 )
      {
        (*sumHead)->start_time.gpsNanoSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 3 )
      {
        (*sumHead)->end_time.gpsSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 4 )
      {
        (*sumHead)->end_time.gpsNanoSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 5 )
      {
	LALSnprintf( (*sumHead)->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR),
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 6 )
      {
	LALSnprintf( (*sumHead)->name, LIGOMETA_SUMMVALUE_NAME_MAX * 
	    sizeof(CHAR), "%s",
	    env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 7 )
      {
	(*sumHead)->value = r4colData;
      }
      else if ( tableDir[j].idx == 8 )
      {
        LALSnprintf( (*sumHead)->comment, LIGOMETA_SUMMVALUE_NAME_MAX * 
	    sizeof(CHAR), "%s",
	    env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else
      {
        LALFree( *sumHead );
        fprintf( stderr, "unknown column while parsing\n" );
        return -1;
      }
    }

    /* increase the count of rows parsed */
    ++nrows;
  }

  if ( mioStatus == -1 )
  {
    fprintf( stderr, "error parsing after row %d\n", i );
    LALFree( *sumHead );
    MetaioClose( env );
    return -1;
  }

  /* we have sucesfully parsed table */
  MetaioClose( env );
  return nrows;  
}
