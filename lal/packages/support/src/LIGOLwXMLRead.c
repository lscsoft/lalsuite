/*----------------------------------------------------------------------- 
 * 
 * File Name: ligolwbank.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

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


void
LALCreateMetaTableDir(
    LALStatus              *status,
    MetaTableDirectory    **tableDir,
    const MetaioParseEnv    env,
    MetadataTableType       table
    )
{
  INT4 i;

  INITSTATUS( status, "LALCreateMetaTableDir", LIGOLWXMLREADC );
  ATTATCHSTATUSPTR (status);

  /* check the inputs */
  ASSERT( *tableDir, status, LIGOLWXMLREADH_ENNUL, LIGOLWXMLREADH_MSGENNUL );

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

void
LALSnglBurstTableFromLIGOLw (
    LALStatus          *status,
    SnglBurstTable    **eventHead,
    CHAR               *fileName
    )
{
  int                                   i, j, nrows;
  int                                   mioStatus;
  SnglBurstTable                    *thisEvent = NULL;
  struct MetaioParseEnvironment         parseEnv;
  const  MetaioParseEnv                 env = &parseEnv;
  MetaTableDirectory                   *tableDir;

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
      fprintf( stderr, "could not allocate inspiral template\n" );
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


int
LALSnglInspiralTableFromLIGOLw (
    SnglInspiralTable **eventHead,
    CHAR               *fileName,
    INT4                startEvent,
    INT4                stopEvent
    )
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
    {"impulse_time",            -1, 5},
    {"impulse_time_ns",         -1, 6},
    {"template_duration",       -1, 7},
    {"event_duration",          -1, 8},
    {"amplitude",               -1, 9},
    {"eff_distance",            -1, 10},
    {"coa_phase",               -1, 11},
    {"mass1",                   -1, 12},
    {"mass2",                   -1, 13},
    {"mchirp",                  -1, 14},
    {"eta",                     -1, 15},
    {"tau0",                    -1, 16},
    {"tau2",                    -1, 17},
    {"tau3",                    -1, 18},
    {"tau4",                    -1, 19},
    {"tau5",                    -1, 20},
    {"ttotal",                  -1, 21},
    {"snr",                     -1, 22},
    {"chisq",                   -1, 23},
    {"chisq_dof",               -1, 24},
    {"sigmasq",                 -1, 25},
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
          thisEvent->impulse_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 6 )
        {
          thisEvent->impulse_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 7 )
        {
          thisEvent->template_duration = r8colData;
        }
        else if ( tableDir[j].idx == 8 )
        {
          thisEvent->event_duration = r8colData;
        }
        else if ( tableDir[j].idx == 9 )
        {
          thisEvent->amplitude = r4colData;
        }
        else if ( tableDir[j].idx == 10 )
        {
          thisEvent->eff_distance = r4colData;
        }
        else if ( tableDir[j].idx == 11 )
        {
          thisEvent->coa_phase = r4colData;
        }
        else if ( tableDir[j].idx == 12 )
        {
          thisEvent->mass1 = r4colData;
        }
        else if ( tableDir[j].idx == 13 )
        {
          thisEvent->mass2 = r4colData;
        }
        else if ( tableDir[j].idx == 14 )
        {
          thisEvent->mchirp = r4colData;
        }
        else if ( tableDir[j].idx == 15 )
        {
          thisEvent->eta = r4colData;
        }
        else if ( tableDir[j].idx == 16 )
        {
          thisEvent->tau0 = r4colData;
        }
        else if ( tableDir[j].idx == 17 )
        {
          thisEvent->tau2 = r4colData;
        }
        else if ( tableDir[j].idx == 18 )
        {
          thisEvent->tau3 = r4colData;
        }
        else if ( tableDir[j].idx == 19 )
        {
          thisEvent->tau4 = r4colData;
        }
        else if ( tableDir[j].idx == 20 )
        {
          thisEvent->tau5 = r4colData;
        }
        else if ( tableDir[j].idx == 21 )
        {
          thisEvent->ttotal = r4colData;
        }
        else if ( tableDir[j].idx == 22 )
        {
          thisEvent->snr = r4colData;
        }
        else if ( tableDir[j].idx == 23 )
        {
          thisEvent->chisq = r4colData;
        }
        else if ( tableDir[j].idx == 24 )
        {
          thisEvent->chisq_dof = i4colData;
        }
        else if ( tableDir[j].idx == 25 )
        {
          thisEvent->sigmasq = r8colData;
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
