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
\idx{LALStochasticTableFromLIGOLw()}
\idx{LALStochSummTableFromLIGOLw()}
    
\subsubsection*{Description}

The routine \verb+LALCreateMetaTableDir+ constructs a \verb+MetaTableDirectory+
for a given LIGOLwXML table.  It determines the location of each column expected
to be present in the XML table and populates the \verb+pos+ field with this
information.  This then allows other routines to parse the contents of an XML
file and correctly interpret the entries.  Currently, this has been implemented
only for the \verb+sngl_burst+ and \verb+sim_burst+ tables.  When reading these
tables, a call is made to \verb+LALCreateMetaTableDir+.  For all other tables,
the directory is constructed internally by the reading code.

The routine \verb+LALSnglBurstTableFromLIGOLw+ reads in a \verb+single_burst+
table from a LIGOLwXML file  specified in \verb+fileName+; \verb+eventHead+
provides a pointer to the head of a linked list of \verb+SnglBurstTable+s
containing the events.  The routine is passed the \verb+fileName+ of an XML file
containing a \verb+sngl_burst+ table.  First, the table is opened using
\verb+MetaioOpenTable+.  Then a directory of the table is generated using
\verb+LALCreateMetaTableDir+.  Rows of the table are read in sequentially from
the file.  Each entry in the row is stored in the appopriate entry of a
\verb+SnglBurstTable+ which is appended to the end of a linked list of such
tables.  When all rows have been read in, the file is closed using
\verb+MetaioClose+.  \verb+eventHead+ is set to point to the head of the linked
list of \verb+SnglBurstTable+s.  

The routine \verb+LALSimBurstTableFromLIGOLw+ reads in a \verb+sim_burst+ table
from the LIGOLwXML file specified in \verb+fileName+; \verb+eventHead+ provides a pointer to the head of
a linked list of \verb+SimBurstTables+ containing the events.  It operates in a
similar manner to \verb+LALSnglBurstTableFromLIGOLw+.  Additionally, a
\verb+startTime+ and \verb+stopTime+ are be specified.  Only simulated events
occuring between these times are returned.  If the \verb+endTime+ is set to
zero, then all events are returned.

The routine \verb+LALSnglInspiralTableFromLIGOLw+ reads in a
\verb+sngl_inspiral+ table from the LIGOLwXML file specified in \verb+fileName+.
It returns the number of triggers read in and \verb+eventHead+ provides a
pointer to the head of a linked list of \verb+SnglInspiralTable+s containing the
events.  It will return all events between the \verb+startEvent+ and
\verb+stopEvent+; if these are set to 0 and -1 respectively, all events are
returned.

The routine \verb+InspiralTmpltBankFromLIGOLw+ reads in a \verb+sngl_inspiral+
table from the LIGOLwXML file specified in \verb+fileName+. It returns the
number of templates read in and \verb+bankHead+ provides a pointer to the head
of a linked list of \verb+InspiralTemplate+s containing the templates read in.
It will return all events between the \verb+startTmplt+ and \verb+stopTmplt+; if
these are set to 0 and -1 respectively, all events are returned.  Although a
\verb+sngl_inspiral+ table is read in, only those entries relevant for an
InspiralTemplate are read in and stored.

The routine \verb+SimInspiralTableFromLIGOLw+ reads in a \verb+sim_inspiral+
table from the LIGOLwXML file specified in \verb+fileName+.  It returns the
number of rows read in and \verb+SimHead+ provides a pointer to the head of a
linked list of \verb+SimInspiralTable+s containing the events.  Additionally, a
\verb+startTime+ and \verb+endTime+ are specified.  Only simulated events
occuring between these times are returned.  If the \verb+endTime+ is set to
zero, then all events are returned.

The routine \verb+SearchSummaryTableFromLIGOLw+ reads in a \verb+search_summary+
table from the LIGOLwXML file specified in \verb+fileName+.  It returns the
number of rows read in and \verb+sumHead+ provides a pointer to the head of a
linked list of \verb+SearchSummaryTable+s.

The routine \verb+SummValueTableFromLIGOLw+ reads in a \verb+summ_value+
table from the LIGOLwXML file specified in \verb+fileName+.  It returns the
number of rows read in and \verb+sumHead+ provides a pointer to the head of a
linked list of \verb+SummValueTable+s.

The routine \verb+LALStochasticTableFromLIGOLw+ reads in a
\verb+stochastic_table+ table from the LIGOLwXML file specified in
\verb+fileName+.  It returns the number of rows read in and
\verb+stochHead+ provides a pointer to the head of a linked list of
\verb+StochasticTable+s.

The routine \verb+LALStochSummTableFromLIGOLw+ reads in a
\verb+stoch_summ_table+ table from the LIGOLwXML file specified in
\verb+fileName+.  It returns the number of rows read in and
\verb+stochSummHead+ provides a pointer to the head of a linked list of
\verb+StochSummTable+s.


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
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
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
          LALSnprintf( thisEvent->waveform, 
              LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s", 
              env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
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
    {"mtotal",                  -1, 16},
    {"eta",                     -1, 17},
    {"tau0",                    -1, 18},
    {"tau2",                    -1, 19},
    {"tau3",                    -1, 20},
    {"tau4",                    -1, 21},
    {"tau5",                    -1, 22},
    {"ttotal",                  -1, 23},
    {"psi0",                    -1, 24},
    {"psi3",                    -1, 25},
    {"alpha",                   -1, 26},
    {"alpha1",                  -1, 27},
    {"alpha2",                  -1, 28},
    {"alpha3",                  -1, 29},
    {"alpha4",                  -1, 30},
    {"alpha5",                  -1, 31},
    {"alpha6",                  -1, 32},
    {"beta",                    -1, 33},
    {"f_final",                 -1, 34}, 
    {"snr",                     -1, 35},
    {"chisq",                   -1, 36},
    {"chisq_dof",               -1, 37},
    {"sigmasq",                 -1, 38},
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
          thisEvent->mtotal = r4colData;
        }
        else if ( tableDir[j].idx == 17 )
        {
          thisEvent->eta = r4colData;
        }
        else if ( tableDir[j].idx == 18 )
        {
          thisEvent->tau0 = r4colData;
        }
        else if ( tableDir[j].idx == 19 )
        {
          thisEvent->tau2 = r4colData;
        }
        else if ( tableDir[j].idx == 20 )
        {
          thisEvent->tau3 = r4colData;
        }
        else if ( tableDir[j].idx == 21 )
        {
          thisEvent->tau4 = r4colData;
        }
        else if ( tableDir[j].idx == 22 )
        {
          thisEvent->tau5 = r4colData;
        }
        else if ( tableDir[j].idx == 23 )
        {
          thisEvent->ttotal = r4colData;
        }
        else if ( tableDir[j].idx == 24 )
        {
          thisEvent->psi0 = r4colData;
        }
        else if ( tableDir[j].idx == 25 )
        {
          thisEvent->psi3 = r4colData;
        }
        else if ( tableDir[j].idx == 26 )
        {
          thisEvent->alpha = r4colData;
        }
        else if ( tableDir[j].idx == 27 )
        {
          thisEvent->alpha1 = r4colData;
        }
        else if ( tableDir[j].idx == 28 )
        {
          thisEvent->alpha2 = r4colData;
        }
        else if ( tableDir[j].idx == 29 )
        {
          thisEvent->alpha3 = r4colData;
        }
        else if ( tableDir[j].idx == 30 )
        {
          thisEvent->alpha4 = r4colData;
        }
        else if ( tableDir[j].idx == 31 )
        {
          thisEvent->alpha5 = r4colData;
        }
        else if ( tableDir[j].idx == 32 )
        {
          thisEvent->alpha6 = r4colData;
        }
        else if ( tableDir[j].idx == 33 )
        {
          thisEvent->beta = r4colData;
        }
        else if ( tableDir[j].idx == 34 )
        {
          thisEvent->f_final = r4colData;
        }
        else if ( tableDir[j].idx == 35 )
        {
          thisEvent->snr = r4colData;
        }
        else if ( tableDir[j].idx == 36 )
        {
          thisEvent->chisq = r4colData;
        }
        else if ( tableDir[j].idx == 37 )
        {
          thisEvent->chisq_dof = i4colData;
        }
        else if ( tableDir[j].idx == 38 )
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
    {"beta",    -1, 12},
    {"f_final", -1, 13},
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
          thisTmplt->beta = colData;
        }
        else if ( tableDir[j].idx == 13 )
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
    {"waveform",            -1, 0},
    {"geocent_end_time",    -1, 1},
    {"geocent_end_time_ns", -1, 2},
    {"h_end_time",          -1, 3},
    {"h_end_time_ns",       -1, 4},
    {"l_end_time",          -1, 5},
    {"l_end_time_ns",       -1, 6},
    {"g_end_time",          -1, 7},
    {"g_end_time_ns",       -1, 8},
    {"t_end_time",          -1, 9},
    {"t_end_time_ns",       -1, 10},
    {"v_end_time",          -1, 11},
    {"v_end_time_ns",       -1, 12},
    {"end_time_gmst",       -1, 13},
    {"source",              -1, 14},
    {"mass1",               -1, 15},
    {"mass2",               -1, 16},
    {"eta",                 -1, 17},
    {"distance",            -1, 18},
    {"longitude",           -1, 19},
    {"latitude",            -1, 20},
    {"inclination",         -1, 21},
    {"coa_phase",           -1, 22},
    {"polarization",        -1, 23},
    {"psi0",                -1, 24},
    {"psi3",                -1, 25},
    {"alpha",               -1, 26},
    {"alpha1",              -1, 27},
    {"alpha2",              -1, 28},
    {"alpha3",              -1, 29},
    {"alpha4",              -1, 30},
    {"alpha5",              -1, 31},
    {"alpha6",              -1, 32},
    {"beta",                -1, 33},
    {"spin1x",              -1, 34},
    {"spin1y",              -1, 35},
    {"spin1z",              -1, 36},
    {"spin2x",              -1, 37},
    {"spin2y",              -1, 38},
    {"spin2z",              -1, 39},
    {"theta0",              -1, 40},
    {"phi0",                -1, 41},
    {"fLower",             -1, 42},
    {"fFinal",             -1, 43},
    {"mchirp",              -1, 44},
    {"eff_dist_h",          -1, 45},
    {"eff_dist_l",          -1, 46},
    {"eff_dist_g",          -1, 47},
    {"eff_dist_t",          -1, 48},
    {"eff_dist_v",          -1, 49},
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
          thisSim->g_end_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 8 )
        {
          thisSim->g_end_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 9 )    
        {
          thisSim->t_end_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 10 )
        {
          thisSim->t_end_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 11 )     
        {
          thisSim->v_end_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 12 )
        {
          thisSim->v_end_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 13 )
        {
          thisSim->end_time_gmst = r8colData;
        }
        else if ( tableDir[j].idx == 14 )
        {
          LALSnprintf(thisSim->source, LIGOMETA_SOURCE_MAX * sizeof(CHAR),
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
        }
        else if ( tableDir[j].idx == 15 )
        {
          thisSim->mass1 = r4colData;
        }
        else if ( tableDir[j].idx == 16 )
        {
          thisSim->mass2 = r4colData;
        }
        else if ( tableDir[j].idx == 17 )
        {
          thisSim->eta = r4colData;
        }
        else if ( tableDir[j].idx == 18 )
        {
          thisSim->distance = r4colData;
        }
        else if ( tableDir[j].idx == 19 )
        {
          thisSim->longitude = r4colData;
        }
        else if ( tableDir[j].idx == 20 )
        {
          thisSim->latitude = r4colData;
        }
        else if ( tableDir[j].idx == 21 )
        {
          thisSim->inclination = r4colData;
        }
        else if ( tableDir[j].idx == 22 )
        {
          thisSim->coa_phase = r4colData;
        } 
	else if ( tableDir[j].idx == 23 )
	  {
	    thisSim->polarization = r4colData;
	  }
	else if ( tableDir[j].idx == 24 )
	  {
	    thisSim->psi0 = r4colData;
	  }
	else if ( tableDir[j].idx == 25 )
	  {
	    thisSim->psi3 = r4colData;
	  }
	else if ( tableDir[j].idx == 26 )
	  {
	    thisSim->alpha = r4colData;
	  }
	else if ( tableDir[j].idx == 27 )
	  {
	    thisSim->alpha1 = r4colData;
	  }
	else if ( tableDir[j].idx == 28 )
	  {
	    thisSim->alpha2 = r4colData;
	  }	
	else if ( tableDir[j].idx == 29 )
	  {
	    thisSim->alpha3 = r4colData;
	  }	
	else if ( tableDir[j].idx == 30 )
	  {
	    thisSim->alpha4 = r4colData;
	  }
	else if ( tableDir[j].idx == 31 )
	  {
	    thisSim->alpha5 = r4colData;
	  }
	else if ( tableDir[j].idx == 32 )
	  {
	    thisSim->alpha6 = r4colData;
	  }
	else if ( tableDir[j].idx == 33 )
	  {
	    thisSim->beta = r4colData;
	  }    
	else if ( tableDir[j].idx == 34 )
	  {
	    thisSim->spin1x = r4colData;
	  }
	else if ( tableDir[j].idx == 35 )
	  {
	    thisSim->spin1y = r4colData;
	  }
	else if ( tableDir[j].idx == 36 )
	  {
	    thisSim->spin1z = r4colData;
	  }
	else if ( tableDir[j].idx == 37 )
	  {
	    thisSim->spin2x = r4colData;
	  }
	else if ( tableDir[j].idx == 38 )
	  {
	    thisSim->spin2y = r4colData;
	  }
	else if ( tableDir[j].idx == 39 )
	  {
	    thisSim->spin2z = r4colData;
	  }
	else if ( tableDir[j].idx == 40 )
	  {
	    thisSim->theta0 = r4colData;
	  }
	else if ( tableDir[j].idx == 41 )
	  {
	    thisSim->phi0 = r4colData;
	  }
	else if ( tableDir[j].idx == 42 )
	  {
	    thisSim->fLower = r4colData;
	  }
	else if ( tableDir[j].idx == 43 )
	  {
	    thisSim->fFinal = r4colData;
	  }
	else if ( tableDir[j].idx == 43 )
	  {
	    thisSim->mchirp = r4colData;
	  }
        else if ( tableDir[j].idx == 44 )
	  {
	    thisSim->eff_dist_h = r4colData;
	  }
        else if ( tableDir[j].idx == 45 )
	  {
	    thisSim->eff_dist_l = r4colData;
	  }
        else if ( tableDir[j].idx == 46 )
	  {
	    thisSim->eff_dist_g = r4colData;
	  }
        else if ( tableDir[j].idx == 47 )
	  {
	    thisSim->eff_dist_t = r4colData;
	  }
        else if ( tableDir[j].idx == 48 )
	  {
	    thisSim->eff_dist_v = r4colData;
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


#define CLOBBER_VAL \
    while ( *sumHead ); \
    { \
      thisValue = *sumHead; \
      *sumHead = (*sumHead)->next; \
      LALFree( thisValue ); \
      thisValue = NULL; \
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
  SummValueTable                       *thisValue = NULL;
  struct MetaioParseEnvironment         parseEnv;
  const  MetaioParseEnv                 env = &parseEnv;
  MetaTableDirectory tableDir[] =
  {
    {"program",               -1, 0},
    {"start_time",            -1, 1},
    {"start_time_ns",         -1, 2},
    {"end_time",              -1, 3},
    {"end_time_ns",           -1, 4},
    {"ifo",                   -1, 5},
    {"name",                  -1, 6},
    {"value",                 -1, 7},
    {"comment",               -1, 8},
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

  /* loop over the rows in the file */
  i = nrows = 0;
  while ( (mioStatus = MetaioGetRow(env)) == 1 ) 
  {
    /* allocate memory for the table */
    if ( ! *sumHead )
    {
      thisValue = *sumHead = (SummValueTable *) 
    LALCalloc( 1, sizeof(SummValueTable) );
    }
    else
    {
      thisValue = thisValue->next = (SummValueTable *) 
          LALCalloc( 1, sizeof(SummValueTable) );
    }
    if ( ! thisValue )
    {
      fprintf( stderr, "could not allocate summ value\n" );
      CLOBBER_VAL;
      MetaioClose( env );
      return -1;
    }

    /* parse the rows into the SummValue structure */
    for ( j = 0; tableDir[j].name; ++j )
    {
      REAL4 r4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
      INT4  i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;

      if ( tableDir[j].idx == 0 )
      {
        LALSnprintf( thisValue->program, LIGOMETA_PROGRAM_MAX * sizeof(CHAR),
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 1 )
      {
        thisValue->start_time.gpsSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 2 )
      {
        thisValue->start_time.gpsNanoSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 3 )
      {
        thisValue->end_time.gpsSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 4 )
      {
        thisValue->end_time.gpsNanoSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 5 )
      {
        LALSnprintf( thisValue->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR),
          "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 6 )
      {
        LALSnprintf( thisValue->name, LIGOMETA_SUMMVALUE_NAME_MAX * 
          sizeof(CHAR), "%s",
          env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 7 )
      {
        thisValue->value = r4colData;
      }
      else if ( tableDir[j].idx == 8 )
      {
        LALSnprintf( thisValue->comment, LIGOMETA_SUMMVALUE_NAME_MAX * 
          sizeof(CHAR), "%s",
          env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else
      {
        CLOBBER_VAL;
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
    CLOBBER_VAL;
    MetaioClose( env );
    return -1;
  }

  /* we have sucesfully parsed table */
  MetaioClose( env );
  return nrows;  
}

#undef CLOBBER_VAL

#define CLOBBER_STOCH_VAL \
    while (*stochHead); \
    { \
      thisValue = *stochHead; \
      *stochHead = (*stochHead)->next; \
      LALFree( thisValue ); \
      thisValue = NULL; \
    }

/* <lalVerbatim file="LIGOLwXMLReadCP"> */
int
LALStochasticTableFromLIGOLw (
    StochasticTable **stochHead,
    CHAR *fileName)
/* </lalVerbatim> */
{
  int i, j, nrows;
  int mioStatus;
  StochasticTable *thisValue = NULL;

  struct MetaioParseEnvironment parseEnv;
  const MetaioParseEnv env = &parseEnv;
  MetaTableDirectory tableDir[] =
  {
    {"ifo_one",       -1,  0},
    {"ifo_two",       -1,  1},
    {"channel_one",   -1,  2},
    {"channel_two",   -1,  3},
    {"start_time",    -1,  4},
    {"start_time_ns", -1,  5},
    {"duration",      -1,  6},
    {"duration_ns",   -1,  7},
    {"f_min",         -1,  8},
    {"f_max",         -1,  9},
    {"cc_stat",       -1, 10},
    {"cc_sigma",      -1, 11},
    {NULL,             0,  0}
  };

  /* check that the table handle and pointer are valid */
  if (!stochHead)
  {
    fprintf(stderr, "null pointer passed as handle to stochastic value\n");
    return -1;
  }
  if (*stochHead)
  {
    fprintf(stderr, "non-null pointer passed as pointer to stochastic value\n");
    return -1;
  }

  /* open the stochastic)table in the file file */
  mioStatus = MetaioOpenTable(env, fileName, "stochastic");
  if (mioStatus)
  {
    fprintf(stderr, "error opening stochastic table from file %s\n", \
        fileName);
    return -1;
  }

  /* figure out the column positions */
  for (i = 0; tableDir[i].name; ++i)
  {
    if ((tableDir[i].pos = MetaioFindColumn(env, tableDir[i].name)) < 0)
    {
      fprintf(stderr, "unable to find column %s\n", tableDir[i].name);
      MetaioClose(env);
      return -1;
    }
  }

  /* loop over the rows in the file */
  i = nrows = 0;
  while ((mioStatus = MetaioGetRow(env)) == 1) 
  {
    /* allocate memory for the table */
    if (!*stochHead)
    {
      thisValue = *stochHead = (StochasticTable *) \
                  LALCalloc(1, sizeof(StochasticTable));
    }
    else
    {
      thisValue = thisValue->next = (StochasticTable *) \
                  LALCalloc( 1, sizeof(StochasticTable) );
    }
    if (!thisValue)
    {
      fprintf(stderr, "could not allocate stochastic table\n");
      CLOBBER_STOCH_VAL;
      MetaioClose(env);
      return -1;
    }

    /* parse the rows into the StochasticTable structure */
    for ( j = 0; tableDir[j].name; ++j )
    {
      REAL8 r8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_8;
      INT4  i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;

      if (tableDir[j].idx == 0)
      {
        LALSnprintf(thisValue->ifo_one, LIGOMETA_IFO_MAX * sizeof(CHAR), \
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if (tableDir[j].idx == 1)
      {
        LALSnprintf(thisValue->ifo_two, LIGOMETA_IFO_MAX * sizeof(CHAR), \
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if ( tableDir[j].idx == 2 )
      {
        LALSnprintf(thisValue->channel_one, LIGOMETA_CHANNEL_MAX * \
            sizeof(CHAR), "%s", \
            env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if ( tableDir[j].idx == 3 )
      {
        LALSnprintf(thisValue->channel_two, LIGOMETA_CHANNEL_MAX * \
            sizeof(CHAR), "%s", \
            env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if ( tableDir[j].idx == 4 )
      {
        thisValue->start_time.gpsSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 5 )
      {
        thisValue->start_time.gpsNanoSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 6 )
      {
        thisValue->duration.gpsSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 7 )
      {
        thisValue->duration.gpsNanoSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 8 )
      {
        thisValue->f_min = r8colData;
      }
      else if ( tableDir[j].idx == 9 )
      {
        thisValue->f_max = r8colData;
      }
      else if ( tableDir[j].idx == 10 )
      {
        thisValue->cc_stat = r8colData;
      }
      else if ( tableDir[j].idx == 11 )
      {
        thisValue->cc_sigma = r8colData;
      }
      else
      {
        CLOBBER_STOCH_VAL;
        fprintf(stderr, "unknown column while parsing\n");
        return -1;
      }
    }

    /* increase the count of rows parsed */
    ++nrows;
  }

  if (mioStatus == -1)
  {
    fprintf(stderr, "error parsing after row %d\n", i);
    CLOBBER_STOCH_VAL;
    MetaioClose(env);
    return -1;
  }

  /* we have sucesfully parsed table */
  MetaioClose(env);
  return nrows;  
}

#undef CLOBBER_STOCH_VAL

#define CLOBBER_STOCH_SUMM_VAL \
    while (*stochSummHead); \
    { \
      thisValue = *stochSummHead; \
      *stochSummHead = (*stochSummHead)->next; \
      LALFree( thisValue ); \
      thisValue = NULL; \
    }

/* <lalVerbatim file="LIGOLwXMLReadCP"> */
int
LALStochSummTableFromLIGOLw (
    StochSummTable **stochSummHead,
    CHAR *fileName)
/* </lalVerbatim> */
{
  int i, j, nrows;
  int mioStatus;
  StochSummTable *thisValue = NULL;

  struct MetaioParseEnvironment parseEnv;
  const MetaioParseEnv env = &parseEnv;
  MetaTableDirectory tableDir[] =
  {
    {"ifo_one",       -1,  0},
    {"ifo_two",       -1,  1},
    {"channel_one",   -1,  2},
    {"channel_two",   -1,  3},
    {"start_time",    -1,  4},
    {"start_time_ns", -1,  5},
    {"end_time",      -1,  6},
    {"end_time_ns",   -1,  7},
    {"f_min",         -1,  8},
    {"f_max",         -1,  9},
    {"y_opt",         -1, 10},
    {"error",         -1, 11},
    {NULL,             0,  0}
  };

  /* check that the table handle and pointer are valid */
  if (!stochSummHead)
  {
    fprintf(stderr, "null pointer passed as handle to stoch_summ value\n");
    return -1;
  }
  if (*stochSummHead)
  {
    fprintf(stderr, "non-null pointer passed as pointer to stoch_summ value\n");
    return -1;
  }

  /* open the stoch_summ_table in the file file */
  mioStatus = MetaioOpenTable(env, fileName, "stoch_summ");
  if (mioStatus)
  {
    fprintf(stderr, "error opening stoch_summ table from file %s\n", \
        fileName);
    return -1;
  }

  /* figure out the column positions */
  for (i = 0; tableDir[i].name; ++i)
  {
    if ((tableDir[i].pos = MetaioFindColumn(env, tableDir[i].name)) < 0)
    {
      fprintf(stderr, "unable to find column %s\n", tableDir[i].name);
      MetaioClose(env);
      return -1;
    }
  }

  /* loop over the rows in the file */
  i = nrows = 0;
  while ((mioStatus = MetaioGetRow(env)) == 1) 
  {
    /* allocate memory for the table */
    if (!*stochSummHead)
    {
      thisValue = *stochSummHead = (StochSummTable *) \
                  LALCalloc(1, sizeof(StochSummTable));
    }
    else
    {
      thisValue = thisValue->next = (StochSummTable *) \
                  LALCalloc( 1, sizeof(StochSummTable) );
    }
    if (!thisValue)
    {
      fprintf(stderr, "could not allocate stoch_summ table\n");
      CLOBBER_STOCH_SUMM_VAL;
      MetaioClose(env);
      return -1;
    }

    /* parse the rows into the StochSummTable structure */
    for ( j = 0; tableDir[j].name; ++j )
    {
      REAL8 r8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_8;
      INT4  i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;

      if (tableDir[j].idx == 0)
      {
        LALSnprintf(thisValue->ifo_one, LIGOMETA_IFO_MAX * sizeof(CHAR), \
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if (tableDir[j].idx == 1)
      {
        LALSnprintf(thisValue->ifo_two, LIGOMETA_IFO_MAX * sizeof(CHAR), \
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if ( tableDir[j].idx == 2 )
      {
        LALSnprintf(thisValue->channel_one, LIGOMETA_CHANNEL_MAX * \
            sizeof(CHAR), "%s", \
            env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if ( tableDir[j].idx == 3 )
      {
        LALSnprintf(thisValue->channel_two, LIGOMETA_CHANNEL_MAX * \
            sizeof(CHAR), "%s", \
            env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if ( tableDir[j].idx == 4 )
      {
        thisValue->start_time.gpsSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 5 )
      {
        thisValue->start_time.gpsNanoSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 6 )
      {
        thisValue->end_time.gpsSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 7 )
      {
        thisValue->end_time.gpsNanoSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 8 )
      {
        thisValue->f_min = r8colData;
      }
      else if ( tableDir[j].idx == 9 )
      {
        thisValue->f_max = r8colData;
      }
      else if ( tableDir[j].idx == 10 )
      {
        thisValue->y_opt = r8colData;
      }
      else if ( tableDir[j].idx == 11 )
      {
        thisValue->error = r8colData;
      }
      else
      {
        CLOBBER_STOCH_SUMM_VAL;
        fprintf(stderr, "unknown column while parsing\n");
        return -1;
      }
    }

    /* increase the count of rows parsed */
    ++nrows;
  }

  if (mioStatus == -1)
  {
    fprintf(stderr, "error parsing after row %d\n", i);
    CLOBBER_STOCH_SUMM_VAL;
    MetaioClose(env);
    return -1;
  }

  /* we have sucesfully parsed table */
  MetaioClose(env);
  return nrows;  
}

#undef CLOBBER_STOCH_SUMM_VAL


#define CLOBBER_EVENTS \
    while ( *eventHead ); \
    { \
      thisEvent = *eventHead; \
      *eventHead = (*eventHead)->next; \
      LALFree( thisEvent ); \
      thisEvent = NULL; \
    }



/* <lalVerbatim file="LIGOLwXMLReadCP"> */
int
LALExtTriggerTableFromLIGOLw (
    ExtTriggerTable   **eventHead,
    CHAR               *fileName,
    INT4                startEvent,
    INT4                stopEvent
    )
/* </lalVerbatim> */
{
  int                                   i, j, nrows;
  int                                   mioStatus;
  ExtTriggerTable                       *thisEvent = NULL;
  struct MetaioParseEnvironment         parseEnv;
  const  MetaioParseEnv                 env = &parseEnv;
  MetaTableDirectory tableDir[] =
  {
    {"det_alts",               -1,  0},
    {"det_band",               -1,  1},
    {"det_fluence",            -1,  2},
    {"det_fluence_int",        -1,  3},
    {"det_name",               -1,  4},
    {"det_peak",               -1,  5},
    {"det_peak_int",           -1,  6},
    {"det_snr",                -1,  7},
    {"email_time",             -1,  8},
    {"event_dec",              -1,  9},
    {"event_dec_err",          -1, 10},
    {"event_epoch",            -1, 11},
    {"event_err_type",         -1, 12},
    {"event_ra",               -1, 13},
    {"event_ra_err",           -1, 14},
    {"start_time",             -1, 15},
    {"start_time_ns",          -1, 16},
    {"event_type",             -1, 17},
    {"event_z",                -1, 18},
    {"event_z_err",            -1, 19},
    {"notice_comments",        -1, 20},
    {"notice_id",              -1, 21},
    {"notice_sequence",        -1, 22},
    {"notice_time",            -1, 23},
    {"notice_type",            -1, 24},
    {"notice_url",             -1, 25},
    {"obs_fov_dec",            -1, 26},
    {"obs_fov_dec_width",      -1, 27},
    {"obs_fov_ra",             -1, 28},
    {"obs_fov_ra_width",       -1, 29},   
    {"obs_loc_ele",            -1, 30},
    {"obs_loc_lat",            -1, 31},
    {"obs_loc_long",           -1, 32},
    {NULL,                      0, 0}
  };


  /* check that the bank handle and pointer are void */
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
  mioStatus = MetaioOpenTable( env, fileName, "external_trigger" );
  if ( mioStatus )
  {
    fprintf( stdout, "no ext_trigger table in file %s\n", fileName );
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
        thisEvent = *eventHead = (ExtTriggerTable*) 
          LALCalloc( 1, sizeof(ExtTriggerTable) );
      }
      else
      {
        thisEvent = thisEvent->next = (ExtTriggerTable *) 
          LALCalloc( 1, sizeof(ExtTriggerTable) );
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
          LALSnprintf( thisEvent->det_alts, LIGOMETA_STD * sizeof(CHAR), 
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 1 )
        {
          LALSnprintf( thisEvent->det_band, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 2 )
        {
          LALSnprintf( thisEvent->det_fluence, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }  
	else if ( tableDir[j].idx == 3 )
        {
          LALSnprintf( thisEvent->det_fluence_int, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }  
	else if ( tableDir[j].idx == 4 )
        {
          LALSnprintf( thisEvent->det_name, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
	else if ( tableDir[j].idx == 5 )
        {
          LALSnprintf( thisEvent->det_peak, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }  
	else if ( tableDir[j].idx == 6 )
        {
          LALSnprintf( thisEvent->det_peak_int, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }  
	else if ( tableDir[j].idx == 7 )
	  {
	    LALSnprintf( thisEvent->det_snr, LIGOMETA_STD * sizeof(CHAR),
			 "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
	  } 
	else if ( tableDir[j].idx == 8 )
	  {
	    thisEvent->email_time = i4colData;
	  }
	else if ( tableDir[j].idx == 9 )
	  {
	    thisEvent->event_dec = r4colData;
	  }
	else if ( tableDir[j].idx == 10 )
	  {
	    thisEvent->event_dec_err = r4colData;
	  }
	else if ( tableDir[j].idx == 11 )
	  {  
	    LALSnprintf( thisEvent->event_epoch, LIGOMETA_STD * sizeof(CHAR),
			 "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
	  }
	else if ( tableDir[j].idx == 12 )
	  {
	    LALSnprintf( thisEvent->event_err_type, LIGOMETA_STD * sizeof(CHAR),
			 "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
	  }
	else if ( tableDir[j].idx == 13 )
	  {
	    thisEvent->event_ra = r4colData;
	  }
	else if ( tableDir[j].idx == 14 )
	  {
	    thisEvent->event_ra_err = r4colData;
	  }
	else if ( tableDir[j].idx == 15 )
	  {
	    thisEvent->start_time = i4colData;
	    /*  printf("start time:%d\n",i4colData); */
	  }
	else if ( tableDir[j].idx == 16 )
	  {
	    thisEvent->start_time_ns = i4colData;
	  }	
	else if ( tableDir[j].idx == 17 )
	  {
	    LALSnprintf( thisEvent->event_type, LIGOMETA_STD * sizeof(CHAR),
			 "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
	  }
	else if ( tableDir[j].idx == 18 )
	  {
	    thisEvent->event_z = r4colData;
	  }
	else if ( tableDir[j].idx == 19 )
	  {
	    thisEvent->event_z_err = r4colData;
	  }
	else if ( tableDir[j].idx == 20 )
	  {
	    LALSnprintf( thisEvent->notice_comments, LIGOMETA_STD * sizeof(CHAR),
			 "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
	  }	
	else if ( tableDir[j].idx == 21 )
	  {
	    LALSnprintf( thisEvent->notice_id, LIGOMETA_STD * sizeof(CHAR),
			 "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
	  }
	else if ( tableDir[j].idx == 22 )
	  {
	    LALSnprintf( thisEvent->notice_sequence, LIGOMETA_STD * sizeof(CHAR),
			 "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
	  }
	else if ( tableDir[j].idx == 23 )
	  {
	    thisEvent->notice_time = i4colData;
	  }
	else if ( tableDir[j].idx == 24 )
	  {
	    LALSnprintf( thisEvent->notice_type, LIGOMETA_STD * sizeof(CHAR),
			 "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
	  }
	else if ( tableDir[j].idx == 25 )
	  {
	    LALSnprintf( thisEvent->notice_url, LIGOMETA_STD * sizeof(CHAR),
			 "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
	  }	
        else if ( tableDir[j].idx == 26 )
	  {
	    thisEvent->obs_fov_dec = r4colData;
	  }
        else if ( tableDir[j].idx == 27 )
	  {
	    thisEvent->obs_fov_dec_width = r4colData;
	  }
        else if ( tableDir[j].idx == 28 )
	  {
	    thisEvent->obs_fov_ra = r4colData;
	  }
        else if ( tableDir[j].idx == 29 )
	  {
	    thisEvent->obs_fov_ra_width = i4colData;
	  }
        else if ( tableDir[j].idx == 30 )
	  {
	    thisEvent->obs_loc_ele = r4colData;
	  }
        else if ( tableDir[j].idx == 31 )
	  {
	    thisEvent->obs_loc_lat = r4colData;
	  }
        else if ( tableDir[j].idx == 32 )
	  {
	    thisEvent->obs_loc_long = r4colData;
	  }
        else
        {
          CLOBBER_EVENTS;
          fprintf( stderr, "unknown column while parsing ext_trigger\n" );
          return -1;
        }
      }

      /* count the number of template parsed */
      nrows++;
    }
  }

  /* must be reduced to avoid stopping psocessing with triggers.xml
     because that file is generated corrupted (by just adding new triggers
     in new lines */
  /*
  if ( mioStatus == -1 )
  {
    fprintf( stderr, "error parsing after row %d\n", i );
    CLOBBER_EVENTS;
    MetaioClose( env );
    return -1;
  }
  */

  /* we have sucesfully parsed temples */
  MetaioClose( env );
  return nrows;  
}


#undef CLOBBER_EVENTS
