/*----------------------------------------------------------------------- 
 * 
 * File Name: LIGOLwXML.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="LIGOLwXMLCV">
Author: Brown, D. A.
$Id$
</lalVerbatim> 
#endif

#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLHeaders.h>

NRCSID( LIGOLWXMLC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{LIGOLwXML.c}}

Routines to write LIGO metadata database structures to LIGO lightweight XML
files.

\subsubsection*{Prototypes}
\input{LIGOLwXMLCP}
\idx{LALOpenLIGOLwXMLFile()}
\idx{LALCloseLIGOLwXMLFile()}
\idx{LALBeginLIGOLwXMLTable()}
\idx{LALEndLIGOLwXMLTable()}
\idx{LALWriteLIGOLwXMLTable()}

\subsubsection*{Description}

The routine \verb+LALOpenLIGOLwXMLFile+ calls the C standard library function
\verb+fopen+ to open a file specified by the \verb+path+ argument. The file is
truncated to zero length if already exists. The standard LIGO lightweight XML
header is then written to the file and the the pointer to the file stream is
retuened in the \verb+xmlfp+ argument.

The routine \verb+LALCloseLIGOLwXMLFile+ prints the standard LIGO lightweight
XML footer and closes the file stream pointed to by \verb+xmlfp+.

The routine \verb+LALBeginLIGOLwXMLTable+

\subsubsection*{Algorithm}

None.

\subsubsection*{Uses}

\verb+fopen()+
\verb+fprintf()+
\verb+fclose()+

\subsubsection*{Notes}
 
%% Any relevant notes.
 
\vfill{\footnotesize\input{LIGOLwXMLCV}}

</lalLaTeX>
#endif

/* <lalVerbatim file="LIGOLwXMLCP"> */
void
LALOpenLIGOLwXMLFile (
    LALStatus          *status,
    LIGOLwXMLStream    *xml,
    const CHAR         *path
    )
/* </lalVerbatim> */
{
  /*  open the file and print the xml header */
  INITSTATUS( status, "LALOpenLIGOLwXMLFile", LIGOLWXMLC );
  ASSERT( xml, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  ASSERT( ! xml->fp, status, LIGOLWXMLH_ENNUL, LIGOLWXMLH_MSGENNUL );
  if ( ! (xml->fp = fopen( path, "w" )) )
  {
    ABORT( status, LIGOLWXMLH_EOPEN, LIGOLWXMLH_MSGEOPEN );
  }
  fprintf( xml->fp, LIGOLW_XML_HEADER );
  xml->table = no_table;
  RETURN( status );
}

/* <lalVerbatim file="LIGOLwXMLCP"> */
void
LALCloseLIGOLwXMLFile (
    LALStatus          *status,
    LIGOLwXMLStream    *xml
    )
/* </lalVerbatim> */
{
  /* print the xml footer and close the file handle */
  INITSTATUS( status, "LALCloseLIGOLwXMLFile", LIGOLWXMLC );
  ASSERT( xml, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  ASSERT( xml->fp, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  if ( xml->table != no_table )
  {
    ABORT( status, LIGOLWXMLH_ECLOS, LIGOLWXMLH_MSGECLOS );
  }
  fprintf( xml->fp, LIGOLW_XML_FOOTER );
  fclose( xml->fp );
  xml->fp = NULL;
  RETURN( status );
}

/* <lalVerbatim file="LIGOLwXMLCP"> */
void
LALBeginLIGOLwXMLTable (
    LALStatus           *status,
    LIGOLwXMLStream     *xml,
    MetadataTableType    table
    )
/* </lalVerbatim> */
{
  /* print the header for the xml table */
  INITSTATUS( status, "LALBeginLIGOLwXMLTable", LIGOLWXMLC );
  ASSERT( xml, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  ASSERT( xml->fp, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  if ( xml->table != no_table )
  {
    ABORT( status, LIGOLWXMLH_EBGNT, LIGOLWXMLH_MSGEBGNT );
  }

  switch( table )
  {
    case no_table:
      ABORT( status, LIGOLWXMLH_ENTAB, LIGOLWXMLH_MSGENTAB );
      break;
    case process_table:
      fprintf( xml->fp, LIGOLW_XML_PROCESS );
      break;
    case process_params_table:
      fprintf( xml->fp, LIGOLW_XML_PROCESS_PARAMS );
      break;
    case search_summary_table:
      fprintf( xml->fp, LIGOLW_XML_SEARCH_SUMMARY );
      break;
    case search_summvars_table:
      fprintf( xml->fp, LIGOLW_XML_SEARCH_SUMMVARS );
      break;
    case sngl_burst_table:
      fprintf( xml->fp, LIGOLW_XML_SNGL_BURST );
      break;
    case sngl_inspiral_table:
      fprintf( xml->fp, LIGOLW_XML_SNGL_INSPIRAL );
      break;
    case summ_value_table:
      fprintf( xml->fp, LIGOLW_XML_SUMM_VALUE );
      break;
    default:
      ABORT( status, LIGOLWXMLH_EUTAB, LIGOLWXMLH_MSGEUTAB );
  }
  xml->first = 1;
  xml->table = table;
  RETURN( status );
}

/* <lalVerbatim file="LIGOLwXMLCP"> */
void
LALEndLIGOLwXMLTable (
    LALStatus           *status,
    LIGOLwXMLStream     *xml
    )
/* </lalVerbatim> */
{
  /* print the header for the xml table */
  INITSTATUS( status, "LALEndLIGOLwXMLTable", LIGOLWXMLC );
  ASSERT( xml, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  ASSERT( xml->fp, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  if ( xml->table == no_table )
  {
    ABORT( status, LIGOLWXMLH_EENDT, LIGOLWXMLH_MSGEENDT );
  }
  fprintf( xml->fp, LIGOLW_XML_TABLE_FOOTER );
  xml->table = no_table;
  RETURN( status );
}

/* macro to print a comma on subsequent table rows */
#define FIRST_TABLE_ROW \
        if ( xml->first ) \
        { \
          xml->first = 0; \
        } else \
        { \
          fprintf( xml->fp, ",\n" ); \
        }

/* <lalVerbatim file="LIGOLwXMLCP"> */
void
LALWriteLIGOLwXMLTable (
    LALStatus           *status,
    LIGOLwXMLStream     *xml,
    MetadataTable        tablePtr,
    MetadataTableType    table
    )
/* </lalVerbatim> */
{
  /* print contents of the database struct into the xml table */
  INITSTATUS( status, "LALWriteLIGOLwXMLTable", LIGOLWXMLC );
  ASSERT( xml, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  ASSERT( xml->fp, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  if ( xml->table == no_table )
  {
    ABORT( status, LIGOLWXMLH_ETNOP, LIGOLWXMLH_MSGETNOP );
  }
  if ( xml->table != table )
  {
    ABORT( status, LIGOLWXMLH_ETMSM, LIGOLWXMLH_MSGETMSM );
  }
  switch( table )
  {
    case no_table:
      ABORT( status, LIGOLWXMLH_ENTAB, LIGOLWXMLH_MSGENTAB );
      break;
    case process_table:
      while( tablePtr.processTable )
      {
        FIRST_TABLE_ROW
        fprintf( xml->fp, PROCESS_ROW,
            tablePtr.processTable->program,
            tablePtr.processTable->version,
            tablePtr.processTable->cvs_repository,
            tablePtr.processTable->cvs_entry_time.gpsSeconds,
            tablePtr.processTable->comment,
            tablePtr.processTable->is_online,
            tablePtr.processTable->node,
            tablePtr.processTable->username,
            tablePtr.processTable->unix_procid,
            tablePtr.processTable->start_time.gpsSeconds,
            tablePtr.processTable->end_time.gpsSeconds,
            tablePtr.processTable->jobid,
            tablePtr.processTable->domain,
            tablePtr.processTable->ifos
            );
        tablePtr.processTable = tablePtr.processTable->next;
      }
      break;
    case process_params_table:
      while( tablePtr.processParamsTable )
      {
        FIRST_TABLE_ROW
        fprintf( xml->fp, PROCESS_PARAMS_ROW,
            tablePtr.processParamsTable->program,
            tablePtr.processParamsTable->param,
            tablePtr.processParamsTable->type,
            tablePtr.processParamsTable->value
            );
        tablePtr.processParamsTable = tablePtr.processParamsTable->next;
      }
      break;
    case search_summary_table:
      while( tablePtr.searchSummaryTable )
      {
        FIRST_TABLE_ROW
        fprintf( xml->fp, SEARCH_SUMMARY_ROW,
            tablePtr.searchSummaryTable->comment,
            tablePtr.searchSummaryTable->in_start_time.gpsSeconds,
            tablePtr.searchSummaryTable->in_start_time.gpsNanoSeconds,
            tablePtr.searchSummaryTable->in_end_time.gpsSeconds,
            tablePtr.searchSummaryTable->in_end_time.gpsNanoSeconds,
            tablePtr.searchSummaryTable->out_start_time.gpsSeconds,
            tablePtr.searchSummaryTable->out_start_time.gpsNanoSeconds,
            tablePtr.searchSummaryTable->out_end_time.gpsSeconds,
            tablePtr.searchSummaryTable->out_end_time.gpsNanoSeconds,
            tablePtr.searchSummaryTable->nevents,
            tablePtr.searchSummaryTable->nnodes
            );
        tablePtr.searchSummaryTable = tablePtr.searchSummaryTable->next;
      }
      break;
    case search_summvars_table:
      while( tablePtr.searchSummvarsTable )
      {
        FIRST_TABLE_ROW
        fprintf( xml->fp, SEARCH_SUMMVARS_ROW,
            tablePtr.searchSummvarsTable->name,
            tablePtr.searchSummvarsTable->string,
            tablePtr.searchSummvarsTable->value
            );
        tablePtr.searchSummvarsTable = tablePtr.searchSummvarsTable->next;
      }
      break;
    case sngl_burst_table:
      while( tablePtr.snglBurstTable )
      {
        FIRST_TABLE_ROW
        fprintf( xml->fp, SNGL_BURST_ROW,
            tablePtr.snglBurstTable->ifo,
            tablePtr.snglBurstTable->search,
            tablePtr.snglBurstTable->channel,
            tablePtr.snglBurstTable->start_time.gpsSeconds,
            tablePtr.snglBurstTable->start_time.gpsNanoSeconds,
            tablePtr.snglBurstTable->duration,
            tablePtr.snglBurstTable->central_freq,
            tablePtr.snglBurstTable->bandwidth,
            tablePtr.snglBurstTable->amplitude,
            tablePtr.snglBurstTable->snr,
            tablePtr.snglBurstTable->confidence
            );
        tablePtr.snglBurstTable = tablePtr.snglBurstTable->next;
      }
      break;
    case sngl_inspiral_table:
      while( tablePtr.snglInspiralTable )
      {
        FIRST_TABLE_ROW
        fprintf( xml->fp, SNGL_INSPIRAL_ROW,
            tablePtr.snglInspiralTable->ifo,
            tablePtr.snglInspiralTable->search,
            tablePtr.snglInspiralTable->channel,
            tablePtr.snglInspiralTable->end_time.gpsSeconds,
            tablePtr.snglInspiralTable->end_time.gpsNanoSeconds,
            tablePtr.snglInspiralTable->impulse_time.gpsSeconds,
            tablePtr.snglInspiralTable->impulse_time.gpsNanoSeconds,
            tablePtr.snglInspiralTable->template_duration,
            tablePtr.snglInspiralTable->event_duration,
            tablePtr.snglInspiralTable->amplitude,
            tablePtr.snglInspiralTable->eff_distance,
            tablePtr.snglInspiralTable->coa_phase,
            tablePtr.snglInspiralTable->mass1,
            tablePtr.snglInspiralTable->mass2,
            tablePtr.snglInspiralTable->mchirp,
            tablePtr.snglInspiralTable->eta,
            tablePtr.snglInspiralTable->tau0,
            tablePtr.snglInspiralTable->tau2,
            tablePtr.snglInspiralTable->tau3,
            tablePtr.snglInspiralTable->tau4,
            tablePtr.snglInspiralTable->tau5,
            tablePtr.snglInspiralTable->ttotal,
            tablePtr.snglInspiralTable->snr,
            tablePtr.snglInspiralTable->chisq,
            tablePtr.snglInspiralTable->chisq_dof,
            tablePtr.snglInspiralTable->sigmasq
            );
        tablePtr.snglInspiralTable = tablePtr.snglInspiralTable->next;
      }
      break;
    case summ_value_table:
      while( tablePtr.summValueTable )
      {
        FIRST_TABLE_ROW
        fprintf( xml->fp, SUMM_VALUE_ROW,
            tablePtr.summValueTable->program,
            tablePtr.summValueTable->start_time.gpsSeconds,
            tablePtr.summValueTable->start_time.gpsNanoSeconds,
            tablePtr.summValueTable->end_time.gpsSeconds,
            tablePtr.summValueTable->end_time.gpsNanoSeconds,
            tablePtr.summValueTable->ifo,
            tablePtr.summValueTable->name,
            tablePtr.summValueTable->value,
            tablePtr.summValueTable->comment
            );
        tablePtr.snglInspiralTable = tablePtr.snglInspiralTable->next;
      }
      break;
    default:
      ABORT( status, LIGOLWXMLH_EUTAB, LIGOLWXMLH_MSGEUTAB );
  }
  RETURN( status );
}

#undef FIRST_TABLE_ROW /* undefine first table row macro */
