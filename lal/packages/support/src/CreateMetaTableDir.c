/*----------------------------------------------------------------------- 
 * 
 * File Name: CreateMetaTableDir.c
 *
 * Author: Brown, D. A., and Brady, P. R.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="CreateMetaTableDirCV">
Author: Brown, D. A., and Brady, P. R.
$Id$
</lalVerbatim> 
#endif

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXMLRead.h>

NRCSID( CREATEMETATABLEDIRC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{CreateMetaTableDir.c}}

Routines to create a directory of a LIGO LW XML file.

\subsubsection*{Prototypes}
\input{CreateMetaTableDirCP}
\idx{XLALCreateMetaTableDir()}
\idx{LALCreateMetaTableDir()}

\subsubsection*{Description}

The routine \verb+LALCreateMetaTableDir+ constructs a \verb+MetaTableDirectory+
for a given LIGOLwXML table.  It determines the location of each column expected
to be present in the XML table and populates the \verb+pos+ field with this
information.  This then allows other routines to parse the contents of an XML
file and correctly interpret the entries.  Currently, this has been implemented
only for the \verb+sngl_burst+ and \verb+sim_burst+ tables.  When reading these
tables, a call is made to \verb+LALCreateMetaTableDir+.  For all other tables,
the directory is constructed internally by the reading code.

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

\vfill{\footnotesize\input{CreateMetaTableDirCV}}

</lalLaTeX>
#endif

/* <lalVerbatim file="CreateMetaTableDirCP"> */
MetaTableDirectory* XLALCreateMetaTableDir(
    const MetaioParseEnv    env,
    MetadataTableType       table
    )
/* </lalVerbatim> */
{
  static const char   *func = "XLALCreateMetaTableDir";
  MetaTableDirectory  *tableDir;
  INT4 i;

  switch( table )
  {
    case no_table:
      XLALPrintError( "XLALError - unable to index type no_table\n" );
      XLAL_ERROR_NULL( func, XLAL_EINVAL );
      break;
    case process_table:
      {
        MetaTableDirectory tmpTableDir[] =
        {
          {"program",                  -1, 0},
          {"version",                  -1, 1},
          {"cvs_repository",           -1, 2},
          {"cvs_entry_time",           -1, 3},
          {"comment",                  -1, 4},
          {"is_online",                -1, 5},
          {"node",                     -1, 6},
          {"username",                 -1, 7},
          {"start_time",               -1, 8},
          {"end_time",                 -1, 9},
          {"jobid",                    -1, 10},
          {"domain",                   -1, 11},
          {"ifos",                     -1, 12},
          {NULL,                        0, 0}
        };
        for ( i=0 ; tmpTableDir[i].name; ++i )
        {
          if ( (tmpTableDir[i].pos = 
                MetaioFindColumn( env, tmpTableDir[i].name )) < 0 )
          {
            XLALPrintError( "XLALError - unable to find column %s\n", 
                tmpTableDir[i].name );
            XLAL_ERROR_NULL( func, XLAL_EFAILED );
          }
        }

        tableDir = (MetaTableDirectory *) LALMalloc( (i+1) * 
            sizeof(MetaTableDirectory)) ;
        memcpy(tableDir, tmpTableDir, (i+1)*sizeof(MetaTableDirectory) );
      }
      break;
    case process_params_table:
      {
        MetaTableDirectory tmpTableDir[] =
        {
          {"program",                  -1, 0},
          {"param",                    -1, 1},
          {"type",                     -1, 2},
          {"value",                    -1, 3},
          {NULL,                        0, 0}
        };
        for ( i=0 ; tmpTableDir[i].name; ++i )
        {
          if ( (tmpTableDir[i].pos = 
                MetaioFindColumn( env, tmpTableDir[i].name )) < 0 )
          {
            XLALPrintError( "XLALError - unable to find column %s\n", 
                tmpTableDir[i].name );
            XLAL_ERROR_NULL( func, XLAL_EFAILED );
          }
        }

        tableDir = (MetaTableDirectory *) LALMalloc( (i+1) * 
            sizeof(MetaTableDirectory)) ;
        memcpy(tableDir, tmpTableDir, (i+1)*sizeof(MetaTableDirectory) );
      }
      break;
    case search_summary_table:
      XLALPrintError( "XLALError - "
          "unable to index type search_summary_table\n" );
      XLAL_ERROR_NULL( func, XLAL_EINVAL );
      break;
    case search_summvars_table:
      XLALPrintError( "XLALError - "
          "unable to index type search_summvars_table\n" );
      XLAL_ERROR_NULL( func, XLAL_EINVAL );
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
          {"tfvolume",                -1, 13},
          {"event_id",                -1, 14},
          {NULL,                       0, 0}
        };
        for ( i=0 ; tmpTableDir[i].name; ++i )
        {
          if ( (tmpTableDir[i].pos = 
                MetaioFindColumn( env, tmpTableDir[i].name )) < 0 )
          {
            XLALPrintError( "XLALError - unable to find column %s\n", 
                tmpTableDir[i].name );
            XLAL_ERROR_NULL( func, XLAL_EFAILED );
          }
        }

        tableDir = (MetaTableDirectory *) LALMalloc( (i+1) * 
            sizeof(MetaTableDirectory)) ;
        memcpy(tableDir, tmpTableDir, (i+1)*sizeof(MetaTableDirectory) );
      }
      break;
    case sngl_inspiral_table:
      XLALPrintError( "XLALError - "
          "unable to index type sngl_inspiral_table\n" );
      XLAL_ERROR_NULL( func, XLAL_EINVAL );
      break;
    case sngl_ringdown_table:
      {
        MetaTableDirectory tmpTableDir[] =
        {
          {"ifo",                     -1, 0},
          {"channel",                 -1, 1},
          {"start_time",              -1, 2},
          {"start_time_ns",           -1, 3},
          {"start_time_gmst",         -1, 4},
          {"frequency",               -1, 5},
          {"quality",                 -1, 6},
          {"phase",                   -1, 7},
          {"mass",                    -1, 8},
          {"spin",                    -1, 9},
          {"epsilon",                 -1, 10},
          {"amplitude",               -1, 11},
          {"snr",                     -1, 12},
          {"eff_dist",                -1, 13},
          {"sigma_sq",                -1, 14},
          {"event_id",                -1, 15},
          {NULL,                       0, 0}
        };
        for ( i=0 ; tmpTableDir[i].name; ++i )
        {
          if ( (tmpTableDir[i].pos = 
                MetaioFindColumn( env, tmpTableDir[i].name )) < 0 )
          {
            XLALPrintError( "XLALError - unable to find column %s\n", 
                tmpTableDir[i].name );
            XLAL_ERROR_NULL( func, XLAL_EFAILED );
          }
        }

        tableDir = (MetaTableDirectory *) LALMalloc( (i+1) * 
            sizeof(MetaTableDirectory)) ;
        memcpy(tableDir, tmpTableDir, (i+1)*sizeof(MetaTableDirectory) );
      }
      break;
    case multi_inspiral_table:
      {
        MetaTableDirectory tmpTableDir[] =
        {
          {"ifos",                    -1, 0},
          {"search",                  -1, 1},
          {"end_time",                -1, 2},
          {"end_time_ns",             -1, 3},
          {"end_time_gmst",           -1, 4},
          {"impulse_time",            -1, 5},
          {"impulse_time_ns",         -1, 6},
          {"amplitude",               -1, 7},
          {"ifo1_eff_distance",       -1, 8},
          {"ifo2_eff_distance",       -1, 9},
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
          {"ifo1_snr",                -1, 22},
          {"ifo2_snr",                -1, 23},
          {"snr",                     -1, 24},
          {"chisq",                   -1, 25},
          {"chisq_dof",               -1, 26},
          {"sigmasq",                 -1, 27},
          {"ligo_axis_ra",            -1, 28},
          {"ligo_axis_dec",           -1, 29},
          {"ligo_angle",              -1, 30},
          {"ligo_angle_sig",          -1, 31},
          {"inclination",             -1, 32},
          {"polarization",            -1, 33},
          {NULL,                       0, 0}
        };
        for ( i=0 ; tmpTableDir[i].name; ++i )
        {
          if ( (tmpTableDir[i].pos = 
                MetaioFindColumn( env, tmpTableDir[i].name )) < 0 )
          {
            XLALPrintError( "XLALError - unable to find column %s\n", 
                tmpTableDir[i].name );
            XLAL_ERROR_NULL( func, XLAL_EFAILED );
          }
        }

        tableDir = (MetaTableDirectory *) LALMalloc( (i+1) * 
            sizeof(MetaTableDirectory)) ;
        memcpy(tableDir, tmpTableDir, (i+1)*sizeof(MetaTableDirectory) );
      }
      break;
    case sim_inspiral_table:
      XLALPrintError( "XLALError - "
          "unable to index type sim_inspiral_table\n" );
      XLAL_ERROR_NULL( func, XLAL_EINVAL );
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
            XLALPrintError( "XLALError - unable to find column %s\n", 
                tmpTableDir[i].name );
            XLAL_ERROR_NULL( func, XLAL_EFAILED );
          }
        }

        tableDir = (MetaTableDirectory *) LALMalloc( (i+1) * 
            sizeof(MetaTableDirectory)) ;
        memcpy(tableDir, tmpTableDir, (i+1)*sizeof(MetaTableDirectory) );
      }
      break;
    case sim_ringdown_table:
      {
        MetaTableDirectory tmpTableDir[] =
        {
          {"waveform",                     -1, 0},
          {"coordinates",                  -1, 1},
          {"geocent_start_time",           -1, 2},
          {"geocent_start_time_ns",        -1, 3},
          {"h_start_time",                 -1, 4},
          {"h_start_time_ns",              -1, 5},
          {"l_start_time",                 -1, 6},
          {"l_start_time_ns",              -1, 7},
          {"start_time_gmst",              -1, 8},
          {"longitude",                    -1, 9},
          {"latitude",                     -1, 10},
          {"distance",                     -1, 11},
          {"inclination",                  -1, 12},
          {"polarization",                 -1, 13},
          {"frequency",                    -1, 14},
          {"quality",                      -1, 15},
          {"phase",                        -1, 16},
          {"mass",                         -1, 17},
          {"spin",                         -1, 18},
          {"epsilon",                      -1, 19},
          {"amplitude",                    -1, 20},
          {"eff_dist_h",                   -1, 21},
          {"eff_dist_l",                   -1, 22},
          {"hrss",                         -1, 23},
          {"hrss_h",                       -1, 24},
          {"hrss_l",                       -1, 25},
          {NULL,                            0, 0}
        };
        for ( i=0 ; tmpTableDir[i].name; ++i )
        {
          if ( (tmpTableDir[i].pos =
                MetaioFindColumn( env, tmpTableDir[i].name )) < 0 )
          {
            XLALPrintError( "XLALError - unable to find column %s\n", 
                tmpTableDir[i].name );
            XLAL_ERROR_NULL( func, XLAL_EFAILED );
          }
        }

        tableDir = (MetaTableDirectory *) LALMalloc( (i+1) *
            sizeof(MetaTableDirectory)) ;
        memcpy(tableDir, tmpTableDir, (i+1)*sizeof(MetaTableDirectory) );
      }
      break;
    case summ_value_table:
      XLALPrintError( "XLALError - "
          "unable to index type summ_value_table\n" );
      XLAL_ERROR_NULL( func, XLAL_EINVAL );
      break;
    default: 
      XLALPrintError( "XLALError - "
          "unable to index table due to unknown table type error\n" );
      XLAL_ERROR_NULL( func, XLAL_EFAILED );
  }

  return tableDir;
}

/* <lalVerbatim file="CreateMetaTableDirCP"> */
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

  INITSTATUS( status, "LALCreateMetaTableDir", CREATEMETATABLEDIRC );
  ATTATCHSTATUSPTR (status);

  /* check the inputs */
  ASSERT( !*tableDir, status, LIGOLWXMLREADH_ENNUL, LIGOLWXMLREADH_MSGENNUL );

  switch( table )
  {
    case no_table:
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
      break;
    case process_table:
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
      break;
    case process_params_table:
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
      break;
    case search_summary_table:
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
      break;
    case search_summvars_table:
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
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
          {"tfvolume",                -1, 13},
          {"event_id",                -1, 14},
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
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
      break;
    case multi_inspiral_table:
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
      break;
    case sim_inspiral_table:
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
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
          {"distance",                     -1, 16},
          {"freq",                         -1, 17},
          {"tau",                          -1, 18},
          {"zm_number",                    -1, 19},
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
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
      break;
    default:
      ABORT( status, LIGOLWXMLREADH_EUTAB, LIGOLWXMLREADH_MSGEUTAB );
  }

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN( status );
}
