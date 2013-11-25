/*
*  Copyright (C) 2007 Duncan Brown, Kipp Cannon, Lisa M. Goggin
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
 * File Name: CreateMetaTableDir.c
 *
 * Author: Brown, D. A., and Brady, P. R.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown, D. A., and Brady, P. R.
 * \file
 * \ingroup lalmetaio
 * \brief Construct a \c MetaTableDirectory for a given LIGOLwXML table.
 *
 * ### Description ###
 *
 * The routine \c LALCreateMetaTableDir constructs a
 * \c MetaTableDirectory for a given LIGOLwXML table.  It determines the
 * location of each column expected to be present in the XML table and
 * populates the \c pos field with this information.  This then allows
 * other routines to parse the contents of an XML file and correctly interpret
 * the entries.  When reading these tables, a call is made to
 * \c LALCreateMetaTableDir.  For all other tables, the directory is
 * constructed internally by the reading code.
 *
 * ### Algorithm ###
 *
 * None.
 *
 * ### Uses ###
 *
 * Functions in the Metaio library:
 * <ul>
 * <li> \c MetaioFindColumn()</li>
 * <li> \c MetaioGetRow()</li>
 * <li> \c MetaioOpenTable()</li>
 * <li> \c MetaioClose()</li>
 * </ul>
 *
 * ### Notes ###
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <metaio.h>

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXMLRead.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

MetaTableDirectory* XLALCreateMetaTableDir(
    struct MetaioParseEnvironment *const env,
    MetadataTableType       table
    )

{
  MetaTableDirectory  *tableDir;
  INT4 i;

  switch( table )
  {
    case no_table:
      XLALPrintError( "XLALError - unable to index type no_table\n" );
      XLAL_ERROR_NULL( XLAL_EINVAL );
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
            XLAL_ERROR_NULL( XLAL_EFAILED );
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
            XLAL_ERROR_NULL( XLAL_EFAILED );
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
      XLAL_ERROR_NULL( XLAL_EINVAL );
      break;
    case search_summvars_table:
      XLALPrintError( "XLALError - "
          "unable to index type search_summvars_table\n" );
      XLAL_ERROR_NULL( XLAL_EINVAL );
      break;
    case sngl_inspiral_table:
      XLALPrintError( "XLALError - "
          "unable to index type sngl_inspiral_table\n" );
      XLAL_ERROR_NULL( XLAL_EINVAL );
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
          {"num_clust_trigs",         -1, 11},
          {"ds2_H1H2",                -1, 12},
          {"ds2_H1L1",                -1, 13},
          {"ds2_H1V1",                -1, 14},
          {"ds2_H2L1",                -1, 15},
          {"ds2_H2V1",                -1, 16},
          {"ds2_L1V1",                -1, 17},
          {"amplitude",               -1, 18},
          {"snr",                     -1, 19},
          {"eff_dist",                -1, 20},
          {"sigma_sq",                -1, 21},
          {"event_id",                -1, 22},
          {NULL,                       0, 0}
        };
        for ( i=0 ; tmpTableDir[i].name; ++i )
        {
          if ( (tmpTableDir[i].pos =
                MetaioFindColumn( env, tmpTableDir[i].name )) < 0 )
          {
            XLALPrintError( "XLALError - unable to find column %s\n",
                tmpTableDir[i].name );
            XLAL_ERROR_NULL( XLAL_EFAILED );
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
          {"event_id",                -1, 34},
          {"null_statistic",          -1, 35},
          {NULL,                       0, 0}
        };
        for ( i=0 ; tmpTableDir[i].name; ++i )
        {
          if ( (tmpTableDir[i].pos =
                MetaioFindColumn( env, tmpTableDir[i].name )) < 0 )
          {
            XLALPrintError( "XLALError - unable to find column %s\n",
                tmpTableDir[i].name );
            XLAL_ERROR_NULL( XLAL_EFAILED );
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
      XLAL_ERROR_NULL( XLAL_EINVAL );
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
          {"v_start_time",                 -1, 8},
          {"v_start_time_ns",              -1, 9},
          {"start_time_gmst",              -1, 10},
          {"longitude",                    -1, 11},
          {"latitude",                     -1, 12},
          {"distance",                     -1, 13},
          {"inclination",                  -1, 14},
          {"polarization",                 -1, 15},
          {"frequency",                    -1, 16},
          {"quality",                      -1, 17},
          {"phase",                        -1, 18},
          {"mass",                         -1, 19},
          {"spin",                         -1, 20},
          {"epsilon",                      -1, 21},
          {"amplitude",                    -1, 22},
          {"eff_dist_h",                   -1, 23},
          {"eff_dist_l",                   -1, 24},
          {"eff_dist_v",                   -1, 25},
          {"hrss",                         -1, 26},
          {"hrss_h",                       -1, 27},
          {"hrss_l",                       -1, 28},
          {"hrss_v",                       -1, 29},
          {NULL,                            0, 0}
        };
        for ( i=0 ; tmpTableDir[i].name; ++i )
        {
          if ( (tmpTableDir[i].pos =
                MetaioFindColumn( env, tmpTableDir[i].name )) < 0 )
          {
            XLALPrintError( "XLALError - unable to find column %s\n",
                tmpTableDir[i].name );
            XLAL_ERROR_NULL( XLAL_EFAILED );
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
      XLAL_ERROR_NULL( XLAL_EINVAL );
      break;
    default:
      XLALPrintError( "XLALError - "
          "unable to index table due to unknown table type error\n" );
      XLAL_ERROR_NULL( XLAL_EFAILED );
  }

  return tableDir;
}


void
LALCreateMetaTableDir(
    LALStatus              *status,
    MetaTableDirectory    **tableDir,
    struct MetaioParseEnvironment *const env UNUSED,
    MetadataTableType       table
    )

{
  INITSTATUS(status);
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
    case sngl_inspiral_table:
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
      break;
    case multi_inspiral_table:
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
      break;
    case sim_inspiral_table:
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
      break;
    case summ_value_table:
      ABORT( status, LIGOLWXMLREADH_EMTAB, LIGOLWXMLREADH_MSGEMTAB );
      break;
    default:
      ABORT( status, LIGOLWXMLREADH_EUTAB, LIGOLWXMLREADH_MSGEUTAB );
  }

  (void)tableDir;

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN( status );
}
