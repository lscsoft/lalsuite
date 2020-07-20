/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Lisa M. Goggin, Alexander Dietz, Kipp Cannon, Patrick Brady, Robert Adam Mercer, Stephen Fairhurst, Thomas Cokelaer
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
 *-----------------------------------------------------------------------
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
#include <lal/LIGOLwXMLRingdownRead.h>
#include <lal/StringInput.h>
#include <lal/XLALError.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * \author Brown, D. A. and Goggin, L. M.
 * \file
 *
 * \brief Routines to read the various ringdown search XML data into LAL structures.
 *
 * ### Description ###
 *
 *
 * ### Algorithm ###
 *
 * None.
 *
 * ### Uses ###
 *
 * Functions in the Metaio library:
 * <ul>
 * <li>MetaioFindColumn()</li>
 * <li>MetaioGetRow()</li>
 * <li>MetaioOpenTable()</li>
 * <li>MetaioClose()</li>
 * </ul>
 *
 * ### Notes ###
 *
 */


typedef struct
tagMetaTableDirectory
{
  const CHAR *name;
  INT4   pos;
  INT4   idx;
}
MetaTableDirectory;


static
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
    default:
      XLALPrintError( "XLALError - "
          "unable to index table due to unknown table type error\n" );
      XLAL_ERROR_NULL( XLAL_EFAILED );
  }

  return tableDir;
}


#define XLAL_CLOBBER_EVENTS \
  while ( eventHead ) \
{ \
  thisEvent = eventHead; \
  eventHead = (eventHead)->next; \
  LALFree( thisEvent ); \
  thisEvent = NULL; \
}


SimRingdownTable* XLALSimRingdownTableFromLIGOLw (
    CHAR               *fileName,
    INT4                startTime,
    INT4                stopTime
    )

{
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
    XLAL_ERROR_NULL( XLAL_EIO );
  }

  /* create table directory to find columns in file */
  tableDir = XLALCreateMetaTableDir(env, sim_ringdown_table);
  if ( ! tableDir )
  {
    XLALPrintError( "XLAL Error - "
        "unable to create sim_ringdown table directory\n" );
    XLAL_ERROR_NULL( XLAL_EIO );
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
      XLALPrintError( "XLAL Error - bad table directory for element %d\n", i );
      XLAL_CLOBBER_EVENTS;
      XLAL_ERROR_NULL( XLAL_EIO );
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
        XLAL_ERROR_NULL( XLAL_ENOMEM );
      }

      /* parse the contents of the row into the SimRingdownTable structure */
      for ( j = 0; tableDir[j].name; ++j )
      {
        if ( tableDir[j].pos < 0 )
        {
          XLALPrintError( "XLAL Error - bad table directory for element %d\n", j );
          XLAL_CLOBBER_EVENTS;
          XLAL_ERROR_NULL( XLAL_EIO );
        }

        /* dereference the data stored in the table */
        i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;
        r4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
        r8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_8;

        if ( tableDir[j].idx == 0 )
        {
          snprintf( thisEvent->waveform,
              LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s",
              env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 1 )
        {
          snprintf( thisEvent->coordinates,
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
          thisEvent->v_start_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 9 )
        {
          thisEvent->v_start_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 10 )
        {
          thisEvent->start_time_gmst = r8colData;
        }
        else if ( tableDir[j].idx == 11 )
        {
          thisEvent->longitude = r4colData;
        }
        else if ( tableDir[j].idx == 12 )
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
        else if ( tableDir[j].idx == 16 )
        {
          thisEvent->frequency = r4colData;
        }
        else if ( tableDir[j].idx == 17 )
        {
          thisEvent->quality = r4colData;
        }
        else if ( tableDir[j].idx == 18 )
        {
          thisEvent->phase = r4colData;
        }
        else if ( tableDir[j].idx == 19 )
        {
          thisEvent->mass = r4colData;
        }
        else if ( tableDir[j].idx == 20 )
        {
          thisEvent->spin = r4colData;
        }
        else if ( tableDir[j].idx == 21 )
        {
          thisEvent->epsilon= r4colData;
        }
        else if ( tableDir[j].idx == 22 )
        {
          thisEvent->amplitude = r4colData;
        }
        else if ( tableDir[j].idx == 23 )
        {
          thisEvent->eff_dist_h = r4colData;
        }
        else if ( tableDir[j].idx == 24 )
        {
          thisEvent->eff_dist_l = r4colData;
        }
        else if ( tableDir[j].idx == 25 )
        {
          thisEvent->eff_dist_v = r4colData;
        }
        else if ( tableDir[j].idx == 26 )
        {
          thisEvent->hrss = r4colData;
        }
        else if ( tableDir[j].idx == 27 )
        {
          thisEvent->hrss_h = r4colData;
        }
        else if ( tableDir[j].idx == 28 )
        {
          thisEvent->hrss_l = r4colData;
        }
        else if ( tableDir[j].idx == 29 )
        {
          thisEvent->hrss_v = r4colData;
        }
        else
        {
          XLALPrintError( "XLAL Error - "
              "table directory index %d out of bounds\n", j );
          XLAL_CLOBBER_EVENTS;
          XLAL_ERROR_NULL( XLAL_EIO );
        }
      }
    }
  }

  if ( mioStatus == -1 )
  {
    XLALPrintError( "XLAL Error - error parsing after row %d\n", i );
    XLAL_CLOBBER_EVENTS;
    MetaioClose( env );
    XLAL_ERROR_NULL( XLAL_EIO);
  }

  /* Normal exit */
  LALFree( tableDir );
  MetaioClose( env );

  return eventHead;
}

#undef XLAL_CLOBBER_EVENTS
