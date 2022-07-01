/*
 * LIGOLwXMLInspiralRead.c
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 */

#include <stdio.h>
#include <metaio.h>

#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>

#include <LALExtTriggerTableFromLIGOLw.h>


/*
 *
 * LAL Functions
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


#define CLOBBER_EVENTS \
  while ( *eventHead ) \
{ \
  thisEvent = *eventHead; \
  *eventHead = (*eventHead)->next; \
  LALFree( thisEvent ); \
  thisEvent = NULL; \
}


int
LALExtTriggerTableFromLIGOLw (
    ExtTriggerTable   **eventHead,
    CHAR               *fileName,
    INT4                startEvent,
    INT4                stopEvent
    )

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
    {"ligo_fave_lho",          -1, 33},
    {"ligo_fave_llo",          -1, 34},
    {"ligo_delay",             -1, 35},
    {"event_number_gcn",       -1, 36},
    {"event_number_grb",       -1, 37},
    {"event_status",           -1, 38},
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
        /* REAL8 r8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_8; */
        INT4  i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;

        if ( tableDir[j].idx == 0 )
        {
          snprintf( thisEvent->det_alts, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 1 )
        {
          snprintf( thisEvent->det_band, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 2 )
        {
          snprintf( thisEvent->det_fluence, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 3 )
        {
          snprintf( thisEvent->det_fluence_int, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 4 )
        {
          snprintf( thisEvent->det_name, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 5 )
        {
          snprintf( thisEvent->det_peak, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 6 )
        {
          snprintf( thisEvent->det_peak_int, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 7 )
        {
          snprintf( thisEvent->det_snr, LIGOMETA_STD * sizeof(CHAR),
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
          snprintf( thisEvent->event_epoch, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 12 )
        {
          snprintf( thisEvent->event_err_type, LIGOMETA_STD * sizeof(CHAR),
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
          snprintf( thisEvent->event_type, LIGOMETA_STD * sizeof(CHAR),
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
          snprintf( thisEvent->notice_comments, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 21 )
        {
          snprintf( thisEvent->notice_id, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 22 )
        {
          snprintf( thisEvent->notice_sequence, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 23 )
        {
          thisEvent->notice_time = i4colData;
        }
        else if ( tableDir[j].idx == 24 )
        {
          snprintf( thisEvent->notice_type, LIGOMETA_STD * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 25 )
        {
          snprintf( thisEvent->notice_url, LIGOMETA_STD * sizeof(CHAR),
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
        else if ( tableDir[j].idx == 33 )
        {
          thisEvent->ligo_fave_lho = r4colData;
        }
        else if ( tableDir[j].idx == 34 )
        {
          thisEvent->ligo_fave_llo = r4colData;
        }
        else if ( tableDir[j].idx == 35 )
        {
          thisEvent->ligo_delay = r4colData;
        }
        else if ( tableDir[j].idx == 36 )
        {
          thisEvent->event_number_gcn= i4colData;
        }
        else if ( tableDir[j].idx == 37 )
        {
          snprintf( thisEvent->event_number_grb, 8 * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 38 )
        {
          thisEvent->event_status = i4colData;
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
