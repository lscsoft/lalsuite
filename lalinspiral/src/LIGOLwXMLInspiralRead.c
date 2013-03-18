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
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXMLInspiralRead.h>


#define XLAL_CLOBBER_EVENTS \
  while ( eventHead ) \
{ \
  thisEvent = eventHead; \
  eventHead = (eventHead)->next; \
  LALFree( thisEvent ); \
  thisEvent = NULL; \
}


MultiInspiralTable    * XLALMultiInspiralTableFromLIGOLw (
    CHAR               *fileName
    )

{
  int                                   i, j, nrows;
  int                                   mioStatus=0;
  MultiInspiralTable                   *thisEvent = NULL;
  MultiInspiralTable                   *eventHead = NULL;
  struct MetaioParseEnvironment         parseEnv;
  const  MetaioParseEnv                 env = &parseEnv;
  MetaTableDirectory                   *tableDir = NULL;

  /* open the multi_inspiral XML file */
  mioStatus = MetaioOpenTable( env, fileName, "multi_inspiral" );
  if ( mioStatus )
  {
    XLAL_ERROR_NULL( XLAL_EIO );
  }

  /* create table directory to find columns in file*/
  tableDir = XLALCreateMetaTableDir(env, multi_inspiral_table);

  /* loop over the rows in the file */
  i = nrows = 0;
  while ( (mioStatus = MetaioGetRow(env)) == 1 )
  {
    /* count the rows in the file */
    i++;
    /* allocate memory for the template we are about to read in */
    if ( ! eventHead )
    {
      thisEvent = eventHead = (MultiInspiralTable *)
        LALCalloc( 1, sizeof(MultiInspiralTable) );
    }
    else
    {
      thisEvent = thisEvent->next = (MultiInspiralTable *)
        LALCalloc( 1, sizeof(MultiInspiralTable) );
    }
    if ( ! thisEvent )
    {
      fprintf( stderr, "could not allocate multi inspiral event\n" );
      XLAL_CLOBBER_EVENTS;
      MetaioClose( env );
      XLAL_ERROR_NULL( XLAL_ENOMEM );
    }


    /* parse the contents of the row into the InspiralTemplate structure */
    for ( j = 0; tableDir[j].name; ++j )
    {
      enum METAIO_Type column_type = env->ligo_lw.table.col[tableDir[j].pos].data_type;
      REAL4 r4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
      REAL8 r8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_8;
      INT4  i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;

      if ( tableDir[j].idx == 0 )
      {
        snprintf( thisEvent->ifos, LIGOMETA_IFOS_MAX * sizeof(CHAR),
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 1 )
      {
        snprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 2 )
      {
        thisEvent->end_time.gpsSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 3 )
      {
        thisEvent->end_time.gpsNanoSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 4 )
      {
        thisEvent->end_time_gmst = r8colData;
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
        thisEvent->amplitude = r4colData;
      }
      else if ( tableDir[j].idx == 8 )
      {
        thisEvent->distance = r4colData;
      }
      else if ( tableDir[j].idx == 9 )
      {
        thisEvent->eff_dist_h1 = r4colData;
      }
      else if ( tableDir[j].idx == 10 )
      {
        thisEvent->eff_dist_h2 = r4colData;
      }
      else if ( tableDir[j].idx == 11 )
      {
        thisEvent->eff_dist_l = r4colData;
      }
      else if ( tableDir[j].idx == 12 )
      {
        thisEvent->eff_dist_g = r4colData;
      }
      else if ( tableDir[j].idx == 13 )
      {
        thisEvent->eff_dist_t = r4colData;
      }
      else if ( tableDir[j].idx == 14 )
      {
        thisEvent->eff_dist_v = r4colData;
      }
      else if ( tableDir[j].idx == 15 )
      {
        thisEvent->eff_dist_h1h2 = r4colData;
      }
      else if ( tableDir[j].idx == 16 )
      {
        thisEvent->chi = r4colData;
      }
      else if ( tableDir[j].idx == 17 )
      {
        thisEvent->kappa = r4colData;
      }
      else if ( tableDir[j].idx == 18 )
      {
        thisEvent->coa_phase = r4colData;
      }
      else if ( tableDir[j].idx == 19 )
      {
        thisEvent->mass1 = r4colData;
      }
      else if ( tableDir[j].idx == 20 )
      {
        thisEvent->mass2 = r4colData;
      }
      else if ( tableDir[j].idx == 21 )
      {
        thisEvent->mchirp = r4colData;
      }
      else if ( tableDir[j].idx == 22 )
      {
        thisEvent->eta = r4colData;
      }
      else if ( tableDir[j].idx == 23 )
      {
        thisEvent->tau0 = r4colData;
      }
      else if ( tableDir[j].idx == 24 )
      {
        thisEvent->tau2 = r4colData;
      }
      else if ( tableDir[j].idx == 25 )
      {
        thisEvent->tau3 = r4colData;
      }
      else if ( tableDir[j].idx == 26 )
      {
        thisEvent->tau4 = r4colData;
      }
      else if ( tableDir[j].idx == 27 )
      {
        thisEvent->tau5 = r4colData;
      }
      else if ( tableDir[j].idx == 28 )
      {
        thisEvent->ttotal = r4colData;
      }
      else if ( tableDir[j].idx == 29 )
      {
        thisEvent->snr = r4colData;
      }
      else if ( tableDir[j].idx == 30 )
      {
        thisEvent->snr_dof = i4colData;
      }
      else if ( tableDir[j].idx == 31 )
      {
        thisEvent->chisq = r4colData;
      }
      else if ( tableDir[j].idx == 32 )
      {
        thisEvent->chisq_dof = i4colData;
      }
      else if ( tableDir[j].idx == 33 )
      {
        thisEvent->bank_chisq = r4colData;
      }
      else if ( tableDir[j].idx == 34 )
      {
        thisEvent->bank_chisq_dof = i4colData;
      }
      else if ( tableDir[j].idx == 35 )
      {
        thisEvent->cont_chisq = r4colData;
      }
      else if ( tableDir[j].idx == 36 )
      {
        thisEvent->cont_chisq_dof = i4colData;
      }
      else if ( tableDir[j].idx == 37 )
      {
        thisEvent->sigmasq_h1 = r8colData;
      }
      else if ( tableDir[j].idx == 38 )
      {
        thisEvent->sigmasq_h2 = r8colData;
      }
      else if ( tableDir[j].idx == 39 )
      {
        thisEvent->sigmasq_l = r8colData;
      }
      else if ( tableDir[j].idx == 40 )
      {
        thisEvent->sigmasq_g = r8colData;
      }
      else if ( tableDir[j].idx == 41 )
      {
        thisEvent->sigmasq_t = r8colData;
      }
      else if ( tableDir[j].idx == 42 )
      {
        thisEvent->sigmasq_v = r8colData;
      }
      else if ( tableDir[j].idx == 43 )
      {
        thisEvent->chisq_h1 = r4colData;
      }
      else if ( tableDir[j].idx == 44 )
      {
        thisEvent->chisq_h2 = r4colData;
      }
      else if ( tableDir[j].idx == 45 )
      {
        thisEvent->chisq_l = r4colData;
      }
      else if ( tableDir[j].idx == 46 )
      {
        thisEvent->chisq_g = r4colData;
      }
      else if ( tableDir[j].idx == 47 )
      {
        thisEvent->chisq_t = r4colData;
      }
      else if ( tableDir[j].idx == 48 )
      {
        thisEvent->chisq_v = r4colData;
      }
      else if ( tableDir[j].idx == 49 )
      {
        thisEvent->ra = r4colData;
      }
      else if ( tableDir[j].idx == 50 )
      {
        thisEvent->dec = r4colData;
      }
      else if ( tableDir[j].idx == 51 )
      {
        thisEvent->ligo_angle = r4colData;
      }
      else if ( tableDir[j].idx == 52 )
      {
        thisEvent->ligo_angle_sig = r4colData;
      }
      else if ( tableDir[j].idx == 53 )
      {
        thisEvent->inclination = r4colData;
      }
      else if ( tableDir[j].idx == 54 )
      {
        thisEvent->polarization = r4colData;
      }
      else if ( tableDir[j].idx == 55 )
      {
        thisEvent->null_statistic = r4colData;
      }
      else if ( tableDir[j].idx == 56 )
      {
        thisEvent->null_stat_h1h2 = r4colData;
      }
      else if ( tableDir[j].idx == 57 )
      {
        thisEvent->null_stat_degen = r4colData;
      }
      else if ( tableDir[j].idx == 58 )
      {
        if ( tableDir[j].pos > 0 )
        {
          INT8 i8colData;
          if ( column_type == METAIO_TYPE_INT_8S )
            i8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_8s;
          else
          {
            i8colData = XLALLIGOLwParseIlwdChar(env, tableDir[j].pos, "multi_inspiral", "event_id");
            if ( i8colData < 0 )
              XLAL_ERROR_NULL( XLAL_EFUNC );
          }
          if ( i8colData )
          {
            thisEvent->event_id = LALCalloc( 1, sizeof(*thisEvent->event_id) );
            thisEvent->event_id->id = i8colData;
            thisEvent->event_id->multiInspiralTable = thisEvent;
          }
        }
      }
      else if ( tableDir[j].idx == 59 )
      {
        thisEvent->h1quad.re = r4colData;
      }
      else if ( tableDir[j].idx == 60 )
      {
        thisEvent->h1quad.im = r4colData;
      }
      else if ( tableDir[j].idx == 61 )
      {
        thisEvent->h2quad.re = r4colData;
      }
      else if ( tableDir[j].idx == 62 )
      {
        thisEvent->h2quad.im = r4colData;
      }
      else if ( tableDir[j].idx == 63 )
      {
        thisEvent->l1quad.re = r4colData;
      }
      else if ( tableDir[j].idx == 64 )
      {
        thisEvent->l1quad.im = r4colData;
      }
      else if ( tableDir[j].idx == 65 )
      {
        thisEvent->g1quad.re = r4colData;
      }
      else if ( tableDir[j].idx == 66 )
      {
        thisEvent->g1quad.im = r4colData;
      }
      else if ( tableDir[j].idx == 67 )
      {
        thisEvent->t1quad.re = r4colData;
      }
      else if ( tableDir[j].idx == 68 )
      {
        thisEvent->t1quad.im = r4colData;
      }
      else if ( tableDir[j].idx == 69 )
      {
        thisEvent->v1quad.re = r4colData;
      }
      else if ( tableDir[j].idx == 70 )
      {
        thisEvent->v1quad.im = r4colData;
      }
      else if ( tableDir[j].idx == 71 )
      {
        thisEvent->coh_snr_h1h2 = r4colData;
      }
      else if ( tableDir[j].idx == 72 )
      {
        thisEvent->cohSnrSqLocal = r4colData;
      }
      else if ( tableDir[j].idx == 73 )
      {
        thisEvent->autoCorrCohSq = r4colData;
      }
      else if ( tableDir[j].idx == 74 )
      {
        thisEvent->crossCorrCohSq = r4colData;
      }
      else if ( tableDir[j].idx == 75 )
      {
        thisEvent->autoCorrNullSq = r4colData;
      }
      else if ( tableDir[j].idx == 76 )
      {
        thisEvent->crossCorrNullSq = r4colData;
      }
      else if ( tableDir[j].idx == 77 )
      {
        thisEvent->ampMetricEigenVal1 = r8colData;
      }
      else if ( tableDir[j].idx == 78 )
      {
        thisEvent->ampMetricEigenVal2 = r8colData;
      }
      else
      {
        XLAL_CLOBBER_EVENTS;
        XLAL_ERROR_NULL( XLAL_EIO);
      }
    }
    /* count the number of triggers parsed */
    nrows++;
  }

  if ( mioStatus == -1 )
  {
    fprintf( stderr, "error parsing after row %d\n", i );
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


/*
 *
 * LAL Functions
 *
 */


#define CLOBBER_EVENTS \
  while ( *eventHead ) \
{ \
  thisEvent = *eventHead; \
  *eventHead = (*eventHead)->next; \
  LALFree( thisEvent ); \
  thisEvent = NULL; \
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
    {"bank_chisq",                   -1, 38},
    {"bank_chisq_dof",               -1, 39},
    {"cont_chisq",                   -1, 40},
    {"cont_chisq_dof",               -1, 41},
    {"sigmasq",                 -1, 42},
    {"rsqveto_duration",        -1, 43},
    {"event_id",                -1, 44},
    {"Gamma0",                  -1, 45},
    {"Gamma1",                  -1, 46},
    {"Gamma2",                  -1, 47},
    {"Gamma3",                  -1, 48},
    {"Gamma4",                  -1, 49},
    {"Gamma5",                  -1, 50},
    {"Gamma6",                  -1, 51},
    {"Gamma7",                  -1, 52},
    {"Gamma8",                  -1, 53},
    {"Gamma9",                  -1, 54},
    {"kappa",                   -1, 55},
    {"chi",                     -1, 56},
    {"spin1x",                  -1, 57},
    {"spin1y",                  -1, 58},
    {"spin1z",                  -1, 59},
    {"spin2x",                  -1, 60},
    {"spin2y",                  -1, 61},
    {"spin2z",                  -1, 62},
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
  mioStatus = MetaioOpenFile( env, fileName );
  if ( mioStatus )
  {
    MetaioAbort(env);
    MetaioClose(env);
    fprintf( stderr, "unable to open file %s\n", fileName );
    return -1;
  }

  /* open the sngl_inspiral table template bank file */
  mioStatus = MetaioOpenTableOnly( env, "sngl_inspiral" );
  if ( mioStatus )
  {
    MetaioAbort(env);
    MetaioClose(env);
    fprintf( stdout, "no sngl_inspiral table in file %s\n", fileName );
    return 0;
  }

  /* figure out the column positions of the template parameters */
  for ( i = 0; tableDir[i].name; ++i )
  {
    if ( (tableDir[i].pos = MetaioFindColumn( env, tableDir[i].name )) < 0
        &&  tableDir[i].idx != 39 )
    {
      fprintf( stderr, "unable to find column %s\n", tableDir[i].name );

      if ( ! strcmp(tableDir[i].name, "event_id") )
      {
        fprintf( stderr,
            "The event_id column is not populated, continuing anyway\n");
      }
      else if ( strstr(tableDir[i].name, "Gamma") )
      {
        fprintf( stderr,
            "The %s column is not populated, continuing anyway\n", tableDir[i].name);
      }
      else
      {
        MetaioClose(env);
        return -1;
      }
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
        enum METAIO_Type column_type = env->ligo_lw.table.col[tableDir[j].pos].data_type;
        REAL4 r4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
        REAL8 r8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_8;
        INT4  i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;

        if ( tableDir[j].pos < 0 ) continue;

        if ( tableDir[j].idx == 0 )
        {
          snprintf( thisEvent->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 1 )
        {
          snprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 2 )
        {
          snprintf( thisEvent->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
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
          thisEvent->bank_chisq = r4colData;
        }
        else if ( tableDir[j].idx == 39 )
        {
          thisEvent->bank_chisq_dof = i4colData;
        }
        else if ( tableDir[j].idx == 40 )
        {
          thisEvent->cont_chisq = r4colData;
        }
        else if ( tableDir[j].idx == 41 )
        {
          thisEvent->cont_chisq_dof = i4colData;
        }
        else if ( tableDir[j].idx == 42 )
        {
          thisEvent->sigmasq = r8colData;
        }
        else if ( tableDir[j].idx == 43 )
        {
          thisEvent->rsqveto_duration = r4colData;
        }
        else if ( tableDir[j].idx == 44 )
        {
          if ( tableDir[j].pos > 0 )
          {
            INT8 i8colData;
            if ( column_type == METAIO_TYPE_INT_8S )
              i8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_8s;
            else
            {
              i8colData = XLALLIGOLwParseIlwdChar(env, tableDir[j].pos, "sngl_inspiral", "event_id");
              if ( i8colData < 0 )
                return -1;
            }
            if ( i8colData )
            {
              thisEvent->event_id = LALCalloc( 1, sizeof(*thisEvent->event_id) );
              thisEvent->event_id->id = i8colData;
              thisEvent->event_id->snglInspiralTable = thisEvent;
            }
          }
        }
        else if ( tableDir[j].idx == 45 )
        {
          thisEvent->Gamma[0] = r4colData;
        }
        else if ( tableDir[j].idx == 46 )
        {
          thisEvent->Gamma[1] = r4colData;
        }
        else if ( tableDir[j].idx == 47 )
        {
          thisEvent->Gamma[2] = r4colData;
        }
        else if ( tableDir[j].idx == 48 )
        {
          thisEvent->Gamma[3] = r4colData;
        }
        else if ( tableDir[j].idx == 49 )
        {
          thisEvent->Gamma[4] = r4colData;
        }
        else if ( tableDir[j].idx == 50 )
        {
          thisEvent->Gamma[5] = r4colData;
        }
        else if ( tableDir[j].idx == 51 )
        {
          thisEvent->Gamma[6] = r4colData;
        }
        else if ( tableDir[j].idx == 52 )
        {
          thisEvent->Gamma[7] = r4colData;
        }
        else if ( tableDir[j].idx == 53 )
        {
          thisEvent->Gamma[8] = r4colData;
        }
        else if ( tableDir[j].idx == 54 )
        {
          thisEvent->Gamma[9] = r4colData;
        }
        else if ( tableDir[j].idx == 55 )
        {
          thisEvent->kappa = r4colData;
        }
        else if ( tableDir[j].idx == 56 )
        {
          thisEvent->chi = r4colData;
        }
        else if ( tableDir[j].idx == 57 )
        {
          thisEvent->spin1x = r4colData;
        }
        else if ( tableDir[j].idx == 58 )
        {
          thisEvent->spin1y = r4colData;
        }
        else if ( tableDir[j].idx == 59 )
        {
          thisEvent->spin1z = r4colData;
        }
        else if ( tableDir[j].idx == 60 )
        {
          thisEvent->spin2x = r4colData;
        }
        else if ( tableDir[j].idx == 61 )
        {
          thisEvent->spin2y = r4colData;
        }
        else if ( tableDir[j].idx == 62 )
        {
          thisEvent->spin2z = r4colData;
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
  while ( *bankHead ) \
{ \
  thisTmplt = *bankHead; \
  *bankHead = (*bankHead)->next; \
  LALFree( thisTmplt ); \
  thisTmplt = NULL; \
}


int
InspiralTmpltBankFromLIGOLw (
    InspiralTemplate  **bankHead,
    const CHAR         *fileName,
    INT4                startTmplt,
    INT4                stopTmplt
    )

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
    {"end_time", -1, 14},
    {"end_time_ns", -1, 15},
    {"event_id", -1, 16},
    {"ifo", -1, 17},
    {"Gamma0", -1, 18},
    {"Gamma1", -1, 19},
    {"Gamma2", -1, 20},
    {"Gamma3", -1, 21},
    {"Gamma4", -1, 22},
    {"Gamma5", -1, 23},
    {"Gamma6", -1, 24},
    {"Gamma7", -1, 25},
    {"Gamma8", -1, 26},
    {"Gamma9", -1, 27},
    {"kappa", -1, 28},
    {"chi", -1, 29},
    {"spin1x", -1, 30},
    {"spin1y", -1, 31},
    {"spin1z", -1, 32},
    {"spin2x", -1, 33},
    {"spin2y", -1, 34},
    {"spin2z", -1, 35},
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

      if ( ! strcmp(tableDir[i].name, "event_id") )
      {
        fprintf( stderr,
            "The event_id column is not populated, continuing anyway\n");
      }
      else if ( strstr(tableDir[i].name, "Gamma") )
      {
        fprintf( stderr,
            "The %s column is not populated, continuing anyway\n", tableDir[i].name);
      }
      else
      {
        MetaioClose(env);
        return -1;
      }
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
        enum METAIO_Type column_type = env->ligo_lw.table.col[tableDir[j].pos].data_type;
        REAL4 colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
        INT4 i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;

        if ( tableDir[j].pos < 0 ) continue;
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
        else if ( tableDir[j].idx == 14 )
        {
          thisTmplt->end_time.gpsSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 15 )
        {
          thisTmplt->end_time.gpsNanoSeconds = i4colData;
        }
        else if ( tableDir[j].idx == 16 )
        {
          if ( tableDir[j].pos > 0 )
          {
            INT8 i8colData;
            if ( column_type == METAIO_TYPE_INT_8S )
              i8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_8s;
            else
            {
              i8colData = XLALLIGOLwParseIlwdChar(env, tableDir[j].pos, "sngl_inspiral", "event_id");
              if ( i8colData < 0 )
                return -1;
            }
            if ( i8colData )
            {
              thisTmplt->event_id = LALCalloc( 1, sizeof(*thisTmplt->event_id) );
              thisTmplt->event_id->id = i8colData;
              thisTmplt->event_id->inspiralTemplate = thisTmplt;
            }
          }
        }
        else if ( tableDir[j].idx == 17 )
        {
          snprintf( thisTmplt->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
        }
        else if ( tableDir[j].idx == 18 )
        {
          thisTmplt->Gamma[0] = colData;
        }
        else if ( tableDir[j].idx == 19 )
        {
          thisTmplt->Gamma[1] = colData;
        }
        else if ( tableDir[j].idx == 20 )
        {
          thisTmplt->Gamma[2] = colData;
        }
        else if ( tableDir[j].idx == 21 )
        {
          thisTmplt->Gamma[3] = colData;
        }
        else if ( tableDir[j].idx == 22 )
        {
          thisTmplt->Gamma[4] = colData;
        }
        else if ( tableDir[j].idx == 23 )
        {
          thisTmplt->Gamma[5] = colData;
        }
        else if ( tableDir[j].idx == 24 )
        {
          thisTmplt->Gamma[6] = colData;
        }
        else if ( tableDir[j].idx == 25 )
        {
          thisTmplt->Gamma[7] = colData;
        }
        else if ( tableDir[j].idx == 26 )
        {
          thisTmplt->Gamma[8] = colData;
        }
        else if ( tableDir[j].idx == 27 )
        {
          thisTmplt->Gamma[9] = colData;
        }
        else if ( tableDir[j].idx == 28 )
        {
          thisTmplt->kappa = colData;
        }
        else if ( tableDir[j].idx == 29 )
        {
          thisTmplt->chi = colData;
        }
        else if ( tableDir[j].idx == 30 )
        {
          thisTmplt->spin1[0] = colData;
        }
        else if ( tableDir[j].idx == 31 )
        {
          thisTmplt->spin1[1] = colData;
        }
        else if ( tableDir[j].idx == 32 )
        {
          thisTmplt->spin1[2] = colData;
        }
        else if ( tableDir[j].idx == 33 )
        {
          thisTmplt->spin2[0] = colData;
        }
        else if ( tableDir[j].idx == 34 )
        {
          thisTmplt->spin2[1] = colData;
        }
        else if ( tableDir[j].idx == 35 )
        {
          thisTmplt->spin2[2] = colData;
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
  while ( *simHead ) \
{ \
  thisSim = *simHead; \
  *simHead = (*simHead)->next; \
  LALFree( thisSim ); \
  thisSim = NULL; \
}


int
SimInspiralTableFromLIGOLw (
    SimInspiralTable   **simHead,
    const CHAR          *fileName,
    INT4                 startTime,
    INT4                 endTime
    )

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
    {"f_lower",             -1, 42},
    {"f_final",             -1, 43},
    {"mchirp",              -1, 44},
    {"eff_dist_h",          -1, 45},
    {"eff_dist_l",          -1, 46},
    {"eff_dist_g",          -1, 47},
    {"eff_dist_t",          -1, 48},
    {"eff_dist_v",          -1, 49},
    {"numrel_mode_min",     -1, 50},
    {"numrel_mode_max",     -1, 51},
    {"numrel_data",         -1, 52},
    {"amp_order",           -1, 53},
    {"taper",               -1, 54},
    {"bandpass",            -1, 55},
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

  /* open the sim_inspiral table file */
  mioStatus = MetaioOpenTable( env, fileName, "sim_inspiral" );
  if ( mioStatus )
  {
    fprintf( stderr, "error opening sim_inspiral table from file %s\n",
        fileName );
    return -1;
  }

  /* figure out the column positions of the simulated parameters */
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
    /* JC: AGAIN... HOPE PARENTHESES ARE RIGHT! */
    if ( ! endTime || ( geo_time > startTime && geo_time < endTime ) )
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
          snprintf(thisSim->waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
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
          snprintf(thisSim->source, LIGOMETA_SOURCE_MAX * sizeof(CHAR),
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
          thisSim->f_lower = r4colData;
        }
        else if ( tableDir[j].idx == 43 )
        {
          thisSim->f_final = r4colData;
        }
        else if ( tableDir[j].idx == 44 )
        {
          thisSim->mchirp = r4colData;
        }
        else if ( tableDir[j].idx == 45 )
        {
          thisSim->eff_dist_h = r4colData;
        }
        else if ( tableDir[j].idx == 46 )
        {
          thisSim->eff_dist_l = r4colData;
        }
        else if ( tableDir[j].idx == 47 )
        {
          thisSim->eff_dist_g = r4colData;
        }
        else if ( tableDir[j].idx == 48 )
        {
          thisSim->eff_dist_t = r4colData;
        }
        else if ( tableDir[j].idx == 49 )
        {
          thisSim->eff_dist_v = r4colData;
        }
	else if ( tableDir[j].idx == 50 )
	{
	  thisSim->numrel_mode_min = i4colData;
	}
	else if ( tableDir[j].idx == 51 )
	{
	  thisSim->numrel_mode_max = i4colData;
	}
	else if ( tableDir[j].idx == 52 )
	{
          snprintf(thisSim->numrel_data, LIGOMETA_STRING_MAX * sizeof(CHAR),
              "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
	}
        else if ( tableDir[j].idx == 53 )
        {
            thisSim->amp_order = i4colData;
        }
        else if ( tableDir[j].idx == 54 )
        {
            snprintf(thisSim->taper, LIGOMETA_INSPIRALTAPER_MAX * sizeof(CHAR),
                    "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
        }
        else if ( tableDir[j].idx == 55 )
        {
            thisSim->bandpass = i4colData;
        }
        else if ( tableDir[j].idx == 56 ) {
        	thisSim->qmParameter1 = r4colData;
        }
        else if ( tableDir[j].idx == 57 ) {
        	thisSim->qmParameter2 = r4colData;
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


#define CLOBBER_VAL \
  while ( *sumHead ) \
{ \
  thisValue = *sumHead; \
  *sumHead = (*sumHead)->next; \
  LALFree( thisValue ); \
  thisValue = NULL; \
}



int
SummValueTableFromLIGOLw (
    SummValueTable **sumHead,
    CHAR           *fileName
    )

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
    MetaioAbort(env);
    MetaioClose(env);
    fprintf( stderr,
        "warning: unable to open summ_value table from file %s\n",
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
        snprintf( thisValue->program, LIGOMETA_PROGRAM_MAX * sizeof(CHAR),
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
        snprintf( thisValue->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR),
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 6 )
      {
        snprintf( thisValue->name, LIGOMETA_SUMMVALUE_NAME_MAX *
            sizeof(CHAR), "%s",
            env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 7 )
      {
        thisValue->value = r4colData;
      }
      else if ( tableDir[j].idx == 8 )
      {
        snprintf( thisValue->comment, LIGOMETA_SUMMVALUE_COMM_MAX *
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


int
XLALReadSummValueFile (
    SummValueTable **summValueList,
    CHAR                  *fileName
    )

{
#if 0
  INT4 numFileTriggers = 0;
#endif
  INT4 haveSummValue = 0;
  SummValueTable  *thisSummValue = NULL;
  SummValueTable  *inputSummValue = NULL;


  /* read in the summ value table and store */
  XLALPrintInfo(
      "XLALReadSummValueFile(): Reading summ_value table\n");

  haveSummValue = SummValueTableFromLIGOLw(&inputSummValue, fileName);

  if ( ! inputSummValue || haveSummValue < 1)
  {
    XLALPrintInfo("No valid summ_value table in %s\n", fileName );
  }
  else
  {
    /* store the summ value table in SummValueList list */
    XLALPrintInfo("XLALReadSummValueFile(): Checking summ_value_table\n");

    /* keep only relevant information (inspiral_effective_distance)*/
    XLALCleanSummValueTable(&inputSummValue);

    if ( ! *summValueList )
    {
      *summValueList = thisSummValue = inputSummValue;
    }
    else
    {
      for ( thisSummValue = *summValueList; thisSummValue->next;
          thisSummValue = thisSummValue->next);
      thisSummValue = thisSummValue->next = inputSummValue;
    }
  }
 return 1;
}


int
XLALReadInspiralTriggerFile (
    SnglInspiralTable    **inspiralEventList,
    SnglInspiralTable    **lastTrigger,
    SearchSummaryTable   **searchSummList,
    SearchSummvarsTable  **inputFileList,
    CHAR                  *fileName
    )

{
  INT4 numFileTriggers = 0;
  /*INT4 haveSummValue = 0;*/
  SnglInspiralTable  *inputData = NULL;
  SearchSummaryTable *inputSummary = NULL;
  SearchSummaryTable *thisSearchSumm = NULL;
  SearchSummvarsTable  *thisInputFile = NULL;
#if 0
  SummValueTable  *thisSummValue = NULL;
  SummValueTable  *inputSummValue = NULL;
#endif


  /* store the file name in search summvars */
  XLALPrintInfo(
      "XLALReadInspiralTriggerFile(): storing input file name %s\n"
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
  snprintf( thisInputFile->name, LIGOMETA_NAME_MAX,
      "input_file" );
  snprintf( thisInputFile->string, LIGOMETA_NAME_MAX,
      "%s", fileName );


  /* read in the search summary and store */
  XLALPrintInfo(
      "XLALReadInspiralTriggerFile(): Reading search_summary table\n");

  inputSummary = XLALSearchSummaryTableFromLIGOLw(fileName);

  if ( ! inputSummary )
  {
    XLALPrintError("No valid search_summary table in %s, exiting\n",
        fileName );
    LALFree(thisInputFile);
    XLAL_ERROR(XLAL_EIO);
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
#if 0
  /* read in the summ value table and store */
  XLALPrintInfo(
      "XLALReadInspiralTriggerFile(): Reading summ_value table\n");

  haveSummValue = SummValueTableFromLIGOLw(&inputSummValue, fileName);

  if ( ! inputSummValue || haveSummValue < 1)
  {
    XLALPrintError("No valid summ_value table in %s, exiting\n",
        fileName );
    LALFree(thisInputFile);
    XLAL_ERROR(XLAL_EIO);
  }
  else
  {
    /* store the summ value table in SummValueList list */
    XLALPrintInfo("XLALReadInspiralTriggerFile(): Checking summ_value_table\n");

    /* keep only relevant information (inspiral_effective_distance)*/
    XLALCleanSummValueTable(&inputSummValue);

    if ( ! *summValueList )
    {
      *summValueList = thisSummValue = inputSummValue;
    }
    else
    {
      for ( thisSummValue = *summValueList; thisSummValue->next;
          thisSummValue = thisSummValue->next);
      thisSummValue = thisSummValue->next = inputSummValue;
    }
  }
#endif
  /* read in the triggers */
  numFileTriggers =
    LALSnglInspiralTableFromLIGOLw( &inputData, fileName, 0, -1 );

  if ( numFileTriggers < 0 )
  {
    XLALPrintError("Unable to read sngl_inspiral table from %s\n",
        fileName );
    LALFree(thisInputFile);
    XLAL_ERROR(XLAL_EIO);
  }
  else if ( numFileTriggers > 0 )
  {

    XLALPrintInfo(
        "XLALReadInspiralTriggerFile(): Got %d sngl_inspiral rows from %s\n",
        numFileTriggers, fileName );

    /* store the triggers */
    if ( ! *inspiralEventList )
    {
      /* store the head of the linked list */
      *inspiralEventList = *lastTrigger = inputData;
    }
    else
    {
      /* append to the end of the linked list and set current    */
      /* trigger to the first trigger of the list being appended */
      *lastTrigger = (*lastTrigger)->next = inputData;
    }

    /* scroll to the end of the linked list of triggers */
    for ( ; (*lastTrigger)->next; *lastTrigger = (*lastTrigger)->next );
  }

  return( numFileTriggers );
}



/* function which reads a summ_value table and remove
   all rows which do not contain name set
   to "effective_inspiral_distance". */
void XLALCleanSummValueTable(SummValueTable **inputSummValue)
{
  SummValueTable *this = NULL;
  SummValueTable *prev = NULL;
  SummValueTable *head = NULL;

  this  = *inputSummValue;
  head = NULL;

  while( this )
  {
    INT4 discard = 0;
    SummValueTable *tmp = this;

    this = this->next;

    if (strncmp( tmp->name, "inspiral_effective_distance",
	LIGOMETA_SUMMVALUE_NAME_MAX) )
    {
      /* not an effective distance -- discard */
      XLALPrintInfo(
  	"XLALReadIspiralTriggerFile(): Removing entry with  \"%s\" name \n",
  	 tmp->name);
      discard = 1;
    }
    else
    {
      discard = 0;
      XLALPrintInfo(
	"XLALReadIspiralTriggerFile(): Got inspiral effective distance of %f for a %s system.\n",
	tmp->value, tmp->comment);
    }

    if (discard)
    {
      LALFree(tmp);
    }
    else
    {
      if (!head)
      {
        head = tmp;
      }
      else
      {
        prev->next = tmp;
      }
      tmp->next = NULL;
      prev = tmp;
    }
  }
  *inputSummValue = head;
}


#define CLOBBER_EVENTS \
  while ( *eventHead ) \
{ \
  thisEvent = *eventHead; \
  *eventHead = (*eventHead)->next; \
  LALFree( thisEvent ); \
  thisEvent = NULL; \
}


int
LALMultiInspiralTableFromLIGOLw (
    MultiInspiralTable **eventHead,
    CHAR                *fileName
    )

{
  int                                   i, j, nrows;
  int                                   mioStatus;
  MultiInspiralTable                   *thisEvent = NULL;
  struct MetaioParseEnvironment         parseEnv;
  const  MetaioParseEnv                 env = &parseEnv;
  MetaTableDirectory tableDir[] =
    {
          {"ifos",                    -1, 0},
          {"search",                  -1, 1},
          {"end_time",                -1, 2},
          {"end_time_ns",             -1, 3},
          {"end_time_gmst",           -1, 4},
          {"impulse_time",            -1, 5},
          {"impulse_time_ns",         -1, 6},
          {"amplitude",               -1, 7},
          {"distance"         ,       -1, 8},
          {"eff_dist_h1",             -1, 9},
          {"eff_dist_h2",             -1, 10},
          {"eff_dist_l",              -1, 11},
          {"eff_dist_g",              -1, 12},
          {"eff_dist_t",              -1, 13},
          {"eff_dist_v",              -1, 14},
          {"eff_dist_h1h2",           -1, 15},
          {"chi",                     -1, 16},
          {"kappa",                   -1, 17},
          {"coa_phase",               -1, 18},
          {"mass1",                   -1, 19},
          {"mass2",                   -1, 20},
          {"mchirp",                  -1, 21},
          {"eta",                     -1, 22},
          {"tau0",                    -1, 23},
          {"tau2",                    -1, 24},
          {"tau3",                    -1, 25},
          {"tau4",                    -1, 26},
          {"tau5",                    -1, 27},
          {"ttotal",                  -1, 28},
          {"snr",                     -1, 29},
          {"snr_dof",                 -1, 30},
          {"chisq",                   -1, 31},
          {"chisq_dof",               -1, 32},
          {"bank_chisq",              -1, 33},
          {"bank_chisq_dof",          -1, 34},
          {"cont_chisq",              -1, 35},
          {"cont_chisq_dof",          -1, 36},
          {"trace_snr",               -1, 37},
          {"snr_h1",                  -1, 38},
          {"snr_h2",                  -1, 39},
          {"snr_l",                   -1, 40},
          {"snr_g",                   -1, 41},
          {"snr_t",                   -1, 42},
          {"snr_v",                   -1, 43},
          {"amp_term_1",              -1, 44},
          {"amp_term_2",              -1, 45},
          {"amp_term_3",              -1, 46},
          {"amp_term_4",              -1, 47},
          {"amp_term_5",              -1, 48},
          {"amp_term_6",              -1, 49},
          {"amp_term_7",              -1, 50},
          {"amp_term_8",              -1, 51},
          {"amp_term_9",              -1, 52},
          {"amp_term_10",             -1, 53},
          {"sigmasq_h1",              -1, 54},
          {"sigmasq_h2",              -1, 55},
          {"sigmasq_l",               -1, 56},
          {"sigmasq_g",               -1, 57},
          {"sigmasq_t",               -1, 58},
          {"sigmasq_v",               -1, 59},
          {"chisq_h1",                -1, 60},
          {"chisq_h2",                -1, 61},
          {"chisq_l",                 -1, 62},
          {"chisq_g",                 -1, 63},
          {"chisq_t",                 -1, 64},
          {"chisq_v",                 -1, 65},
          {"sngl_chisq_dof",          -1, 66},
          {"bank_chisq_h1",           -1, 67},
          {"bank_chisq_h2",           -1, 68},
          {"bank_chisq_l",            -1, 69},
          {"bank_chisq_g",            -1, 70},
          {"bank_chisq_t",            -1, 71},
          {"bank_chisq_v",            -1, 72},
          {"sngl_bank_chisq_dof",     -1, 73},
          {"cont_chisq_h1",           -1, 74},
          {"cont_chisq_h2",           -1, 75},
          {"cont_chisq_l",            -1, 76},
          {"cont_chisq_g",            -1, 77},
          {"cont_chisq_t",            -1, 78},
          {"cont_chisq_v",            -1, 79},
          {"sngl_cont_chisq_dof",     -1, 80},
          {"ra",                      -1, 81},
          {"dec",                     -1, 82},
          {"ligo_angle",              -1, 83},
          {"ligo_angle_sig",          -1, 84},
          {"inclination",             -1, 85},
          {"polarization",            -1, 86},
          {"null_statistic",          -1, 87},
	  {"null_stat_h1h2",          -1, 88},
	  {"null_stat_degen",         -1, 89},
          {"event_id",                -1, 90},
          {"h1quad_re",               -1, 91},
          {"h1quad_im",               -1, 92},
          {"h2quad_re",               -1, 93},
          {"h2quad_im",               -1, 94},
          {"l1quad_re",               -1, 95},
          {"l1quad_im",               -1, 96},
          {"g1quad_re",               -1, 97},
          {"g1quad_im",               -1, 98},
          {"t1quad_re",               -1, 99},
          {"t1quad_im",               -1, 100},
          {"v1quad_re",               -1, 101},
          {"v1quad_im",               -1, 102},
          {"coh_snr_h1h2",            -1, 103},
          {"cohSnrSqLocal",           -1, 104},
          {"autoCorrCohSq",           -1, 105},
          {"crossCorrCohSq",          -1, 106},
          {"autoCorrNullSq",          -1, 107},
          {"crossCorrNullSq",         -1, 108},
          {"ampMetricEigenVal1",      -1, 109},
          {"ampMetricEigenVal2",      -1, 110},
          {"time_slide_id",           -1, 111},
          {NULL,                       0, 0}
    };

  /* check that the bank handle and pointer are valid */
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

  /* open the multi_inspiral XML file */
  mioStatus = MetaioOpenFile( env, fileName );
  if ( mioStatus )
  {
    fprintf( stderr, "unable to open file %s\n", fileName );
    return -1;
  }

  /* open the multi_inspiral table template bank file */
  mioStatus = MetaioOpenTableOnly( env, "multi_inspiral" );
  if ( mioStatus )
  {
    fprintf( stdout, "no multi_inspiral table in file %s\n", fileName );
    return 0;
  }

  /* figure out the column positions of the template parameters */
  for ( i = 0; tableDir[i].name; ++i )
  {
    if ( (tableDir[i].pos = MetaioFindColumn( env, tableDir[i].name )) < 0 )
    {
      fprintf( stderr, "unable to find column %s\n", tableDir[i].name );

      if ( ! strcmp(tableDir[i].name, "event_id") )
      {
        fprintf( stderr,
            "The event_id column is not populated, continuing anyway\n");
      }
      else
      {
        MetaioClose(env);
        return -1;
      }
    }
  }

  /* loop over the rows in the file */
  i = nrows = 0;
  while ( (mioStatus = MetaioGetRow(env)) == 1 )
  {
    /* count the rows in the file */
    i++;
    /* allocate memory for the template we are about to read in */
    if ( ! *eventHead )
    {
      thisEvent = *eventHead = (MultiInspiralTable *)
        LALCalloc( 1, sizeof(MultiInspiralTable) );
    }
    else
    {
      thisEvent = thisEvent->next = (MultiInspiralTable *)
        LALCalloc( 1, sizeof(MultiInspiralTable) );
    }
    if ( ! thisEvent )
    {
      fprintf( stderr, "could not allocate multi inspiral event\n" );
        CLOBBER_EVENTS;
        MetaioClose( env );
        return -1;
    }

    /* parse the contents of the row into the MultiInspiralTable structure */
    for ( j = 0; tableDir[j].name; ++j )
    {
      enum METAIO_Type column_type = env->ligo_lw.table.col[tableDir[j].pos].data_type;
      REAL4 r4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_4;
      REAL8 r8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.real_8;
      INT4  i4colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_4s;

      if ( tableDir[j].pos < 0 ) continue;

      /* dereference the data stored in the table */
      if ( tableDir[j].idx == 0 )
      {
        snprintf( thisEvent->ifos, LIGOMETA_IFOS_MAX * sizeof(CHAR),
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 1 )
      {
        snprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data );
      }
      else if ( tableDir[j].idx == 2 )
      {
        thisEvent->end_time.gpsSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 3 )
      {
        thisEvent->end_time.gpsNanoSeconds = i4colData;
      }
      else if ( tableDir[j].idx == 4 )
      {
        thisEvent->end_time_gmst = r8colData;
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
        thisEvent->amplitude = r4colData;
      }
      else if ( tableDir[j].idx == 8 )
      {
        thisEvent->distance = r4colData;
      }
      else if ( tableDir[j].idx == 9 )
      {
        thisEvent->eff_dist_h1 = r4colData;
      }
      else if ( tableDir[j].idx == 10 )
      {
        thisEvent->eff_dist_h2 = r4colData;
      }
      else if ( tableDir[j].idx == 11 )
      {
        thisEvent->eff_dist_l = r4colData;
      }
      else if ( tableDir[j].idx == 12 )
      {
        thisEvent->eff_dist_g = r4colData;
      }
      else if ( tableDir[j].idx == 13 )
      {
        thisEvent->eff_dist_t = r4colData;
      }
      else if ( tableDir[j].idx == 14 )
      {
        thisEvent->eff_dist_v = r4colData;
      }
      else if ( tableDir[j].idx == 15 )
      {
        thisEvent->eff_dist_h1h2 = r4colData;
      }
      else if ( tableDir[j].idx == 16 )
      {
        thisEvent->chi = r4colData;
      }
      else if ( tableDir[j].idx == 17 )
      {
        thisEvent->kappa = r4colData;
      }
      else if ( tableDir[j].idx == 18 )
      {
        thisEvent->coa_phase = r4colData;
      }
      else if ( tableDir[j].idx == 19 )
      {
        thisEvent->mass1 = r4colData;
      }
      else if ( tableDir[j].idx == 20 )
      {
        thisEvent->mass2 = r4colData;
      }
      else if ( tableDir[j].idx == 21 )
      {
        thisEvent->mchirp = r4colData;
      }
      else if ( tableDir[j].idx == 22 )
      {
        thisEvent->eta = r4colData;
      }
      else if ( tableDir[j].idx == 23 )
      {
        thisEvent->tau0 = r4colData;
      }
      else if ( tableDir[j].idx == 24 )
      {
        thisEvent->tau2 = r4colData;
      }
      else if ( tableDir[j].idx == 25 )
      {
        thisEvent->tau3 = r4colData;
      }
      else if ( tableDir[j].idx == 26 )
      {
        thisEvent->tau4 = r4colData;
      }
      else if ( tableDir[j].idx == 27 )
      {
        thisEvent->tau5 = r4colData;
      }
      else if ( tableDir[j].idx == 28 )
      {
        thisEvent->ttotal = r4colData;
      }
      else if ( tableDir[j].idx == 29 )
      {
        thisEvent->snr = r4colData;
      }
      else if ( tableDir[j].idx == 30 )
      {
        thisEvent->snr_dof = i4colData;
      }
      else if ( tableDir[j].idx == 31 )
      {
        thisEvent->chisq = r4colData;
      }
      else if ( tableDir[j].idx == 32 )
      {
        thisEvent->chisq_dof = i4colData;
      }
      else if ( tableDir[j].idx == 33 )
      {
        thisEvent->bank_chisq = r4colData;
      }
      else if ( tableDir[j].idx == 34 )
      {
        thisEvent->bank_chisq_dof = i4colData;
      }
      else if ( tableDir[j].idx == 35 )
      {
        thisEvent->cont_chisq = r4colData;
      }
      else if ( tableDir[j].idx == 36 )
      {
        thisEvent->cont_chisq_dof = i4colData;
      }
      else if ( tableDir[j].idx == 37 )
      {
        thisEvent->trace_snr = r4colData;
      }
      else if ( tableDir[j].idx == 38 )
      {
        thisEvent->snr_h1 = r4colData;
      }
      else if ( tableDir[j].idx == 39 )
      {
        thisEvent->snr_h2 = r4colData;
      }
      else if ( tableDir[j].idx == 40 )
      {
        thisEvent->snr_l = r4colData;
      }
      else if ( tableDir[j].idx == 41 )
      {
        thisEvent->snr_g = r4colData;
      }
      else if ( tableDir[j].idx == 42 )
      {
        thisEvent->snr_t = r4colData;
      }
      else if ( tableDir[j].idx == 43 )
      {
        thisEvent->snr_v = r4colData;
      }
      else if ( tableDir[j].idx == 44 )
      {
        thisEvent->amp_term_1 = r4colData;
      }
      else if ( tableDir[j].idx == 45 )
      {
        thisEvent->amp_term_2 = r4colData;
      }
      else if ( tableDir[j].idx == 46 )
      {
        thisEvent->amp_term_3 = r4colData;
      }
      else if ( tableDir[j].idx == 47 )
      {
        thisEvent->amp_term_4 = r4colData;
      }
      else if ( tableDir[j].idx == 48 )
      {
        thisEvent->amp_term_5 = r4colData;
      }
      else if ( tableDir[j].idx == 49 )
      {
        thisEvent->amp_term_6 = r4colData;
      }
      else if ( tableDir[j].idx == 50 )
      {
        thisEvent->amp_term_7 = r4colData;
      }
      else if ( tableDir[j].idx == 51 )
      {
        thisEvent->amp_term_8 = r4colData;
      }
      else if ( tableDir[j].idx == 52 )
      {
        thisEvent->amp_term_9 = r4colData;
      }
      else if ( tableDir[j].idx == 53 )
      {
        thisEvent->amp_term_10 = r4colData;
      }
      else if ( tableDir[j].idx == 54 )
      {
        thisEvent->sigmasq_h1 = r8colData;
      }
      else if ( tableDir[j].idx == 55 )
      {
        thisEvent->sigmasq_h2 = r8colData;
      }
      else if ( tableDir[j].idx == 56 )
      {
        thisEvent->sigmasq_l = r8colData;
      }
      else if ( tableDir[j].idx == 57 )
      {
        thisEvent->sigmasq_g = r8colData;
      }
      else if ( tableDir[j].idx == 58 )
      {
        thisEvent->sigmasq_t = r8colData;
      }
      else if ( tableDir[j].idx == 59 )
      {
        thisEvent->sigmasq_v = r8colData;
      }
      else if ( tableDir[j].idx == 60 )
      {
        thisEvent->chisq_h1 = r4colData;
      }
      else if ( tableDir[j].idx == 61 )
      {
        thisEvent->chisq_h2 = r4colData;
      }
      else if ( tableDir[j].idx == 62 )
      {
        thisEvent->chisq_l = r4colData;
      }
      else if ( tableDir[j].idx == 63 )
      {
        thisEvent->chisq_g = r4colData;
      }
      else if ( tableDir[j].idx == 64 )
      {
        thisEvent->chisq_t = r4colData;
      }
      else if ( tableDir[j].idx == 65 )
      {
        thisEvent->chisq_v = r4colData;
      }
      else if ( tableDir[j].idx == 66 )
      {
        thisEvent->sngl_chisq_dof = i4colData;
      }
      else if ( tableDir[j].idx == 67 )
      {
        thisEvent->bank_chisq_h1 = r4colData;
      }
      else if ( tableDir[j].idx == 68 )
      {
        thisEvent->bank_chisq_h2 = r4colData;
      }
      else if ( tableDir[j].idx == 69 )
      {
        thisEvent->bank_chisq_l = r4colData;
      }
      else if ( tableDir[j].idx == 70 )
      {
        thisEvent->bank_chisq_g = r4colData;
      }
      else if ( tableDir[j].idx == 71 )
      {
        thisEvent->bank_chisq_t = r4colData;
      }
      else if ( tableDir[j].idx == 72 )
      {
        thisEvent->bank_chisq_v = r4colData;
      }
      else if ( tableDir[j].idx == 73 )
      {
        thisEvent->sngl_bank_chisq_dof = i4colData;
      }
      else if ( tableDir[j].idx == 74 )
      {
        thisEvent->cont_chisq_h1 = r4colData;
      } 
      else if ( tableDir[j].idx == 75 )
      {
        thisEvent->cont_chisq_h2 = r4colData;
      } 
      else if ( tableDir[j].idx == 76 )
      {
        thisEvent->cont_chisq_l = r4colData;
      } 
      else if ( tableDir[j].idx == 77 )
      {
        thisEvent->cont_chisq_g = r4colData;
      }
      else if ( tableDir[j].idx == 78 )
      {
        thisEvent->cont_chisq_t = r4colData;
      }
      else if ( tableDir[j].idx == 79 )
      {
        thisEvent->cont_chisq_v = r4colData;
      }
      else if ( tableDir[j].idx == 80 )
      {
        thisEvent->sngl_cont_chisq_dof = i4colData;
      }
      else if ( tableDir[j].idx == 81 )
      {
        thisEvent->ra = r4colData;
      }
      else if ( tableDir[j].idx == 82 )
      {
        thisEvent->dec = r4colData;
      }
      else if ( tableDir[j].idx == 83 )
      {
        thisEvent->ligo_angle = r4colData;
      }
      else if ( tableDir[j].idx == 84 )
      {
        thisEvent->ligo_angle_sig = r4colData;
      }
      else if ( tableDir[j].idx == 85 )
      {
        thisEvent->inclination = r4colData;
      }
      else if ( tableDir[j].idx == 86 )
      {
        thisEvent->polarization = r4colData;
      }
      else if ( tableDir[j].idx == 87 )
      {
        thisEvent->null_statistic = r4colData;
      }
      else if ( tableDir[j].idx == 88 )
      {
        thisEvent->null_stat_h1h2 = r4colData;
      }
      else if ( tableDir[j].idx == 89 )
      {
        thisEvent->null_stat_degen = r4colData;
      }
      else if ( tableDir[j].idx == 90 )
      {
        if ( tableDir[j].pos > 0 )
        {
          INT8 i8colData;
          if ( column_type == METAIO_TYPE_INT_8S )
            i8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_8s;
          else
          {
            i8colData = XLALLIGOLwParseIlwdChar(env, tableDir[j].pos, "multi_inspiral", "event_id");
            if ( i8colData < 0 )
              return -1;
          }
          if ( i8colData >= 0 )
          {
            thisEvent->event_id = LALCalloc( 1, sizeof(*thisEvent->event_id) );
            thisEvent->event_id->id = i8colData;
            thisEvent->event_id->multiInspiralTable = thisEvent;
          }
        }
      }
      else if ( tableDir[j].idx == 91 )
      {
        thisEvent->h1quad.re = r4colData;
      }
      else if ( tableDir[j].idx == 92 )
      {
        thisEvent->h1quad.im = r4colData;
      }
      else if ( tableDir[j].idx == 93 )
      {
        thisEvent->h2quad.re = r4colData;
      }
      else if ( tableDir[j].idx == 94 )
      {
        thisEvent->h2quad.im = r4colData;
      }
      else if ( tableDir[j].idx == 95 )
      {
        thisEvent->l1quad.re = r4colData;
      }
      else if ( tableDir[j].idx == 96 )
      {
        thisEvent->l1quad.im = r4colData;
      }
      else if ( tableDir[j].idx == 97 )
      {
        thisEvent->g1quad.re = r4colData;
      }
      else if ( tableDir[j].idx == 98 )
      {
        thisEvent->g1quad.im = r4colData;
      }
      else if ( tableDir[j].idx == 99 )
      {
        thisEvent->t1quad.re = r4colData;
      }
      else if ( tableDir[j].idx == 100 )
      {
        thisEvent->t1quad.im = r4colData;
      }
      else if ( tableDir[j].idx == 101 )
      {
        thisEvent->v1quad.re = r4colData;
      }
      else if ( tableDir[j].idx == 102 )
      {
        thisEvent->v1quad.im = r4colData;
      }
      else if ( tableDir[j].idx == 103 )
      {
        thisEvent->coh_snr_h1h2 = r4colData;
      }
      else if ( tableDir[j].idx == 104 )
      {
        thisEvent->cohSnrSqLocal = r4colData;
      }
      else if ( tableDir[j].idx == 105 )
      {
        thisEvent->autoCorrCohSq = r4colData;
      }
      else if ( tableDir[j].idx == 106 )
      {
        thisEvent->crossCorrCohSq = r4colData;
      }
      else if ( tableDir[j].idx == 107 )
      {
        thisEvent->autoCorrNullSq = r4colData;
      }
      else if ( tableDir[j].idx == 108 )
      {
        thisEvent->crossCorrNullSq = r4colData;
      }
      else if ( tableDir[j].idx == 109 )
      {
        thisEvent->ampMetricEigenVal1 = r8colData;
      }
      else if ( tableDir[j].idx == 110 )
      {
        thisEvent->ampMetricEigenVal2 = r8colData;
      }
      else if ( tableDir[j].idx == 111 )
      {
        if ( tableDir[j].pos > 0 )
        {
          INT8 i8colData;
          if ( column_type == METAIO_TYPE_INT_8S )
            i8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_8s;
          else
          {
            i8colData = XLALLIGOLwParseIlwdChar(env, tableDir[j].pos, "multi_inspiral", "time_slide_id");
            if ( i8colData < 0 )
              return -1;
          }
          if ( i8colData )
          {
            thisEvent->time_slide_id = LALCalloc( 1,\
                                        sizeof(*thisEvent->time_slide_id) );
            thisEvent->time_slide_id->id = i8colData;
            thisEvent->time_slide_id->multiInspiralTable = thisEvent;
          }
        }
      }
      else
      {
        CLOBBER_EVENTS;
        fprintf( stderr, "unknown column while parsing multi_inspiral\n" );
        return -1;
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
    return -1;
  }

  /* Normal exit */
  MetaioClose( env );
  return nrows;
}

#undef CLOBBER_EVENTS


int
XLALReadMultiInspiralTriggerFile (
    MultiInspiralTable    **inspiralEventList,
    MultiInspiralTable    **lastTrigger,
    SearchSummaryTable    **searchSummList,
    SearchSummvarsTable   **inputFileList,
    CHAR                   *fileName
    )

{
  INT4 numFileTriggers = 0;
  /*INT4 haveSummValue = 0;*/
  MultiInspiralTable   *inputData = NULL;
  SearchSummaryTable   *inputSummary = NULL;
  SearchSummaryTable   *thisSearchSumm = NULL;
  SearchSummvarsTable  *thisInputFile = NULL;
#if 0
  SummValueTable  *thisSummValue = NULL;
  SummValueTable  *inputSummValue = NULL;
#endif


  /* store the file name in search summvars */
  XLALPrintInfo(
      "XLALReadMultiInspiralTriggerFile(): storing input file name %s\n"
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
  snprintf( thisInputFile->name, LIGOMETA_NAME_MAX,
      "input_file" );
  snprintf( thisInputFile->string, LIGOMETA_NAME_MAX,
      "%s", fileName );


  /* read in the search summary and store */
  XLALPrintInfo(
      "XLALReadMultiInspiralTriggerFile(): Reading search_summary table\n");

  inputSummary = XLALSearchSummaryTableFromLIGOLw(fileName);

  if ( ! inputSummary )
  {
    XLALPrintError("No valid search_summary table in %s, exiting\n",
        fileName );
    LALFree(thisInputFile);
    XLAL_ERROR(XLAL_EIO);
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
#if 0
  /* read in the summ value table and store */
  XLALPrintInfo(
      "XLALReadMultiInspiralTriggerFile(): Reading summ_value table\n");

  haveSummValue = SummValueTableFromLIGOLw(&inputSummValue, fileName);

  if ( ! inputSummValue || haveSummValue < 1)
  {
    XLALPrintError("No valid summ_value table in %s, exiting\n",
        fileName );
    LALFree(thisInputFile);
    XLAL_ERROR(XLAL_EIO);
  }
  else
  {
    /* store the summ value table in SummValueList list */
    XLALPrintInfo("XLALReadMultiInspiralTriggerFile(): Checking summ_value_table\n");

    /* keep only relevant information (inspiral_effective_distance)*/
    XLALCleanSummValueTable(&inputSummValue);

    if ( ! *summValueList )
    {
      *summValueList = thisSummValue = inputSummValue;
    }
    else
    {
      for ( thisSummValue = *summValueList; thisSummValue->next;
          thisSummValue = thisSummValue->next);
      thisSummValue = thisSummValue->next = inputSummValue;
    }
  }
#endif
  /* read in the triggers */
  numFileTriggers =
    LALMultiInspiralTableFromLIGOLw( &inputData, fileName);

  if ( numFileTriggers < 0 )
  {
    XLALPrintError("Unable to read multi_inspiral table from %s\n",
        fileName );
    LALFree(thisInputFile);
    XLAL_ERROR(XLAL_EIO);
  }
  else if ( numFileTriggers > 0 )
  {

    XLALPrintInfo(
        "XLALReadMultiInspiralTriggerFile(): Got %d multi_inspiral rows from %s\n",
        numFileTriggers, fileName );

    /* store the triggers */
    if ( ! *inspiralEventList )
    {
      /* store the head of the linked list */
      *inspiralEventList = *lastTrigger = inputData;
    }
    else
    {
      /* append to the end of the linked list and set current    */
      /* trigger to the first trigger of the list being appended */
      *lastTrigger = (*lastTrigger)->next = inputData;
    }

    /* scroll to the end of the linked list of triggers */
    for ( ; (*lastTrigger)->next; *lastTrigger = (*lastTrigger)->next );
  }

  return( numFileTriggers );
}
