/*
 * LIGOLwXMLStochasticRead.c
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <metaio.h>

#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXMLStochasticRead.h>

#define CLOBBER_STOCH_VAL \
  while (*stochHead) \
{ \
  thisValue = *stochHead; \
  *stochHead = (*stochHead)->next; \
  LALFree( thisValue ); \
  thisValue = NULL; \
}

int
LALStochasticTableFromLIGOLw (
    StochasticTable **stochHead,
    CHAR *fileName)
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
        snprintf(thisValue->ifo_one, LIGOMETA_IFO_MAX * sizeof(CHAR), \
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if (tableDir[j].idx == 1)
      {
        snprintf(thisValue->ifo_two, LIGOMETA_IFO_MAX * sizeof(CHAR), \
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if ( tableDir[j].idx == 2 )
      {
        snprintf(thisValue->channel_one, LIGOMETA_CHANNEL_MAX * \
            sizeof(CHAR), "%s", \
            env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if ( tableDir[j].idx == 3 )
      {
        snprintf(thisValue->channel_two, LIGOMETA_CHANNEL_MAX * \
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
  while (*stochSummHead) \
{ \
  thisValue = *stochSummHead; \
  *stochSummHead = (*stochSummHead)->next; \
  LALFree( thisValue ); \
  thisValue = NULL; \
}


int
LALStochSummTableFromLIGOLw (
    StochSummTable **stochSummHead,
    CHAR *fileName)

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
        snprintf(thisValue->ifo_one, LIGOMETA_IFO_MAX * sizeof(CHAR), \
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if (tableDir[j].idx == 1)
      {
        snprintf(thisValue->ifo_two, LIGOMETA_IFO_MAX * sizeof(CHAR), \
            "%s", env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if ( tableDir[j].idx == 2 )
      {
        snprintf(thisValue->channel_one, LIGOMETA_CHANNEL_MAX * \
            sizeof(CHAR), "%s", \
            env->ligo_lw.table.elt[tableDir[j].pos].data.lstring.data);
      }
      else if ( tableDir[j].idx == 3 )
      {
        snprintf(thisValue->channel_two, LIGOMETA_CHANNEL_MAX * \
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
