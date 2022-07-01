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
#include <lal/LALInspiral.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>
#include <InspiralTmpltBankFromLIGOLw.h>


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

      if ( strstr(tableDir[i].name, "Gamma") )
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
          INT8 i8colData;
          if ( column_type == METAIO_TYPE_INT_8S )
            i8colData = env->ligo_lw.table.elt[tableDir[j].pos].data.int_8s;
          else
          {
            i8colData = XLALLIGOLwParseIlwdChar(env, tableDir[j].pos, "sngl_inspiral", "event_id");
            if ( i8colData < 0 )
              return -1;
          }
          thisTmplt->event_id = i8colData;
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
#undef CLOBBER_EVENTS
