/*----------------------------------------------------------------------- 
 * 
 * File Name: ligolw_tmpltbank.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include "ligolwbank.h"

#define CLOBBER_BANK \
    while ( *bankHead ); \
    { \
      thisTmplt = *bankHead; \
      *bankHead = (*bankHead)->next; \
      LALFree( thisTmplt ); \
      thisTmplt = NULL; \
    }

struct MetaTableDirectory
{
  char *name;
  int   pos;
  int   idx;
};

int
InspiralTmpltBankFromLIGOLw (
    InspiralTemplate  **bankHead,
    CHAR               *fileName,
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
  struct MetaTableDirectory tableDir[] =
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
    MetaioClose( env );
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
          "--minimum-match" ) )
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
    fprintf( stderr, "error opening sngl_inspiral table from file %s\n", 
        fileName );
    MetaioClose( env );
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

        switch ( tableDir[j].idx )
        {
          case 0:
            thisTmplt->mass1 = colData;
            break;
          case 1:
            thisTmplt->mass2 = colData;
            break;
          case 2:
            thisTmplt->chirpMass = colData;
            break;
          case 3:
            thisTmplt->eta = colData;
            break;
          case 4:
            thisTmplt->t0 = colData;
            break;
          case 5:
            thisTmplt->t2 = colData;
            break;
          case 6:
            thisTmplt->t3 = colData;
            break;
          case 7:
            thisTmplt->t4 = colData;
            break;
          case 8:
            thisTmplt->t5 = colData;
            break;
          case 9:
            thisTmplt->tC = colData;
            break;
          default:
            CLOBBER_BANK;
            fprintf( stderr, "unknown column while parsing\n" );
            return -1;
        }
      }

      /* compute derived mass parameters */
      thisTmplt->totalMass = thisTmplt->mass1 + thisTmplt->mass2;
      thisTmplt->mu = thisTmplt->mass1 * thisTmplt->mass2 / 
        thisTmplt->totalMass;
      
      /* set the match determined from the bank generation process params */
      thisTmplt->minMatch = minMatch;

      /* increase the count of rows parsed */
      ++nrows;
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

#if 0
int
InspiralTmpltBankToLIGOLw (
    InspiralTemplate   **bankHead,
    CHAR                *fileName,
    INT4                startTmplt,
    INT4                stopTmplt
    )
{
  *bankHead = NULL;
  fileName = NULL;
}
#endif

#undef CLOBBER_BANK
