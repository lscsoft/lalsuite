/************************************ <lalVerbatim file="ConfigFileCV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{ConfigFile}
\label{ss:ConfigFile.c}

Some general-purpose routines for config-file reading

\subsubsection*{Prototypes}
\input{ConfigFileCP}
\idx{NextSkyPosition()}

\subsubsection*{Description}

allow reading of "variable = value" type config-files, that may also contain 
comments

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{ConfigFileCV}}

******************************************************* </lalLaTeX> */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALStdio.h>
#include <lal/StreamInput.h>

#include "ConfigFile.h"

NRCSID( CONFIGFILEC, "$Id$" );

extern INT4 lalDebugLevel;


#define ERR -1
#define OK  0

#define FMT_STRING "string"    /* reading in quoted strings needs some special treatment */
#define WHITESPACE " \t"
#define COMMENT_CHARS	"#%"

/* local prototypes */
static void cleanConfig (CHARSequence *text);

/*----------------------------------------------------------------------
 * parse a config-file stream into a token-list
 * 
 * gets rid of comments, empty lines, and does line-continuation of '\'
 *
 *----------------------------------------------------------------------*/
void
LALLoadConfigFile (LALStatus *stat, LALConfigData_t *cfgdata, FILE *instream)
{
  CHARSequence *rawdata = NULL;
  UINT4 i;

  INITSTATUS( stat, "LALLoadConfigFile", CONFIGFILEC );
  ATTATCHSTATUSPTR (stat);

  ASSERT (cfgdata->lines == NULL, stat, CONFIGFILEH_ENONULL, CONFIGFILEH_MSGENONULL);
  ASSERT (cfgdata->wasRead == NULL, stat, CONFIGFILEH_ENONULL, CONFIGFILEH_MSGENONULL);
  ASSERT (instream != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);

  TRY (LALCHARReadSequence (stat->statusPtr, &rawdata, instream), stat);

  if (rawdata == NULL) {
    ABORT (stat, CONFIGFILEH_EFILE, CONFIGFILEH_MSGEFILE);
  }

  /* get rid of comments and do line-continuation */
  cleanConfig (rawdata);

  /* parse this into individual lines */
  TRY (LALCreateTokenList (stat->statusPtr, &(cfgdata->lines), rawdata->data, "\n"), stat);

  LALFree (rawdata->data);
  LALFree (rawdata);

  if (lalDebugLevel >= 3)
    {
      LALPrintError ("ConfigFile DEBUG: parsed config-file contents:\n");
      for (i=0; i < cfgdata->lines->nTokens; i++)
	printf ( "%d: '%s'\n", i, cfgdata->lines->tokens[i] );
    }


  /* initialize the 'wasRead' flags for the lines */
  cfgdata->wasRead = LALCalloc (1, cfgdata->lines->nTokens * sizeof(cfgdata->wasRead[0]));

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* LALLoadConfigFile() */

/*----------------------------------------------------------------------
 * free memory associated with a LALConfigData_t structure
 *----------------------------------------------------------------------*/
void
LALDestroyConfigData (LALStatus *stat, LALConfigData_t *cfgdata)
{
  INITSTATUS( stat, "LALDestroyConfigData", CONFIGFILEC );
  ATTATCHSTATUSPTR (stat);

  ASSERT (cfgdata, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT (cfgdata->lines, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT (cfgdata->wasRead, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);

  TRY ( LALDestroyTokenList (stat->statusPtr, &(cfgdata->lines)), stat);
  cfgdata->lines = NULL;
  LALFree (cfgdata->wasRead);
  cfgdata->wasRead = NULL;

  DETATCHSTATUSPTR (stat);
  RETURN (stat);
} /* LALDestroyConfigData() */

/*----------------------------------------------------------------------
 * specialization to BOOLEAN variables
 *----------------------------------------------------------------------*/
void
LALReadConfigBOOLVariable (LALStatus *stat, BOOLEAN *varp, LALConfigData_t *cfgdata, CHAR *varName)
{
  static LALConfigVar_t param;
  INT2 tmp;

  param.varName = varName;
  param.fmt = "%" LAL_INT2_FORMAT;
  param.strictness = CONFIGFILE_WARN;

  LALReadConfigVariable (stat, (void*) &tmp, cfgdata, &param);
  
  *varp = (BOOLEAN) tmp;

  return;
}


/*----------------------------------------------------------------------
 * specialization to INT2 variables
 *----------------------------------------------------------------------*/
void
LALReadConfigINT2Variable (LALStatus *stat, INT2 *varp, LALConfigData_t *cfgdata, CHAR *varName)
{
  static LALConfigVar_t param;

  param.varName = varName;
  param.fmt = "%" LAL_INT2_FORMAT;
  param.strictness = CONFIGFILE_WARN;

  LALReadConfigVariable (stat, (void*) varp, cfgdata, &param);
  
  return;
}

/*----------------------------------------------------------------------
 * specialization to INT4 variables
 *----------------------------------------------------------------------*/
void
LALReadConfigINT4Variable (LALStatus *stat, INT4 *varp, LALConfigData_t *cfgdata, CHAR *varName)
{
  static LALConfigVar_t param;

  param.varName = varName;
  param.fmt = "%" LAL_INT4_FORMAT;
  param.strictness = CONFIGFILE_WARN;

  LALReadConfigVariable (stat, (void*) varp, cfgdata, &param);
  
  return;
}

/*----------------------------------------------------------------------
 * specialization to REAL4 variables
 *----------------------------------------------------------------------*/
void
LALReadConfigREAL4Variable (LALStatus *stat, REAL4 *varp, LALConfigData_t *cfgdata, CHAR *varName)
{
  static LALConfigVar_t param;

  param.varName = varName;
  param.fmt = "%" LAL_REAL4_FORMAT;
  param.strictness = CONFIGFILE_WARN;

  LALReadConfigVariable (stat, (void*) varp, cfgdata, &param);
  
  return;
}

/*----------------------------------------------------------------------
 * specialization to REAL8 variables
 *----------------------------------------------------------------------*/
void
LALReadConfigREAL8Variable (LALStatus *stat, REAL8 *varp, LALConfigData_t *cfgdata, CHAR *varName)
{
  static LALConfigVar_t param;

  param.varName = varName;
  param.fmt = "%" LAL_REAL8_FORMAT;
  param.strictness = CONFIGFILE_WARN;

  LALReadConfigVariable (stat, (void*) varp, cfgdata, &param);
  
  return;
}

/*----------------------------------------------------------------------
 * specialization to STRING variables 
 * NOTE: this means the rest of the line after the variable, and NOT "%s" ! 
 * here we need the pointer to the char-pointer
 *----------------------------------------------------------------------*/
void
LALReadConfigSTRINGVariable (LALStatus *stat, CHAR **varp, LALConfigData_t *cfgdata, CHAR *varName)
{
  static LALConfigVar_t param;

  param.varName = varName;
  param.fmt = FMT_STRING;
  param.strictness = CONFIGFILE_WARN;

  LALReadConfigVariable (stat, (void*) varp, cfgdata, &param);
  
  return;
}

/*----------------------------------------------------------------------
 *  parser for config-file: can read config-variables of the form
 *	VARIABLE [=:] VALUE
 * input is a TokenList containing the 'logical' lines of the cleaned config-file
 *
 * param->varName is the name of the config-variable to read
 * param->fmt    is the format string to use for reading
 *  
 * NOTE1: a special format-string is FMT_STRING, which means read the whole remaining line 
 *   which is different from "%s"! (reads only one word)
 *   In this case, this also does the memory-allocation!
 *
 * ----------------------------------------------------------------------*/
/* <lalVerbatim file="ConfigFileCP"> */
void
LALReadConfigVariable (LALStatus *stat, void *varp, LALConfigData_t *cfgdata, LALConfigVar_t *param)
{ /* </lalVerbatim> */
  CHAR *found = NULL;
  INT2 ret = 0;

  UINT4 i;
  INT4 linefound = -1;
  size_t len;
  size_t searchlen = strlen (param->varName);

  INITSTATUS( stat, "LALReadConfigVariable", CONFIGFILEC );

  /* This traps coding errors in the calling routine. */
  ASSERT( cfgdata != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( cfgdata->lines != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( cfgdata->wasRead != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( varp != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( param->varName != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL );  
  ASSERT( param->fmt != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL );  
  ASSERT( param->strictness >= 0 && param->strictness < CONFIGFILE_LAST, stat, CONFIGFILEH_ESTRICT, CONFIGFILEH_MSGESTRICT );

  /* let's look for the variable-name in the token-list (has to at beginning of line!) */
  for (i=0; i<cfgdata->lines->nTokens; i++)
    {
      len = strcspn (cfgdata->lines->tokens[i], WHITESPACE "=:"); /* get length of variable-name */
      if (len == 0) { /* malformed token-list */
	ABORT (stat, CONFIGFILEH_ETOKENS, CONFIGFILEH_MSGETOKENS);
      }
      /* pre-select based on length of variable-name */
      if ( len != searchlen )
	continue;

      /* same len, but are they identical ? */
      if ( strncmp (param->varName, cfgdata->lines->tokens[i], len) == 0)
	{
	  found = cfgdata->lines->tokens[i] + len;
	  found += strspn (found, WHITESPACE "=:");  /* skip all whitespace and define-chars */
	  linefound = i;
	  break; /* ok, we've found it */
	}

    } /* for lines */
  

  /* FIXME: make abort or not dependent on some 'strictness' parameter! */
  if (!found)
    {
      switch (param->strictness) 
	{
	case CONFIGFILE_ERROR:
	  LALPrintError ("Error: Config-file variable %s was not found!\n", param->varName);
	  ABORT (stat, CONFIGFILEH_EVAR, CONFIGFILEH_MSGEVAR );
	  break;
	case CONFIGFILE_WARN:
	  if (lalDebugLevel & LALWARNING)
	    LALPrintError ("Warning: Config-file variable '%s' was not found!\n", param->varName);
	  RETURN (stat);
	  break;
	case CONFIGFILE_IGNORE:
	  RETURN (stat);
	  break;
	default: 
	  ABORT( stat, CONFIGFILEH_ESTRICT, CONFIGFILEH_MSGESTRICT );
	  break;
	} /* switch (strictness) */

    } /* if not found */

  /* now read the value into the variable */
  
  /* reading a quoted string needs some special treatment: */
  if ( !strcmp(param->fmt, FMT_STRING) )
    {
      /* NOTE: varp here is supposed to be a pointer to CHAR* !! */
      CHAR **cstr = (CHAR**) varp;

      ASSERT ( *cstr == NULL, stat, CONFIGFILEH_ENONULL, CONFIGFILEH_MSGENONULL);

      (*cstr) = (CHAR*) LALMalloc( strlen (found) + 1); 
      strcpy ( (*cstr), found);
      ret = 1;
    }
  else  /* but the default case is just sscanf... */
    ret = sscanf (found, param->fmt, varp);

  if ( (ret == 0) || (ret == EOF) )
    {
      LALPrintError("ERROR: Config-file variable %s was not readable using the format %s\n", param->varName, param->fmt);
      ABORT( stat, CONFIGFILEH_EFMT, CONFIGFILEH_MSGEFMT );
    }

  /* ok, we have successfully read in the config-variable: let's make a note of it */
  cfgdata->wasRead[linefound] = 1;

  RETURN (stat);

} /* LALReadConfigVariable() */

/*----------------------------------------------------------------------
 * check if all lines of config-file have been successfully read in 
 * and issue a warning or error (depending on strictness) if not
 *----------------------------------------------------------------------*/
void
LALCheckCfgReadComplete (LALStatus *stat, LALConfigData_t *cfgdata, INT2 strictness)
{
  UINT4 i;

  INITSTATUS( stat, "LALCheckCfgReadComplete", CONFIGFILEC );  

  ASSERT (cfgdata != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT (cfgdata->lines != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT (cfgdata->wasRead != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( strictness >= 0 && strictness < CONFIGFILE_LAST, stat, CONFIGFILEH_ESTRICT, CONFIGFILEH_MSGESTRICT );

  for (i=0; i < cfgdata->lines->nTokens; i++)
    if (cfgdata->wasRead[i] == 0)
      break;

  if (i != cfgdata->lines->nTokens)
    {
      switch (strictness)
	{
	case CONFIGFILE_IGNORE:
	  RETURN (stat);
	  break;
	case CONFIGFILE_WARN:
	  if (lalDebugLevel & LALWARNING)
	    {
	      LALPrintError ("Warning: config-file entry #%d has not been read!\n", i);
	      LALPrintError ("Line was: '%s'\n", cfgdata->lines->tokens[i]);
	    }
	  RETURN(stat);
	  break;
	case CONFIGFILE_ERROR:
	  LALPrintError ("ERROR: config-file entry #%d has not been read!\n", i);
	  LALPrintError ("Line was: '%s'\n", cfgdata->lines->tokens[i]);
	  ABORT (stat, CONFIGFILEH_EUNKNOWN, CONFIGFILEH_MSGEUNKNOWN);
	  break;
	default:
	  ABORT ( stat, CONFIGFILEH_ESTRICT, CONFIGFILEH_MSGESTRICT );
	  break;
	  
	} /* switch strictness */
    } /* if some line not read */
  
  RETURN (stat);

} /* LALCheckCfgReadComplete() */




/* ---------------------------------------------------------------------- 
 *   INTERNAL FUNCTIONS FOLLOW HERE
 *----------------------------------------------------------------------*/

/* ----------------------------------------------------------------------
 * cleanConfig(): do some preprocessing on the config-file, namely 'erase' 
 * all comments by '\n', and glue '\'-continued lines
 *----------------------------------------------------------------------*/
void
cleanConfig (CHARSequence *text)
{
  size_t len;  /* comment length */
  CHAR *ptr, *ptr2;

  /* clean out comments, by replacing them by '\n' */
  ptr = text->data;

  while ( (ptr += strcspn(ptr, COMMENT_CHARS)) < (text->data + text->length-1) )
    {
      len = strcspn (ptr, "\n"); 
      memset ( (void*)ptr, '\n', len); 	
    } /* while comment found */
  
  /* do line-gluing when '\' is found at end-of-line */
  ptr = text->data;
  while ( (ptr = strchr(ptr, '\\')) != NULL )
    if ( ptr[1] == '\n' )
      ptr[0] = ptr[1] = ' ';

  /* let's get rid of tabs just in case.. */
  ptr = text->data;
  while ( (ptr = strchr(ptr, '\t')) != NULL )
    *ptr = ' ';


  /* finally lets get rid of initial and trailing whitespace (we replace it by '\n') */
  ptr = text->data;

  while (ptr < (text->data + text->length -1) )
    {
      len = strspn (ptr, WHITESPACE); 
      if (len) memset ( (void*)ptr, '\n', len);
      ptr = strchr (ptr, '\n') + 1; /* point to next line-start */
      /* now clean away all trailing whitespace */
      ptr2 = ptr - 2; 
      while ( *ptr2 == ' ' )
	*ptr2-- = '\n';
    }

  return;

} /* cleanConfig() */

void
testConfigFile(void)
{
  /*--testbed --------------------------------------------------------------------*/
  FILE *fp;
  static LALConfigData_t cfgdata;
  BOOLEAN relat;
  INT2 freq;
  INT4 mermax;
  REAL4 precis;
  REAL8 relax;
  CHAR *string = NULL;
  static LALStatus status;

  fp = fopen ("settings.par", "r");
  
  LALLoadConfigFile (&status, &cfgdata, fp);

  LALReadConfigBOOLVariable (&status, &relat, &cfgdata, "relat");
  LALReadConfigINT2Variable (&status, &freq, &cfgdata, "freq_si");
  LALReadConfigINT4Variable (&status, &mermax, &cfgdata, "mer_max");
  LALReadConfigREAL4Variable (&status, &precis, &cfgdata, "precis");
  LALReadConfigREAL8Variable (&status, &relax, &cfgdata, "relax");
  LALReadConfigSTRINGVariable (&status, &string, &cfgdata, "stringy");
  
  fclose (fp);

  printf ("\nrelat = %d\n", relat);
  printf ("freq_si = %d\n", freq);
  printf ("mermax = %d\n", mermax);
  printf ("precis = %e\n", precis);
  printf ("relax = %lf\n", relax);
  printf ("stringy = '%s'\n", string);


  LALFree (string);
  REPORTSTATUS (&status);

  /*----------------------------------------------------------------------*/

  return;
}
