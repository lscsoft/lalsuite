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

#define WHITESPACE " \t"
#define COMMENT_CHARS	"#%/"


static CHAR *va(CHAR *format, ...);
static void cleanConfig (CHARSequence *text);

/*----------------------------------------------------------------------
 * parse a config-file stream into a token-list
 * 
 * gets rid of comments, empty lines, and does line-continuation of '\'
 *
 *----------------------------------------------------------------------*/
void
LALParseConfigFile (LALStatus *stat, TokenList **lines, FILE *instream)
{
  CHARSequence *rawdata = NULL;
  UINT4 i;

  INITSTATUS( stat, "LALReadConfigVariable", CONFIGFILEC );
  ATTATCHSTATUSPTR (stat);

  TRY (LALCHARReadSequence (stat->statusPtr, &rawdata, instream), stat);

  if (rawdata == NULL) {
    ABORT (stat, CONFIGFILEH_EFILE, CONFIGFILEH_MSGEFILE);
  }

  /* get rid of comments and do line-continuation */
  cleanConfig (rawdata);

  printf (rawdata->data);

  /* parse this into individual lines */
  TRY (LALCreateTokenList (stat->statusPtr, lines, rawdata->data, "\n"), stat);

  LALFree (rawdata);

  for (i=0; i < (*lines)->nTokens; i++)
    printf ( "%d: '%s'\n", i, (*lines)->tokens[i] );

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* LALParseConfigFile() */


/*----------------------------------------------------------------------
 * specialization to BOOLEAN variables
 *----------------------------------------------------------------------*/
void
LALReadConfigBOOLVariable (LALStatus *stat, BOOLEAN *varp, TokenList *lines, CHAR *varName)
{
  static LALConfigVar_t param;
  INT2 tmp;

  param.varName = varName;
  param.fmt = "%" LAL_INT2_FORMAT;

  LALReadConfigVariable (stat, (void*) &tmp, lines, &param);
  
  *varp = (BOOLEAN) tmp;

  return;
}


/*----------------------------------------------------------------------
 * specialization to INT2 variables
 *----------------------------------------------------------------------*/
void
LALReadConfigINT2Variable (LALStatus *stat, INT2 *varp, TokenList *lines, CHAR *varName)
{
  static LALConfigVar_t param;

  param.varName = varName;
  param.fmt = "%" LAL_INT2_FORMAT;

  LALReadConfigVariable (stat, (void*) varp, lines, &param);
  
  *varp = (INT2) tmp;

  return;
}

/*----------------------------------------------------------------------
 * specialization to INT4 variables
 *----------------------------------------------------------------------*/
void
LALReadConfigINT4Variable (LALStatus *stat, INT2 *varp, TokenList *lines, CHAR *varName)
{
  static LALConfigVar_t param;

  param.varName = varName;
  param.fmt = "%" LAL_INT2_FORMAT;

  LALReadConfigVariable (stat, (void*) varp, lines, &param);
  
  *varp = (INT2) tmp;

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
LALReadConfigVariable (LALStatus *stat, void *varp, TokenList *lines, LALConfigVar_t *param)
{ /* </lalVerbatim> */
  CHAR *found = NULL;
  INT2 ret = 0;

  UINT4 i;
  size_t len;
  size_t searchlen = strlen (param->varName);

  INITSTATUS( stat, "LALReadConfigVariable", CONFIGFILEC );

  /* This traps coding errors in the calling routine. */
  ASSERT( lines != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( varp != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( param->varName != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL );  
  ASSERT( param->fmt != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL );  

  /* let's look for the variable-name in the token-list (has to at beginning of line!) */
  for (i=0; i<lines->nTokens; i++)
    {
      len = strcspn (lines->tokens[i], WHITESPACE "=:"); /* get length of variable-name */
      if (len == 0) { /* malformed token-list */
	ABORT (stat, CONFIGFILEH_ETOKENS, CONFIGFILEH_MSGETOKENS);
      }
      /* pre-select based on length of variable-name */
      if ( len != searchlen )
	continue;

      /* same len, but are they identical ? */
      if ( strncmp (param->varName, lines->tokens[i], len) == 0)
	{
	  found = lines->tokens[i] + len;
	  found += strspn (found, WHITESPACE "=:");  /* skip all whitespace and define-chars */
	  break; /* ok, we've found it */
	}

    } /* for lines */
  

  if (!found)
    {
      LALWarning (stat, param->varName);
      ABORT( stat, CONFIGFILEH_EVAR, CONFIGFILEH_MSGEVAR );
    }

  /* now read the value into the variable */
  
  /* reading a quoted string needs some special treatment: */
  if ( !strcmp(param->fmt, FMT_STRING) )
    {
      /* NOTE: varp here is supposed to be a pointer to CHAR* !! */
      CHAR **cstr = (CHAR**) varp;
      (*cstr) = (CHAR*) LALMalloc( strlen (found) + 1); 
      strcpy ( (*cstr), found);
      ret = 1;
    }
  else  /* but the default case is just sscanf... */
    ret = sscanf (found, param->fmt, varp);

  if ( (ret == 0) || (ret == EOF) )
    {
      LALError( stat, va("Variable %s was not readable using the format %s\n", param->varName, param->fmt) );
      ABORT( stat, CONFIGFILEH_EFMT, CONFIGFILEH_MSGEFMT );
    }

  RETURN (stat);

} /* LALReadConfigVariable() */

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
  CHAR *ptr;

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


  /* finally lets get rid of initial whitespace (we replace it by '\n') */
  ptr = text->data;

  while (ptr < (text->data + text->length -1) )
    {
      len = strspn (ptr, WHITESPACE); 
      if (len) memset ( (void*)ptr, '\n', len);
      ptr = strchr (ptr, '\n') + 1; /* point to next line-start */
    }

  return;

} /* cleanConfig() */

/*============
va()

does a varargs printf into a temp buffer, so I don't need to have
varargs versions of all text functions.
FIXME: make this buffer size safe someday
============*/
char	*va(char *format, ...)
{
	va_list		argptr;
	static char		string[1024];
	
	va_start (argptr, format);
	LALVsnprintf (string, 1024, format,argptr);
	va_end (argptr);

	string[1023] = '\0';

	return (string);	
} /* va() */
