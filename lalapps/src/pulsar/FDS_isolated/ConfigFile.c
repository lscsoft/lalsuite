/************************************ <lalVerbatim file="ConfigFileCV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{ConfigFile.c}}
\label{ss:ConfigFile.c}

Some general-purpose routines for config-file reading

\subsubsection*{Prototypes}
\input{ConfigFileCP}
\idx{LALLoadConfigFile()}
\idx{LALDestroyConfigData()}
\idx{LALReadConfigVariable()}
\idx{LALReadConfigBOOLVariable()}
\idx{LALReadConfigINT2Variable()}
\idx{LALReadConfigINT4Variable()}
\idx{LALReadConfigREAL4Variable()}
\idx{LALReadConfigREAL8Variable()}
\idx{LALReadConfigSTRINGVariable()}
\idx{LALReadConfigSTRINGNVariable()}
\idx{LALCheckConfigReadComplete()}

\subsubsection*{Description}

allow reading of "variable = value" type config-files, which may contain 
comments and line-continuation.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\verb+LALReadConfigSTRINGVariable()+ is a special case, as it reads
the rest of the logical line (excluding comments) as a string (different
from \verb+"%s"+). The required memory is allocated and has to be
freed by the caller. 

\verb+LALReadConfigSTRINGVariable()+ is a somewhat different version
of the same, but reads only up to a maximum of $N$ bytes, and the
corresponding memory has to be provided by the user.
This is achieved by using a \verb+CHARVector+, where the
\verb+length+-entry specifies the maximal length of the string.

\vfill{\footnotesize\input{ConfigFileCV}}

******************************************************* </lalLaTeX> */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <getopt.h>

#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALStdio.h>
#include <lal/StreamInput.h>

#include "ConfigFile.h"

NRCSID( CONFIGFILEC, "$Id$" );

extern INT4 lalDebugLevel;

/* needed for command-line parsing */
extern char *optarg;
extern int optind, opterr, optopt;

#define ERR -1
#define OK  0

#define FMT_STRING "string"    /* reading in quoted strings needs some special treatment */
#define WHITESPACE " \t"

/* local prototypes */
static void cleanConfig (CHARSequence *text);


/*----------------------------------------------------------------------
 * parse a config-file stream into a token-list
 * 
 * gets rid of comments, empty lines, and does line-continuation of '\'
 *
 *----------------------------------------------------------------------*/
void
LALLoadConfigFile (LALStatus *stat, LALConfigData **cfgdata, const CHAR *fname)
{
  CHARSequence *rawdata = NULL;
  UINT4 i;
  FILE *fp;

  INITSTATUS( stat, "LALLoadConfigFile", CONFIGFILEC );
  ATTATCHSTATUSPTR (stat);

  ASSERT (*cfgdata == NULL, stat, CONFIGFILEH_ENONULL, CONFIGFILEH_MSGENONULL);
  ASSERT (fname != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);


  if ( (fp = fopen(fname, "r")) == NULL) {
    ABORT (stat, CONFIGFILEH_EFILE, CONFIGFILEH_MSGEFILE);
  }

  LALCHARReadSequence (stat->statusPtr, &rawdata, fp);
  fclose (fp);
  CHECKSTATUSPTR (stat);

  if (rawdata == NULL) {
    ABORT (stat, CONFIGFILEH_EFILE, CONFIGFILEH_MSGEFILE);
  }

  /* get rid of comments and do line-continuation */
  cleanConfig (rawdata);

  if ( (*cfgdata = LALCalloc (1, sizeof(LALConfigData))) == NULL) {
    ABORT (stat, CONFIGFILEH_EMEM, CONFIGFILEH_MSGEMEM);
  }
  
  /* parse this into individual lines */
  LALCreateTokenList (stat->statusPtr, &((*cfgdata)->lines), rawdata->data, "\n");
  LALFree (rawdata->data);
  LALFree (rawdata);
  
  BEGINFAIL (stat)
    LALFree (*cfgdata);
  ENDFAIL (stat);

  if (lalDebugLevel >= 3)
    {
      LALPrintError ("ConfigFile DEBUG: parsed config-file contents:\n");
      for (i=0; i < (*cfgdata)->lines->nTokens; i++)
	printf ( "%d: '%s'\n", i, (*cfgdata)->lines->tokens[i] );
    }


  /* initialize the 'wasRead' flags for the lines */
  if ( ((*cfgdata)->wasRead = LALCalloc (1, (*cfgdata)->lines->nTokens * sizeof( (*cfgdata)->wasRead[0]))) == NULL) {
    LALFree ((*cfgdata)->lines);
    ABORT (stat, CONFIGFILEH_EMEM, CONFIGFILEH_MSGEMEM);
  }

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* LALLoadConfigFile() */

/*----------------------------------------------------------------------
 * free memory associated with a LALConfigData structure
 *----------------------------------------------------------------------*/
void
LALDestroyConfigData (LALStatus *stat, LALConfigData **cfgdata)
{
  INITSTATUS( stat, "LALDestroyConfigData", CONFIGFILEC );
  ATTATCHSTATUSPTR (stat);

  ASSERT (cfgdata != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT (*cfgdata != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT ((*cfgdata)->lines != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT ( (*cfgdata)->wasRead != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);

  TRY ( LALDestroyTokenList (stat->statusPtr, &((*cfgdata)->lines)), stat);
  LALFree ( (*cfgdata)->wasRead);
  LALFree ( *cfgdata );
  
  *cfgdata = NULL;

  DETATCHSTATUSPTR (stat);
  RETURN (stat);
} /* LALDestroyConfigData() */

/*----------------------------------------------------------------------
 * specialization to BOOLEAN variables
 *----------------------------------------------------------------------*/
void
LALReadConfigBOOLVariable (LALStatus *stat, BOOLEAN *varp, LALConfigData *cfgdata, const CHAR *varName)
{
  static LALConfigVar param;
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
LALReadConfigINT2Variable (LALStatus *stat, INT2 *varp, LALConfigData *cfgdata, const CHAR *varName)
{
  static LALConfigVar param;

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
LALReadConfigINT4Variable (LALStatus *stat, INT4 *varp, LALConfigData *cfgdata, const CHAR *varName)
{
  static LALConfigVar param;

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
LALReadConfigREAL4Variable (LALStatus *stat, REAL4 *varp, LALConfigData *cfgdata, const CHAR *varName)
{
  static LALConfigVar param;

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
LALReadConfigREAL8Variable (LALStatus *stat, REAL8 *varp, LALConfigData *cfgdata, const CHAR *varName)
{
  static LALConfigVar param;

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
LALReadConfigSTRINGVariable (LALStatus *stat, CHAR **varp, LALConfigData *cfgdata, const CHAR *varName)
{
  static LALConfigVar param;

  param.varName = varName;
  param.fmt = FMT_STRING;
  param.strictness = CONFIGFILE_WARN;

  LALReadConfigVariable (stat, (void*) varp, cfgdata, &param);
  
  return;
}

/*----------------------------------------------------------------------
 * READING OF FIXED LENGTH STRINGS:
 * another variant of string-reading:similar to ReadConfigSTRING, but
 * here a fixed-size CHAR-array is used as input, no memory is allocated
 * NOTE: you have to provide the length of your string-array as input!
 *      in varp->length
 *
 * (this is basically a wrapper for ReadConfigSTRINGVariable())
 *
 * NOTE2: the behaviour is similar to strncpy, i.e. we silently clip the
 *       string to the right length, BUT we also 0-terminate it properly.
 *       No error or warning is generated when clipping occurs!
 *
 * NOTE3: at return, the value varp->length is set to the length of the
 *        string copied
 *
 *----------------------------------------------------------------------*/
void
LALReadConfigSTRINGNVariable (LALStatus *stat, CHARVector *varp, LALConfigData *cfgdata, const CHAR *varName)
{
  CHAR *tmp = NULL;

  INITSTATUS( stat, "LALReadConfigVariable", CONFIGFILEC );
  ATTATCHSTATUSPTR (stat);
  
  /* This traps coding errors in the calling routine. */
  ASSERT( varp != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( varp->data != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( varp->length != 0, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);

  TRY (LALReadConfigSTRINGVariable (stat->statusPtr, &tmp, cfgdata, varName), stat);

  if (tmp != NULL)
    {
      strncpy (varp->data, tmp, varp->length - 1);
      varp->data[varp->length-1] = '\0';
      LALFree (tmp);
      varp->length = strlen (varp->data);
    }
    
  DETATCHSTATUSPTR (stat);
  RETURN (stat);  
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
LALReadConfigVariable (LALStatus *stat, void *varp, LALConfigData *cfgdata, LALConfigVar *param)
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
LALCheckConfigReadComplete (LALStatus *stat, LALConfigData *cfgdata, INT4 strictness)
{
  UINT4 i;

  INITSTATUS( stat, "LALCheckConfigReadComplete", CONFIGFILEC );  

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

} /* LALCheckConfigReadComplete() */


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
  CHAR *ptr, *ptr2, *eol;
  BOOLEAN inQuotes = 0;

  /* clean out comments, by replacing them by '\n' */
  ptr = text->data;

  while ( *ptr )
    {
      if ( (*ptr) == '\"' )
	inQuotes = !inQuotes;

      if ( ((*ptr) == '#') || ( (*ptr) == ';') )
	if ( !inQuotes )	/* only consider as comments if not quoted */
	  {
	    len = strcspn (ptr, "\n"); 
	    memset ( (void*)ptr, '\n', len); 	
	  }
	
      ptr ++;

    } /* while *ptr */

  
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
      eol = strchr (ptr, '\n'); /* point to end-of-line */
      if (eol != NULL)
	ptr = eol;
      else
	ptr = strchr (ptr, '\0'); /* or end of file */

      /* clean away all trailing whitespace of last line*/
      ptr2 = ptr - 1; 
      while ( *ptr2 == ' ' )
	*ptr2-- = '\n';

      /* step to next line */
      ptr += 1;
    }

  return;

} /* cleanConfig() */

void
testConfigFile(void)
{
  /*--testbed --------------------------------------------------------------------*/
  static LALConfigData *cfgdata;
  BOOLEAN relat;
  INT2 freq;
  INT4 mermax;
  REAL4 precis;
  REAL8 relax;
  CHAR *string = NULL;
  static LALStatus status;

  LALLoadConfigFile (&status, &cfgdata, "settings.par");

  LALReadConfigBOOLVariable (&status, &relat, cfgdata, "relat");
  LALReadConfigINT2Variable (&status, &freq, cfgdata, "freq_si");
  LALReadConfigINT4Variable (&status, &mermax, cfgdata, "mer_max");
  LALReadConfigREAL4Variable (&status, &precis, cfgdata, "precis");
  LALReadConfigREAL8Variable (&status, &relax, cfgdata, "relax");
  LALReadConfigSTRINGVariable (&status, &string, cfgdata, "stringy");
  
  LALDestroyConfigData (&status, &cfgdata);

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


/*----------------------------------------------------------------------
 * parse command-line into UserVariable array
 *----------------------------------------------------------------------*/
void
ReadCmdlineInput (LALStatus *stat, int argc, char *argv[], UserVariable *uvars)
{
  INT4 c;
  UINT4 index, numvars;
  UINT4 i, pos;
  char optstring[512] = "\0";	/* should be easily enough */
  struct option *long_options;

  INITSTATUS( stat, "ReadCmdlineInput", CONFIGFILEC );

  ASSERT (argv, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT (uvars, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);

  /* count number of user-variables */
  numvars = 0;
  while (uvars[numvars].name != NULL)
    numvars ++;

  /* build optstring */
  for (i=0; i< numvars; i++)
    {
      pos = strlen(optstring);
      optstring[pos++] = uvars[i].optchar;
      if (uvars[i].type != UVAR_BOOL)	/* everything but bool takes an argument */
	optstring[pos++] = ':';
      optstring[pos] = '\0';
    } /* for i < numvars */

  long_options = LALMalloc ( (numvars+1) * sizeof(struct option));
  
  for (i=0; i < numvars; i++)
    {
      long_options[i].name = uvars[i].name;
      long_options[i].has_arg = (uvars[i].type == UVAR_BOOL) ? no_argument : required_argument;
      long_options[i].flag = 0;
      long_options[i].val = uvars[i].optchar;
    }
  long_options[numvars].name = 0;
  long_options[numvars].has_arg = 0;
  long_options[numvars].flag = 0;
  long_options[numvars].val = 0;

  while ( (c = getopt_long(argc, argv, optstring, long_options, NULL)) != -1 )
    {
      index = 0;
      while ( long_options[index].name && ( c != long_options[index].val) )
	index ++;
      if (index >= numvars) {
	ABORT (stat, CONFIGFILEH_EOPT, CONFIGFILEH_MSGEOPT);
      }

      switch (uvars[index].type)
	{
	case UVAR_BOOL:
	  *(BOOLEAN*)(uvars[index].varp)  = 1;
	  break;
	case UVAR_INT4:
	  *(INT4*)(uvars[index].varp) = (INT4) atoi (optarg);
	  break;
	case UVAR_REAL8:
	  *(REAL8*)(uvars[index].varp) = (REAL8) atof (optarg);
	  break;
	case UVAR_CHAR:
	  if (!optarg) break;
	  if ( *(CHAR**)uvars[index].varp != NULL)	 /* something allocated here before? */
	    LALFree ( *(CHAR**)uvars[index].varp );

	  *(CHAR**)(uvars[index].varp) = LALMalloc (strlen(optarg) + 1);
	  strcpy ( *(CHAR**)uvars[index].varp, optarg);
	  break;
	} /* switch index */

    } /* while getopt() */

  LALFree (long_options);

  RETURN (stat);

} /* ReadCmdlineInput() */

/*----------------------------------------------------------------------
 * Read config-variables from cfgfile and parse into input-structure
 *
 * an error is reported if the config-file reading fails, but the 
 * individual variable-reads are treated as optional
 *----------------------------------------------------------------------*/
void
ReadCfgfileInput (LALStatus *stat, const CHAR *cfgfile, UserVariable *uvars)
{
  LALConfigData *cfg = NULL;
  UINT4 i, numvars;
  CHAR *stringbuf;
		  

  INITSTATUS( stat, "ReadCfgFile", CONFIGFILEC );
  ATTATCHSTATUSPTR (stat);

  TRY (LALLoadConfigFile (stat->statusPtr, &cfg, cfgfile), stat);


  /* count number of user-variables */
  numvars = 0;
  while (uvars[numvars].name != NULL)
    numvars ++;

  for (i=0; i< numvars; i++)
    {
      switch (uvars[i].type)
	{
	case UVAR_BOOL:
	  TRY( LALReadConfigBOOLVariable (stat->statusPtr, uvars[i].varp, cfg, uvars[i].name ), stat);
	  break;
	case UVAR_INT4:
	  TRY( LALReadConfigINT4Variable (stat->statusPtr, uvars[i].varp, cfg, uvars[i].name ), stat);
	  break;
	case UVAR_REAL8:
	  TRY( LALReadConfigREAL8Variable (stat->statusPtr, uvars[i].varp, cfg, uvars[i].name ), stat);
	  break;
	case UVAR_CHAR:
	  stringbuf = NULL;
	  TRY( LALReadConfigSTRINGVariable (stat->statusPtr, &stringbuf, cfg, uvars[i].name ), stat);
	  if (stringbuf)
	    {
	      if ( *(CHAR**)uvars[i].varp != NULL)	 /* something allocated here before? */
		LALFree ( *(CHAR**)uvars[i].varp );
	      *(CHAR**)uvars[i].varp = stringbuf;
	    }
	  break;
	} /* switch type */

    } /* for i < numvars */


  /* ok, that should be it: check if there were more definitions we did not read */
  TRY (LALCheckConfigReadComplete (stat->statusPtr, cfg, CONFIGFILE_ERROR), stat);	

  TRY( LALDestroyConfigData (stat->statusPtr, &cfg), stat);


  DETATCHSTATUSPTR(stat);
  RETURN (stat);

} /* ReadCfgfileInput() */

/*--------------------------------------------------------------------------------
 * free memory associated with user-variable array (basically free all the strings 
 *--------------------------------------------------------------------------------*/
void 
FreeUserVars (LALStatus *stat, UserVariable *uvars)
{
  UINT4 i;

  INITSTATUS( stat, "FreeUserVars", CONFIGFILEC );

  ASSERT (uvars != NULL,  stat, CONFIGFILEH_ENULL,  CONFIGFILEH_MSGENULL);

  i = 0;
  while ( uvars[i].name != NULL )
    {
      if ( (uvars[i].type == UVAR_CHAR) && (uvars[i].varp != NULL) && (*(CHAR**)(uvars[i].varp) != NULL) )
	{
	  LALFree ( *((CHAR**)uvars[i].varp) );
	  *(CHAR**)uvars[i].varp = NULL; 	/* important */
	}
      
      i++;
    } /* while uvars != NULL */
  
  RETURN(stat);
} /* FreeUserVars() */

/*----------------------------------------------------------------------
 * assemble all help-info from uvars into a help-string
 *----------------------------------------------------------------------*/
void
GetUvarHelpString (LALStatus *stat, CHAR **helpstring, UserVariable *uvars)
{
  UINT4 numvars, i, mem;
  CHAR strbuf[512];	/* should be enough for one line...*/
  CHAR defaultstr[100];
  /* we need strings for UVAR_BOOL, UVAR_INT4, UVAR_REAL8, UVAR_CHAR: */
  const CHAR *typestr[] = {"BOOL", "INT", "FLOAT", "STRING"}; 

  INITSTATUS (stat, "GetUvarHelpString", CONFIGFILEC);

  ASSERT (helpstring != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT ( *helpstring == NULL, stat, CONFIGFILEH_ENONULL, CONFIGFILEH_MSGENONULL);

  /* count number of user-variables  and estimate lower-limit on memory requirements*/
  numvars = 0;
  mem = 0;
  while (uvars[numvars].name != NULL)
    {
      mem += strlen (uvars[numvars].name) + 10;
      if (uvars[numvars].help)
	mem += strlen (uvars[numvars].help) + 10;

      numvars ++;
    }/* while uvars[numvars] */
  
  if ( (*helpstring = LALCalloc (1, mem)) == NULL ) {
    ABORT (stat,  CONFIGFILEH_EMEM,  CONFIGFILEH_MSGEMEM);
  }
  
  for (i=0; i < numvars; i++)
    {
      switch (uvars[i].type)
	{
	case UVAR_BOOL:
	  sprintf (defaultstr, *(BOOLEAN*)uvars[i].varp ? "True" : "False");
	  break;
	case UVAR_INT4:
	  sprintf (defaultstr, "%d", *(INT4*)uvars[i].varp );
	  break;
	case UVAR_REAL8:
	  sprintf (defaultstr, "%.3g", *(REAL8*)uvars[i].varp);
	  break;
	case UVAR_CHAR:
	  if ( *(CHAR**)uvars[i].varp )
	    strcpy (defaultstr, *(CHAR**)uvars[i].varp );
	  else
	    defaultstr[0] = '\0';
	  break;
	} /* switch type */

      sprintf (strbuf, "   --%-14s (-%c) \t%s\t%s  (Default: %s)\n", 
	       uvars[i].name, uvars[i].optchar, typestr[uvars[i].type], uvars[i].help, defaultstr);

      strcat (*helpstring, strbuf);

    } /* for i < numvars */


  RETURN(stat);
} /* GetUvarHelpString() */


/*----------------------------------------------------------------------
 * This function puts all the above pieces together, and basically does
 * everything: get config-filename from cmd-line (if found),
 * then interpret config-file and then the command-line
 *----------------------------------------------------------------------*/
void
ReadUserInput (LALStatus *stat, int argc, char *argv[], UserVariable *uvars)
{
  INT2 i;
#define MAXFILENAMELENGTH 512
  CHAR fname[MAXFILENAMELENGTH] = "";

  INITSTATUS( stat, "ReadUserInput", CONFIGFILEC);
  ATTATCHSTATUSPTR (stat);

  /*----------------------------------------------------------------------
   * pre-process command-line: we want config-file and debug-level 
   */
  for (i=1; i < argc; i++)
    {
      /* was a config-file specified ? */
      if ( argv[i][0] == '@' )
	{
	  strncpy (fname, argv[i]+1, MAXFILENAMELENGTH);
	  fname[MAXFILENAMELENGTH-1] = '\0'; /* terminate in case it was truncated */
	}

    } /* for cmd-line */
  /*----------------------------------------------------------------------*/

  /* if config-file specified, read from that first */
  if (fname[0] != '\0') {
    TRY (ReadCfgfileInput (stat->statusPtr, fname, uvars ), stat);
  }

  /* now do proper cmdline parsing: overloads config-file settings */
  TRY (ReadCmdlineInput (stat->statusPtr, argc, argv, uvars), stat);

  DETATCHSTATUSPTR (stat);
  RETURN (stat);
  
} /* ReadUserInput() */

