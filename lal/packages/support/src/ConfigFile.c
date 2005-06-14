/*
 * Copyright (C) 2004, 2005 Reinhard Prix
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

/************************************ <lalVerbatim file="ConfigFileCV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{ConfigFile.c}}
\label{ss:ConfigFile.c}

Some general-purpose routines for config-file reading

\subsubsection*{Prototypes}
\idx{LALParseDataFile()}
\idx{LALDestroyParsedDataFile()}
\idx{LALReadConfigVariable()}
\idx{LALReadConfigBOOLVariable()}
\idx{LALReadConfigINT2Variable()}
\idx{LALReadConfigINT4Variable()}
\idx{LALReadConfigREAL4Variable()}
\idx{LALReadConfigREAL8Variable()}
\idx{LALReadConfigSTRINGVariable()}
\idx{LALReadConfigSTRINGNVariable()}
\idx{LALCheckConfigReadComplete()}

\input{ConfigFileCP}

\subsubsection*{Description}

This module provides routines for reading formatted
config-files containing definitions of the form \mbox{\texttt{variable = value}}.
The general syntax is somewhat similar to the one provided by the
perl-module \texttt{ConfigParser} (cf. 
\verb+http://www.python.org/doc/current/lib/module-ConfigParser.html+ )
but (currently) without the possibility of "chapters".
Comments are allowed using either '\#' or ';'. You can also use
standard line-continuation  using a '\verb+\+' at the end of the line.
Also note that '\#' or ';' within double-quotes '\"' are \emph{not}
treated as comment-characters.  The general syntax is best illustrated
using a simple example: 
\begin{verbatim}
# comment line
var1 = 1.0    ; you can also comment using semi-colons
somevar = some text.\
        You can also use\
        line-continuation	
   var3 = 4      # whatever that means
note = "this is also possible, and # here does nothing"
a_switch = true	 #possible values: 0,1,true,false,yes,no, case insensitive
# etc etc.
\end{verbatim}

Note that TABS generally get replaced by a single space, which can be
useful in the case of line-continuation (see example). All leading and
trailing spaces in are ignore (except within double-quotes).

The general approach of reading from such a config-file, is to first
call\\
\verb+LALLoadConfigFile(stat, LALConfigData *cfg, FILE *fp)+,
which loads and pre-parses the contents of the config-file into the
structure \verb+LALConfigData+. Then one can then read in
config-variables either using one of the custom-wrappers:\\
\verb+LALReadConfig<TYPE>Variable(stat, <TYPE> *cvar, LALConfigData *cfg, CHAR *varname)+
or the general-purpose reading function:\\
\verb+LALReadConfigVariable(stat, void *cvar, LALConfigData *cfg, LALConfigVar *var)+


A boolean variable read by \verb+LALReadConfigBOOLVariable()+ can have any of the values 
\verb+{1, 0, yes, no, true, false}+, where the comparison is done \emph{case-insensitively}, 
i.e. you also use "True" or "FALSE"....


If one wishes a ``tight'' sytnax for the config-file, one can check
that there are no "illegal" entries in the config-file. This is done
by checking at the end that all config-file entries have been
successfully parsed, using: \\
\verb+LALCheckConfigReadComplete (stat, LALConfigData *cfg, INT2 strictness)+,
where \verb+strictness+ is either \verb+CONFIGFILE_WARN+ or \verb+CONFIGFILE_ERROR+. 
In the first case only a warning is issued, while in the second it is
treated as a LAL-error if some config-file entries have not been
read-in. (The use of this function is optional).


The configfile-data should be freed at the end using\\
\verb+void LALDestroyConfigData (LALStatus *stat, LALConfigData *cfg)+.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALCHARReadSequence()
LALCreateTokenList()       LALDestroyTokenList()
LALCalloc()                LALMalloc()             LALFree()  
LALPrintError()            LALOpenDataFile()                 fclose()
\end{verbatim}

\subsubsection*{Notes}

\verb+LALReadConfigSTRINGVariable()+ and
\verb+LALReadConfigSTRINGVariable()+ are not the same as using
\verb+"%s"+ as a format string, as they read the \emph{rest} of the
logical line (excluding comments) as a string.


In the case of \verb+LALReadConfigSTRINGVariable()+, the required
memory is allocated and has to be freed by the caller, while for   
\verb+LALReadConfigSTRINGVariable()+ the caller has to provide a 
\verb+CHARVector+ of length $N$, which defines the maximum length of
string to be read.


\textbf{Note:} instead of using these functions directly, it might be
more convenient to use the \verb+UserInput+ infrastructure
(cf.~\ref{s:UserInput.h}).

\vfill{\footnotesize\input{ConfigFileCV}}

******************************************************* </lalLaTeX> */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
/* #include <ctype.h> */  /* don't use this, as it binds us to GLIBC_2.3 symbols!! */

#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/StreamInput.h>

#include <lal/ConfigFile.h>

NRCSID( CONFIGFILEC, "$Id$");

extern INT4 lalDebugLevel;

#define FMT_STRING "string"    /* reading in quoted strings needs some special treatment */
#define WHITESPACE " \t"

#define TRUE   (1==1)
#define FALSE  (1==0)

/* local prototypes */
static void cleanConfig (CHARSequence *text);
CHAR my_tolower (CHAR in);
/* ctype replacements w/o locale */
static int TOLOWER(int c);
static const char upper_chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
static const char lower_chars[] = "abcdefghijklmnopqrstuvwxyz";


/* static const char digit_chars[] = "0123456789"; */

#define ISDIGIT(c) \
    ((c) && strchr(digit_chars, (c)) != NULL)

/* static const char space_chars[] = " \f\n\r\t\v"; */

#define ISSPACE(c) \
    ((c) && strchr(space_chars, (c)) != NULL)



/*----------------------------------------------------------------------
 * parse an ASCII data-file into a pre-cleaned array of lines
 * 
 * the cleaning gets rid of comments ('#', ';'), empty lines, 
 * and performs line-continuation if '\' is found at EOL
 *
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="ConfigFileCP"> */
void
LALParseDataFile (LALStatus *stat, 
		  LALParsedDataFile **cfgdata, 	/* output: pre-parsed data-file lines */
		  const CHAR *fname)		/* name of config-file to be read */
{ /* </lalVerbatim> */

  CHARSequence *rawdata = NULL;
  FILE *fp;

  INITSTATUS( stat, "LALParseDataFile", CONFIGFILEC );
  ATTATCHSTATUSPTR (stat);

  ASSERT (*cfgdata == NULL, stat, CONFIGFILEH_ENONULL, CONFIGFILEH_MSGENONULL);
  ASSERT (fname != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);

  if ( (fp = LALOpenDataFile (fname)) == NULL) {
    LALPrintError ("Could not open data-file: `%s`\n\n", fname);
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

  if ( (*cfgdata = LALCalloc (1, sizeof(LALParsedDataFile))) == NULL) {
    ABORT (stat, CONFIGFILEH_EMEM, CONFIGFILEH_MSGEMEM);
  }
  
  /* parse this into individual lines */
  LALCreateTokenList (stat->statusPtr, &((*cfgdata)->lines), rawdata->data, "\n");
  LALFree (rawdata->data);
  LALFree (rawdata);
  
  BEGINFAIL (stat)
    LALFree (*cfgdata);
  ENDFAIL (stat);

  /* initialize the 'wasRead' flags for the lines */
  if ( ((*cfgdata)->wasRead = LALCalloc (1, (*cfgdata)->lines->nTokens * sizeof( (*cfgdata)->wasRead[0]))) == NULL) {
    LALFree ((*cfgdata)->lines);
    ABORT (stat, CONFIGFILEH_EMEM, CONFIGFILEH_MSGEMEM);
  }

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* LALLoadConfigFile() */

/*----------------------------------------------------------------------
 * free memory associated with a LALParsedDataFile structure
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="ConfigFileCP"> */
void
LALDestroyParsedDataFile (LALStatus *stat, LALParsedDataFile **cfgdata)
{ /* </lalVerbatim> */
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
LALReadConfigVariable (LALStatus *stat, 
		       void *varp, 			/* output: result gets written here! */
		       const LALParsedDataFile *cfgdata, 	/* input: pre-parsed config-data */
		       const LALConfigVar *param,	/* var-name, fmt-string, strictness */
		       BOOLEAN *wasRead)		/* output: did we succeed in reading? */
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

  *wasRead = FALSE;

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
  
  if (!found)
    {
      switch (param->strictness) 
	{
	case CONFIGFILE_IGNORE:
	  RETURN (stat);
	  break;
	case CONFIGFILE_WARN:
	  if (lalDebugLevel & LALWARNING)
	    LALPrintError ("\nWarning: Config-file variable '%s' was not found!\n", param->varName);
	  RETURN (stat);
	  break;

	case CONFIGFILE_ERROR:
	default: 
	  LALPrintError ("\nError: Config-file variable %s was not found!\n", param->varName);
	  ABORT (stat, CONFIGFILEH_EVAR, CONFIGFILEH_MSGEVAR );
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
      LALPrintError("\nERROR: Config-file variable %s was not readable using the format %s\n\n", param->varName, param->fmt);
      ABORT( stat, CONFIGFILEH_EFMT, CONFIGFILEH_MSGEFMT );
    }

  /* ok, we have successfully read in the config-variable: let's make a note of it */
  cfgdata->wasRead[linefound] = 1;
  
  *wasRead = TRUE;

  RETURN (stat);

} /* LALReadConfigVariable() */



/*----------------------------------------------------------------------
 * specialization to BOOLEAN variables
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="ConfigFileCP"> */
void
LALReadConfigBOOLVariable (LALStatus *stat, 
			   BOOLEAN *varp, 		 /* output: variable to store result */
			   const LALParsedDataFile *cfgdata, /* input: pre-parsed config-data */
			   const CHAR *varName,		 /* input: variable-name to read */
			   BOOLEAN *wasRead)		 /* output: did we succeed in reading? */
{ /* </lalVerbatim> */

  CHAR *tmp = NULL;
  INT2 ret = -1;	/* -1 means no legal value has been parsed */

  INITSTATUS( stat, "LALReadConfigBOOLVariable", CONFIGFILEC );
  ATTATCHSTATUSPTR (stat);

  *wasRead = FALSE;

  /* first read the value as a string */
  TRY (LALReadConfigSTRINGVariable (stat->statusPtr, &tmp, cfgdata, varName, wasRead), stat);

  if (*wasRead && tmp) /* if we read anything at all... */
    {
      /* get rid of case ambiguities */
      TRY (LALLowerCaseString (stat->statusPtr, tmp), stat);

      /* try to parse it as a bool */
      if (      !strcmp(tmp, "yes") || !strcmp(tmp, "true") || !strcmp(tmp,"1"))
	ret = 1;
      else if ( !strcmp (tmp, "no") || !strcmp(tmp,"false") || !strcmp(tmp,"0"))
	ret = 0;
      else
	{
	  LALPrintError ( "illegal bool-value `%s`\n", tmp);
	  LALFree (tmp);
	  ABORT (stat, CONFIGFILEH_EBOOL, CONFIGFILEH_MSGEBOOL);
	}
      LALFree (tmp);
      
      if (ret != -1)	/* only set value of something has been found */
	{
	  *varp = (BOOLEAN) ret;
	  *wasRead = TRUE;
	}

    } /* if wasRead && tmp */

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* LALReadConfigBOOLVariable() */

/*----------------------------------------------------------------------
 * specialization to INT4 variables
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="ConfigFileCP"> */
void
LALReadConfigINT4Variable (LALStatus *stat, 
			   INT4 *varp, 
			   const LALParsedDataFile *cfgdata, 
			   const CHAR *varName,
			   BOOLEAN *wasRead)
{ /* </lalVerbatim> */
  LALConfigVar param = {0,0,0};

  INITSTATUS( stat, "LALReadConfigINT4Variable", CONFIGFILEC );

  param.varName = varName;
  param.fmt = "%" LAL_INT4_FORMAT;
  param.strictness = CONFIGFILE_IGNORE;

  LALReadConfigVariable (stat, (void*) varp, cfgdata, &param, wasRead);
  
  RETURN (stat);

} /* LALReadConfigINT4Variable() */

/*----------------------------------------------------------------------
 * specialization to REAL8 variables
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="ConfigFileCP"> */
void
LALReadConfigREAL8Variable (LALStatus *stat, 
			    REAL8 *varp, 
			    const LALParsedDataFile *cfgdata, 
			    const CHAR *varName,
			    BOOLEAN *wasRead)
{ /* </lalVerbatim> */
  LALConfigVar param = {0,0,0};

  INITSTATUS( stat, "LALReadConfigREAL8Variable", CONFIGFILEC );

  param.varName = varName;
  param.fmt = "%" LAL_REAL8_FORMAT;
  param.strictness = CONFIGFILE_IGNORE;

  LALReadConfigVariable (stat, (void*) varp, cfgdata, &param, wasRead);
  
  RETURN (stat);

} /* LALReadConfigREAL8Variable() */

/*----------------------------------------------------------------------
 * specialization to STRING variables 
 * NOTE: this means the rest of the line, NOT "%s"! (but excluding comments of course), 
 * NOTE2: if string is quoted by \", everything within quotes is read,
 *       and the quotes are removed here 
 *
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="ConfigFileCP"> */
void
LALReadConfigSTRINGVariable (LALStatus *stat, 
			     CHAR **varp, 		/* output: string, allocated here! */
			     const LALParsedDataFile *cfgdata, /* pre-parsed config-data */
			     const CHAR *varName,	/* variable-name to be read */
			     BOOLEAN *wasRead)			     
{ /* </lalVerbatim> */
  LALConfigVar param = {0,0,0};
  CHAR *str = NULL;
  CHAR *ret = NULL;

  INITSTATUS( stat, "LALReadConfigSTRINGVariable", CONFIGFILEC );

  ASSERT ( *varp == NULL, stat, CONFIGFILEH_ENONULL, CONFIGFILEH_MSGENONULL);

  param.varName = varName;
  param.fmt = FMT_STRING;
  param.strictness = CONFIGFILE_IGNORE;

  LALReadConfigVariable (stat, (void*) &str, cfgdata, &param, wasRead);

  if ( *wasRead && (str!=NULL) )
    {
      INT2 numQuotes = 0;
      CHAR *ptr = str;
      /* count number of quotation marks */
      while ( (ptr = strchr(ptr, '"')) )
	{
	  numQuotes ++;
	  ptr ++;
	} /* while quotes found */

      /* check balanced quotes (don't allow escaping for now) */
      if ( (numQuotes !=0) && (numQuotes != 2) ) {
	ABORT (stat, CONFIGFILEH_ESTRING, CONFIGFILEH_MSGESTRING);
      }
      if ( numQuotes==2 )  
	{
	  /* allowed only at end and beginning */
	  if ( (str[0] != '"') || (str[strlen(str)-1] != '"') ) {
	    ABORT (stat, CONFIGFILEH_ESTRING, CONFIGFILEH_MSGESTRING);
	  }
	  /* quotes ok, now remove them */
	  if ( (ret = LALMalloc( strlen(str) -2 + 1)) == NULL ) {
	    ABORT (stat, CONFIGFILEH_EMEM, CONFIGFILEH_MSGEMEM);
	  }
	  str[strlen(str)-1] = 0;
	  strcpy (ret, str+1);
	  LALFree (str);
	} /* if 2 quotation marks */
      else
	ret = str;	/* no quotes, just return string */
      
      *varp = ret;

    } /* if wasRead */
  else
    *varp = NULL;

  RETURN (stat);

} /* LALReadConfigSTRINGVariable() */



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
/* <lalVerbatim file="ConfigFileCP"> */
void
LALReadConfigSTRINGNVariable (LALStatus *stat, 
			      CHARVector *varp, 	/* output: must be allocated! */
			      const LALParsedDataFile *cfgdata, /* pre-parsed config-data */
			      const CHAR *varName,	/* variable-name */
			      BOOLEAN *wasRead)
{ /* </lalVerbatim> */
  CHAR *tmp = NULL;

  INITSTATUS( stat, "LALReadSTRINGNVariable", CONFIGFILEC );
  ATTATCHSTATUSPTR (stat);
  
  /* This traps coding errors in the calling routine. */
  ASSERT( varp != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( varp->data != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( varp->length != 0, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  
  TRY (LALReadConfigSTRINGVariable (stat->statusPtr, &tmp, cfgdata, varName, wasRead), stat);

  if (*wasRead && tmp)
    {
      strncpy (varp->data, tmp, varp->length - 1);
      varp->data[varp->length-1] = '\0';
      LALFree (tmp);
      varp->length = strlen (varp->data);
      *wasRead = TRUE;
    }
  else
    *wasRead = FALSE;
    
  DETATCHSTATUSPTR (stat);
  RETURN (stat);  

} /* LALReadConfigSTRINGNVariable() */



/*----------------------------------------------------------------------
 * check if all lines of config-file have been successfully read in 
 * and issue a warning or error (depending on strictness) if not
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="ConfigFileCP"> */
void
LALCheckConfigReadComplete (LALStatus *stat, 
			    const LALParsedDataFile *cfgdata, /* config-file data */
			    ConfigStrictness strict)  /* what to do if unparsed lines */
{ /* </lalVerbatim> */
  UINT4 i;

  INITSTATUS( stat, "LALCheckConfigReadComplete", CONFIGFILEC );  

  ASSERT (cfgdata != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT (cfgdata->lines != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT (cfgdata->wasRead != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);

  for (i=0; i < cfgdata->lines->nTokens; i++)
    if (cfgdata->wasRead[i] == 0)
      break;

  if (i != cfgdata->lines->nTokens)
    {
      switch (strict)
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
	default:
	  LALPrintError ("ERROR: config-file entry #%d has not been read!\n", i);
	  LALPrintError ("Line was: '%s'\n", cfgdata->lines->tokens[i]);
	  ABORT (stat, CONFIGFILEH_EUNKNOWN, CONFIGFILEH_MSGEUNKNOWN);
	  break;
	} /* switch strict */
    } /* if some line not read */
  
  RETURN (stat);

} /* LALCheckConfigReadComplete() */



/* ---------------------------------------------------------------------- 
 *   INTERNAL FUNCTIONS FOLLOW HERE
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
 * little helper function:  turn a string into lowercase
 *----------------------------------------------------------------------*/
void
LALLowerCaseString (LALStatus *stat, CHAR *string)
{
  UINT4 i;

  INITSTATUS( stat, "LALLowerCaseString", CONFIGFILEC );  

  ASSERT (string != NULL, stat, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);

  for (i=0; i < strlen (string); i++)
    string[i] = TOLOWER( string[i] );
  
  RETURN (stat);

} /* LALLowerCaseString() */

/*----------------------------------------------------------------------
 * tolower() replacement w/o locale
 *----------------------------------------------------------------------*/
static int
TOLOWER(int c)
{
    if (c) {
        char *p = strchr(upper_chars, c);

        if (p) {
            c = lower_chars[p - upper_chars];
        }
    }
    return c;
} /* TOLOWER() */

/*----------------------------------------------------------------------*/


/* ----------------------------------------------------------------------
 * cleanConfig(): do some preprocessing on the config-file, namely 'erase' 
 * all comments by '\n', and glue '\'-continued lines
 *----------------------------------------------------------------------*/
void
cleanConfig (CHARSequence *text)
{
  size_t len;  
  CHAR *ptr, *ptr2, *eol;
  BOOLEAN inQuotes = 0;

  /*----------------------------------------------------------------------
   * RUN 1: clean out comments, by replacing them by '\n' 
   */
  ptr = text->data;
  while ( *ptr )
    {
      if ( (*ptr) == '\"' )
	inQuotes = !inQuotes;	/* flip state */

      if ( ((*ptr) == '#') || ( (*ptr) == ';') )
	if ( !inQuotes )	/* only consider as comments if not quoted */
	  {
	    len = strcspn (ptr, "\n"); 
	    memset ( (void*)ptr, '\n', len); 	
	  }
	
      ptr ++;

    } /* while *ptr */

  /*----------------------------------------------------------------------
   * RUN 2: do line-gluing when '\' is found at end-of-line 
   */
  ptr = text->data;
  while ( (ptr = strchr(ptr, '\\')) != NULL )
    {
      if ( ptr[1] == '\n' ) 
	{	
	  /* ok, now it gets a bit tricky: to avoid getting spurious spaces from
	   * the line-continuation, we shift the rest of the file forward by 2 positions 
	   * to nicely fit to the previous line... 
	   */
	  len = strlen (ptr+2);
	  memmove(ptr, ptr+2, len+1);	/* move the whole rest (add +1 for '\0') */
	}
    } /* while '\' found in text */

  /*----------------------------------------------------------------------
   * RUN 3: turn all tabs into single spaces.. 
   */
  ptr = text->data;
  while ( (ptr = strchr(ptr, '\t')) != NULL )
    *ptr = ' ';

  /*---------------------------------------------------------------------- 
   * RUN 4: get rid of initial and trailing whitespace (replace it by '\n') 
   */
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

