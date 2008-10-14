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

/** \file
 * \ingroup UserInput
 * \author Reinhard Prix
 *
 * \brief General-purpose routines for config-file reading.
 *
 *
\par Description

This module provides routines for reading formatted
config-files containing definitions of the form <tt>variable = value</tt>.
The general syntax is somewhat similar to the one provided by the
perl-module <tt>ConfigParser</tt> (cf.
http://www.python.org/doc/current/lib/module-ConfigParser.html)
but (currently) without the possibility of "chapters".
Comments are allowed using either '<tt>\#</tt>', '<tt>;</tt>' or <tt>%</tt>.
You can also use line-continuation  using a '<tt>\\</tt>' at the end of the line.
Also note that comment-signs '<tt>\#;%</tt>' within double-quotes '<tt>"..."</tt>'
are <em>not</em> treated as comment-characters.  The general syntax is best illustrated
using a simple example:
\code
# comment line
var1 = 1.0    ; you can also comment using semi-colons
somevar = some text.\
        You can also use\
        line-continuation
   var3 = 4      # whatever that means
note = "this is also possible, and # here does nothing"
a_switch = true	 #possible values: 0,1,true,false,yes,no, case insensitive
%% etc etc.
\endcode

Note that TABS generally get replaced by a single space, which can be
useful in the case of line-continuation (see example). All leading and
trailing spaces in are ignore (except within double-quotes).

The general approach of reading from such a config-file, is to first
call LALParseDataFile() which loads and pre-parses the contents of the
config-file into the structure LALParsedDataFile. Then one can read in
config-variables either using one of the type-strict custom-wrappers
<tt>LALReadConfig<TYPE>Variable()</tt> or the general-purpose reading function
LALReadConfigVariable().


A boolean variable read by LALReadConfigBOOLVariable() can have any of the values
<tt>{"1", "0", "yes", "no", "true", "false"}</tt>, where the comparison is done
<em>case-insensitively</em>, i.e. you can also use "True" or "FALSE"....


If one wishes a tight sytnax for the config-file, one can check
that there are no illegal entries in the config-file. This is done
by checking at the end that all config-file entries have been
successfully parsed, using:
LALCheckConfigReadComplete(), where \a strictness is either
CONFIGFILE_WARN or CONFIGFILE_ERROR.
In the first case only a warning is issued, while in the second it is
treated as a LAL-error if some config-file entries have not been
read-in. (The use of this function is optional).


The configfile-data should be freed at the end using
LALDestroyParsedDataFile().

\par Uses
\code
LALCHARReadSequence()
LALCreateTokenList()       LALDestroyTokenList()
LALCalloc()                LALMalloc()             LALFree()
LALPrintError()            LALOpenDataFile()                 fclose()
\endcode

\par Notes

LALReadConfigSTRINGVariable() and LALReadConfigSTRINGVariable() are not
the same as using <tt>"%s"</tt> as a format string, as they read the
<em>rest</em> of the logical line (excluding comments) as a string.


In the case of LALReadConfigSTRINGVariable(), the required
memory is allocated and has to be freed by the caller, while for
LALReadConfigSTRINGVariable() the caller has to provide a
CHARVector of length \f$N\f$, which defines the maximum length of
string to be read.


\note instead of using these functions directly, it might be
more convenient to use the UserInput.c infrastructure

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

\subsubsection*{Description}

This module provides routines for reading formatted
config-files containing definitions of the form \mbox{\texttt{variable = value}}.
The general syntax is somewhat similar to the one provided by the
perl-module \texttt{ConfigParser} (cf.
\verb+http://www.python.org/doc/current/lib/module-ConfigParser.html+ )
but (currently) without the possibility of "chapters".
Comments are allowed using either '\#', ';' or '%'. You can also use
standard line-continuation  using a '\verb+\+' at the end of the line.
Also note that '\#', ';' or '%' within double-quotes '\"' are \emph{not}
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
\verb+LALLoadConfigFile(status, LALConfigData *cfg, FILE *fp)+,
which loads and pre-parses the contents of the config-file into the
structure \verb+LALConfigData+. Then one can then read in
config-variables either using one of the custom-wrappers:\\
\verb+LALReadConfig<TYPE>Variable(status, <TYPE> *cvar, LALConfigData *cfg, CHAR *varname)+
or the general-purpose reading function:\\
\verb+LALReadConfigVariable(status, void *cvar, LALConfigData *cfg, LALConfigVar *var)+


A boolean variable read by \verb+LALReadConfigBOOLVariable()+ can have any of the values
\verb+{1, 0, yes, no, true, false}+, where the comparison is done \emph{case-insensitively},
i.e. you also use "True" or "FALSE"....


If one wishes a ``tight'' sytnax for the config-file, one can check
that there are no "illegal" entries in the config-file. This is done
by checking at the end that all config-file entries have been
successfully parsed, using: \\
\verb+LALCheckConfigReadComplete (status, LALConfigData *cfg, INT2 strictness)+,
where \verb+strictness+ is either \verb+CONFIGFILE_WARN+ or \verb+CONFIGFILE_ERROR+.
In the first case only a warning is issued, while in the second it is
treated as a LAL-error if some config-file entries have not been
read-in. (The use of this function is optional).


The configfile-data should be freed at the end using\\
\verb+void LALDestroyParsedDataFile(LALStatus *status, LALConfigData *cfg)+.

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

#include <errno.h>

#if HAVE_SYS_STAT_H
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

/* #include <ctype.h> */  /* don't use this, as it binds us to GLIBC_2.3 symbols!! */

#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/StreamInput.h>
#include <lal/LogPrintf.h>

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



/** Parse an ASCII data-file into a pre-cleaned array of lines.
 *
 * The cleaning gets rid of comments ('#', ';'), empty lines,
 * and performs line-continuation if '\' is found at EOL
 *
 */
void
LALParseDataFile (LALStatus *status,
		  LALParsedDataFile **cfgdata, 	/**< [out] pre-parsed data-file lines */
		  const CHAR *fname)		/**< [in] name of config-file to be read */
{

  CHARSequence *rawdata = NULL;
  FILE *fp;
#if HAVE_STAT
  struct stat stat_out;
#endif

  INITSTATUS( status, "LALParseDataFile", CONFIGFILEC );
  ATTATCHSTATUSPTR (status);

  ASSERT (*cfgdata == NULL, status, CONFIGFILEH_ENONULL, CONFIGFILEH_MSGENONULL);
  ASSERT (fname != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);

#if HAVE_STAT
  if (  stat ( fname, &stat_out ) )
    {
      LogPrintf ( LOG_CRITICAL, "Could not stat data-file: `%s` : \n\n", fname, strerror(errno) );
      ABORT (status, CONFIGFILEH_EFILE, CONFIGFILEH_MSGEFILE);
    }

  if ( ! S_ISREG(stat_out.st_mode)  )
    {
      LogPrintf ( LOG_CRITICAL, "'%s' does not seem to be a regular file!\n");
      ABORT (status, CONFIGFILEH_EFILE, CONFIGFILEH_MSGEFILE);
    }
#endif

  if ( (fp = LALOpenDataFile (fname)) == NULL) {
    LogPrintf ( LOG_CRITICAL, "Could not open data-file: `%s`\n\n", fname);
    ABORT (status, CONFIGFILEH_EFILE, CONFIGFILEH_MSGEFILE);
  }

  LALCHARReadSequence (status->statusPtr, &rawdata, fp);
  fclose (fp);
  CHECKSTATUSPTR (status);

  if (rawdata == NULL) {
    ABORT (status, CONFIGFILEH_EFILE, CONFIGFILEH_MSGEFILE);
  }

  /* get rid of comments and do line-continuation */
  cleanConfig (rawdata);

  if ( (*cfgdata = LALCalloc (1, sizeof(LALParsedDataFile))) == NULL) {
    ABORT (status, CONFIGFILEH_EMEM, CONFIGFILEH_MSGEMEM);
  }

  /* parse this into individual lines */
  LALCreateTokenList (status->statusPtr, &((*cfgdata)->lines), rawdata->data, "\n");
  LALFree (rawdata->data);
  LALFree (rawdata);

  BEGINFAIL (status)
    LALFree (*cfgdata);
  ENDFAIL (status);

  /* initialize the 'wasRead' flags for the lines */
  if ( (*cfgdata)->lines->nTokens )
    {
      if ( ((*cfgdata)->wasRead =
	    LALCalloc(1,(*cfgdata)->lines->nTokens * sizeof( (*cfgdata)->wasRead[0]))) == NULL)
	{
	  LALFree ((*cfgdata)->lines);
	  ABORT (status, CONFIGFILEH_EMEM, CONFIGFILEH_MSGEMEM);
	}
    }
  else
    (*cfgdata)->wasRead = NULL;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALLoadConfigFile() */



/** Parse an ASCII data-file into a pre-cleaned array of lines.
 *
 * The cleaning gets rid of comments ('#', ';'), empty lines,
 * and performs line-continuation if '\' is found at EOL
 */
int
XLALParseDataFile (LALParsedDataFile **cfgdata, /**< [out] pre-parsed data-file lines */
		   const CHAR *fname)		/**< [in] name of config-file to be read */
{

  CHARSequence *rawdata = NULL;
  FILE *fp;
  int err = 0;  /* error code */
#if HAVE_STAT
  struct stat stat_out;
#endif

  if (*cfgdata != NULL) {
      fprintf(stderr, CONFIGFILEH_MSGENONULL);
      return CONFIGFILEH_ENONULL;
  }
  if (fname == NULL) {
      fprintf(stderr, CONFIGFILEH_MSGENULL);
      return CONFIGFILEH_ENULL;
  }

#if HAVE_STAT
  if (  stat ( fname, &stat_out ) )
    {
      LogPrintf ( LOG_CRITICAL, "Could not stat data-file: `%s` : \n\n", fname, strerror(errno) );
      fprintf(stderr, CONFIGFILEH_MSGEFILE);
      return CONFIGFILEH_EFILE;
    }

  if ( ! S_ISREG(stat_out.st_mode)  )
    {
      LogPrintf ( LOG_CRITICAL, "'%s' does not seem to be a regular file!\n");
      fprintf(stderr, CONFIGFILEH_MSGEFILE);
      return CONFIGFILEH_EFILE;
    }
#endif

  if ( (fp = LALOpenDataFile (fname)) == NULL) {
    LogPrintf ( LOG_CRITICAL, "Could not open data-file: `%s`\n\n", fname);
    fprintf(stderr, CONFIGFILEH_MSGEFILE);
    return CONFIGFILEH_EFILE;
  }

  err = XLALCHARReadSequence (&rawdata, fp);
  fclose (fp);
  if (err)
      return err;

  if (rawdata == NULL) {
    fprintf(stderr, CONFIGFILEH_MSGEFILE);
    return CONFIGFILEH_EFILE;
  }

  /* get rid of comments and do line-continuation */
  cleanConfig (rawdata);

  if ( (*cfgdata = LALCalloc (1, sizeof(LALParsedDataFile))) == NULL) {
    fprintf(stderr, CONFIGFILEH_MSGEMEM);
    return CONFIGFILEH_EMEM;
  }

  /* parse this into individual lines */
  err = XLALCreateTokenList (&((*cfgdata)->lines), rawdata->data, "\n");
  LALFree (rawdata->data);
  LALFree (rawdata);

  if (err) {
    LALFree (*cfgdata);
    return err;
  }

  /* initialize the 'wasRead' flags for the lines */
  if ( (*cfgdata)->lines->nTokens )
    {
      if ( ((*cfgdata)->wasRead =
	    LALCalloc(1,(*cfgdata)->lines->nTokens * sizeof( (*cfgdata)->wasRead[0]))) == NULL)
	{
	  LALFree ((*cfgdata)->lines);
	  fprintf(stderr, CONFIGFILEH_MSGEMEM);
	  return CONFIGFILEH_EMEM;
	}
    }
  else
    (*cfgdata)->wasRead = NULL;

  return 0;
}



/** Free memory associated with a LALParsedDataFile structure.
 */
void
LALDestroyParsedDataFile (LALStatus *status,
			  LALParsedDataFile **cfgdata)	/**< [in/out] config-file data */
{
  INITSTATUS( status, "LALDestroyConfigData", CONFIGFILEC );
  ATTATCHSTATUSPTR (status);

  ASSERT (cfgdata != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT (*cfgdata != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT ((*cfgdata)->lines != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);

  TRY ( LALDestroyTokenList (status->statusPtr, &((*cfgdata)->lines)), status);

  if ( (*cfgdata)->wasRead )
    LALFree ( (*cfgdata)->wasRead);

  LALFree ( *cfgdata );

  *cfgdata = NULL;

  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* LALDestroyConfigData() */




/** Free memory associated with a LALParsedDataFile structure.
 */
int
XLALDestroyParsedDataFile (LALParsedDataFile **cfgdata)	/**< [in/out] config-file data */
{
  int err = 0;  /* error code */

  if ( ! cfgdata || ! *cfgdata || ! (*cfgdata)->lines ) {
    fprintf(stderr, CONFIGFILEH_MSGENULL);
    return CONFIGFILEH_ENULL;
  }

  err = XLALDestroyTokenList ( &((*cfgdata)->lines) );
  if (err) {
    fprintf(stderr, "Error destroying token list.\n");
    return err;
  }
  
  if ( (*cfgdata)->wasRead )
    LALFree ( (*cfgdata)->wasRead );

  LALFree ( *cfgdata );

  *cfgdata = NULL;

  return 0;
}



/** Parser for config-file: can read config-variables of the form
 *	VARIABLE [=:] VALUE.
 * Input is a TokenList containing the 'logical' lines of the cleaned config-file
 *
 * - <tt>param->varName</tt> is the name of the config-variable to read
 * - <tt>param->fmt</tt>     is the format string to use for reading
 *
 * \note a special format-string is FMT_STRING, which means read the whole remaining line
 *   which is different from "%s"! (reads only one word)
 *   In this case, this also does the memory-allocation!
 *
 */
void
LALReadConfigVariable (LALStatus *status,
		       void *varp, 			/**< [out] result gets written here! */
		       const LALParsedDataFile *cfgdata,/**< [in] pre-parsed config-data */
		       const LALConfigVar *param,	/**< [in]  var-name, fmt-string, strictness */
		       BOOLEAN *wasRead)		/**< [out] did we succeed in reading? */
{
  CHAR *found = NULL;
  INT2 ret = 0;

  UINT4 i;
  INT4 linefound = -1;
  size_t len;
  size_t searchlen = strlen (param->varName);

  INITSTATUS( status, "LALReadConfigVariable", CONFIGFILEC );

  /* This traps coding errors in the calling routine. */
  ASSERT( cfgdata != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( cfgdata->lines != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( cfgdata->wasRead != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( varp != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( param->varName != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL );
  ASSERT( param->fmt != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL );

  *wasRead = FALSE;

  /* let's look for the variable-name in the token-list (has to at beginning of line!) */
  for (i=0; i<cfgdata->lines->nTokens; i++)
    {
      len = strcspn (cfgdata->lines->tokens[i], WHITESPACE "=:"); /* get length of variable-name */
      if (len == 0) { /* malformed token-list */
	ABORT (status, CONFIGFILEH_ETOKENS, CONFIGFILEH_MSGETOKENS);
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
	  RETURN (status);
	  break;
	case CONFIGFILE_WARN:
	  if (lalDebugLevel & LALWARNING)
	    LogPrintf ( LOG_CRITICAL, "Warning: Config-file variable '%s' was not found!\n", param->varName);
	  RETURN (status);
	  break;

	case CONFIGFILE_ERROR:
	default:
	  LogPrintf ( LOG_CRITICAL, "Error: Config-file variable %s was not found!\n", param->varName);
	  ABORT (status, CONFIGFILEH_EVAR, CONFIGFILEH_MSGEVAR );
	  break;
	} /* switch (strictness) */

    } /* if not found */

  /* now read the value into the variable */

  /* reading a quoted string needs some special treatment: */
  if ( !strcmp(param->fmt, FMT_STRING) )
    {
      /* NOTE: varp here is supposed to be a pointer to CHAR* !! */
      CHAR **cstr = (CHAR**) varp;

      ASSERT ( *cstr == NULL, status, CONFIGFILEH_ENONULL, CONFIGFILEH_MSGENONULL);

      (*cstr) = (CHAR*) LALMalloc( strlen (found) + 1);
      strcpy ( (*cstr), found);
      ret = 1;
    }
  else  /* but the default case is just sscanf... */
    ret = sscanf (found, param->fmt, varp);

  if ( (ret == 0) || (ret == EOF) )
    {
      LogPrintf (LOG_CRITICAL, "ERROR: Config-file variable %s was not readable using the format %s\n\n", param->varName, param->fmt);
      ABORT( status, CONFIGFILEH_EFMT, CONFIGFILEH_MSGEFMT );
    }

  /* ok, we have successfully read in the config-variable: let's make a note of it */
  cfgdata->wasRead[linefound] = 1;

  *wasRead = TRUE;

  RETURN (status);

} /* LALReadConfigVariable() */



/** Type-specialization of generic reading-function LALReadConfigVariable() to BOOLEAN variables.
 */
void
LALReadConfigBOOLVariable (LALStatus *status,
			   BOOLEAN *varp, 		 /**< [out] variable to store result */
			   const LALParsedDataFile *cfgdata,/**< [in] pre-parsed config-data */
			   const CHAR *varName,		 /**< [in] variable-name to read */
			   BOOLEAN *wasRead)		 /**< [out] did we succeed in reading? */
{
  CHAR *tmp = NULL;
  INT2 ret = -1;	/* -1 means no legal value has been parsed */

  INITSTATUS( status, "LALReadConfigBOOLVariable", CONFIGFILEC );
  ATTATCHSTATUSPTR (status);

  *wasRead = FALSE;

  /* first read the value as a string */
  TRY (LALReadConfigSTRINGVariable (status->statusPtr, &tmp, cfgdata, varName, wasRead), status);

  if (*wasRead && tmp) /* if we read anything at all... */
    {
      /* get rid of case ambiguities */
      TRY (LALLowerCaseString (status->statusPtr, tmp), status);

      /* try to parse it as a bool */
      if (      !strcmp(tmp, "yes") || !strcmp(tmp, "true") || !strcmp(tmp,"1"))
	ret = 1;
      else if ( !strcmp (tmp, "no") || !strcmp(tmp,"false") || !strcmp(tmp,"0"))
	ret = 0;
      else
	{
	  LogPrintf ( LOG_CRITICAL,  "illegal bool-value `%s`\n", tmp);
	  LALFree (tmp);
	  ABORT (status, CONFIGFILEH_EBOOL, CONFIGFILEH_MSGEBOOL);
	}
      LALFree (tmp);

      if (ret != -1)	/* only set value of something has been found */
	{
	  *varp = (BOOLEAN) ret;
	  *wasRead = TRUE;
	}

    } /* if wasRead && tmp */

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadConfigBOOLVariable() */


/** Type-specialization of generic reading-function LALReadConfigVariable() to INT4 variables.
 */
void
LALReadConfigINT4Variable (LALStatus *status,
			   INT4 *varp,
			   const LALParsedDataFile *cfgdata,
			   const CHAR *varName,
			   BOOLEAN *wasRead)
{
  LALConfigVar param = {0,0,0};

  INITSTATUS( status, "LALReadConfigINT4Variable", CONFIGFILEC );

  param.varName = varName;
  param.fmt = "%" LAL_INT4_FORMAT;
  param.strictness = CONFIGFILE_IGNORE;

  LALReadConfigVariable (status, (void*) varp, cfgdata, &param, wasRead);

  RETURN (status);

} /* LALReadConfigINT4Variable() */

/** Type-specialization of generic reading-function LALReadConfigVariable() to REAL8 variables.
 */
void
LALReadConfigREAL8Variable (LALStatus *status,
			    REAL8 *varp,
			    const LALParsedDataFile *cfgdata,
			    const CHAR *varName,
			    BOOLEAN *wasRead)
{
  LALConfigVar param = {0,0,0};

  INITSTATUS( status, "LALReadConfigREAL8Variable", CONFIGFILEC );

  param.varName = varName;
  param.fmt = "%" LAL_REAL8_FORMAT;
  param.strictness = CONFIGFILE_IGNORE;

  LALReadConfigVariable (status, (void*) varp, cfgdata, &param, wasRead);

  RETURN (status);

} /* LALReadConfigREAL8Variable() */


/** Type-specialization of generic reading-function LALReadConfigVariable() to
 * STRING variables
 * \note this means the rest of the line, NOT "%s"! (but excluding comments of course),
 * \par Note2: if string is quoted by ", everything within quotes is read,
 *       and the quotes are removed here
 *
 */
void
LALReadConfigSTRINGVariable (LALStatus *status,
			     CHAR **varp, 		/**< [out] string, allocated here! */
			     const LALParsedDataFile *cfgdata, /**< [in] pre-parsed config-data */
			     const CHAR *varName,	/**< [in] variable-name to be read */
			     BOOLEAN *wasRead)		/**< [out] did we succeed in reading? */
{
  LALConfigVar param = {0,0,0};
  CHAR *str = NULL;
  CHAR *ret = NULL;

  INITSTATUS( status, "LALReadConfigSTRINGVariable", CONFIGFILEC );

  ASSERT ( *varp == NULL, status, CONFIGFILEH_ENONULL, CONFIGFILEH_MSGENONULL);

  param.varName = varName;
  param.fmt = FMT_STRING;
  param.strictness = CONFIGFILE_IGNORE;

  LALReadConfigVariable (status, (void*) &str, cfgdata, &param, wasRead);

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
	ABORT (status, CONFIGFILEH_ESTRING, CONFIGFILEH_MSGESTRING);
      }
      if ( numQuotes==2 )
	{
	  /* allowed only at end and beginning */
	  if ( (str[0] != '"') || (str[strlen(str)-1] != '"') ) {
	    ABORT (status, CONFIGFILEH_ESTRING, CONFIGFILEH_MSGESTRING);
	  }
	  /* quotes ok, now remove them */
	  if ( (ret = LALMalloc( strlen(str) -2 + 1)) == NULL ) {
	    ABORT (status, CONFIGFILEH_EMEM, CONFIGFILEH_MSGEMEM);
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

  RETURN (status);

} /* LALReadConfigSTRINGVariable() */




/** Type-specialization of generic reading-function LALReadConfigVariable() to
 * reading of <em>fixed-length</em> strings.
 * Another variant of string-reading:similar to ReadConfigSTRINGVariable(), but
 * here a fixed-size CHAR-array is used as input, no memory is allocated by
 * the function.
 * \note you have to provide the length of your string-array as input in <tt>varp->length</tt>
 * (this is basically a wrapper for ReadConfigSTRINGVariable())
 *
 * \par Note 2: the behaviour is similar to strncpy, i.e. we silently clip the
 *       string to the right length, BUT we also 0-terminate it properly.
 *       No error or warning is generated when clipping occurs!
 *
 * \par Note 3: at return, the value <tt>varp->length</tt> is set to the length of the
 *        string copied
 *
 */
void
LALReadConfigSTRINGNVariable (LALStatus *status,
			      CHARVector *varp, 	/**< [out] must be allocated! */
			      const LALParsedDataFile *cfgdata, /**< [in] pre-parsed config-data */
			      const CHAR *varName,	/**< [in] variable-name */
			      BOOLEAN *wasRead)		/**< [out] did we succeed in reading? */
{
  CHAR *tmp = NULL;

  INITSTATUS( status, "LALReadSTRINGNVariable", CONFIGFILEC );
  ATTATCHSTATUSPTR (status);

  /* This traps coding errors in the calling routine. */
  ASSERT( varp != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( varp->data != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT( varp->length != 0, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);

  TRY (LALReadConfigSTRINGVariable (status->statusPtr, &tmp, cfgdata, varName, wasRead), status);

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

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadConfigSTRINGNVariable() */



/** Check if all lines of config-file have been successfully read in
 * and issue a warning or error (depending on strictness) if not.
 */
void
LALCheckConfigReadComplete (LALStatus *status,
			    const LALParsedDataFile *cfgdata, /**< [in] config-file data */
			    ConfigStrictness strict)  	/**< [in] what to do if unparsed lines */
{
  UINT4 i;

  INITSTATUS( status, "LALCheckConfigReadComplete", CONFIGFILEC );

  ASSERT (cfgdata != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT (cfgdata->lines != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);
  ASSERT (cfgdata->wasRead != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);

  for (i=0; i < cfgdata->lines->nTokens; i++)
    {
      if (cfgdata->wasRead[i] == 0)
	{
	  switch (strict)
	    {
	    case CONFIGFILE_IGNORE:
	      continue;
	    case CONFIGFILE_WARN:
	      LogPrintf ( LOG_CRITICAL, "Warning: Ignoring unknown config-file entry '%s'.\n",
			     cfgdata->lines->tokens[i] );
	      continue;
	    case CONFIGFILE_ERROR:
	    default:
	      LogPrintf ( LOG_CRITICAL, "ERROR: config-file entry #%d has not been read!\n", i);
	      LogPrintf ( LOG_CRITICAL, "Line was: '%s'\n", cfgdata->lines->tokens[i]);
	      ABORT (status, CONFIGFILEH_EUNKNOWN, CONFIGFILEH_MSGEUNKNOWN);
	      break;
	    } /* switch strict */
	} /* if some line not read */

    } /* for i < lines */

  RETURN (status);

} /* LALCheckConfigReadComplete() */



/* ----------------------------------------------------------------------
 *   INTERNAL FUNCTIONS FOLLOW HERE
 *----------------------------------------------------------------------*/

/** Helper function:  turn a string into lowercase without using locale-functions.
 */
void
LALLowerCaseString (LALStatus *status,
		    CHAR *string)	/**< [in/out] string to convert */
{
  UINT4 i;

  INITSTATUS( status, "LALLowerCaseString", CONFIGFILEC );

  ASSERT (string != NULL, status, CONFIGFILEH_ENULL, CONFIGFILEH_MSGENULL);

  for (i=0; i < strlen (string); i++)
    string[i] = TOLOWER( string[i] );

  RETURN (status);

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

      if ( ((*ptr) == '#') || ( (*ptr) == ';') || ( (*ptr) == '%') )
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
      eol = strchr (ptr, '\n'); /* point to end-of-line */

      len = strspn (ptr, WHITESPACE);
      if (len) memset ( (void*)ptr, '\n', len);

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

