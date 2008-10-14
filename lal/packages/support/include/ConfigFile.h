/*
 * Copyright (C) 2004 Reinhard Prix
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
 * \date $Date$
 * \brief Header file defining the API for ConfigFile.c.
 */

/************************************ <lalVerbatim file="ConfigFileHV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/**************************************************** <lalLaTeX>
\section{Header \texttt{ConfigFile.h}}
\label{s:ConfigFile.h}

Routines for general config-file reading.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/ConfigFile.h>
\end{verbatim}

***************************************************** </lalLaTeX> */

#ifndef _CONFIGFILE_H  /* Double-include protection. */
#define _CONFIGFILE_H

#include <lal/LALDatatypes.h>
#include <lal/StringInput.h>

/* C++ protection. */
#ifdef  __cplusplus   
extern "C" {
#endif

NRCSID( CONFIGFILEH, "$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
***************************************************** <lalErrTable> */
/** \name Error codes */
/*@{*/
#define CONFIGFILEH_ENULL 		1
#define CONFIGFILEH_EFILE		2
#define CONFIGFILEH_EVAR		3
#define CONFIGFILEH_EFMT		4
#define CONFIGFILEH_ETOKENS		5
#define CONFIGFILEH_ENONULL		6
#define CONFIGFILEH_EUNKNOWN		8
#define CONFIGFILEH_EMEM		9
#define CONFIGFILEH_EBOOL		10
#define CONFIGFILEH_ESTRING		11

#define CONFIGFILEH_MSGENULL 		"Arguments contained an unexpected null pointer."
#define CONFIGFILEH_MSGEFILE		"File error."
#define CONFIGFILEH_MSGEVAR		"Config variable not found."
#define CONFIGFILEH_MSGEFMT		"Config variable not readable using given format-string."
#define CONFIGFILEH_MSGETOKENS		"The input ConfigData seems corrupted."
#define CONFIGFILEH_MSGENONULL		"Output pointer is not NULL"
#define CONFIGFILEH_MSGEUNKNOWN		"Unknown config-file entry found"
#define CONFIGFILEH_MSGEMEM		"Out of memory"
#define CONFIGFILEH_MSGEBOOL		"Illegal BOOLEAN entry"
#define CONFIGFILEH_MSGESTRING		"Malformed quoted string"
/*@}*/
/*************************************************** </lalErrTable> */

/** Levels of strictness for config-file parsing. */
typedef enum {
  CONFIGFILE_IGNORE = 0,	/**< ignore missing config-variable or unparsed config-entries */
  CONFIGFILE_WARN,		/**< issue a warning but don't report an error. */
  CONFIGFILE_ERROR,		/**< issue an error-message and report a LAL-error */
  CONFIGFILE_LAST
} ConfigStrictness;


/** This structure defines a config-variable to be read in using the
 * general-purpose reading function LALReadConfigVariable(). */
typedef struct {
  const CHAR *varName;		/**< Variable-name to be read in the config-file */
  const CHAR *fmt;		/**< Format string for reading (<tt>sscanf()</tt>-style) */
  ConfigStrictness strictness;	/**< what to do if variable not found: ignore, warn, error */
} LALConfigVar;


/** This structure is returned by LALParseDataFile() and holds the contents of an 
 * ASCII data-file in a pre-parsed form, namely stripped from all comments ('#', ';'+), 
 * spurious whitespaces, and separated into lines (taking into account line-continuation 
 * by '\' at the end of lines).
 * This is used as the input structure in the config-variable reading routines.
 */
typedef struct {
  TokenList *lines;	/**< list of pre-parsed data-file lines */
  BOOLEAN *wasRead;	/**< keep track of successfully read lines for strictness-checking */
} LALParsedDataFile;

/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{ConfigFileHV}}
\newpage\input{ConfigFileC}
\newpage\input{ConfigFileTestC}
******************************************************* </lalLaTeX> */

/* Function prototypes */
void LALParseDataFile (LALStatus *, LALParsedDataFile **cfgdata, const CHAR *fname);
void LALDestroyParsedDataFile (LALStatus *, LALParsedDataFile **cfgdata);

int XLALParseDataFile (LALParsedDataFile **cfgdata, const CHAR *fname);
int XLALDestroyParsedDataFile (LALParsedDataFile **cfgdata);

void 
LALReadConfigBOOLVariable (LALStatus *, 
			  BOOLEAN *varp, 
			  const LALParsedDataFile *cfgdata, 
			  const CHAR *varName, 
			  BOOLEAN *wasRead);

void
LALReadConfigINT4Variable (LALStatus *,
			   INT4 *varp, 
			   const LALParsedDataFile *cfgdata, 
			   const CHAR *varName, 
			   BOOLEAN *wasRead);

void
LALReadConfigREAL8Variable (LALStatus *, 
			    REAL8 *varp, 
			    const LALParsedDataFile *cfgdata, 
			    const CHAR *varName, 
			    BOOLEAN *wasRead);

void 
LALReadConfigSTRINGVariable (LALStatus *, 
			     CHAR **varp, 
			     const LALParsedDataFile *cfgdata, 
			     const CHAR *varName,
			     BOOLEAN *wasRead);

void
LALReadConfigSTRINGNVariable (LALStatus *, 
			      CHARVector *varp,
			      const LALParsedDataFile *cfgdata, 
			      const CHAR *varName,
			      BOOLEAN *wasRead);

void
LALReadConfigVariable (LALStatus *, 
		       void *varp,
		       const LALParsedDataFile *cfgdata,
		       const LALConfigVar *param,
		       BOOLEAN *wasRead);

void LALCheckConfigReadComplete (LALStatus *, const LALParsedDataFile *cfgdata, ConfigStrictness strict);

void LALLowerCaseString (LALStatus *, CHAR *string);

/* C++ protection. */
#ifdef  __cplusplus
}
#endif  

#endif  /* Double-include protection. */



