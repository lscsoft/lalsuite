/************************************ <lalVerbatim file="ConfigFileHV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\section{Header \texttt{ConfigFile.h}}
\label{s:ConfigFile.h}

Header file for ConfigFile

\subsection*{Synopsis}
\begin{verbatim}
#include "ConfigFile.h"
\end{verbatim}

\noindent This header provides two trivial functions to divide real
numbers.  It exists primarily to demonstrate documentation and coding
standards for LAL headers.

******************************************************* </lalLaTeX> */

#ifndef _CONFIGFILE_H  /* Double-include protection. */
#define _CONFIGFILE_H

#include <lal/LALDatatypes.h>
#include <lal/StringInput.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID( CONFIGFILEH, "$Id" );

/********************************************************** <lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
***************************************************** <lalErrTable> */
#define CONFIGFILEH_ENULL 		1
#define CONFIGFILEH_EFILE		2
#define CONFIGFILEH_EVAR		3
#define CONFIGFILEH_EFMT		4
#define CONFIGFILEH_ETOKENS		5

#define CONFIGFILEH_MSGENULL 		"Arguments contained an unexpected null pointer."
#define CONFIGFILEH_MSGEFILE		"File error."
#define CONFIGFILEH_MSGEVAR		"Config variable not found."
#define CONFIGFILEH_MSGEFMT		"Config variable not readable using given format-string."
#define CONFIGFILEH_MSGETOKENS		"The input TokenList is malformed."

/*************************************************** </lalErrTable> */

typedef struct {
  CHAR *varName;
  CHAR *fmt;
} LALConfigVar_t;

#define FMT_STRING "string"    /* reading in quoted strings needs some special treatment */

/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{ConfigFileHV}}
\newpage\input{ConfigFileC}
******************************************************* </lalLaTeX> */

/* Function prototypes */
void LALParseConfigFile (LALStatus *stat, TokenList **lines, FILE *instream);
void LALReadConfigVariable (LALStatus *stat, void *varp, TokenList *lines, LALConfigVar_t *param);
void LALReadConfigBOOLVariable (LALStatus *stat, BOOLEAN *varp, TokenList *lines, CHAR *varName);
  
/********************************************************** <lalLaTeX>
\newpage\input{LALSampleTestC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */



