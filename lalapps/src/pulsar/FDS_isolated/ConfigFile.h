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
#define CONFIGFILEH_ENONULL		6
#define CONFIGFILEH_ESTRICT		7
#define CONFIGFILEH_EUNKNOWN		8

#define CONFIGFILEH_MSGENULL 		"Arguments contained an unexpected null pointer."
#define CONFIGFILEH_MSGEFILE		"File error."
#define CONFIGFILEH_MSGEVAR		"Config variable not found."
#define CONFIGFILEH_MSGEFMT		"Config variable not readable using given format-string."
#define CONFIGFILEH_MSGETOKENS		"The input TokenList is malformed."
#define CONFIGFILEH_MSGENONULL		"Input pointer is not NULL"
#define CONFIGFILEH_MSGESTRICT		"Strictness parameter out of range"
#define CONFIGFILEH_MSGEUNKNOWN		"Unknown config-file entry found"

/*************************************************** </lalErrTable> */
enum {
  CONFIGFILE_IGNORE = 0,
  CONFIGFILE_WARN,
  CONFIGFILE_ERROR,
  CONFIGFILE_LAST
};


typedef struct {
  CHAR *varName;
  CHAR *fmt;
  INT2 strictness;	/* what to do if variable not found: ignore, warn, error */
} LALConfigVar_t;


typedef struct {
  TokenList *lines;	/* the pre-cleaned contents of the config-file */
  BOOLEAN *wasRead;	/* keep track of successfully read lines for checking */
} LALConfigData_t;

/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{ConfigFileHV}}
\newpage\input{ConfigFileC}
******************************************************* </lalLaTeX> */

/* Function prototypes */
void LALLoadConfigFile (LALStatus *stat, LALConfigData_t *cfgdata, FILE *instream);
void LALDestroyConfigData (LALStatus *stat, LALConfigData_t *cfgdata);

void LALReadConfigVariable (LALStatus *stat, void *varp, LALConfigData_t *cfgdata, LALConfigVar_t *param);

void LALReadConfigBOOLVariable (LALStatus *stat, BOOLEAN *varp, LALConfigData_t *cfgdata, CHAR *varName);
void LALReadConfigINT2Variable (LALStatus *stat, INT2 *varp, LALConfigData_t *cfgdata, CHAR *varName);
void LALReadConfigINT4Variable (LALStatus *stat, INT4 *varp, LALConfigData_t *cfgdata, CHAR *varName);
void LALReadConfigREAL4Variable (LALStatus *stat, REAL4 *varp, LALConfigData_t *cfgdata, CHAR *varName);
void LALReadConfigREAL8Variable (LALStatus *stat, REAL8 *varp, LALConfigData_t *cfgdata, CHAR *varName);
void LALReadConfigSTRINGVariable (LALStatus *stat, CHAR **varp, LALConfigData_t *cfgdata, CHAR *varName);

void LALCheckCfgReadComplete (LALStatus *stat, LALConfigData_t *cfgdata, INT2 strictness);

void testConfigFile(void);
  
/********************************************************** <lalLaTeX>
\newpage\input{LALSampleTestC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */



