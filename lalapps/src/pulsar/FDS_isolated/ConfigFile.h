/************************************ <lalVerbatim file="ConfigFileHV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/**************************************************** <lalLaTeX>
\section{Header \texttt{ConfigFile.h}}
\label{s:ConfigFile.h}

Header file for ConfigFile reading.

\subsection*{Synopsis}
\begin{verbatim}
#include "ConfigFile.h"
\end{verbatim}

\noindent This module provides routines for reading formatted
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
somevar = some text. \
       You can also use \
       line-continuation	
   var3 = 4      # whatever that means
note : "this is also possible, and # here does nothing"
# etc etc.
\end{verbatim}

The general approach of reading from such a config-file, is to first
call \verb+LALLoadConfigFile(stat, LALConfigData *, FILE *)+,\\
which loads and pre-parses the contents of the config-file into the
structure \verb+LALConfigData+. Then one can then read in
config-variables either using one of the custom-wrappers:\\
\verb+LALReadConfig<TYPE>Variable(stat, <TYPE> *, LALConfigData *, CHAR *)+

or the general-purpose reading function:\\
\verb+LALReadConfigVariable(stat, void *, LALConfigData *, LALConfigVar *)+

If one wishes a ``tight'' sytnax for the config-file, one can check at
the end that all config-file entries have been successfully read-in, using:\\
\verb+LALCheckConfigReadComplete (stat, LALConfigData *, INT2 strictness)+

where \verb+strictness+ is either \verb+CONFIGFILE_WARN+ or \verb+CONFIGFILE_ERROR+.
In the first case only a warning is issued, while in the second it is
treated as a LAL-error if some config-file entries have not been
read-in. (The use of this function is obviously optional).

The configfile-data should be freed at the end using\\
\verb+void LALDestroyConfigData (LALStatus *stat, LALConfigData *)+.

******************************************************* </lalLaTeX> */

#ifndef _CONFIGFILE_H  /* Double-include protection. */
#define _CONFIGFILE_H

#include <lal/LALDatatypes.h>
#include <lal/StringInput.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID( CONFIGFILEH, "$Id$" );

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
#define CONFIGFILEH_EMEM		9
#define CONFIGFILEH_EOPT		10

#define CONFIGFILEH_MSGENULL 		"Arguments contained an unexpected null pointer."
#define CONFIGFILEH_MSGEFILE		"File error."
#define CONFIGFILEH_MSGEVAR		"Config variable not found."
#define CONFIGFILEH_MSGEFMT		"Config variable not readable using given format-string."
#define CONFIGFILEH_MSGETOKENS		"The input ConfigData seems corrupted."
#define CONFIGFILEH_MSGENONULL		"Output pointer is not NULL"
#define CONFIGFILEH_MSGESTRICT		"Strictness parameter out of range"
#define CONFIGFILEH_MSGEUNKNOWN		"Unknown config-file entry found"
#define CONFIGFILEH_MSGEMEM		"Out of memory"
#define CONFIGFILEH_MSGEOPT		"Unknown command-line option encountered"

/*************************************************** </lalErrTable> */



/***************************************************** <lalLaTeX>

\subsection*{Constants}
%% \idx[Constant]{LAL\_INT2\_FORMAT}
The following are constants used to define various \texttt{strictness}-levels
\begin{description}
\item[\texttt{CONFIGFILE\_IGNORE}] don't care if a config-variable was
  not found, or if config-file entries have remained unread.

\item[\texttt{CONFIGFILE\_WARN}] issue a warning but don't report an
  error.

\item[\texttt{CONFIGFILE\_ERROR}] issue an error-message and report a LAL-error.
\end{description}

The following constants are custom \texttt{fmt}-strings for reading
extended variable formats:
\begin{description}
\item[\texttt{CONFIGFILE\_FMT\_STRING = "string"}] signals reading of
  the whole remaining logical line (excluding comments) as a string. 
\end{description}

******************************************************* </lalLaTeX> */

enum {
  CONFIGFILE_IGNORE = 0,
  CONFIGFILE_WARN,
  CONFIGFILE_ERROR,
  CONFIGFILE_LAST		/* unused: just marks the last entry */
};


/****************************************************** <lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{LALConfigVar}}
\idx[Type]{LALConfigVar}

This structure defines a config-variable to be read in using the
general-purpose reading function \verb+LALReadConfigVariable()+:

\begin{description}
\item[\texttt{CHAR *varName}] Variable-name in the config-file to be read.
\item[\texttt{CHAR *fmt}] Format string used for
  reading. \texttt{sscanf()}-style formats can be used, but also
  custom-defined ones like \texttt{CONFIGFILE\_FMT\_STRING}.
\item[\texttt{INT2 strictness}] Defines wether a variable not found in
  the config-file should be ignored (\texttt{CONFIGFILE\_IGNORE}),
  issue a warning only (\texttt{CONFIGFILE\_WARN}) or also report an
  error (\texttt{CONFIGFILE\_ERROR}).
\end{description}

\subsubsection*{Structure \texttt{LALConfigData}}
\idx[Type]{LALConfigData}

This structure holds the actual config-file contents in a pre-parsed
form, namely stripped from all comments and ignored whitespaces, and
separated into lines (taking into account line-continuation).
This is used as the input structure for the config-variable reading routines.

\begin{description}
\item[\texttt{TokenList *lines}] list of pre-parsed config-file lines
\item[\texttt{BOOLEAN *wasRead}] list of flags signalling if a given
  line has been read successfully into a C-variable 
(used by \verb+LALCheckConfigReadComplete()+). 
\end{description}

******************************************************* </lalLaTeX> */

typedef struct {
  const CHAR *varName;
  const CHAR *fmt;
  INT4 strictness;	/* what to do if variable not found: ignore, warn, error */
} LALConfigVar;


typedef struct {
  TokenList *lines;	/* the pre-cleaned contents of the config-file */
  BOOLEAN *wasRead;	/* keep track of successfully read lines for checking */
} LALConfigData;


typedef enum {
  UVAR_BOOL,
  UVAR_INT4,
  UVAR_REAL8,
  UVAR_CHAR
} UserVarType;

typedef struct {
  const CHAR *name;	/* full name */
  UserVarType type;	/* bool, int, float or char */
  CHAR optchar;		/* cmd-line character */
  const CHAR *help;
  void *varp;		/* pointer to the actual C-variable */
} UserVariable;


#define regUserVar(name,type,option,help) { #name, type, option, help, &(uvar_ ## name) }



/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{ConfigFileHV}}
\newpage\input{ConfigFileC}
******************************************************* </lalLaTeX> */

/* Function prototypes */
void LALLoadConfigFile (LALStatus *stat, LALConfigData **cfgdata, const CHAR *fname);
void LALDestroyConfigData (LALStatus *stat, LALConfigData **cfgdata);

void LALReadConfigVariable (LALStatus *stat, void *varp, LALConfigData *cfgdata, LALConfigVar *param);
void LALReadConfigBOOLVariable (LALStatus *stat, BOOLEAN *varp, LALConfigData *cfgdata, const CHAR *varName);
void LALReadConfigINT2Variable (LALStatus *stat, INT2 *varp, LALConfigData *cfgdata, const CHAR *varName);
void LALReadConfigINT4Variable (LALStatus *stat, INT4 *varp, LALConfigData *cfgdata, const CHAR *varName);
void LALReadConfigREAL4Variable (LALStatus *stat, REAL4 *varp, LALConfigData *cfgdata, const CHAR *varName);
void LALReadConfigREAL8Variable (LALStatus *stat, REAL8 *varp, LALConfigData *cfgdata, const CHAR *varName);
void LALReadConfigSTRINGVariable (LALStatus *stat, CHAR **varp, LALConfigData *cfgdata, const CHAR *varName);
void LALReadConfigSTRINGNVariable (LALStatus *stat, CHARVector *varp, LALConfigData *cfgdata, const CHAR *varName);
void LALCheckConfigReadComplete (LALStatus *stat, LALConfigData *cfgdata, INT4 strictness);


void ReadCmdlineInput (LALStatus *stat, int argc, char *argv[], UserVariable *uvars);
void ReadCfgfileInput (LALStatus *stat, const CHAR *cfgfile, UserVariable *uvars );
void FreeUserVars (LALStatus *stat, UserVariable *uvars);
void GetUvarHelpString (LALStatus *stat, CHAR **helpstring, UserVariable *uvars);
void ReadUserInput (LALStatus *stat, int argc, char *argv[], UserVariable *uvars);

#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */



